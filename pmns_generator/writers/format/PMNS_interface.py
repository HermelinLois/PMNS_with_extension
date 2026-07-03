import pickle
from pathlib import Path
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Tuple
from sage.all import ceil, block_matrix, matrix, Integer, Rational

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent.parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

from ..format import format_element as format
from .format_element import _get_dimension, _ensure_list, _pad_coefficients
import config
from pmns_factory.core.operations.convertions_gestion import compute_conversion_tables, compute_nb_internal_reductions, gen_transition_matrix
from pmns_factory.core.operations.reductions.babai_reduction import gen_params_for_babai
from pmns_factory.core.operations.reductions.montgomery_reduction import gen_mn_reduction_matrix, search_polynomial_m
from pmns_factory.core.parameters.matrix_gestion import gen_overflow_matrix, gen_toeplitz_representation
from pmns_factory.core.parameters.params_gestion import search_memory_overhead
from config import STRUCT_SPARSE, STRUCT_GENERIC, ETYPE0, ETYPE1

# =====================================================================
#                           DATA STRUCTURES
# =====================================================================

@dataclass
class PMNSData:
    """
    PMNS data structure to hold all relevant parameters and matrices
    for a given PMNS configuration.
    """
    etype: int
    structure: str
    n: int
    n_limbs: int
    
    toeplitz_is_usable: bool
    elements_are_in_gamma_basis: bool
    lattice_is_double_sparse: bool
    babai_is_usable: bool
    
    pmns_params: Dict[str, Any] = field(default_factory=dict)
    matrices: Dict[str, Any] = field(default_factory=dict)
    reductions: Dict[str, Any] = field(default_factory=dict)
    conversions: Dict[str, Any] = field(default_factory=dict)


# =====================================================================
#                    PMNS CONTAINER STRUCTURE
# =====================================================================

class PMNSContainer:
    """
    PMNSContainer provides a unified interface to access PMNS parameters and matrices.
    It allows retrieval of parameters in their raw form or formatted for specific use cases.
    And it also provides methods to save and load the container state.
    """

    type0 = ETYPE0
    type1 = ETYPE1
    
    def __init__(self, data: PMNSData) -> None:
        self.data = data


    def _search_in_spaces(self, name: str) -> Any:
        """Search for a parameter in the PMNS data"""

        for space in [self.data.matrices, self.data.reductions, self.data.conversions, self.data.pmns_params, self.data.__dict__]:
            if name in space:
                return space[name]
            
        raise KeyError(f"Unknown parameter '{name}'.")
    

    def get(self, name: str, to_format: bool = False) -> Any:
        val = self._search_in_spaces(name)
                
        if not to_format:
            return val

        n_limbs = self.data.n_limbs
        dim = _get_dimension(val)
        
        # Convert to MPN format if the name is in the list of known MPN parameters
        mpn_names = ['p', 'gamma', 'T_mat', 'gamma_pows']
        if name in mpn_names:
            if dim >= 2: 
                return format.format_matrix_to_mpn(val, n_limbs)
            if hasattr(val, '_vector_'):
                return format.format_extension_field_to_mpn(val, n_limbs)
            return format.format_element_to_mpn(int(val), n_limbs)

        # Convert to int64 format if the name is in the list of known int64 parameters
        if dim == 3:
            return format.format_tensor_to_int64(val)
        if dim == 2:
            return format.format_matrix_to_int64(_ensure_list(val))
        if dim == 1:
            return format.format_list_to_int64(val)
        
        # If dim == 0, try to convert to int64, otherwise return as string
        try:
            return format.format_to_int64(int(val))
        except (TypeError, ValueError):
            return str(val)

    def save(self) -> None:
        """
        Save the PMNSContainer to a pickle file for later retrieval.
        The filename is constructed based on the PMNS parameters to ensure uniqueness.
        """
        path = config.OUTPUT_DIR_SAVES
        path.mkdir(parents=True, exist_ok=True)

        filename = f"pmns_p{self.get('p').nbits()}_k{self.get('k')}_E{self.data.etype}_{self.data.structure}.pkl"
        with open(path / filename, 'wb') as f:
            pickle.dump(self, f)
            

    @classmethod
    def load(cls, pbits: int, k: int, Etype: int, structure: str = "generic") -> "PMNSContainer":
        """
        Load a PMNSContainer from a pickle file based on the provided parameters.
        Args:
            pbits (int): The bit size of the prime p.
            k (int): The degree of the polynomial.
            Etype (int): The type of the extension field.
            structure (str): The structure type, either "generic", "specific" or "sparse".
        Returns:
            PMNSContainer: The loaded PMNSContainer object.
        """
        path = config.OUTPUT_DIR_SAVES / f"pmns_p{pbits}_k{k}_E{Etype}_{structure}.pkl"
        if not path.exists():
            raise FileNotFoundError(f"Error: The PMNSContainer file '{path}' does not exist. Please ensure the parameters are correct and the file has been generated.")
        
        with open(path, 'rb') as f:
            return pickle.load(f)


# =====================================================================
#                 PMNS PARAMETERS BUILDER STRUCTURE
# =====================================================================

class PMNSBuilder:
    """
    Compute and build PMNS parameters and matrices based on the provided PMNS data.
    This class encapsulates the logic for computing matrices, reduction parameters, and conversion parameters.
    """

    def __init__(self, Etype: int, pmns: dict, structure: str=STRUCT_GENERIC) -> None:
        self.etype = Etype
        self.pmns = pmns
        self.structure = structure

    def build(self) -> PMNSContainer:
        """
        Build the PMNSContainer by computing all necessary parameters and matrices.
        Returns:
            PMNSContainer: The container with all computed parameters and matrices.
        """
        n = self.pmns['E'].degree()
        n_limbs = ceil(self.pmns['p'].nbits() / self.pmns['phi_pow'])
        lattice_is_double_sparse = (self.structure == STRUCT_SPARSE)
        
        # Initialize the PMNSData object to hold all relevant parameters and matrices
        data = PMNSData(
            etype=self.etype,
            structure=self.structure,
            n=n,
            n_limbs=n_limbs,
            toeplitz_is_usable=(self.etype == 0),
            elements_are_in_gamma_basis=(self.pmns['gamma'] == self.pmns['gamma'].parent().gen()),
            lattice_is_double_sparse=lattice_is_double_sparse,
            babai_is_usable=not lattice_is_double_sparse,
            pmns_params=self.pmns
        )

        container = PMNSContainer(data)
        
        self._compute_matrix_parameters(container)
        self._compute_reduction_parameters(container)
        self._compute_conversions_parameters(container)
        
        return container

    def _get_reduction_bounds(self, container: PMNSContainer, theta_pow: int) -> Tuple[int, int, int]:
        """
        Compute the number of internal reductions required for exact, pseudo-fast, and fast conversions based on the PMNS parameters.
        Args:
            container (PMNSContainer): The PMNS container with necessary parameters.
            theta_pow (int): The power of theta used in the conversion calculations.
        Returns:
            Tuple[int, int, int]: A tuple containing the number of exact, pseudo-fast, and fast internal reductions required.
        """
        L = container.get('L')
        pol_e = container.get('E')
        phi = 2 ** container.get('phi_pow')
        rho = container.get('rho')
        k = container.get('k')
        n = container.data.n

        w = search_memory_overhead(pol_e)
        # Compute a factor of the coefficients bound to accurately estimate the number of internal reductions required for pseudo-fast conversions
        factor = (1 / 2 * L.norm(1)) ** (1 if container.data.lattice_is_double_sparse else 2)
        
        n_red_pseudo = compute_nb_internal_reductions(w * k * 2**theta_pow * factor, phi, rho, L)
        n_red_fast = compute_nb_internal_reductions(k * 2**theta_pow * 1 / 2 * L.norm(1), phi, rho, L)
        n_red_exact = compute_nb_internal_reductions((2 * rho) ** (n / k), phi, rho, L)

        # Ensure that the computed number of reductions meet the expected constraints for the PMNS implementation
        assert n_red_fast == 1, f"Program suppose n_red_fast == 1, got {n_red_fast}."
        assert n_red_pseudo <= 2, f"Program suppose n_red_pseudo <= 2, got {n_red_pseudo}." 

        return n_red_exact, n_red_pseudo, n_red_fast

    def _compute_matrix_parameters(self, container: PMNSContainer) -> None:
        """
        Compute the matrix parameters needed for the PMNS implementation.
        Args:
            container (PMNSContainer): The PMNS container with necessary parameters.
        Returns:
            None : Updates the container's matrices space with computed matrices.
        """
        E, L = container.get('E'), container.get('L')
        p, k = container.get('p'), container.get('k')
        phi_pow = container.get('phi_pow')
        gamma = container.get('gamma')
        n = container.data.n
        phi = 2 ** phi_pow
        
        m_space = container.data.matrices
        # compute the inverse of L modulo phi for use in reduction operations
        L_inv = -L.inverse() % phi
        m_space['L_inv'] = L_inv

        # compute the external reduction matrix for use in reduction operations
        # and store it as a square matrix of size n where the last row is a row 
        # of zeros
        external_mat_red = block_matrix([[gen_overflow_matrix(E)], [matrix(1, n)]])
        m_space['ext_red_mat'] = external_mat_red
        
        # compute the internal reduction matrices M and N for use in reduction operations
        # and various representations of these matrices for use in the PMNS implementation
        M = search_polynomial_m(L, k, p, gamma, E, sparse=container.data.lattice_is_double_sparse)
        M_mat, N_mat = gen_mn_reduction_matrix(M, E, phi)
        m_space['M_mat'] =  M_mat
        m_space['N_mat'] = N_mat

        if container.data.toeplitz_is_usable:
            m_space['toeplitz_mat_m'] = gen_toeplitz_representation(M_mat)
            m_space['toeplitz_mat_n'] = gen_toeplitz_representation(N_mat)

        T_mat = gen_transition_matrix(gamma, k)
        m_space['T_mat'] = T_mat

    def _compute_reduction_parameters(self, container: PMNSContainer) -> None:
        """
        Compute the reduction parameters needed for the PMNS implementation.
        Args:
            container (PMNSContainer): The PMNS container with necessary parameters.
        Returns:
            None: Updates the container's reductions space with computed parameters.
        """
        E = container.get('E')
        phi_pow = container.get('phi_pow')
        r_space = container.data.reductions
        
        # Compute different parameters for the PMNS implementation based on the structure of the PMNS
        if container.data.lattice_is_double_sparse:
            k = container.get('k')
            phi = 2 ** phi_pow
            gamma = container.get('gamma')

            # compute specific parameters for the double sparse structure to optimize the reduction operations
            #  and store them in minimal form 
            lamb = -E[0]
            lamb_inv_mod = pow(lamb, -1, phi)
            lambda_inv_mod = lamb_inv_mod if lamb_inv_mod <= phi // 2 else lamb_inv_mod - phi

            gam_pow = int(pow(gamma, k)) % phi
            gamma_pow_k_mod = gam_pow if gam_pow <= phi // 2 else gam_pow - phi

            r_space['lambda'] = lamb
            r_space['lambda_inv_mod'] = lambda_inv_mod
            r_space['gamma_pow_k_mod'] = gamma_pow_k_mod

            prod = (gamma_pow_k_mod * lambda_inv_mod) % phi
            r_space['gamma_pow_n_lambdda_mod'] = prod if prod <= phi // 2 else prod - phi
        
        if container.data.babai_is_usable:
            L = container.get('L')
            rho = container.get('rho')

            # compute parameters for the Babai reduction method to optimize the reduction operations
            # if h1 or h2 is nul during computation, it means that the Babai reduction method is not usable for this PMNS configuration
            result = gen_params_for_babai(L, phi_pow, rho, E)
            if result is None:
                container.data.babai_is_usable = False
                return
            
            h1, h2, L_inv_babai = result
            r_space['h1'] = h1
            r_space['h2'] = h2
            r_space['L_inv_babai'] = L_inv_babai

    def _compute_conversions_parameters(self, container: PMNSContainer) -> None:
        """ 
        Compute the conversion parameters needed for the PMNS implementation.
        Args:
            container (PMNSContainer): The PMNS container with necessary parameters.
        Returns:
            None: Updates the container's conversions space with computed parameters.
        """
        c_space = container.data.conversions
        n = container.data.n
        p, k = container.get('p'), container.get('k')
        
        theta_pow = ceil(p.nbits() * k / n)
        n_pols = ceil(n / k)

        n_red_exact, n_red_pseudo, n_red_fast = self._get_reduction_bounds(container, theta_pow)

        c_space['theta_pow'] = theta_pow     
        c_space['n_pols'] = n_pols
        c_space['n_red_exact'] = n_red_exact
        c_space['n_red_pseudo'] = n_red_pseudo
        c_space['n_red_fast'] = n_red_fast
        
        # compute the conversion tables for fast, pseudo fast and exact conversions
        c_space['fast_pols'] = _pad_coefficients(compute_conversion_tables(container, theta_pow, 1, n_pols, over_field=True), n)

        i_pols_full = _pad_coefficients(compute_conversion_tables(container, theta_pow, n_red_pseudo, n_pols, over_field=False), n)
        c_space['int_pols'] = i_pols_full[0]

        z_pols_full = _pad_coefficients(compute_conversion_tables(container, 0, 0, k, over_field=True), n)
        c_space['z_pols'] = [z_pols_full[i][0] for i in range(len(z_pols_full))]

        # compute the powers of gamma for use to convert element from PMNS representation to the extension field representation
        gamma = container.get('gamma')       
        gamma_pows = []
        for i in range(n):
            gamma_pows.append(list((gamma ** i)._vector_()))
        c_space['gamma_pows'] = gamma_pows