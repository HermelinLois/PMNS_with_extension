import pickle
from pathlib import Path
import sys
from sage.all import ceil, block_matrix, matrix, Integer

# Gestion des imports relatifs/système
CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent.parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

from ..format import format_element as format
import config
from pmns_factory.core.operations.convertions_gestion import compute_conversion_tables, compute_nb_internal_reductions, gen_transition_matrix
from pmns_factory.core.operations.reductions.babai_reduction import gen_params_for_babai
from pmns_factory.core.operations.reductions.montgomery_reduction import gen_mn_reduction_matrix, search_polynomial_m
from pmns_factory.core.parameters.matrix_gestion import gen_overflow_matrix, gen_toeplitz_representation
from pmns_factory.core.parameters.params_gestion import search_memory_overhead
from config import STRUCT_SPARSE, STRUCT_GENERIC, STRUCT_SPECIFIC

def _ensure_list(obj):
        if hasattr(obj, 'rows'):
            return [list(row) for row in obj.rows()]
        return obj

def _pad_coefficients(poly, n):
        if hasattr(poly, 'list'):
            poly = list(poly)

        if isinstance(poly, list) and len(poly) > 0 and isinstance(poly[0], list):
            return [_pad_coefficients(item, n) for item in poly]
        
        if isinstance(poly, list):
            return poly + [0] * (n - len(poly))

        return poly

    
class PMNSContainer:
    struct_generic = STRUCT_GENERIC
    struct_specific = STRUCT_SPECIFIC
    struct_sparse = STRUCT_SPARSE
    
    def __init__(self, Etype:int, pmns: dict, structure:str=struct_generic):
        self.params = {**pmns}
        self.params['type'] = Etype
        self.params['struct'] = structure

        self.params['n'] = pmns['E'].degree()
        self.params['n_limbs'] = ceil(pmns['p'].nbits() / pmns['phi_pow'])

        self.params['is_toeplitz_usable'] = (Etype == 0)
        self.params['is_elements_in_gamma_basis'] = pmns['gamma'] == pmns['gamma'].parent().gen()
        self.params['is_double_sparse'] = structure == self.struct_sparse
        self.params['is_babai_usable'] = not self.params['is_double_sparse']        
        
        self._build()


    def _build(self):
        self.update_matrix_parameters()
        self.update_reduction_parameters()
        self.update_conversions_parameters()

    def _get_reduction_parameters(self, theta_pow):
        # without external reduction polynomial product is bounded by k * theta * (1/2 * L.norm(1))**2) or k * theta * 1/2 * L.norm(1) if 
        # used structure is sparse with our implementation
        # this is due to the fact that precomputed polynomial have coefficients bounded by (1/2 * L.norm(1))**2) following F. PALMA thesis
        # by construction, we add at most a coefficient of size theta = 2**theta_pow
        # k factor is comming from the polynomials addition (cf. montgomery_pseudo_fast_conversion)
        # as we realise an external reduction to ensure a representation bounded by rho, w parameters appears in the bound
        
        L = self.params['L_origin']
        pol_e = self.params['E']
        phi = 2**self.params['phi_pow']
        rho = self.params['rho']
        k = self.params['k']
        n = self.params['n']

        w = search_memory_overhead(pol_e)
            
        factor = (1/2 * L.norm(1))**(1 if self.params['struct'] == PMNSContainer.struct_sparse else 2)
        n_red_pseudo = compute_nb_internal_reductions(w * k * 2**theta_pow * factor , phi, rho, L)
        n_red_fast = compute_nb_internal_reductions(k * 2**theta_pow * 1/2 * L.norm(1), phi, rho, L)
        n_red_extact = compute_nb_internal_reductions((2*rho)**(n/k), phi, rho, L)

        assert n_red_fast == 1, f"Program suppose n_red_fast == 1, got {n_red_fast}."
        assert n_red_pseudo <= 2, f"Program suppose n_red_pseudo <= 2, got {n_red_pseudo}. Change function use in conversion to switch from mpn use to mpz use." 

        return n_red_extact, n_red_pseudo, n_red_fast

        

    def update_matrix_parameters(self):
        E, L = self.params['E'], self.params['L']
        p, k = self.params['p'], self.params['k']
        phi_pow = self.params['phi_pow']
        gamma = self.params['gamma']
        n = self.params['n']
        phi = 2**phi_pow
        
        self.params['L'] = _ensure_list(L)
        self.params['L_origin'] = L

        L_inv = -L.inverse() % phi
        self.params['L_inv'] = _ensure_list(L_inv)
        self.params['L_inv_origin'] = L_inv

        external_mat_red = block_matrix([[gen_overflow_matrix(E)], [matrix(1, n)]])
        self.params['ext_red_mat'] = _ensure_list(external_mat_red)
        
        M = search_polynomial_m(L, k, p, gamma, E, sparse=self.params['struct']==STRUCT_SPARSE)
        M_mat, N_mat = gen_mn_reduction_matrix(M, E, phi)
        self.params['M_mat'] = _ensure_list(M_mat)
        self.params['M_mat_origin'] = M_mat
        self.params['N_mat'] = _ensure_list(N_mat)
        self.params['N_mat_origin'] = N_mat
        self.params['toeplitz_mat_m'] = _ensure_list(gen_toeplitz_representation(M_mat))
        self.params['toeplitz_mat_n'] = _ensure_list(gen_toeplitz_representation(N_mat))

        T_mat = gen_transition_matrix(gamma, k)
        self.params['T_mat'] = _ensure_list(T_mat)
        self.params['T_mat_origin'] = T_mat



    def update_reduction_parameters(self):
        E = self.params['E']
        phi_pow = self.params['phi_pow']
        
        if self.params['is_double_sparse']:
            # compute parameters for sparse reduction with External reduction of form X^n - lambda with lambda = -E[0]
            # and with gamma such that gamma^k is an integer
            k = self.params['k']
            phi = 2**phi_pow
            gamma = self.params['gamma']

            lamb = - E[0]
            lamb_inv_mod = pow(lamb, -1, phi)
            lambda_inv_mod = lamb_inv_mod if lamb_inv_mod <= phi//2 else lamb_inv_mod - phi

            gam_pow = int(pow(gamma,k)) % phi
            gamma_pow_k_mod = gam_pow if gam_pow <= phi//2 else gam_pow - phi

            self.params['lambda'] = lamb
            self.params['lambda_inv_mod'] = lambda_inv_mod
            self.params['gamma_pow_k_mod'] = gamma_pow_k_mod

            prod = (gamma_pow_k_mod * lambda_inv_mod) % phi
            self.params['gamma_pow_n_lambdda_mod'] = prod if prod <= phi//2 else prod - phi
        
        if self.params['is_babai_usable']:
            # compute parameters for Babai reduction
            L = self.params['L_origin']
            rho = self.params['rho']

            result = gen_params_for_babai(L, phi_pow, rho, E)
            if result is None:
                self.params['is_babai_usable'] = False
                return
            
            h1, h2, L_inv_babai = result
            self.params['h1'] = h1
            self.params['h2'] = h2
            self.params['L_inv_babai'] = _ensure_list(L_inv_babai)
            self.params['L_inv_babai_origin'] = L_inv_babai


    def update_conversions_parameters(self):
        tpow = "theta_pow"
        nb_pols = 'n_pols'
        nb_exact = 'n_red_exact'
        nb_pseudo = 'n_red_pseudo'
        nb_fast = 'n_red_fast'
        fpols = 'fast_pols'
        ipols = 'int_pols'
        zpols = 'z_pols'
        
        required_keys = [tpow, nb_pols, nb_exact, nb_pseudo, nb_fast, fpols, ipols, zpols]
        if all(key in self.params for key in required_keys):
            return

        n = self.params['n']
        p, k = self.params['p'], self.params['k']
        
        theta_pow = ceil(p.nbits() * k / n)
        n_pols = ceil(n/k)

        n_red_exact, n_red_pseudo, n_red_fast = self._get_reduction_parameters(theta_pow)

        self.params[tpow] = theta_pow     
        self.params[nb_pols] = n_pols
        self.params[nb_exact] = n_red_exact
        self.params[nb_pseudo] = n_red_pseudo
        self.params[nb_fast] = n_red_fast
        
        self.params[fpols] = _pad_coefficients(compute_conversion_tables(self, theta_pow, 1, n_pols, over_field=True), n)

        i_pols_full = _pad_coefficients(compute_conversion_tables(self, theta_pow, n_red_pseudo, n_pols, over_field=False), n)
        self.params[ipols] = i_pols_full[0]

        z_pols_full = _pad_coefficients(compute_conversion_tables(self, 0, 0, k, over_field=True), n)
        self.params[zpols] = [z_pols_full[i][0] for i in range(len(z_pols_full))]

        gamma = self.params['gamma']       
        gamma_pows = []
        for i in range(n):
            vec = list((gamma**i)._vector_())
            gamma_pows.append([int(c) for c in vec])
        self.params['gamma_pows'] = gamma_pows


    def save(self):
        path = config.OUTPUT_DIR_SAVES
        path.mkdir(parents=True, exist_ok=True)
        
        filename = f"pmns_p{self.params['p'].nbits()}_k{self.params['k']}_E{self.params['type']}_{self.params['struct']}.pkl"
        full_path = path / filename
        with open(full_path, 'wb') as f:
            pickle.dump(self, f)
            
            
    @classmethod
    def load(cls, pbits:int, k:int, Etype:int, structure:str="generic"):
        path = config.OUTPUT_DIR_SAVES / f"pmns_p{pbits}_k{k}_E{Etype}_{structure}.pkl"
        
        if not path.exists():
            raise FileNotFoundError(f"Path to file {path} doesn't seems to exist.")

        with open(path, 'rb') as f:
            return pickle.load(f)
        
        
    def get(self, name: str, to_format: bool = False):
        val = self.params.get(name)
        
        if val is None:
            raise KeyError(f"Unknow parameter '{name}'.")

        if not to_format:
            return val

        n_limbs = self.params.get('n_limbs')
        
        is_list = isinstance(val, (list, tuple))
        is_nested = is_list and len(val) > 0 and isinstance(val[0], (list, tuple))
        is_tensor = is_nested and len(val[0]) > 0 and isinstance(val[0][0], (list, tuple))
        
        mpn_names = ['p', 'gamma', 'T_mat', 'gamma_pows']
        if name in mpn_names:
            if is_nested: 
                return format.format_matrix_to_mpn(val, n_limbs)
            if hasattr(val, '_vector_'):
                return format.format_extension_field_to_mpn(val, n_limbs)
            return format.format_element_to_mpn(int(val), n_limbs)

        if is_tensor:
            return format.format_tensor_to_int64(val)
        if is_nested:
            return format.format_matrix_to_int64(val)
        if is_list:
            return format.format_list_to_int64(val)
            
        try:
            return format.format_to_int64(int(val))
        except (TypeError, ValueError):
            return str(val)