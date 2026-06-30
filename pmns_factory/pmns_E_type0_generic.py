# ==================================================
# pmns_E_type0.py
# Functions to construct a PMNS in extension field GF(p^k)
# with E = X^n - lamb as the external reduction polynomial
# ==================================================

from sage.all import PolynomialRing, ZZ, random_prime, GF
from pmns_factory.core.parameters.params_gestion import search_minimal_degree as SMD, search_base_rho_and_gamma, search_memory_overhead, cast_polynomial_to_minimal_representation
from pmns_factory.core.parameters.roots_gestion import search_roots

PR = PolynomialRing(ZZ, "X")
X = PR("X")

INIT_LAMB = 2

def gen_pol_e(n:int, lamb:int):
    """
    Function that create the polynomial for external reduction in PMNS

    Args:
        n (int): degree of pol_e
        lamb (int): value of the coefficient of degree 0

    Returns:
        Polynomial: return the external reduction polynomial
    """
    return X**n - lamb

def increase_parameters(pol_e, p:int, k:int, phi:int) -> tuple:
    """
    Function that increament parameters of the PMNS with condition to increase lamb or n

    Args:
        pol_e (Polynomial): external reduction use in the PMNS
        p (int): prime use to contruct extension field
        k (int): extension degree
        phi (int): word size interger of the architecture (ie 2**(word size))

    Returns:
        lamb (int): either lamb +1 or INIT_LAMB=2
        n (int): either n or n+k
    """
    n = pol_e.degree()
    lamb = - pol_e[0]
    n_lamb = lamb + 1

    E = gen_pol_e(n, n_lamb)

    # compute overhead based on: "PMNS for efficient arithmetic and small memory cost" 
    # (J. Robert, P. Véron, F. Dosso, 2022)
    w = search_memory_overhead(E)
    if 2 * w * p**(k/n) >= phi:
        return INIT_LAMB, n + 1
    return n_lamb, n


def search_minimal_degree(p: int, k: int, phi_pow: int) -> int:
    """
    Function that search minimal degree from wich we can possibly construct a PMNS from our initial coefficients

    Args:
        p (int): prime use to construct extension field
        k (int): extension degree
        phi_pow (int): word size of the architecture

    Returns:
        int: degree from wich we can possibly construct a PMNS
    """

    init_polynomial = lambda n : gen_pol_e(n, INIT_LAMB)
    n = SMD(p, k, phi_pow, init_polynomial)
    
    return n


def gen_parameters(psize:int, k:int, phi_pow:int=64, n:int=None, name:str="z") -> dict:
    """
    Function use to generate PMNS parameters given the prime size psize, the extension degree and the word size parameter of the architecture

    Args:
        psize (int): prime size use to generate prime
        k (int): extension degree
        name (str): name given to element to extension field element
        phi_pow (int, optional): word size use by the arcitecture (usully 2**word size). Defaults to 64.

    Returns:
        p (Integer): prime used to construct pmns
        rho (int): upper bound of the coefficient for polynomial is PMNS
        gamma (extension field element): root of the external reduction polynomial and such that gamma^k 
            is an interger and consecutiv power of gamma are a base of the polynomial extension
        phi_pow (int): word size architectur from wich we construct the PMNS
        L (matrix): reduce lattice matrix of null polynomial over gamma
        E (Polynomial): external reduction polynomial use by the PMNS
        mod (Polynomial): Polynomial used to construct the extension field
    """
    assert psize >= phi_pow, f"construction only works if the number of bits in prime (here {p=}) is greater or equal to {phi_pow=}"
    
    
    # initailisation of elements
    p = random_prime(2**psize, lbound=2**(psize-1))
    if n is None:
        n = search_minimal_degree(p, k, phi_pow)

    lamb = INIT_LAMB
    phi = 2**phi_pow

    # extension field creation
    K = GF(p**k, name)
    
    parameters_not_found = True
    result = None
    iteration = 0
    while parameters_not_found:
        iteration += 1
        # construction of the polynomial pol_e
        pol_e = gen_pol_e(n, lamb)
        print(pol_e)
        roots = search_roots(K, pol_e)

        if roots:
            # search suitable element to construct a PMNS using found roots
            result = search_base_rho_and_gamma(roots, k, p, phi, pol_e)
            parameters_not_found = (result is None)

        if parameters_not_found:
            lamb, n = increase_parameters(pol_e, p, k, phi)

    
    L, rho, gamma = result
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': cast_polynomial_to_minimal_representation(K.modulus(), p), 'p': p, 'k':k, 'it':iteration}