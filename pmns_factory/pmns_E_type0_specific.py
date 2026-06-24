# ==================================================
# Objectif is to construct E = X^n - lambda
# with fast search of root by constructing a 
# specific structure of the extension field
# ==================================================

from sage.all import PolynomialRing, ZZ, Integer, GF, random_prime, factor, gcd, primitive_root
from pmns_factory.core.parameters.params_gestion import search_minimal_degree as SMD, search_base_rho_and_gamma, search_memory_overhead, cast_polynomial_to_minimal_representation
from pmns_factory.core.parameters.roots_gestion import is_gamma_feasible, search_roots

PR = PolynomialRing(ZZ, "X")
X = PR.gen()
INIT_LAMB = 2

def gen_pol_e(n:int, lamb:int):
    return X**n - lamb

def construct_irreducible_polynomial(k, p):
    # check if k = 0 mod(4) and in that case check 
    # if p = 1 mod(4) 
    if k%4 == 0 and p%4 != 1:
        return None

    primes_k = [q for q, _ in factor(k)]

    # fast verification to see if it's impossible to 
    # construct an irreducible polynomial with p and k
    for q in primes_k:
        if (p - 1) % q != 0:
            return None

    # construct irreducible polynomial using generator of Fp
    b = primitive_root(p)
    return X**k - b


def increase_parameters(pol_e, p:int, k:int, phi:int) -> tuple:
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
    Search the minimal degree n of the reduction polynomial E = X^n - lambda, 
    such that we can construct a PMNS with p, k and phi_pow.

    Args:
        p (int): the prime
        k (int): the extension degree of the field
        phi_pow (int): the power of phi used for generating the sparse reduction polynomial
    Returns:
        int: the minimal degree n of the reduction polynomial
    """

    init_polynomial = lambda n : gen_pol_e(n, INIT_LAMB)
    n = SMD(p, k, phi_pow, init_polynomial)
    return n



def gen_parameters(psize:int, k:int, phi_pow:int=64, n:int=None, name:str="z") -> dict: 
    """
    Function use to generate PMNS parameters given the prime size, the extension degree and the word size parameter of the architecture
    This implementation use a specific construction of irreducible polynomial to construct extension field.
    here, this polynomial as X^k - a with a in Fp

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
    
    p = random_prime(2**psize, lbound=2**(psize-1))
    p = Integer(p)
    
    # this condition permit to know if we can construct irreducible polynomial to construct extension field
    assert gcd(k, p-1) == k, f"impossible to construct an irreducible polynomial over Z/pZ with {p=} and {k=}"

    if n is None:
        n = search_minimal_degree(p, k, phi_pow)
    lamb = INIT_LAMB
    phi = 2**phi_pow

    # extension field creation
    mod = construct_irreducible_polynomial(k, p)
    assert mod is not None, f"Impossible to construct specific extension field with {p= } and {k= }"
    K = GF(p**k, name=name, modulus=mod)

    parameters_not_found = True
    result = None
    iteration = 0
    while parameters_not_found:
        iteration +=1
        # construction of the polynomial pol_e
        pol_e = gen_pol_e(n, lamb)
        roots = search_roots(p, k, pol_e, K)

        if roots:
            # search suitable element to construct a PMNS using found roots
            result = search_base_rho_and_gamma(roots, k, p, phi, pol_e)
            parameters_not_found = (result is None)

        if parameters_not_found:
            lamb, n = increase_parameters(pol_e, p, k, phi)

    L, rho, gamma = result
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': cast_polynomial_to_minimal_representation(K.modulus(), p), 'p': p, 'k':k, 'it':iteration}
