# ==================================================
# Objectif is to construct E = X^n - lambda
# with fast search of root by constructing a 
# specific structure of the extension field
# and sparse subLattice
# ==================================================
from pmns_factory.core.parameters.matrix_gestion import gen_null_base
from pmns_factory.core.parameters.params_gestion import cast_polynomial_to_minimal_representation, is_irreducible_over_Fp
from pmns_factory.core.parameters.params_gestion import search_minimal_degree as SMD
from sage.all import ceil, log, sqrt, next_probable_prime, PolynomialRing, ZZ, GF, Integer

PR = PolynomialRing(ZZ, "X")
X = PR.gen()


def gen_pol(n:int, lamb:int):
    """
    Generate the reduction polynomial E = X^n - lambda. 
    """
    return X**n - lamb


def compute_gammak_bounds(psize: int, k: int, n: int, phi_pow: int) -> tuple:
    """ 
    Using F. palma thesis, compute generation bound for gamma.
    Here we suppose that gamma^k is an integer.

    Args:        
        psize (int): bit-size of the prime p
        k (int): extension degree of the field
        n (int): degree of the reduction polynomial
        phi_pow (int): word size of the architecture

    Returns:        
        tuple: (low_bound, high_bound) for gamma^k
    """

    # As we search a prime of psize bits and we assume that E = X^n - lambda and E(gamma) = p,
    # we can bound gamma^k as follows:
    low_bound = ceil(2**((psize-1)* k/n))
    high_bound = int(2**(psize*k/n))

    # To have sparse element, we search gamma^k such that gamma^k is
    # a multiple of phi**1/2. Round the lower bound down and the upper bound
    # up to avoid collapsing both into the same multiple.
    divide_pow = phi_pow // 2
    min_pow_multiple = 2**divide_pow
    low_bound = (low_bound // min_pow_multiple) * min_pow_multiple
    if low_bound == 0:
        low_bound = min_pow_multiple

    high_bound = ceil(high_bound / min_pow_multiple)*min_pow_multiple

    assert low_bound < high_bound, f"No valid gamma^k found in the computed bounds.\nHere, the bounds are: {low_bound= }, {high_bound= }"
    return low_bound, high_bound



def compute_lambda_max(gammak:int, n:int, phi_pow:int) -> int:
    """
    Compute the maximum value of lambda for a given gamma^k, n and phi_pow.
    Args:
        gammak (int): the value of gamma^k
        n (int): the degree of the reduction polynomial
        phi_pow (int): the power of phi used for generating the sparse reduction polynomial
    Returns:
        int: the maximum value of lambda
    """

    # As we set that E = X^n - lambda and E(gamma) = p, we can compute maximum value of lambda.
    # this value depend on the basis. Here we construct a basis inspired by Bouvier and Imbert approach.
    phi = 2**phi_pow
    part = sqrt(1 + (2*phi - 2*gammak + 2)/((n-1) * (gammak - 2)**2) + 1/((n-1) * (gammak - 2))**2) - 1
    max_lambd = ceil((gammak -2)/2 * part - 1/(2*n - 2))

    return max_lambd


def increase_gammak(gammak:None|int, low_bound:int, high_bound:int, phi_pow:int) -> int:
    """
    Increase gammak by a multiple of phi**1/2. If gammak is None, return low_bound.
    Args:
        gammak (None|int): the current value of gamma^k
        low_bound (int): the lower bound for gamma^k
        high_bound (int): the upper bound for gamma^k
        phi_pow (int): the power of phi used for generating the sparse reduction polynomial
    Returns:
        int: the new value of gamma^k
    """

    if gammak is None:
        return low_bound

    divide_pow = phi_pow // 2
    min_pow_multiple = 2**divide_pow

    # we increase by a multiple of phi**1/2 to make gammak a multiple of phi**1/2 and construct
    # a sparse reduction matrix
    n_gammak =  gammak + min_pow_multiple if gammak is not None else min_pow_multiple

    assert n_gammak < high_bound, f"No valid gamma^k found in the computed bounds.\nAll possible gamma^k have been tested."
    return n_gammak
    
    
def construct_specific_field(p:int, k:int, gammak:int, name:str='z'):
    """
    Construct an extension field using X^k - gamma^k as the irreducible polynomial.

    Since gamma satisfying a PMNS construction is a field generator,
    we can create a specific field using the irreducible polynomial
    X^k - gamma^k. By doing so, a trivial root is the generator 'z',
    which helps reduce the cost of exact and pseudo-fast conversions
    to PMNS.

    Args:
        p (Integer): Prime used to construct the extension field.
        k (int): Extension degree of the field.
        gammak (int): Represents gamma^k.
        name (str, optional): Name used to represent the root of the
            irreducible polynomial in the extension field.

    Returns:
        An extension field of degree k using X^k - gamma^k as the
        irreducible polynomial.
    """
    
    irred_pol = gen_pol(k, gammak)
    return GF(p**k, name=name, modulus=irred_pol)
    


def search_minimal_degree(psize: int, k: int, phi_pow: int) -> int:
    """
    Function that search minimal degree from wich we can possibly construct a PMNS from our initial coefficients

    Args:
        p (int): prime use to construct extension field
        k (int): extension degree
        phi_pow (int): word size of the architecture

    Returns:
        int: degree from wich we can possibly construct a PMNS
    """

    init_polynomial = lambda n : gen_pol(n, 1)
    size_holder = Integer(2**psize - 1)
    n = SMD(size_holder, k, phi_pow, init_polynomial)
    
    return k * ceil(n/k)


def gen_parameters(psize:int, k:int, max_wanted_lamb:int=None, phi_pow:int=64, n:int=None, name:str="z"):
    """
    Generate PMNS parameters of type 0 with sparse reduction polynomial.
    Here gammak represent gamma^k and we search for a gamma such that gamma^k 
    is a multiple of 2**32 to have a sparse reduction polynomial.
    
    Args:
        psize (int): bit-size of the prime p
        k (int): extension degree of the field
        n (int): degree of the reduction polynomial
        max_wanted_lambd (int): the maximum value of lambda to consider
        phi_pow (int): power of phi to use for generating the sparse reduction polynomial
        
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
    if max_wanted_lamb is None:
        max_wanted_lamb = 2**(phi_pow - 1) - 1  
          
    if n is None:
        n = search_minimal_degree(psize, k, phi_pow)
    low_bound, high_bound = compute_gammak_bounds(psize, k, n, phi_pow)

    # test if construction works
    assert compute_lambda_max(low_bound, n, phi_pow) > 0, f"Construction does not work for the computed bounds. Here, with {low_bound= }, {n= }, {phi_pow= } :\n{compute_lambda_max(low_bound, n, phi_pow)= }"

    # start research
    gammak = None
    not_found = True
    iteration = 0
    while not_found :
        iteration += 1
        gammak = increase_gammak(gammak, low_bound, high_bound, phi_pow)
        lamb_max = min(compute_lambda_max(gammak, n, phi_pow), max_wanted_lamb)

        if lamb_max < 1:
            continue

        # delta = compute_delta(n, lamb_max, gammak, phi_pow)
        gamma_n = gammak**(n//k)
        p = next_probable_prime(gamma_n - lamb_max)
        curr_lamb = gamma_n - p
        
        if 0 < abs(curr_lamb) <= lamb_max and p.nbits() == psize and is_irreducible_over_Fp(gen_pol(k, gammak), p):
            not_found = False

    pol_e = gen_pol(n, curr_lamb)
    
    K = construct_specific_field(p, k, gammak, name=name)
    gamma = K.gen()
        
    L = gen_null_base(k, p, pol_e, gamma, sparse=True)
    rho = Integer(abs(gammak) + abs(curr_lamb) - 1)
    return {'rho': rho, 'gamma': gamma, 'phi_pow': phi_pow, 'L': L, 'E': pol_e, 'mod': cast_polynomial_to_minimal_representation(K.modulus(), p), 'p': p, 'k':k, 'it':iteration}