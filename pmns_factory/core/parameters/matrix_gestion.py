# ==================================================
# matrix_gestion.py
# Generic functions around matrix used for PMNS construction
# and C implementation
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing, Integer
from pmns_factory.core.math_utils import square_and_multiply

PR = PolynomialRing(ZZ, "X")
X = PR.gen()

def gen_null_base_general(k:int, p:int, pol_e, gamma):
    """
    Generate a reduced base of null polynomial for gamma.
    
    Args:
        p (int): prime used to construct the extension field
        k (int): extension degree of the field
        pol_e (Polynomial): external reduction polynomial
        gamma: element used to construct PMNS
    
    Returns:
        matrix (LLL) : reduced base where each row evaluates to zero at gamma
    """

    n = pol_e.degree()
    base = matrix(ZZ, n, n, 0)
    pol = PR(gamma.minpoly())

    # fill the diagonal
    for i in range(k):
        base[i, i] = p
        
    for i in range(k, n):
        vect = (pol * X**(i-k)).list()
        base[i] = vect + [0] * (n - len(vect))
        
    # LLL reduction
    return base


def gen_null_base_sparse(k:int, pol_e, gamma):
    """
    Generate a reduced base of null polynomial for gamma.
    Approach based on: An Alternative Approach for SIDH Arithmetic (C.Bouvier & L.Imbert)
    Here, delta is gamma^k and is suposed to be an integer.
    
    Args:
        p (int): prime used to construct the extension field
        k (int): extension degree of the field
        n (int): degree of pol_e used for polynomial reduction in PMNS
            gamma: element used to construct PMNS (gamma^k is integer)
    
    Returns:
        matrix (LLL) : reduced base where each row evaluates to zero at gamma
    """
    n = pol_e.degree()
    lambd = pol_e[0]

    delta = Integer(square_and_multiply(gamma, k))
    base = matrix(ZZ, n, n, 0)

    Gamma = X**k - delta
    M = delta * X**(n-k) + lambd

    for i in range(k):
        vect = (M * X**i).list()
        base[i] = vect + [0] * (n - len(vect))
        
    for i in range(n-k):
        vect = (Gamma * X**i).list()
        base[i+k] = vect + [0] * (n - len(vect))

    return base

def gen_null_base(k:int, p:int, pol_e, gamma, sparse:bool=False):
    if sparse:
        return gen_null_base_sparse(k, pol_e, gamma)
    return gen_null_base_general(k, p, pol_e, gamma)


def gen_overflow_matrix(pol_e):
    """
    Create a matrix used to compute the coefficients after reduction modulo pol_e 
    for powers of X greater or equal than the degree of polynomial pol_e.

    Args:
        pol_e (Polynomial): polynomial used for external reduction in PMNS

    Returns:
        matrix (ZZ): matrix representing the reduction of X^(n+i) modulo pol_e
    """
    n = pol_e.degree()
    X = PR("X")

    matrix_coefficients = []
    for i in range(n-1):
        poly_mod = X**(n + i) % pol_e
        # Pad coefficients to length n
        coeffs = list(poly_mod) + [0] * (n - poly_mod.degree()-1)
        matrix_coefficients.append(coeffs)
    return matrix(ZZ, matrix_coefficients)


def gen_toeplitz_representation(M):
    """
    Generate the Toeplitz representation of a matrix M.
    The Toeplitz representation is a list which reprsenrta the first column and the first row of the matrix M.

    Args:
        M (matrix): the input matrix
    Returns:
        list: the Toeplitz representation of the matrix M
    """
    first_col = M.column(0)
    first_row = M.row(0)[1:]
    return list(first_col)[::-1] + list(first_row)