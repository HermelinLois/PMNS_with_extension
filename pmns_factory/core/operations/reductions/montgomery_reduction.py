# ==================================================
# montgomery_reduction.py
# functions wich permit to apply Montgomery reduction 
# to a polynomial 
# ==================================================

from sage.all import matrix, ZZ, PolynomialRing, Integer, vector, GF, gcd, xgcd
from ...math_utils import square_and_multiply
from ...parameters.matrix_gestion import gen_null_base

PR = PolynomialRing(ZZ, "X")
PR2 = PolynomialRing(GF(2), "X")


def search_m_with_even_degs(base, pol_e):
    """
    Optimized search for an invertible polynomial modulo pol_e
    by looking for elements whose degree-0 coefficient is odd.

    Args:
        base (matrix): matrix of null polynomials evaluated at gamma
        pol_e (Polynomial): polynomial used for external PMNS reduction

    Returns:
        Polynomial: polynomial invertible modulo pol_e and null over gamma
    """
    for row in base:
        if row[0] & 1:
            return PR(list(row)) % pol_e
    return None


def search_m_with_odd_deg(k: int, p: int, gamma, pol_e, sparse:bool =False):
    """
    General search for an invertible polynomial modulo pol_e.

    Args:
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma: root of the external reduction polynomial
        pol_e (Polynomial): external reduction polynomial
        sparse (bool, OPtional): is construct PMNS is sparse

    Returns:
        Polynomial: polynomial invertible modulo pol_e
    """
    n = pol_e.degree()
    base = gen_null_base(k, p, pol_e, gamma, sparse=sparse)

    for i in range(n):
        base[i] = list(
            map(
                lambda enum: int(enum[1]) + p * (int(enum[1]) & 1) * (enum[0] != i),
                enumerate(base[i]),
            )
        )

    reduced_base = base.LLL()
    reference_polynomial = PR2(pol_e)

    for linear_combination in range(1, 2 ** n):
        result = sum(reduced_base[i] for i in range(n) if (linear_combination >> i) & 1)
        polynomial = PR(list(result)) % pol_e
        bin_polynomial = PR2(polynomial)

        if bin_polynomial != 1 and gcd(bin_polynomial, reference_polynomial) == 1:
            return polynomial

    return None


def search_polynomial_m(base, k: int, p: int, gamma, pol_e, sparse:bool =False):
    """
    Search for an invertible element modulo pol_e with optimization
    when the coefficients of pol_e are all even in GF(2).

    Args:
        base (matrix): matrix of null polynomials when evaluated over gamma
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma: root of the external reduction polynomial
        pol_e (Polynomial): external reduction polynomial
        sparse (bool, OPtional): is construct PMNS is sparse

    Returns:
        Polynomial: polynomial invertible modulo pol_e
    """
    no_optimisation = any(coef & 1 for coef in pol_e)
    if no_optimisation:
        return search_m_with_odd_deg(k, p, gamma, pol_e, sparse=sparse)
    return search_m_with_even_degs(base, pol_e)


def search_m_and_n(k: int, p: int, gamma, base, pol_e, sparse:bool=False, phi: int = 2 ** 64):
    """
    Retrieve M invertible modulo pol_e and N = -M^(-1) mod phi.

    Args:
        k (int): extension degree
        p (int): prime used to construct the extension field
        gamma (extension field element): root of E suitable for PMNS construction
        base (matrix): reduced base of null polynomial over gamma
        pol_e (Polynomial): polynomial used for external reduction in PMNS
        phi (int, Optional): word size bound. Equal to 2**64 by default
        sparse (bool, OPtional): is construct PMNS is sparse

    Returns:
        Polynomial: M, a polynomial null over gamma and invertible modulo pol_e
        Polynomial: N, a polynomial such that N = -M^(-1) mod phi
    """
    M = search_polynomial_m(base, k, p, gamma, pol_e, sparse=sparse)

    assert M(gamma) == 0, "problem occuring with the base. Polynomial M isn't inverssible"

    d, u, _ = xgcd(M, pol_e)
    d = int(d)
    d_inv = pow(d, -1, phi)
    N = PR((-d_inv * u) % phi)

    return M, N

def gen_mn_reduction_matrix(M, E, phi: int):
    """
    Generate matrices M and N for external reduction in PMNS.

    Args:
        M (Polynomial): polynomial null in the chosen root of E
        E (Polynomial): polynomial used for external reduction
        phi (int): word size in bits (used for modulo)

    Returns:
        mat_m (matrix): matrix representing M for Montgomery reduction
        mat_n (matrix): matrix representing N = -M^(-1) modulo phi
    """
    n = E.degree()
    X = PR("X")

    matrix_coefficients = []
    for i in range(n):
        # Compute (M * X^i) % E
        poly_mod = (M * square_and_multiply(X,i)) % E
        coeffs = list(poly_mod) + [0] * (n - len(list(poly_mod)))
        matrix_coefficients.append(coeffs)

    mat_m = matrix(ZZ, matrix_coefficients)
    mat_n = -mat_m.inverse() % phi

    return mat_m, mat_n


def montgomery_reduction_polynomial(pol_p, M, N, E, gamma, phi:int =2**64):
    """
    Reduction of a polynomial with Montgomery reduction

    Args:
        pol_p (Polynomial | matrix): elment wich need to be reduced 
        M (Polynomial | matrix): element wich represent M such that M(gamma)=0
        N (Polynomial | matrix): element wich represent N = -M^-1
        E (Polynomial): polynomial for external reduction
        phi (int, opt): represent bit word size use by the architecture. usualy a power of 2 and is 2^64 by default
        gamma (extension field element): root of E

    Returns:
        Polynomial : reduction of pol_p which still represent the same element
    """
    phi = Integer(phi)

    Q = ((pol_p * N) % E) % phi
    Q = PR([elmt - phi * (elmt > phi//2) for elmt in Q])

    T = (Q * M) % E 
    reduction = (pol_p + T) // phi

    assert reduction(gamma) * phi == pol_p(gamma), f"Error occurring during reduction. Please check parameters"
    
    return reduction


def montgomery_reduction_lattice(pol_p, mat_sublattice, mat_sublattice_inv, phi:int = 2**64):
    """
    Fast reduction of a polynomial with Montgomery reduction using precomputed matrices of the sublattice.
    Args:
        pol_p (Polynomial): element to be reduced
        mat_sublattice (matrix): matrix representing the sublattice for reduction
        mat_sublattice_inv (matrix): minus the inverse of the sublattice matrix modulo phi
        phi (int, opt): represent bit word size use by the architecture. usualy a power of 2 and is 2^64 by default
    Returns:
        Polynomial: reduction of pol_p which still represent the same element
    """
    n = mat_sublattice.nrows()
    phi = Integer(phi)

    v_p = vector(ZZ, pol_p.list() + [0] * (n - len(pol_p.list())))
    q_raw = (v_p * mat_sublattice_inv) % phi
    q_centered = vector(ZZ, [c - phi *(c > phi // 2) for c in q_raw])
    T = PR((q_centered * mat_sublattice).list())

    return (pol_p + T)//phi