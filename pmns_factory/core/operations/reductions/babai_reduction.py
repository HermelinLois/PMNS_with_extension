# ==================================================
# babai_reduction.py
# Functions to apply Babai lattice reduction for PMNS.
# Provides three variants:
#   > nearest plane reduction
#   > rounding reduction with unlimited precision
#   > rounding reduction with limited precision
# ==================================================

from sage.all import vector, ZZ, floor, matrix, RR, Integer, PolynomialRing, log
from ...parameters.params_gestion import search_memory_overhead

PR = PolynomialRing(ZZ, "X")

def babai_nearest_plane_reduction(L, pol_p, gamma, rho):
    """
    Reduce a polynomial using Babai's nearest plane algorithm.
    Code adapted from N. Meloni's Sage implementation.

    Args:
        pol_p (Polynomial): polynomial to reduce.
        L (matrix): lattice basis of null polynomials for the PMNS.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """

    n = L.nrows()

    # Gram-Schmidt orthogonalization of the lattice basis
    gram, _ = L.gram_schmidt()

    # normalized Gram vectors
    G = [gram[i] / (gram[i].norm()**2) for i in range(n)]

    s = vector(pol_p.list() + [0]*(n - pol_p.degree() - 1))

    # reduction
    for idx in range(n-1, -1, -1):
        coef = round(s.dot_product(G[idx]))
        s -= coef * L[idx]

    reduction =  PR(list(s))

    assert reduction(gamma) == pol_p(gamma), "Error during reduction. Reduction doesn't represent the same element"
    assert all(abs(c) < rho for c in reduction), "Error during reduction. Reduction coefficient are greater than rho"

    return reduction


def babai_rounding_unlimited_reduction(L, pol_p, gamma, rho):
    """
    Reduce a polynomial using Babai rounding with unlimited precision.
    Code adapted from N. Meloni's Sage implementation.

    Args:
        pol_p (Polynomial): polynomial to reduce.
        L (matrix): lattice basis of null polynomials for the PMNS.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """

    n = L.nrows()
    s = vector(pol_p.list() + [0]*(n - pol_p.degree() - 1))
    l_inv = L.inverse()
    coefs = vector(map(round, s * l_inv))

    reduction = PR(list(s - coefs * L))

    assert reduction(gamma) == pol_p(gamma), "Error during reduction. Reduction doesn't represent the same element"
    assert all(abs(c) < rho for c in reduction), "Error during reduction. Reduction coefficient are greater than rho"

    return reduction

def babai_rounding_limited_reduction(pol_p, container):
    """
    Reduce a polynomial using Babai rounding with limited precision.
    Code adapted from N. Meloni's Sage implementation.

    This version is designed for architectures with fixed precision
    (e.g., C implementations).

    Args:
        pol_p (Polynomial): polynomial to reduce.
        container (PMNSContainer) : object with:
            h1 (int): scaling parameter for Babai rounding.
            h2 (int): secondary scaling parameter for coefficient truncation.
            L (matrix): lattice basis of null polynomials for the PMNS.
            L_inv_babai (matrix): precomputed scaled inverse of the lattice basis.

    Returns:
        Polynomial: reduced polynomial representing the same element.
    """
    n = container.get('n')
    h1 = container.get('h1')
    h2 = container.get('h2')
    L = container.get('L_origin')
    L_inv_babai = container.get('L_inv_babai_origin')
    gamma = container.get('gamma')
    
    resize_vect = pol_p.list() + [0]*(n - len(pol_p.list()))
    v = vector([floor(x / 2**h2) for x in resize_vect])

    s = v * L_inv_babai
    s = vector([floor(x / 2**(h1 - h2)) for x in s])

    reduction = PR(list(vector(resize_vect) - s * L))
    
    assert reduction(gamma) == pol_p(gamma)

    return reduction


def gen_params_for_babai(L, phi_pow:int, rho:int, pol_e):
    """
    Generate parameters for babi limited implementation.

    Args:
        L (matrix): lattice base of null polynomial in gamma
        phi_pow (int): represent bit word size use by the architecture
        rho (int): limit for absolute coefficients values
        pol_e (Polynomial): polynomial use in external reduction

    Returns:
        h1 (int): scaling parameter for Babai rounding.
        h2 (int): secondary scaling parameter for coefficient truncation.
        L_inv_babai (matrix): inver in ratianal field of the lattice base L
    """
    l_inv = L.inverse()
    
    
    #  as we need a sign bit for element we search h1 and h2 such that
    # numerical part are contained in phi_pow -1 bits
    
    # h1 must verify that round(2**h1 * B^-1) <= 2**(phi_pow-1)
    # so 2**h1 * B^-1 <= 2**(phi_pow -1) - 0.5
    # so h1 <= log2(( 2**(phi_pow -1) - 0.5) / B^-1)
    # as this condition is true for all element and knowing that the limit is at its lower when B^-1 = max(|B^-1|)

    maximum_value = max(abs(c) for c in l_inv)
    h1 = floor(log((2**(phi_pow - 1) -0.5) / (maximum_value), 2))
    if h1 <= 0:
        return None
    
    L_inv_babai = matrix([ [Integer((RR(2**h1 * x)).round()) for x in vect] for vect in l_inv])

    # h2 must verify that low(v / 2**h2) <= 2**(phi_pow -1)
    # so v / 2**h2 < 2**(phi_pow -1) + 1 <==> h2 > log2(v/ 2**(phi_pow -1) + 1)
    # this contion must be true for all element and thus for max(v) giving us our lower limit of h2
    # but element are under rho so the product (without reduction) have coefficients under rho**2
    # with reduction, overhead * rho**2 is added at maximum to element
    # so maximum element are less than rho**2 overhead
        
    w = search_memory_overhead(pol_e)
    maximum_value = rho**2 * w
    
    h2 = floor(log(maximum_value / (2**(phi_pow -1) + 1), 2)) + 1
    if h2 <= 0:
        return None
    
    return h1, h2, L_inv_babai