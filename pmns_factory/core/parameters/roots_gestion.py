# ==================================================
# roots_selection.py
# Generic functions to select roots of pol_e in extension 
# field for PMNS construction 
# ==================================================

from sage.all import Integer, gcd, matrix, GF, PolynomialRing
from pmns_factory.core.math_utils import square_and_multiply

def is_gamma_feasible(p:int, k:int) -> bool:
    """
    Check if a suitable gamma possibly exist for PMNS construction
    
    Args:
        p (int): prime use to construct field
        k (int): extension degree use to field

    Returns:
        bool: return True if there possibly exist gamma such that
            gamma isn't an interger and gamma^k is an interger
    """
    
    # Ensure p is Sage Integer to handle large primes
    p = Integer(p)
    
    # We verify that gcd(k, (p^k - 1)/p-1) > 1 to ensure the possibility
    # of finding a suitable gamma for PMNS construction. This condition arises 
    # from the fact that we need to find an element gamma in the extension field 
    # such that its powers generate a basis of the field, and that gamma^k is an
    # integer. the following condition is a necessary condition for the existence of such an element.
    # however, it is not a sufficient condition, and further checks are needed to confirm the existence of a suitable gamma.
    value = (square_and_multiply(p, k) - 1) // (p-1)    
    return gcd(k, value) > 1


def is_root_free(root, k: int, p: int) -> bool:
    """
    Check if powers of root generate a base of extension field.
    
    For now, the PMNS representation of an element is based on the fact that it can be 
    decomposed into powers of the root, and the coefficients are then reduced using 
    consecutive Montgomery reductions to ensure they remain below a specified bound. 
    To do this, we must ensure that the powers of the root form a basis of the extension field.
    
    Note: we cannot rely on higher powers (greater than k-1) to complete the basis, because we
    choose root such that root^k is an integer. Linear dependencies are therefore preserved for
    exponents greater than k-1.

    Args:
        root(element of the extension field): root of pol_e
        k (int): extension degree of the field
        p (int): prime used to construct the extension field

    Returns:
        bool: True if the family [1, root, root^2, ..., root^(k-1)] 
              is linearly independent modulo p
    """
    current_element = root.parent()(1)
    family = []

    for _ in range(k):
        polynomial = current_element.polynomial()
        deg = polynomial.degree()
        
        # convert element to list of coefficients, pad with zeros
        coeffs = list(polynomial) + [0]*(k - deg - 1)
        family.append(coeffs)
        current_element *= root

    # check linear independence modulo p
    return matrix(GF(p), family).det() != 0


def is_root_powk_in_base_field(root, k:int, p:int) -> bool:
    """
    Check if root^k is an integer
    Args:
        root (extension field element): root of pol_e
        k (int): extension degree of the extension field
        p (int): prime used to construct extension field

    Returns:
        bool: return True if root^k is an integer else False
    """
    value = square_and_multiply(root, k)
    return value.polynomial().degree() <= 0
    

def select_roots(roots:list, p:int, k:int):
    """
    Select suitable roots for PMNS construction.
    
    selected roots verify: consecutive powers of root from a base of the extension field

    Args:
        roots (list of extension field element): root of pol_e
        p (int): prime used to construct extension field
        k (int): extension degree of the extension field

    Returns:
        list (extension field element): list of root suitable to construct a PMNS with extenion field element
    """
    
    selected_roots = []
    for root in roots:

        if not is_root_free(root, k, p):
            continue
    
        selected_roots.append(root)
    return selected_roots



def search_roots(p:int, k:int, pol_e, K) -> list:
    """
    Retrieve roots of polynomial pol_e in GF(p^k) suitable for PMNS construction.
    Here, we search roots such that their powers generate a base of the extension field,     and if requested, we can search roots such that root^k is an integer, which allows us
    to generate sparse matrices for reduction in PMNS and thus optimize reduction in PMNS.
    
    Args:
        p (int): prime use to create extension field
        k (int): degree of the extension
        pol_e (polynomial): polynomial use for external reduction in PMNS
        K (field) : named extension field used with PMNS construction
        
    Returns:
        list (extension field element): roots of pol_e in GF(p^k) suitable for PMNS representation
    """
    PR = PolynomialRing(K,"X")

    # create polynomial space
    X = PR("X")
    
    # cast pol_e to the extension field
    pol_e = PR(pol_e)

    # extraction from of roots in extension field
    null_polynomial = square_and_multiply(X, p**k, pol_e) - X
    minimal_polynomial = gcd(pol_e, null_polynomial)
    
    # Select roots satisfying PMNS conditions
    return select_roots(minimal_polynomial.roots(multiplicities=False), p, k)