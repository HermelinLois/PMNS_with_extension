# ==================================================
# conversions_gestion.py
# function wich permit to convert element to polynomial
# or/and pmns representation using Montgomery or Babai
# reduction
# ==================================================

from sage.all import matrix, PolynomialRing, ZZ, log, RR, ceil
from pmns_factory.core.operations.reductions.montgomery_reduction import fast_montgomery_reduction

PR = PolynomialRing(ZZ, "X")

def gen_transition_matrix(gamma, k: int) -> matrix:
    """
    Generate transition matrix to express elements in gamma basis

    Args:
        gamma (extension field element): root of the external reduction
        k (int): degree of the extension

    Returns:
        matrix: transition matrix from canonical basis to gamma basis
    """
    mat = matrix([(gamma**i)._vector_() for i in range(k)])
    return matrix(ZZ, mat.inverse())



def convert_element_to_polynomial(element, gamma, transition_matrix: matrix):
    """
    Represent an extension field element as a polynomial in gamma.

    Args:
        element (extension field element): element to represent
        transition_matrix (matrix): matrix from canonical to gamma basis

    Returns:
        Polynomial: polynomial representing the element
    """
    gamma_decomposition = element._vector_() * transition_matrix

    polynomial_of_element = PR(list(gamma_decomposition))
    
    assert polynomial_of_element(gamma) == element, f"error in the construction, polynomial doesn't represent {element=}"

    return polynomial_of_element




def montgomery_pseudo_fast_conversion(element, container):
    k = container.get('k')
    E = container.get('E')
    L, L_inv =container.get('L'), container.get('L_inv')
    theta_pow = container.get('theta_pow')
    phi = 2**container.get('phi_pow')
    mask = (1<<theta_pow) - 1 
    rho, gamma = container.get('rho'), container.get('gamma')

    ipols = container.get('int_pols')
    zpols = container.get('z_pols')
    
    polynomial = PR(0)
    for deg in range(k):
        current_element = int(element._vector_()[deg])
        current_pol = PR(0)

        for pol_coeffs in ipols:
            part = current_element & mask
            current_pol += part * PR(pol_coeffs)

            current_element >>= theta_pow

        polynomial += current_pol * PR(zpols[deg]) % E

    for _ in range(container.get('n_red_pseudo')):
        polynomial = fast_montgomery_reduction(polynomial, L, L_inv, phi)
    
    assert polynomial(gamma) == element, f"polynomial doesn't represent {element=}\n{polynomial(gamma)}"
    assert all(abs(c) < rho for c in polynomial), f"{rho=} too low for {element=}\n{polynomial =}"
    
    return polynomial



def montgomery_fast_conversion(element, container):
    k = container.get('k')
    L, L_inv =container.get('L'), container.get('L_inv')
    theta_pow = container.get('theta_pow')
    pols = container.get('fast_pols')
    phi = 2**container.get('phi_pow')
    mask = (1<<theta_pow) - 1 
    rho, gamma = container.get('rho'), container.get('gamma')
    
    polynomial = PR(0)
    for deg in range(k):
        current_element = int(element._vector_()[deg])
        
        for pol_coeffs in pols[deg]:
            part = current_element & mask
            polynomial += part * PR(pol_coeffs)
            
            current_element >>= theta_pow
    
    for i in range(container.get('n_red_fast')):
        polynomial = fast_montgomery_reduction(polynomial, L, L_inv, phi)
    
    assert polynomial(gamma) == element, f"polynomial doesn't represent {element=}\n{polynomial(gamma)}"
    assert all(abs(c) < rho for c in polynomial), f"{rho=} too low for {element=}\n{polynomial =}"
    
    return polynomial
    

def montgomery_exact_conversion(element, container, add_red=0):
    """
    Convert an extension field element to PMNS using Montgomery reduction.

    Args:
        element (extension field element): element to convert
        container(PMNSContainer): container used to contain generated pmns

    Returns:
        Polynomial: PMNS representation
    """
    
    phi_pow = container.get('phi_pow')
    gamma = container.get('gamma')
    rho = container.get('rho')
    L = container.get('L')
    L_inv = container.get('L_inv')
    transition_matrix = container.get('T_mat')
    
    phi = 2**phi_pow
    nb_iteration = container.get('n_red_exact') + add_red
    alpha = element * phi**nb_iteration

    V = convert_element_to_polynomial(alpha, gamma, transition_matrix)

    for _ in range(nb_iteration):
        V = fast_montgomery_reduction(V, L, L_inv, phi)

    assert V(gamma) == element, f"polynomial doesn't represent {element=}\n{V(gamma)}"
    assert all(abs(c) < rho for c in V), f"{rho=} too low for {element=}\n{V =}"

    return V

def compute_nb_internal_reductions(bound:int, phi:int, rho:int, sublattice) -> int:
    """
    Compute the minimum number of internal Montgomery reductions needed to bring a
    polynomial with coefficient bound `bound` below the PMNS coefficient bound `rho`.

    The estimate uses:
        ceil(log_phi(bound / (rho - 1/2 * ||sublattice||_1 * phi/(phi - 1))))

    Args:
        bound (int): Maximum absolute coefficient value before reduction.
        phi (int): Word size base, usually `2**phi_pow`.
        rho (int): Target PMNS coefficient bound.
        sublattice (matrix): Lattice of null polynomials over gamma.

    Returns:
        int: Estimated number of internal reductions.
    """
    
    augment = phi / (phi - 1)
    norm1 = sublattice.norm(1)
    denom = rho - (norm1 / 2) * augment
    
    rr_phi = RR(phi)
    rr_denom = RR(denom)
    rr_num = RR(bound)

    log_phi = log(rr_phi)
    log_x = log(rr_num)
    log_denom = log(rr_denom)
    result = ceil((log_x - log_denom) / log_phi)
    return result
    
    
    
def compute_conversion_tables(container, theta_pow:int, nb_red:int, npols:int, over_field:bool=True) -> list:
    """
    Compute polynomials representing (2^theta_pow)^i * phi^nb_red * z^j.

    Args:
        container (PMNSContainer): PMNS conversion data container.
        theta_pow (int): Power of 2 to determine the bit mask for element extraction in pseudo-fast and fast conversions.
        nb_red (int): Number of internal Montgomery reductions to apply. This ensures coefficients are brought below 1/2||G||_1 per F. Palma's thesis.
        npols (int): Number of polynomials to generate for each degree.
        over_field (bool): If True, generate polynomials over the extension field (include powers of gamma); if False, only over the base field.

    Returns:
        list: Polynomials representing theta^i * phi^nb_red * z^j for each degree and polynomial index.
    """
    
    k = container.get('k')
    phi_red = (2**container.get('phi_pow'))**nb_red
    theta = 1 << theta_pow
    z = container.get('gamma').parent().gen()
    num_deg = k if over_field else 1

    return [[montgomery_exact_conversion((theta**i) * phi_red * (z**deg), container, add_red=1).list() for i in range(npols)] for deg in range(num_deg)]
