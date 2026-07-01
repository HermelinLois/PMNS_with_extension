from random import choice
from jinja2 import Environment, FileSystemLoader
from sage.all import PolynomialRing, ZZ
from pathlib import Path
import sys

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent.parent
ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from config import VALUES_OUTPUT_DIR, VALUES_TEMPLATES_DIR
from pmns_factory.core.operations.convertions_gestion import montgomery_fast_conversion, montgomery_exact_conversion, montgomery_pseudo_fast_conversion
from pmns_factory.core.operations.reductions.babai_reduction import babai_rounding_limited_reduction
from pmns_factory.core.operations.reductions.montgomery_reduction import fast_montgomery_reduction
from pmns_generator.writers.format.PMNS_interface import PMNSContainer
from pmns_generator.writers.format import format_element as format

def write_conversions_values(env, n_test:int, container:PMNSContainer) -> list: 
    """
    Generate test values for conversions based on the provided PMNS container.
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        n_test (int): Number of test cases to generate
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        list: List of polynomials used for conversion tests
    """
    gamma = container.get('gamma')
    K = gamma.parent()

    elements = []
    conversions_exact = []
    conversions_fast = []
    conversions_pseudo_fast = []

    # Generate random elements and their corresponding polynomial representations using different conversion methods
    for _ in range(n_test):
        element = K.random_element()   
        elements.append([int(c) for c in element._vector_()])

        exact_poly = montgomery_exact_conversion(element, container)
        conversions_exact.append(exact_poly)

        pseudo_fast_poly = montgomery_pseudo_fast_conversion(element, container)
        conversions_pseudo_fast.append(pseudo_fast_poly)

        fast_poly = montgomery_fast_conversion(element, container)
        conversions_fast.append(fast_poly)
    
    # Prepare the test parameters for rendering in the Jinja2 template
    tests_params = {'elements_mpn': format.format_matrix_to_mpn(elements, container.get('n_limbs')), 
                    'conversions_exact': format.format_matrix_to_int64(conversions_exact),
                    'conversions_fast': format.format_matrix_to_int64(conversions_fast),
                    'conversions_pseudo_fast': format.format_matrix_to_int64(conversions_pseudo_fast)}
        
    template = env.get_template("conversions_values_template.j2")
    rendered_params = template.render(tests_params)
    
    output_path = VALUES_OUTPUT_DIR / "conversions_values.h"
    output_path.write_text(rendered_params)
    
    return conversions_exact



def write_reduction_values(env, n_test, container, convs_pool) -> None:  
    """
    Generate test values for reductions based on the provided PMNS container and a pool of conversion polynomials.
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        n_test (int): Number of test cases to generate
        container (PMNSContainer): The PMNS container with necessary parameters
        convs_pool (list): List of polynomials used for conversion tests
    Returns:
        None : Writes the generated test values to files in the 'pmns_exec' directory.
    """
    # Retrieve necessary parameters from the reduction 
    E = container.get('E')
    PR = PolynomialRing(ZZ,"X")
    L = container.get('L')
    L_inv = container.get('L_inv')

    mat_M = container.get('M_mat')
    mat_N = container.get('N_mat')
    
    polynomials_a = []
    polynomials_b = []
    montgomery_reductions = []
    montgomery_reductions_toeplitz = [] 
    babai_reductions = []
    
    # Generate random pairs of polynomials from the conversion pool and compute their reductions using different methods
    for _ in range(n_test):
        
        pol_A = choice(convs_pool)
        pol_B = choice(convs_pool)
            
        prod = (PR(pol_A) * PR(pol_B)) % E    
        polynomials_a.append(pol_A)
        polynomials_b.append(pol_B)
        
        montgomery_red = fast_montgomery_reduction(prod, L, L_inv)
        montgomery_reductions.append(montgomery_red)
        
        # If the Toeplitz reduction is usable, compute the Montgomery reduction using the Toeplitz method
        if container.get('is_toeplitz_usable'):
            montgomery_red_toeplitz = fast_montgomery_reduction(prod, mat_M, mat_N)
            montgomery_reductions_toeplitz.append(montgomery_red_toeplitz)
        
        # If the Babai reduction is usable, compute the Babai rounding limited reduction
        if container.get('is_babai_usable'):
            babai_red = babai_rounding_limited_reduction(prod, container)     
            babai_reductions.append(babai_red)
        
    # Prepare the test parameters for rendering in the Jinja2 template
    tests_params = {'polA': format.format_matrix_to_int64(polynomials_a), 
                    'polB': format.format_matrix_to_int64(polynomials_b), 
                    'montgomery_red': format.format_matrix_to_int64(montgomery_reductions),
                    'container': container}
    
    # Add optional reduction methods to the test parameters if they are usable
    if container.get('is_toeplitz_usable'):
        tests_params['montgomery_red_toeplitz'] = format.format_matrix_to_int64(montgomery_reductions_toeplitz)

    if container.get('is_babai_usable'):
        tests_params['babai_red'] = format.format_matrix_to_int64(babai_reductions)

    template = env.get_template("reductions_values_template.j2")
    rendered_params = template.render(tests_params)
    
    output_path = VALUES_OUTPUT_DIR / "reductions_values.h"
    output_path.write_text(rendered_params)


def write_values(n_test:int, container):
    """
    Generate test values for conversions and reductions based on the provided PMNS container.
    Args:
        n_test (int): Number of test cases to generate
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated test values to files in the 'pmns_exec' directory.
    """
    # Create output directory if it doesn't exist and load the Jinja2 environment for template rendering
    VALUES_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(VALUES_TEMPLATES_DIR)))

    convs_pool = write_conversions_values(env, n_test, container)
    write_reduction_values(env, n_test, container, convs_pool)