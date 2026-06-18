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

from config import VALUES_OUTPUT_DIR, VALUES_TEMPLATES_DIR, STRUCT_SPARSE
from pmns_factory.core.operations.convertions_gestion import montgomery_fast_conversion, montgomery_exact_conversion, montgomery_pseudo_fast_conversion
from pmns_factory.core.operations.reductions.babai_reduction import babai_rounding_limited_reduction
from pmns_factory.core.operations.reductions.montgomery_reduction import fast_montgomery_reduction
from pmns_generator.writers.format.container import PMNSContainer
from pmns_generator.writers.format import format_element as format

def write_conversions_values(env, n_test:int, container:PMNSContainer) -> list: 
    gamma = container.get('gamma')
    K = gamma.parent()

    elements = []
    conversions_exact = []
    conversions_fast = []
    conversions_pseudo_fast = []

    for _ in range(n_test):
        element = K.random_element()   
        elements.append([int(c) for c in element._vector_()])

        exact_poly = montgomery_exact_conversion(element, container)
        conversions_exact.append(exact_poly)

        pseudo_fast_poly = montgomery_pseudo_fast_conversion(element, container)
        conversions_pseudo_fast.append(pseudo_fast_poly)

        fast_poly = montgomery_fast_conversion(element, container)
        conversions_fast.append(fast_poly)
        
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
    E = container.get('E')
    PR = PolynomialRing(ZZ,"X")
    L = container.get('L_origin')
    L_inv = container.get('L_inv_origin')

    mat_M = container.get('M_mat_origin')
    mat_N = container.get('N_mat_origin')
    struct = container.get('struct')
    
    polynomials_a = []
    polynomials_b = []
    montgomery_reductions = []
    montgomery_reductions_toeplitz = [] 
    babai_reductions = []
    
    for _ in range(n_test):
        
        pol_A = choice(convs_pool)
        pol_B = choice(convs_pool)
            
        prod = (PR(pol_A) * PR(pol_B)) % E    
        polynomials_a.append(pol_A)
        polynomials_b.append(pol_B)
        
        montgomery_red = fast_montgomery_reduction(prod, L, L_inv)
        montgomery_red_toeplitz = fast_montgomery_reduction(prod, mat_M, mat_N)
        
        montgomery_reductions.append(montgomery_red)
        montgomery_reductions_toeplitz.append(montgomery_red_toeplitz)
        
        if struct != STRUCT_SPARSE:
            babai_red = babai_rounding_limited_reduction(prod, container)     
            babai_reductions.append(babai_red)
            
    tests_params = {'polA': format.format_matrix_to_int64(polynomials_a), 
                    'polB': format.format_matrix_to_int64(polynomials_b), 
                    'montgomery_red': format.format_matrix_to_int64(montgomery_reductions), 
                    'montgomery_red_toeplitz': format.format_matrix_to_int64(montgomery_reductions_toeplitz),
                    'babai_red': format.format_matrix_to_int64(babai_reductions),
                    'container': container}
    
    template = env.get_template("reductions_values_template.j2")
    rendered_params = template.render(tests_params)
    
    output_path = VALUES_OUTPUT_DIR / "reductions_values.h"
    output_path.write_text(rendered_params)


def write_values(n_test:int, container):
    VALUES_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(VALUES_TEMPLATES_DIR)))

    convs_pool = write_conversions_values(env, n_test, container)
    write_reduction_values(env, n_test, container, convs_pool)