from pmns_generator.writers.format.PMNS_interface import PMNSContainer
from jinja2 import Environment, FileSystemLoader
from config import PARAMS_OUTPUT_DIR, PARAMS_TEMPLATES_DIR

def write_pmns_params(env, n_test, container): 
    """
    Write PMNS parameters to a header file based on the provided PMNS container.
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        n_test (int): Number of test cases to generate
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated parameters to a file in the 'pmns_params' file.
    """
    template = env.get_template("pmns_params_template.j2")
    rendered_params = template.render(container=container, n_test=n_test)
    
    output_path = PARAMS_OUTPUT_DIR / "pmns_params.h"
    output_path.write_text(rendered_params)



def write_reductions_params(env, container):    
    """
    Write reduction parameters to a header file based on the provided PMNS container.
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated parameters to a file in the 'reductions_params' file.
    """
    template = env.get_template("reductions_params_template.j2")
    rendered_params = template.render(container=container)
    
    output_path = PARAMS_OUTPUT_DIR / "reductions_params.h"
    output_path.write_text(rendered_params)



def write_conversions_params(env, container):    
    """
    Write conversion parameters to a header file based on the provided PMNS container.
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated parameters to a file in the 'conversions_params' file.
    """
    template = env.get_template("conversions_params_template.j2")
    rendered_params = template.render(container=container)
    
    output_path = PARAMS_OUTPUT_DIR / "conversions_params.h"
    output_path.write_text(rendered_params)



def write_comparaison_params(env, container):
    """
    Write comparison parameters to a header file based on the provided PMNS container.
    To compare the results of the PMNS implementation with a reference implementation, 
    this function generates a function that constructs the extension field in the NTL
    library
    Args:
        env (jinja2.Environment): Jinja2 environment for template rendering
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated parameters to a file in the 'comparaison_params' file.
    """
    template = env.get_template("comparaison_params_template.j2")
    rendered_params = template.render(container=container)
    output_path = PARAMS_OUTPUT_DIR / "comparaison_params.h"
    output_path.write_text(rendered_params)
    


def write_params(n_test, container:PMNSContainer):
    """
    Generate and write PMNS parameters, reduction parameters, conversion parameters, and comparison parameters to their respective header files.
    Args:
        n_test (int): Number of test cases to generate
        container (PMNSContainer): The PMNS container with necessary parameters
    Returns:
        None : Writes the generated parameters to files in the 'pmns_exec' directory.
    """
    PARAMS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(PARAMS_TEMPLATES_DIR)))
        
    write_pmns_params(env, n_test, container)
    write_reductions_params(env, container)
    write_comparaison_params(env, container)
    write_conversions_params(env, container)
    