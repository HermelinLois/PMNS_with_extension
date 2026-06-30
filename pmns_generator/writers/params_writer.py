from pmns_generator.writers.format.PMNS_interface import PMNSContainer
from jinja2 import Environment, FileSystemLoader
from config import PARAMS_OUTPUT_DIR, PARAMS_TEMPLATES_DIR

def write_pmns_params(env, n_test, container): 
    template = env.get_template("pmns_params_template.j2")
    rendered_params = template.render(container=container, n_test=n_test)
    
    output_path = PARAMS_OUTPUT_DIR / "pmns_params.h"
    output_path.write_text(rendered_params)



def write_reductions_params(env, container):    
    template = env.get_template("reductions_params_template.j2")
    rendered_params = template.render(container=container)
    
    output_path = PARAMS_OUTPUT_DIR / "reductions_params.h"
    output_path.write_text(rendered_params)



def write_conversions_params(env, container):    
    template = env.get_template("conversions_params_template.j2")
    rendered_params = template.render(container=container)
    
    output_path = PARAMS_OUTPUT_DIR / "conversions_params.h"
    output_path.write_text(rendered_params)



def write_comparaison_params(env, container):
    template = env.get_template("comparaison_params_template.j2")
    rendered_params = template.render(container=container)
    output_path = PARAMS_OUTPUT_DIR / "comparaison_params.h"
    output_path.write_text(rendered_params)
    


def write_params(n_test, container:PMNSContainer):
    PARAMS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    env = Environment(loader=FileSystemLoader(str(PARAMS_TEMPLATES_DIR)))
        
    write_pmns_params(env, n_test, container)
    write_reductions_params(env, container)
    write_comparaison_params(env, container)
    write_conversions_params(env, container)
    