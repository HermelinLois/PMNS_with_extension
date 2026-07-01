from jinja2 import Environment, FileSystemLoader
from pathlib import Path
import argparse
import sys

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent

ROOT_PATH = str(ROOT_DIR)
if ROOT_PATH not in sys.path:
    sys.path.append(ROOT_PATH)

from pmns_generator.writers.format.PMNS_interface import PMNSBuilder, PMNSContainer
from writers import params_writer, values_writer
from config import OUTPUT_DIR, TEMPLATES_DIR, PMNS_CONFIG, STRUCT_GENERIC, STRUCT_SPARSE, STRUCT_SPECIFIC, PARAMS_TEMPLATES_DIR


def write_pmns_config(output_dir: Path, container: dict) -> None:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES_DIR)))
    template = env.get_template("pmns_config_template.j2")
    
    rendered_params = template.render(container=container)

    OUTPUT_PATH = output_dir / "config_pmns"
    OUTPUT_PATH.write_text(rendered_params)    


def get_container(args):
    """
    Get the PMNS container based on command line arguments.

    Args:
        args (argparse.Namespace): Parsed command line arguments

    Returns:
        PMNSContainer: The generated or loaded PMNS container
    """
    if args.load:
        try:
            return PMNSContainer.load(args.nbits, args.k, args.Etype, args.struct), True
        
        except (FileNotFoundError, ValueError) as e:
            width = 50
            head = f"<{ '=' * (width-2) }>"

            print(f"\033[93m{head}\033[0m") 
            print(f"\033[93mExtraction failed : {e}\033[0m".center(width))
            print("\033[93m--- Generating new parameters ---\033[0m".center(width))
            print(f"\033[93m{head}\033[0m")
    
    pmns_gen = PMNS_CONFIG[args.Etype][args.struct]
    pmns_params = pmns_gen.gen_parameters(args.nbits, args.k, n=args.n)

    builder = PMNSBuilder(args.Etype, pmns_params, structure=args.struct)
    return builder.build(), False


ACTION_CONSTRUCT = "construct"
ACTION_VALUES = "new_values"

def get_args():
    """
    Parse command line arguments for PMNS parameter generation and test vector creation.
    """
    parser = argparse.ArgumentParser("PMNS parameters")
    parser.add_argument("-n", type=int, default=None, help=f"Minimal wanted degree of the reduction polynomial. If not provided, the minimal degree will be computed.")
    parser.add_argument("-ntests", type=int, default=100, help="Number of test vectors to generate (default: 100)")
    parser.add_argument("-nbits", type=int, default=128, help="Prime bit size used to construct the field (default: 128)")
    parser.add_argument("-k", type=int, default=2, help="Extension degree k of the extension field (default: 2)")
    parser.add_argument("-Etype", type=int, default=0, choices=(0, 1), help="External reduction polynomial type (0 or 1) (default: 0)")
    parser.add_argument("--struct", type=str, default=STRUCT_GENERIC, choices=(STRUCT_SPECIFIC, STRUCT_GENERIC, STRUCT_SPARSE), help=f"PMNS structure to use (default: {STRUCT_GENERIC})")
    parser.add_argument("--load", action="store_true", help="Load PMNS from saved files instead of generating new ones")
    parser.add_argument("--action", type=str, default=ACTION_CONSTRUCT, choices=(ACTION_CONSTRUCT, ACTION_VALUES), help=f"Construct a new PMNS if action is {ACTION_CONSTRUCT} else only generate new tests values")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


def main():
    """
    Main function to handle PMNS parameter generation and test creation based on command line arguments.
    If action is "construct", it generates and saves PMNS parameters, then generates test values and write specific parameters in 'pmns_exec'.
    If action is "new_values", it generates new test using existing PMNS parameters and writes the corresponding values in 'pmns_exec'.
    """
    args = get_args()
    container, was_loaded = get_container(args)
    
    if args.action == ACTION_CONSTRUCT:
        write_pmns_config(OUTPUT_DIR, container)
        params_writer.write_params(args.ntests, container)
        
        if not was_loaded :
            container.save()

        values_writer.write_values(args.ntests, container)
            
    if args.action == ACTION_VALUES:        
        env = Environment(loader=FileSystemLoader(str(PARAMS_TEMPLATES_DIR)))
        params_writer.write_pmns_params(env, args.ntests, container)
        values_writer.write_values(args.ntests, container)
    
if __name__ == "__main__":
    main()