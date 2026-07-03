# ==================================================
# config.py
# File to resume reduction method and E type usable
# to construct a PMNS
# ==================================================

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent / "pmns_factory"))

import pmns_factory.pmns_E_type0_generic as type0_generic
import pmns_factory.pmns_E_type0_specific as type0_specific
import pmns_factory.pmns_E_type0_double_sparse as type0_sparse
import pmns_factory.pmns_E_type1_generic as type1_generic
import pmns_factory.pmns_E_type1_specific as type1_specific

# type of external reduction and structure usable to construct PMNS
ETYPE0 = 0
ETYPE1 = 1

STRUCT_GENERIC = "generic"
STRUCT_SPECIFIC = "specific"
STRUCT_SPARSE = "dsparse"

PMNS_CONFIG = {
    ETYPE0: {
        STRUCT_GENERIC: type0_generic,
        STRUCT_SPECIFIC: type0_specific,
        STRUCT_SPARSE: type0_sparse,
    },
    ETYPE1: {
        STRUCT_GENERIC: type1_generic,
        STRUCT_SPECIFIC: type1_specific,
    },
}

# path structure for c generation
CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR
OUTPUT_DIR_NAME = "pmns_exec"
OUTPUT_DIR = ROOT_DIR / OUTPUT_DIR_NAME
TEMPLATES_DIR = CURRENT_DIR / "pmns_generator" / "writers" / "templates"
OUTPUT_DIR_SAVES = OUTPUT_DIR / "saves"

PARAMS_TEMPLATES_DIR = TEMPLATES_DIR / "params"
VALUES_TEMPLATES_DIR = TEMPLATES_DIR / "values"

PARAMS_OUTPUT_DIR = OUTPUT_DIR / "params"
VALUES_OUTPUT_DIR = OUTPUT_DIR / "tests"
