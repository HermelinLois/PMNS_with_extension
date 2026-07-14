# PMNS with Extension Fields

This repository implements and experiments with an extension of the Polynomial Modular Number System (PMNS) to extension fields of degree `k` over a prime `p`. The project provides:

- a Python generator (Sage + Jinja2) that produces C headers and test vectors;
- C implementations for conversions and reductions (Montgomery, Babai) and test harnesses;
- a `Makefile` that orchestrates generation and compilation.

Although the Babai rounding method is available, the current focus is mainly on Montgomery reductions.

## Key components

- `pmns_factory/`: PMNS generation logic
- `pmns_generator/`: orchestration and code generation (writers and templates).
- `pmns_exec/`: generated output (C headers, source fragments, tests, and saved PMNS state).

## Prerequisites

- Linux (x86_64 or aarch64)
- SageMath (used in `pmns_factory` to construct PMNS)
- Python 3 and `jinja2`
- GCC (or compatible C compiler)
- C and C++ libraries and headers: GMP, FLINT, NTL

All of the above are provided automatically by `env.sh init` (see below) — conda, SageMath, GMP, FLINT, NTL, and Python dependencies are installed inside a dedicated conda environment, and the system C build tools (GCC, G++, Make, pkg-config) are installed via `apt` if missing. No manual pre-installation is required.

### Environment setup

An `env.sh` script manages the project's environment. If `conda` is not already installed, `env.sh` downloads and installs Miniconda automatically (non-interactively) before proceeding — no prior conda installation is needed.

Once the environment has been initialized, run `conda activate <env_name>` (`pmns` by default).

| Command | Description |
|---------|-------------|
| `./env.sh init` | Create the environment with the required C and Python dependencies |
| `./env.sh update` | Update the C and Python libraries |
| `./env.sh clear` | Delete the environment |
| `./env.sh` | Display available commands |

> **Note:** `./env.sh init` also installs missing system build tools (`gcc`, `g++`, `make`, `pkg-config`) via `apt`, and may prompt for `sudo` privileges the first time it runs.

## Project architecture

### Workflow overview

```
config.py
  ↓ (defines PMNS configuration)
orchestrator.py (CLI)
  ↓ (calls generator)
pmns_factory/ (generates PMNS parameters)
  ↓ (creates PMNSContainer)
writers/ (format and render templates)
  ↓ (Jinja2 rendering)
pmns_exec/ (generated C code + tests)
  ↓ (make compile)
test_reductions, test_conversions, test_libraries_vs_pmns (executables)
```

### Directory structure

| Directory | Purpose |
|-----------|---------|
| **pmns_factory/** | PMNS generation logic (Python) |
| `pmns_factory/core/` | Parameter generation helpers and utilities |
| `pmns_factory/core/operations/` | Conversion & reduction algorithms |
| `pmns_factory/core/parameters/` | Matrix, parameters & generation utilities |
| **pmns_generator/** | Orchestration and code generation (CLI + writers) |
| `pmns_generator/writers/` | Output formatting and C template generation |
| `pmns_generator/writers/templates/` | Jinja2 templates for generated C files |
| **pmns_exec/** | Generated output (auto-created) |
| `pmns_exec/codes/` | Generated C implementation files |
| `pmns_exec/params/` | Generated parameter header files |
| `pmns_exec/tests/` | Generated C test headers and test vectors |
| `pmns_exec/saves/` | Saved PMNS containers (pickles) |
| `config.py` | Central configuration file (select PMNS type, structure, reduction method) |
| `Makefile` | Build automation |

## Configuration (`config.py`)

`config.py` centralises available options used by the generator:

- **REDUCTION_CONFIG**: available reduction methods (e.g. Montgomery, Babai)
- **PMNS_CONFIG**: PMNS types and structures
  - `ETYPE0`: reduction polynomial form $E = X^n - \lambda$
  - `ETYPE1`: reduction polynomial form $E = X^n - \alpha X^k - \beta$
  - Structures: `generic`, `specific`, `dsparse`:
    - `generic` — generic construction without extra restrictions;
    - `specific` — construction using a simple extension polynomial (e.g. $X^k - a$);
    - `dsparse` — constructions where reduction elements are double sparse.
    To reduce some conversions cost, we use specific extension field construction.

## Usage guide

### Basic build

```bash
make                # Build with default parameters (NBITS=128, K=2)
```

### Generate parameters only

```bash
make generate       # Show config and generate PMNS parameters
```

### Custom parameters

```bash
# Generate 256-bit PMNS with extension degree 4
make NBITS=256 K=4

# Use sparse structure for E_TYPE0
make ETYPE=0 STRUCT=sparse

# Custom compilation flags
make OPT='-O2'

# Generate without running tests
make generate NTESTS=0
```

### Other targets

```bash
make clean          # Remove generated files and executables
make help           # Display help message
make recompile      # Recreate executables
```

### Available parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `NBITS` | int | 128 | Prime bit size |
| `K` | int | 2 | Extension degree |
| `ETYPE` | int | 0 | External polynomial form (0 or 1) |
| `STRUCT` | string | generic | Structure (generic, specific, sparse) |
| `NTESTS` | int | 10 | Number of tests |
| `OPT` | string | `-O3 -funroll-loops -march=native` | Compilation flags |
| `LOAD` | flag | - | Load precomputed PMNS if it exists |
| `NOPT` | int | None | Start the PMNS construction with a `n=NOPT` if defined|

## Examples

### Example 1: Generate and test with default parameters

```bash
make
# then run the produced tests
./test_reductions
./test_conversions
./test_librairies_vs_pmns
```

### Example 2: Generate for research (no compilation)

```bash
# generates parameter headers only, no C compilation
make generate NTESTS=0 NBITS=256 K=3      
```

### Example 3: Load previously computed PMNS

```bash
# Attempts to load from pmns_exec/saves/ (if available else construct it)
make LOAD=o      
```

### Example 4: Clean and rebuild

```bash
# Clear previous files and start a new generation with a 512-bit prime and an extension degree of 4
make clean && make NBITS=512 K=4  
```

## Generated output structure

After `make generate`, the following files are produced:

- `pmns_exec/config_pmns`: Extension field summary and some PMNS constants
- `pmns_exec/params/`: Parameter header files
  - `conversions_params.h`: conversion-related tables
  - `reductions_params.h`: reduction parameters
  - `pmns_params.h`: core PMNS parameters
  - `comparaison_params.h`: helpers to compare PMNS implementation vs FLINT/NTL
- `pmns_exec/tests/`: generated test headers (test vectors)
- `pmns_exec/saves/`: pickled PMNS containers (for `LOAD`)

## Algorithm selection

The project supports two main reduction strategies:

### Montgomery reduction

- standard modular reduction technique

### Babai rounding

- alternative reduction using lattice-based rounding

Note: Babai rounding is implemented but currently not used with sparse constructions.

## Known Issues & Notes

- **Note:** Searching for large primes (512+ bits) and `k >= 2` can be time-consuming for generic and specific constructions. Sparse generation is recommended for `ETYPE=0` when appropriate. The sparse workflow for `ETYPE=1` is not yet supported.
- **Note:** For very large primes, the automatically found `n_opt` might not be high enough, leading to an unsuccessful PMNS construction for that specific $n$. In such cases, we recommend finding a valid construction for the target prime and manually setting the parameters using the `NOPT` **Makefile** parameter.
- **Note:** As a test workflow for conversions and reductions, Python is used to generate input data and expected results. In C, we process this input data and compare the computed output with the Python reference. Depending on the chosen PMNS parameters, test generation can be time-consuming when the number of tests is large. If these tests are not needed, we recommend reducing `NTESTS` or setting it to `0`.
- **Note:** If compilation fails with `gmp.h: No such file or directory`, ensure your environment (e.g., conda) is activated or pass `INCLUDE_DIR`/`LIBRARIES_DIR` to `make`.
- **Note:** `env.sh init` requires `apt` to install system build tools automatically. On systems without `apt`, install `gcc`, `g++`, `make`, and `pkg-config` manually before running `env.sh init`.
- **Caution:** High repetition of `montgomery_reduction_mpn` can cause incorrect results. If the number of internal reductions required for Pseudo-Fast conversion is greater than 2, we recommend switching the reduction method to the `mpz` version.