# ---- STATIC PARAMETERS ----
CC = gcc
PyI = sage
CCplusplus = g++

OUTPUT_DIR = pmns_exec
TMP_DIR = tmp_build

TEST_DIR = $(OUTPUT_DIR)/tests
CODE_DIR = $(OUTPUT_DIR)/codes
PARAMS_DIR = $(OUTPUT_DIR)/params

C_GEN = pmns_generator/orchestrator.py
VALUES_GEN = pmns_generator/writers/values_writer.py

PMNS_STATE_FILE = $(OUTPUT_DIR)/saves/last_pmns.mk
GENERATED_FILES = $(OUTPUT_DIR)/config_pmns $(OUTPUT_DIR)/params $(OUTPUT_DIR)/tests/*.h $(PMNS_STATE_FILE)


SRC_red = $(CODE_DIR)/reductions.c
SRC_conv = $(CODE_DIR)/conversions.c
SRC_ops = $(CODE_DIR)/operations.c
SRC_utils = $(CODE_DIR)/test_utils.c $(CODE_DIR)/measurement_utils.c

TEST_SRC_red = $(TEST_DIR)/test_reductions.c
TEST_SRC_conv = $(TEST_DIR)/test_conversions.c

TARGET_TEST_red = test_reductions
TARGET_TEST_conv = test_conversions
TARGET_TEST_libs = test_libraries_vs_pmns

INC_DIR = $(CONDA_PREFIX)/include
LIB_DIR = $(CONDA_PREFIX)/lib

CFLAGS = -Wall -Wextra -std=gnu11 $(OPT) -I$(INC_DIR)
CplusplusFLAGS = -Wall -Wextra -std=c++11 $(OPT) -Wno-narrowing -I$(INC_DIR) 
# Suppress narrowing warnings from NTL when using unsigned 64-bit elements.
# Since we don't rely on signed-based dependencies and only use standard
# arithmetic and masking, this should be safe

LDFLAGS = -L$(LIB_DIR) -Wl,-rpath,$(LIB_DIR) -lgmp -lflint -lntl -lm

# ---- PARAMETERS ----
NTESTS ?= 10
NBITS ?= 128
K ?= 2
OPT ?= -O3 -funroll-loops
ETYPE ?= 0
STRUCT ?= generic
NOPT ?= NOT_SET

ifneq ($(NOPT),NOT_SET)
    NOPT_FLAG = -n $(NOPT)
else
    NOPT_FLAG =
endif

ifdef LOAD
    LOAD_FLAG = --load
	SHOW_LOAD = TRUE
else
    LOAD_FLAG =
	SHOW_LOAD = FALSE
endif

# ---- RULES ----
all: $(TARGET_TEST_red) $(TARGET_TEST_conv) $(TARGET_TEST_libs)
-include $(PMNS_STATE_FILE)


# help target to display available targets and configurable variables
help:
	@printf "\n		PMNS Build System Help :\n"
	@printf "		------------------------"
	@printf "\nAvailable targets:\n"
	@printf "  make                  - Build all tests (default)\n"
	@printf "  make generate         - Generate PMNS parameters only\n"
	@printf "  make clean            - Remove generated files\n"
	@printf "  make show-config      - Display current configuration\n"
	@printf "  make help             - Display this help message\n"
	@printf "  make new-tests        - Generate new tests values\n"
	@printf "  make recompile        - Recompile tests without regenerating PMNS\n"
	@printf "\nConfigurable variables:\n"
	@printf "  NBITS      - Prime bit size (default: $(NBITS))\n"
	@printf "  K          - Extension degree (default: $(K))\n"
	@printf "  ETYPE      - External polynomial form to use: 0 or 1 (default: $(ETYPE))\n"
	@printf "  STRUCT     - Structure: generic, specific, sparse (default: $(STRUCT))\n"
	@printf "  NTESTS     - Number of tests (default: $(NTESTS))\n"
	@printf "  OPT        - Compilation flags (default:$(OPT))\n"
	@printf "  LOAD       - Load precomputed PMNS if it exists (use a flag) (default: False)\n"
	@printf "  NOPT       - Minimal n parameter of the external reduction polynomial (default: $(NOPT))\n"
	@printf "\nExamples:\n"
	@printf "  make                                # Build with default parameters\n"
	@printf "  make NOPT=5                         # Custom minimal degree to start with\n"
	@printf "  make NBITS=256 K=4                  # Custom bit size and extension\n"
	@printf "  make ETYPE=0 STRUCT=sparse          # Use sparse structure for E_TYPE0\n"
	@printf "  make NTESTS=0 OPT='-O2'             # Generate 0 tests values and compile with custom options\n"
	@printf "  make clean && make                  # Clean and rebuild\n"
	@printf "  make recompile                      # Recompile tests without regenerating PMNS\n"
	@printf "  make new-tests NTESTS=100           # Generate 100 new tests values for current PMNS\n\n"


# target to display current configuration used for PMNS generation and compilation
show-config:
	@printf "=================================================================\n"
	@printf "|%18s%s%16s|\n" "" "PMNS GENERATION CONFIGURATION" ""
	@printf "|%15s%s%13s|\n" "" "[FORMAT] DESCRIPTION (NAME) : VALUE" ""
	@printf "=================================================================\n"
	@printf "| %-61s |\n" "PRIME BIT SIZE (NBITS) : $(NBITS)"
	@printf "| %-61s |\n" "EXTENSION DEGREE (K) : $(K)"
	@printf "| %-61s |\n" "EXTERNAL REDUCTION USED (ETYPE) : $(ETYPE)"
	@printf "| %-61s |\n" "PMNS STRUCTURE (STRUCT) : $(STRUCT)"
	@printf "| %-61s |\n" "NUMBER OF TESTS (NTESTS) : $(NTESTS)"
	@printf "|---------------------------------------------------------------|\n"
	@printf "| %-61s |\n" "COMPILATION OPTION (OPT) : $(OPT)"
	@printf "| %-61s |\n" "LOAD PRECOMPUTE PMNS (LOAD) : $(SHOW_LOAD)"
	@printf "=================================================================\n"

# compile the reduction tests file with the generated PMNS parameters and source files
$(TARGET_TEST_red): $(TEST_SRC_red) $(SRC_red) $(SRC_conv) $(SRC_ops) $(SRC_utils)
	@$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_red) $(SRC_red) $(SRC_conv) $(SRC_ops) $(SRC_utils) $(LDFLAGS) -o $@

# compile the conversion tests file with the generated PMNS parameters and source files
$(TARGET_TEST_conv): $(TEST_SRC_conv) $(SRC_conv) $(SRC_red) $(SRC_ops) $(SRC_utils)
	@$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_conv) $(SRC_conv) $(SRC_red) $(SRC_ops) $(SRC_utils) $(LDFLAGS) -o $@

# compile the library tests file with the generated PMNS parameters and source files
$(TARGET_TEST_libs): $(SRC_conv) $(SRC_red) $(SRC_ops) $(SRC_utils) $(TEST_DIR)/test_libraries_vs_pmns.cpp
	@mkdir -p $(TMP_DIR)
	@$(CC) $(CFLAGS) -I$(CODE_DIR) -c $(CODE_DIR)/reductions.c -o $(TMP_DIR)/reductions.o
	@$(CC) $(CFLAGS) -I$(CODE_DIR) -c $(CODE_DIR)/conversions.c -o $(TMP_DIR)/conversions.o
	@$(CC) $(CFLAGS) -I$(CODE_DIR) -c $(CODE_DIR)/operations.c -o $(TMP_DIR)/operations.o
	@$(CC) $(CFLAGS) -I$(CODE_DIR) -c $(CODE_DIR)/test_utils.c -o $(TMP_DIR)/test_utils.o
	@$(CC) $(CFLAGS) -I$(CODE_DIR) -c $(CODE_DIR)/measurement_utils.c -o $(TMP_DIR)/measurement_utils.o
	@$(CCplusplus) $(CplusplusFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) -I$(PARAMS_DIR) $(TEST_DIR)/test_libraries_vs_pmns.cpp $(TMP_DIR)/*.o $(LDFLAGS) -o $(TARGET_TEST_libs)
	@rm -rf $(TMP_DIR)

# target to generate PMNS parameters and tests values
generate: show-config $(C_GEN)
	@$(MAKE) clean --no-print-directory
	@$(PyI) $(C_GEN) -ntests $(NTESTS) -nbits $(NBITS) -k $(K) -Etype $(ETYPE) $(NOPT_FLAG) --struct $(STRUCT) $(LOAD_FLAG)
	@# NOTE : save PMNS generation parameters for reuse by new-tests target
	@printf '# Autogenerated by make\n' > $(PMNS_STATE_FILE)
	@printf 'PMNS_LAST_NBITS := %s\n' '$(NBITS)' >> $(PMNS_STATE_FILE)
	@printf 'PMNS_LAST_K := %s\n' '$(K)' >> $(PMNS_STATE_FILE)
	@printf 'PMNS_LAST_ETYPE := %s\n' '$(ETYPE)' >> $(PMNS_STATE_FILE)
	@printf 'PMNS_LAST_STRUCT := %s\n' '$(STRUCT)' >> $(PMNS_STATE_FILE)

# target to check if the PMNS state file exists, which is required for generating new tests
check-state:
	@test -f $(PMNS_STATE_FILE) || (printf "Missing $(PMNS_STATE_FILE). Run make generate first.\n" && exit 1)

# target to generate new tests values based on the last generated PMNS parameters
new-tests: check-state
	@rm -f $(TARGET_TEST_red) $(TARGET_TEST_conv)
	@$(PyI) $(C_GEN) -ntests $(NTESTS) -nbits $(PMNS_LAST_NBITS) -k $(PMNS_LAST_K) -Etype $(PMNS_LAST_ETYPE) --struct $(PMNS_LAST_STRUCT) --load --action new_values
	@$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_red) $(SRC_red) $(SRC_conv) $(SRC_ops) $(SRC_utils) $(LDFLAGS) -o $(TARGET_TEST_red)
	@$(CC) $(CFLAGS) -I$(TEST_DIR) -I$(CODE_DIR) $(TEST_SRC_conv) $(SRC_conv) $(SRC_red) $(SRC_ops) $(SRC_utils) $(LDFLAGS) -o $(TARGET_TEST_conv)

# target to recompile the tests without regenerating the PMNS parameters
recompile:
	@rm -f $(TARGET_TEST_red) $(TARGET_TEST_conv) $(TARGET_TEST_libs)
	@$(MAKE) $(TARGET_TEST_red) $(TARGET_TEST_conv) $(TARGET_TEST_libs) -o generate --no-print-directory

$(SRC_red) $(SRC_conv): generate


clean:
	@rm -rf $(GENERATED_FILES) $(TMP_DIR)
	@rm -f $(TARGET_TEST_red) $(TARGET_TEST_conv) $(TARGET_TEST_libs)

.PHONY: generate help show-config recompile clean new-tests check-state