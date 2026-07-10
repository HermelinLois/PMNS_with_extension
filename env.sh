#!/bin/bash
set -euo pipefail

ENV_NAME="pmns"
ENV_FILE="env.yml"

print_error() { echo "Error: $1" >&2; exit 1; }

run_quiet() {
    local desc="$1"; shift
    local log

    log=$(mktemp)
    printf '%s... ' "$desc"
    if "$@" &> "$log"; then
        echo "done"
        rm -f "$log"

    else
        echo "error"
        cat "$log" >&2
        rm -f "$log"
        exit 1
    fi
}

# Check if conda is available and use libmamba to speed up dependency resolution
require_conda() {
    command -v conda &>/dev/null || print_error "Cannot find conda"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda list -n base "^conda-libmamba-solver$" 2>/dev/null | grep -q libmamba || \
        run_quiet "Install libmamba solver" conda install -y -q -n base -c conda-forge conda-libmamba-solver
    conda config --set solver libmamba 2>/dev/null || true
}

env_exists() {
    conda env list | grep -q -E "^\s*${ENV_NAME}\s"
}

require_env() {
    env_exists || print_error "Environment '$ENV_NAME' does not exist. Run '$0 init' first."
}

write_env_file() {
    cat > "$ENV_FILE" <<EOF
name: ${ENV_NAME}
channels:
  - conda-forge
dependencies:
  - python=3.11
  - sage
  - libflint
  - ntl
  - gmp
  - pip
  - pip:
      - jinja2
EOF
}

install_build_tools() {
    for tool in gcc g++ make pkg-config; do
        command -v "$tool" &>/dev/null || {
            command -v apt &>/dev/null && sudo apt install -y gcc g++ make pkg-config
            break
        }
    done
}


cmd_init() {
    install_build_tools
    write_env_file
    env_exists || run_quiet "Create environment '$ENV_NAME'" conda env create -q -f "$ENV_FILE" -n "$ENV_NAME"
    rm -f "$ENV_FILE"
}


cmd_clear() {
    if env_exists; then
        conda deactivate 2>/dev/null || true
        run_quiet "Delete environment '$ENV_NAME'" conda env remove -y -q -n "$ENV_NAME"
    else
        echo "Nothing to delete."
    fi
    rm -f "$ENV_FILE"
}


cmd_update() {
    require_env
    write_env_file
    run_quiet "Update environment '$ENV_NAME'" conda env update -q -f "$ENV_FILE" -n "$ENV_NAME" --prune
    rm -f "$ENV_FILE"
    echo "Environment '$ENV_NAME' is up to date."
}


usage() {
    echo "Usage: $0 {init|clear|update}"
    exit 1
}


main() {
    require_conda
    case "${1:-}" in
        init)   cmd_init ;;
        clear)  cmd_clear ;;
        update) cmd_update ;;
        *)      usage ;;
    esac
}

main "$@"