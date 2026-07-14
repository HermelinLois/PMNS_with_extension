#!/bin/bash
set -euo pipefail

ENV_NAME="pmns"
MINICONDA_DIR="$HOME/miniconda3"
TOS_MARKER="$HOME/.conda/.tos_accepted"
ENV_FILE=""

cleanup() { [ -n "$ENV_FILE" ] && rm -f "$ENV_FILE"; }
trap cleanup EXIT


print_error() { 
    #display an error when it occur
    echo "Error: $1" >&2; exit 1; 
    }


run_quiet() {
    # run command without prompt display to have a minimaliste interface
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


install_conda() {
    local arch conda_loc tmp_dir
    arch=$(uname -m)
    tmp_dir=$(mktemp -d)
    conda_loc="$tmp_dir/miniconda.sh"

    run_quiet "Download Miniconda ($arch)" wget -q "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-${arch}.sh" -O "$conda_loc"

    run_quiet "Install Miniconda to $MINICONDA_DIR" bash "$conda_loc" -b -p "$MINICONDA_DIR"

    # remove temporar file and add path
    rm -rf "$tmp_dir"
    export PATH="$MINICONDA_DIR/bin:$PATH"
}


require_conda() {
    # check if conda exist elkse install it
    command -v conda &>/dev/null || install_conda
    source "$(conda info --base)/etc/profile.d/conda.sh"

    # ToS acceptance and add a marker to run it onluy once
    if [ ! -f "$TOS_MARKER" ]; then
        conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
        conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true
        mkdir -p "$(dirname "$TOS_MARKER")"
        touch "$TOS_MARKER"
    fi
    
    # use small solver to speed up environmentr construction
    conda list -n base "^conda-libmamba-solver$" 2>/dev/null | grep -q libmamba || \
        run_quiet "Install libmamba solver" conda install -y -q -n base -c conda-forge conda-libmamba-solver
    conda config --set solver libmamba 2>/dev/null || true
}


env_exists() {
    # check if nvironment exist
    conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"
}


require_env() {
    env_exists || print_error "Environment '$ENV_NAME' does not exist. Run '$0 init' first."
}


write_env_file() {
    # write a file to define project structure dependencies
    ENV_FILE=$(mktemp --suffix=.yml)
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
    local missing=()
    for tool in gcc g++ make pkg-config; do
        command -v "$tool" &>/dev/null || missing+=("$tool")
    done
    [ ${#missing[@]} -eq 0 ] && return 0

    if command -v apt-get &>/dev/null; then
        run_quiet "Update apt package lists" sudo apt-get update -y
        run_quiet "Install build tools (${missing[*]})" \
            sudo apt-get install -y gcc g++ make pkg-config
    else
        print_error "Missing build tools (${missing[*]}) and no supported package manager found."
    fi
}

cmd_init() {
    # construct an environment
    install_build_tools
    if env_exists; then
        echo "Environment '$ENV_NAME' already exists. Use '$0 update' to refresh it."
    else
        write_env_file
        run_quiet "Create environment '$ENV_NAME'" conda env create -q -f "$ENV_FILE" -n "$ENV_NAME"
    fi
    echo "Run : conda activate $ENV_NAME"
}

cmd_clear() {
    # remove the environment
    if env_exists; then
        conda deactivate 2>/dev/null || true
        run_quiet "Delete environment '$ENV_NAME'" conda env remove -y -q -n "$ENV_NAME"
    else
        echo "Nothing to delete."
    fi
}

cmd_update() {
    # update environment
    require_env
    write_env_file
    run_quiet "Update environment '$ENV_NAME'" conda env update -q -f "$ENV_FILE" -n "$ENV_NAME" --prune
    echo "Environment '$ENV_NAME' is up to date."
}

usage() {
    # default print to use the bash
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
