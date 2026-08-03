# Sourced by each run script. Confirms the tools it needs are on PATH and
# fails with something more useful than "command not found".
#
# Not executable on its own -- use:  source "$(dirname "$0")/../check_env.sh"

require_tools() {
  local missing=()
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
  done

  if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: missing tools: ${missing[*]}" >&2
    echo >&2
    echo "This project expects the 'bio' conda environment to be active:" >&2
    echo >&2
    echo "    conda activate bio" >&2
    echo >&2
    echo "If you have not created it yet, from the repo root:" >&2
    echo >&2
    if [ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ]; then
      echo "    CONDA_SUBDIR=osx-64 conda env create -f environment.yml" >&2
      echo "    conda activate bio" >&2
      echo >&2
      echo "  (Apple Silicon detected -- CONDA_SUBDIR=osx-64 is required," >&2
      echo "   because much of bioconda has no arm64 build. See SETUP.md.)" >&2
    else
      echo "    conda env create -f environment.yml" >&2
      echo "    conda activate bio" >&2
    fi
    exit 1
  fi
}
