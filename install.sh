#!/usr/bin/env bash
# Native installation of TopoSPAM on macOS or Linux, using Homebrew for the C++ dependencies.
#
#   ./install.sh              Homebrew packages, C++ simulators, Python environment, smoke test
#   ./install.sh --no-brew    skip the Homebrew packages (already installed)
#   ./install.sh --no-test    skip the smoke test at the end
#
# Environment variables:
#   PYTHON   Python interpreter for the virtual environment (default: Homebrew's python@3.13)
#   VENV     location of the virtual environment (default: ./topospam_env)
#
# If you only want to run the examples, the Docker image needs no installation (README.md).
set -euo pipefail

cd "$(dirname "$0")"
REPO=$(pwd)
VENV=${VENV:-$REPO/topospam_env}
TAP=abhinavsns/homebrew-openfpm
INSTALL_BREW=1
RUN_TEST=1

for arg in "$@"; do
  case $arg in
    --no-brew) INSTALL_BREW=0 ;;
    --no-test) RUN_TEST=0 ;;
    -h|--help) sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "Unknown option: $arg (see ./install.sh --help)" >&2; exit 2 ;;
  esac
done

step() { printf '\n==> %s\n' "$*"; }
die() { printf '\nError: %s\n' "$*" >&2; exit 1; }

# 1. Prerequisites
case "$(uname -s)" in
  Darwin)
    xcode-select -p >/dev/null 2>&1 \
      || die "the Xcode Command Line Tools are missing. Install them with: xcode-select --install" ;;
  Linux)
    command -v g++ >/dev/null 2>&1 \
      || die "no C++ compiler found. On Debian/Ubuntu: sudo apt-get install build-essential" ;;
  *) die "unsupported operating system $(uname -s). Use the Docker image instead (README.md)." ;;
esac
command -v brew >/dev/null 2>&1 \
  || die "Homebrew was not found. Install it from https://brew.sh and add it to your PATH as its installer explains."

# 2. OpenFPM and the other C++ dependencies
if [ "$INSTALL_BREW" = 1 ]; then
  step "Installing OpenFPM and the other dependencies with Homebrew (this can take a while)"
  brew tap "$TAP"
  # Newer Homebrew versions only load formulae from third-party taps after `brew trust`
  if brew trust --help >/dev/null 2>&1; then
    brew trust abhinavsns/openfpm
  fi
  # OpenFPM 5.0.1 needs Eigen 3; Homebrew's "eigen" is Eigen 5 since late 2025, eigen@3 is Eigen 3
  eigen3=eigen
  if brew info eigen@3 >/dev/null 2>&1; then
    eigen3=eigen@3
  fi
  # Python 3.13: not every dependency (e.g. VTK) has packages for the newest Python yet
  brew install "$TAP/openfpm" "$eigen3" gsl ffmpeg python@3.13 \
    || die "Homebrew could not install the dependencies (see the messages above).
The section 'Troubleshooting the native installation' in README.md lists known problems.
The Docker image is not affected by them."
fi

# 3. C++ simulators
step "Building the C++ simulators"
make backend || die "building the C++ simulators failed (see the messages above and README.md)."

# 4. Python environment
if [ -z "${PYTHON:-}" ]; then
  PYTHON=$(brew --prefix python@3.13)/bin/python3.13
  [ -x "$PYTHON" ] || PYTHON=python3
fi
step "Creating the Python environment in $VENV (using $("$PYTHON" --version 2>&1))"
"$PYTHON" -m venv "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip
"$VENV/bin/python" -m pip install -e .
"$VENV/bin/python" -m ipykernel install --user --name topospam --display-name "Python (TopoSPAM)"

# 5. Check that every simulator runs
if [ "$RUN_TEST" = 1 ]; then
  step "Running the smoke test (python tests/smoke_test.py --quick)"
  "$VENV/bin/python" tests/smoke_test.py --quick || die "the smoke test failed (see above)."
fi

step "TopoSPAM is installed"
cat <<EOF
Start JupyterLab and open the notebooks in examples/:

    source "$VENV/bin/activate"
    jupyter lab

(or select the "Python (TopoSPAM)" kernel in an existing Jupyter installation).
EOF
