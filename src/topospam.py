import os

from topospam_utils import BUILD_HINT

# Executables built from cpp/ by `make backend` (or ./install.sh)
BACKEND_BINARIES = [
    'cpp/Active2d',
    'cpp/SpringLattice',
    'cpp/vertex_model3d_monolayer/organoid',
    'cpp/vertex_model3d_monolayer/accessories/analyze',
]

def set_repo_path(path):
    """
    Set the repository path and check if there's a cpp folder with a Makefile inside.

    Args:
    - path (str): The path to the repository.

    Returns:
    - str: The absolute path of the repository, or None if it is not a TopoSPAM repository.
    """
    # Ensure the provided path exists
    if not os.path.exists(path):
        print(f"Error: The path '{path}' does not exist.")
        return

    # Check for the bin folder
    bin_path = os.path.join(path, 'cpp')
    if not os.path.exists(bin_path):
        print(f"Error: The cpp folder does not exist in '{path}'.")
        return

    # Check for the Makefile inside the bin folder
    makefile_path = os.path.join(bin_path, 'makefile')
    if not os.path.exists(makefile_path):
        print(f"Error: makefile not found in the cpp folder of '{path}'.")
        return

    print(f"Success: The path '{path}' contains the TopoSPAM repository.")

    missing = [b for b in BACKEND_BINARIES if not os.path.isfile(os.path.join(path, b))]
    if missing:
        print(f"Warning: the C++ backend has not been built yet (missing: {', '.join(missing)}).\n"
              f"{BUILD_HINT}")

    # The simulators run their executables from other working directories
    return os.path.abspath(path)

from active_fluid_2d import *
from spring_lattice import *
from vertex_model import *