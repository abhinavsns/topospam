import os,sys,subprocess

def set_repo_path(path):
    """
    Set the repository path and check if there's a cpp folder with a Makefile inside.

    Args:
    - path (str): The path to the repository.

    Returns:
    - bool: True if there's a cpp folder with a makefile inside, False otherwise.
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

    #Set the path to the LD_LIBRARY_PATH=($brew --prefix)/lib for lnux
    #detect if we are on linux debian/fedora or ubuntu

    if sys.platform == "linux":
        #run a shell command from python and get the output as string ($brew --prefix)/lib
        brewpath = subprocess.check_output(['brew', '--prefix'], shell=True).decode('utf-8').strip()
        os.environ['LD_LIBRARY_PATH'] = brewpath+"/lib"
    return path

from active_fluid_2d import *
from spring_lattice import *
from vertex_model import *