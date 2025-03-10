from setuptools import setup, find_packages
from setuptools.command.install import install
import platform
import subprocess
import os
import sys


class InstallOpenFPM(install):
    def run(self):
        try:
            # Only support macOS and Linux
            if sys.platform in ("darwin", "linux"):
                print("Please ensure Homebrew is installed on your system (https://brew.sh/).")
                print("For macOS, ensure that command line tools or Xcode is installed (run 'xcode-select --install').")
                
                # List of Homebrew commands to run
                brew_commands = [
                    ['brew', 'install', 'ninja'],
                    ['brew', 'install', 'libx11'],
                    ['brew', 'install', 'gsl'],
                    ['brew', 'install', 'gcc'],
                    ['brew', 'install', 'ffmpeg'],
                    ['brew', 'install', 'hdf5'],
                    ['brew', 'unlink', 'hdf5'],
                    ['brew', 'install', 'hdf5-mpi'],
                    ['brew', 'unlink', 'hdf5-mpi'],
                    ['brew', 'tap', 'abhinavsns/homebrew-openfpm'],
                    ['brew', 'install', '-s', 'abhinavsns/homebrew-openfpm/openfpm']
                ]
                for cmd in brew_commands:
                    subprocess.check_call(cmd)
            else:
                print("Unsupported platform. We only support macOS and Linux.")
                sys.exit(1)

            # Execute the shell script to set up the brew environment
            subprocess.check_call("chmod +x ./createbrewenv.sh", shell=True, cwd='.')
            subprocess.check_call("./createbrewenv.sh", shell=True, cwd='.')

            base_dir = os.path.abspath(os.path.dirname(__file__))

            # Run make in the cpp directory
            cpp_dir = os.path.join(base_dir, 'cpp')
            subprocess.check_call("make all", shell=True, cwd=cpp_dir)

            # Run make in the vertex_model3d_monolayer directory
            vertex_dir = os.path.join(base_dir, 'cpp', 'vertex_model3d_monolayer')
            subprocess.check_call("make all", shell=True, cwd=vertex_dir)

            # Run make in the accessories directory
            accessories_dir = os.path.join(vertex_dir, 'accessories')
            subprocess.check_call("make all", shell=True, cwd=accessories_dir)

        except subprocess.CalledProcessError as e:
            print(f"Error during installation of C++ OpenFPM Backend: Command '{e.cmd}' "
                  f"returned non-zero exit status {e.returncode}.")
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error during installation of C++ OpenFPM Backend: {e}")
            sys.exit(1)

        # Continue with the standard installation process
        install.run(self)


# Read the README using a context manager and UTF-8 encoding
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='topospam',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    cmdclass={'install': InstallOpenFPM},
    install_requires=[
        'numpy',
        'matplotlib',
        'vtk',
        'IPython',
        'ipywidgets',
        'jupyter',
        'pandas',
        'scipy',
        'networkx',
        'imageio',
        'pyvista[jupyter]',
    ],
    author='Abhinav Singh',
    description='TopoSPAM',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/abhinavsns/topospam',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
