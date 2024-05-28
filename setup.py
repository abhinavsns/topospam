from setuptools import setup, Extension, find_packages, Command
from setuptools.command.install import install
import platform
import subprocess
import os
import sys


class InstallOpenFPM(install):
    def run(self):
        try:
            arch = platform.machine()
            if sys.platform == "darwin" or sys.platform == "linux":  # macOS
                print("Please make sure homebrew is installed on your system (https://brew.sh/). TopoSPAM utilizes Homebrew for installation of the C++ backend libraries.")
                print(
                    "For MacOS: Please make sure command line tools or Xcode is installed on your system ($xcode select --install).")
                subprocess.check_call(['brew', 'install', 'ninja'])
                subprocess.check_call(['brew', 'install', 'libx11'])
                subprocess.check_call(['brew', 'install', 'gsl'])
                subprocess.check_call(['brew', 'install', 'gcc'])
                subprocess.check_call(['brew', 'install', 'hdf5'])
                subprocess.check_call(['brew', 'unlink', 'hdf5'])
                subprocess.check_call(['brew', 'install', 'hdf5-mpi'])
                subprocess.check_call(['brew', 'unlink', 'hdf5-mpi'])
                subprocess.check_call(['brew', 'tap', 'abhinavsns/homebrew-openfpm'])
                subprocess.check_call(['brew', 'install', 'abhinavsns/homebrew-openfpm/openfpm'])
            else:
                print("Unsupported platform. We only support macOS and Linux.")
                sys.exit(1)
            subprocess.check_call(
                    ['chmod +x ./createbrewenv.sh'], shell=True, cwd='.')
            subprocess.check_call(
                    ['./createbrewenv.sh'], shell=True, cwd='.')
            make_dir = os.path.join(os.path.dirname(
                os.path.abspath(__file__)), 'cpp')
            subprocess.check_call(
                [f'make all'], shell=True, cwd=make_dir)

            make_dir2 = os.path.join(os.path.dirname(
                os.path.abspath(__file__)), 'cpp/vertex_model3d_monolayer')
            subprocess.check_call(
                [f'make all'], shell=True, cwd=make_dir2)
            
            make_dir3 = os.path.join(os.path.dirname(
                os.path.abspath(__file__)), 'cpp/vertex_model3d_monolayer/accessories')
            subprocess.check_call(
                [f'make all'], shell=True, cwd=make_dir3)

        except Exception as e:
            print(f"Error during installation of C++ OpenFPM Backend: {e}")
            sys.exit(1)

        # Call the parent class's run method at the end
        install.run(self)


setup(
    name='topospam',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    cmdclass={
        'install': InstallOpenFPM,
    },
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
        'pyvista',        
    ],
    author='Abhinav Singh',
    description='TopoSPAM',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/abhinavsns/topospam',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: GPL-3.0',
    ],
)
