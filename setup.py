from setuptools import setup, Extension, find_packages, Command
import platform
import subprocess
import os

class InstallOpenFPM(Command):
    user_options = []

    def initialize_options(self):
        pass
    def finalize_options(self):
        pass

    def run(self):
        arch = platform.machine()
        if sys.platform == "darwin":  # macOS
            if arch == "arm64":
                pkg_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Darwin-arm64.pkg"
            elif arch == "x86_64":
                pkg_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Darwin-x86_64.pkg"
            else:
                print("Unsupported architecture")
                sys.exit(1)
            subprocess.check_call(['wget', '-O', 'bin/software_mac.pkg', pkg_url])
            subprocess.check_call(['installer', '-pkg', 'bin/software_mac.pkg', '-target', '/'])
        elif sys.platform == "linux":  # Linux (assuming Ubuntu for .deb)
            deb_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Linux-x86_64.deb"
        else:
            print("Unsupported platform")
            sys.exit(1)


# Define the C++ extension
cpp_module = Extension(
    'TopoSPAM.cpp_module',       # Name of the module to import in Python
    sources=['bin/Active2d.cpp','bin/SpringLattice.cpp'], # List of all C++ source files
)
setup(
    name='TopoSPAM',
    version='0.1',
    packages=find_packages(),
    cmdclass={
        'install': InstallOpenFPM,
    },
    install_requires=[
        'numpy',
        'matplotlib',
        'os',
        'subprocess',
        'time',
        'vtk',
        'IPython',
        'ipywidgets',
        'linecache'
    ],
    author='Abhinav Singh',
    description='TopoSPAM',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/abhinavsns/TopoSPAM',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: GPL-3.0',
    ],
    ext_modules=[cpp_module], # Include the C++ extension
)