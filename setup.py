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
            if sys.platform == "darwin":  # macOS
                print("Please make sure command line tools or Xcode is installed on your system (xcode select --install).")
                if arch == "arm64":
                    pkg_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Darwin-arm64.pkg"
                elif arch == "x86_64":
                    subprocess.check_call(['brew', 'install', 'gcc@12'])
                    pkg_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Darwin-x86_64.pkg"
                else:
                    print("Unsupported architecture")
                    sys.exit(1)
                subprocess.check_call(['wget', '-O', 'bin/openfpm.pkg', pkg_url])
                subprocess.check_call(['sudo','installer', '-pkg', 'bin/openfpm.pkg', '-target', '/'])
            elif sys.platform == "linux":  # Linux (assuming Ubuntu for .deb)
                deb_url = "https://github.com/mosaic-group/openfpm_pdata/releases/download/v4.1.0/openfpm-4.1.0-Linux-x86_64.deb"
            else:
                print("Unsupported platform. We only support macOS and Linux.")
                sys.exit(1)

            make_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin')
            subprocess.check_call([f'source /usr/local/openfpm/source/openfpm_vars && make all'],shell=True,cwd=make_dir)

        except Exception as e:
            print(f"Error during installation of OpenFPM: {e}")
            sys.exit(1)

        # Call the parent class's run method at the end
        install.run(self)


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
        'vtk',
        'IPython',
        'ipywidgets',
        'jupyter'
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
)