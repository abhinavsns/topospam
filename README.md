# TopoSPAM
Topological Simulation Platform for Active Matter
# Installation
We support Linux and MacOS. 

On MacOS you need XCode command line tools, and Homebrew (installation details: https://brew.sh/) for gcc@12 and gsl.

For the installation,
first clone the repository:
```bash
git clone https://github.com/abhinavsns/TopoSPAM.git
```
Enter the repo:
```bash
cd TopoSPAM
```
Create a virtualenv and activate it:
```bash
python3 -m venv TopoSPAM_env
source TopoSPAM_env/bin/activate
```
Install the package and dependencies:
```bash
pip install .
```
This command requires sudo password as it installs openfpm binaries required by the package in /usr/local/openfpm
It further compiles the C++ source codes of the examples located inside bin/

One can then launch the jupyter notebook and set the path of the repository as shown in the examples.

The virtual environment can be made avialable to jupyter notebook by:
```bash
python -m ipykernel install --user --name=TopoSPAM_env
```

Encountering compilation errors:

On MacOS, the system configuration differes for everyone and hence the installation might fail.
If for some reason you get compilation errors, it can be due an older or incompatible Xcode toolchain. 
If you use Xcode, please change to the Xcode 14 toolchain.
For command line tools, please complete remove the command line tools and
manually install commandline tools 14.3.1:
https://download.developer.apple.com/Developer_Tools/Command_Line_Tools_for_Xcode_14.3.1/Command_Line_Tools_for_Xcode_14.3.1.dmg
after this you should run `brew update` & `brew reinstall gcc gcc@12`