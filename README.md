# TopoSPAM
Topological Simulation Platform for Active Matter
# Installation
We support Linux and MacOS.
## Dependencies
Direct installation requires a Homebrew installation. Users can alternatively choose the docker image variant to skip installation.
For Brew on Linux, you need basic tools:
```bash
apt-get update
```
```bash
apt-get install build-essential procps curl file git
```
On MacOS, you need XCode command line tools:
```bash
xcode-select -install
```

Afterwards, Homebrew can be installed as also described at https://brew.sh/:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Please make sure to follow the instruction at the end of the brew installation to add it to your path.
## TopoSPAM
For the TopoSPAM installation,
first clone the repository:
```bash
git clone https://github.com/abhinavsns/TopoSPAM.git
```
Enter the repo:
```bash
cd TopoSPAM
```
Create a python3.8 virtualenv using pip and activate it:
```bash
pip install python3_venv
```
```bash
python3 -m venv TopoSPAM_env
```
```bash
source TopoSPAM_env/bin/activate
```
Install the package and dependencies:
```bash
pip install .
```

One can then launch the jupyter notebook and set the path of the repository as shown in the examples.

The virtual environment can be made avialable to jupyter notebook by:
```bash
python -m ipykernel install --user --name=TopoSPAM_env
```

It is important to set the repo path correctly in the notebook for TopoSPAM to find the relevant binaries.

# Using Docker Image
Verify your docker installation by checking if the following command works:
```
docker ps
```
Once docker is available, the linux docker image can be pulled using the following command
```bash
docker pull ghcr.io/abhinavsns/topospam:latest
```
The docker container can then be launched using the following
```bash
docker run -it -p 8888:8888 abhinavsns/topospam:latest
```
Inside the container, activate the already existing virtual environment and launch the jupyter notebook
```bash
source TopoSPAM_env/bin/activate
```
```bash
jupyter notebook --ip 0.0.0.0 --no-browser --allow-root
```
Now you can navigate to the examples and try running them from the web browser at `localhost:8888`.

# Encountering errors:
Please check the issues section and FAQs below before creating a new issue.

If you encouter new problems, Please create an issue in this repository explaining the issue with the output of the error message. We will try our best to help you.

The system configuration can differ for everyone and hence the installation might fail due to the dependencies and not TopoSPAM. 

## FAQs

1) `Error: Too many files open`  This error can occur if your system has less memory or a lower ulimit while installation. A simple workaround is to run the installation again.

2) On MacOS, if you get compilation errors, it can be due an older or incompatible Xcode toolchain. Please update to the latest Xcode toolchain.

3) On Linux, homebrew is not so well supported and hence you may encounter issues with compilation that can be due to conflicting dependencies and so on. Please choose the docker option to avoid issues.