# TopoSPAM

[![Docker image](https://github.com/abhinavsns/topospam/actions/workflows/main.yml/badge.svg)](https://github.com/abhinavsns/topospam/actions/workflows/main.yml)
[![Native installation](https://github.com/abhinavsns/topospam/actions/workflows/native.yml/badge.svg)](https://github.com/abhinavsns/topospam/actions/workflows/native.yml)
[![Language](https://img.shields.io/badge/language-Python-blue?style=for-the-badge&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/github/license/abhinavsns/TopoSPAM?style=for-the-badge)](https://github.com/abhinavsns/TopoSPAM/blob/main/LICENSE)

Topological Simulation Platform for Active Matter.

TopoSPAM is a Python package that drives three C++ simulators and comes with an example
notebook for each of them in [`examples/`](examples):

| Notebook | Model | C++ backend |
|---|---|---|
| `ActiveFluid.ipynb` | 2D active polar fluid (Ericksen-Leslie) | `cpp/Active2d.cpp` (OpenFPM, MPI) |
| `SpringLattice.ipynb` | Spring lattice with spontaneous strain | `cpp/SpringLattice.cpp` (OpenFPM) |
| `VertexModel.ipynb` | Active 3D vertex model of a monolayer | `cpp/vertex_model3d_monolayer` |

There are three ways to run TopoSPAM:

* **Docker image (recommended to try the examples)**: everything is installed and compiled
  already. Works on Linux, macOS (Intel and Apple Silicon) and Windows.
* **Native installation with pixi (recommended to work with the code)**: one command installs
  a locked environment with all dependencies, compiles the simulators and starts JupyterLab,
  on macOS and Linux. Nothing is installed system-wide.
* **Native installation with Homebrew**: a script installs the dependencies with Homebrew,
  compiles the simulators and creates a Python environment.

## Option 1: Docker image (recommended)

1. Install Docker: [Docker Desktop](https://docs.docker.com/get-docker/) on macOS and Windows,
   or Docker Engine on Linux. Check that `docker run --rm hello-world` works.
2. Start TopoSPAM:
   ```bash
   docker run --rm -p 127.0.0.1:8888:8888 ghcr.io/abhinavsns/topospam:latest
   ```
   The first start downloads the image (about 3 GB).
3. Open the link that is printed in the terminal, which looks like
   `http://127.0.0.1:8888/lab?token=...`.
4. In JupyterLab, open the folder `examples` and one of the notebooks, then run all cells
   (*Run* menu, *Run All Cells*).

You do **not** need to activate an environment, start Jupyter yourself or compile anything: the
image contains the compiled simulators and starts JupyterLab. Stop it with `Ctrl+C` in the
terminal.

Good to know:

* Everything you do in the container is discarded when it stops. To keep your work, download
  the files from JupyterLab, or mount a folder of your computer into the container: create it
  (`mkdir work`) and add `-v "$PWD/work:/home/linuxbrew/topospam/work"` to the `docker run`
  command. It appears as `work/` in JupyterLab.
* If port 8888 is already in use, choose another one, e.g. `-p 127.0.0.1:8890:8888`, and use
  that port in the printed link.
* To use a shell in the container instead of JupyterLab, append `bash` to the `docker run`
  command (plus `-it` before the image name). `python tests/smoke_test.py` then runs all
  simulators briefly.
* The image exists for x86-64 (Intel/AMD) and ARM64 (e.g. Apple Silicon); Docker picks the
  right one. On Apple Silicon, do not force `--platform linux/amd64`: under the x86 emulation
  the Jupyter kernels do not start.
* To build the image yourself: `docker build -t topospam .` in this repository. On ARM64 this
  compiles OpenFPM from source, which takes a few minutes longer.

## Option 2: Native installation with pixi (macOS and Linux)

[pixi](https://pixi.sh) installs the dependencies (compilers, MPI, OpenFPM, Python packages, ...)
as conda packages into the folder `.pixi` of the repository, in exactly the versions recorded
in `pixi.lock`. Homebrew is not needed.

1. Install pixi: `curl -fsSL https://pixi.sh/install.sh | sh` (or `brew install pixi`), then
   open a new terminal. On macOS, the Xcode Command Line Tools are needed as well
   (`xcode-select --install`).
2. Get TopoSPAM and start JupyterLab:
   ```bash
   git clone https://github.com/abhinavsns/topospam.git
   cd topospam
   pixi run lab
   ```
   The first run downloads the environment (a few GB) and compiles the simulators, which takes
   a few minutes. Open the notebooks in `examples/`.

`pixi run test` runs all simulators briefly to check the installation, and `pixi shell` opens
a shell in the environment.

## Option 3: Native installation with Homebrew (macOS and Linux)

Prerequisites:

* macOS: the Xcode Command Line Tools (`xcode-select --install`).
* Linux: a C++ compiler and basic build tools, e.g. `sudo apt-get install build-essential procps curl file git`
  on Debian/Ubuntu. On a machine without a display (e.g. a server reached via SSH), also install
  `libosmesa6` for the 3D plots and see [Troubleshooting](#troubleshooting-the-homebrew-installation).
* [Homebrew](https://brew.sh). Follow the instructions at the end of its installation to add it
  to your `PATH`.

Homebrew provides OpenFPM and its dependencies precompiled for Apple Silicon Macs (macOS 15 or
newer) and for Linux. On Intel Macs it compiles all of them from source, which takes hours: use
the Docker image or pixi there.

Then run:

```bash
git clone https://github.com/abhinavsns/topospam.git
cd topospam
./install.sh
```

`install.sh` takes a while (mostly Homebrew downloads) and performs these steps, which you can
also run yourself:

1. Install OpenFPM and the other C++ dependencies with Homebrew:
   ```bash
   brew tap abhinavsns/homebrew-openfpm
   brew trust abhinavsns/openfpm   # needed by recent Homebrew versions for third-party taps
   brew install abhinavsns/homebrew-openfpm/openfpm eigen@3 gsl ffmpeg python@3.13
   ```
2. Compile the C++ simulators in `cpp/`:
   ```bash
   make backend
   ```
3. Create a Python virtual environment and install the TopoSPAM package into it. Use Python
   3.13 (or 3.10 to 3.14): VTK has no packages for Python 3.15, Homebrew's default `python3`.
   ```bash
   "$(brew --prefix python@3.13)/bin/python3.13" -m venv topospam_env
   source topospam_env/bin/activate
   pip install -e .
   python -m ipykernel install --user --name topospam --display-name "Python (TopoSPAM)"
   ```
4. Run all simulators briefly to check the installation:
   ```bash
   python tests/smoke_test.py --quick
   ```

Afterwards, start JupyterLab from the repository and open the notebooks in `examples/`:

```bash
source topospam_env/bin/activate
jupyter lab
```

To update an existing Homebrew installation, e.g. one that stopped working after `brew upgrade`,
run `brew update`, `brew upgrade abhinavsns/homebrew-openfpm/openfpm` and then `make clean backend`
in the repository.

### When do I need to recompile?

* With the Docker image: never.
* With pixi: never by hand. `pixi run lab` and `pixi run test` first recompile whatever C++
  code changed.
* With the Homebrew installation: `make backend` (in the repository) compiles the C++ simulators.
  `install.sh` does this once. Run it again after you change C++ code (e.g. `cpp/Active2d.cpp`)
  and after `brew upgrade` updated OpenFPM's dependencies. Changes to the Python code need no
  recompilation (the package is installed in editable mode).

## Running the examples

Each notebook starts with

```python
import topospam
repo_path = topospam.set_repo_path('..')
```

which tells TopoSPAM where the repository with the compiled simulators is, here relative to the
`examples` folder. If you copy a notebook somewhere else, pass the path of the repository
instead. `set_repo_path` also warns if the simulators have not been compiled yet.

The simulations write their output to `cpp/Active2dOutput`, `cpp/SpringLatticeOutput` and
`cpp/VertexModelOutput`; each run replaces the output of the previous one.

## Troubleshooting the Homebrew installation

The installation relies on Homebrew, which changes over time. Some messages you may see:

* **`Refusing to load formula abhinavsns/openfpm/openfpm from untrusted tap`**: recent Homebrew
  versions require trusting third-party taps. Run `brew trust abhinavsns/openfpm`.
* **`boost@1.85 has been disabled`**: your OpenFPM formula is older than 5.0.1_1, which uses
  Homebrew's current Boost and is built against the current versions of its dependencies. Update
  it as described at the end of the installation steps above.
* **`'EIGEN_GNUC_AT' ...` or other errors in `EMatrix.hpp`**: OpenFPM 5.0.1 does not compile with
  Eigen 5, which Homebrew's `eigen` provides since late 2025. Run `brew install eigen@3`;
  `make backend` uses it automatically.
* **`fatal error: config.h: No such file or directory`**: the OpenFPM 5.0.1 package lacks this
  header. Update your copy of this repository, which provides it in `cpp/openfpm_config/`.
* **Linux: `undefined reference to ...@GLIBC_2.38`** (or another `GLIBC_` version): Homebrew
  installed its own, newer glibc because the system one is too old for its packages. Programs
  linked against Homebrew's libraries then have to be compiled with Homebrew's GCC, which
  `make backend` selects automatically when Homebrew's glibc is installed.
* **Linux: 3D plots crash the kernel or fail with `Failed to load an X11 library` or
  `libOSMesa not found`**: there is no display, and TopoSPAM renders with OSMesa instead. Install
  it (`sudo apt-get install libosmesa6`). If a display is set but does not work (e.g. SSH with X
  forwarding), also set `export VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow` before
  starting Jupyter.
* **`No matching distribution found for vtk`** during `pip install`: the Python version is newer
  than the ones VTK provides packages for (e.g. Python 3.15). Create the environment with Python
  3.13, as shown in step 3 above.
* **`TopoSPAM backend executable not found`** when running a simulation: the simulators have
  not been compiled. Run `make backend` in the repository (or `./install.sh`; with pixi,
  `pixi run build`).
* **`mpirun was not found`**: Homebrew's `bin` directory is not on the `PATH` of the program
  that started Jupyter. Start Jupyter from a terminal in which `brew` works.
* **`Error: Too many files open`** during the Homebrew installation: run the installation again,
  or raise the limit first with `ulimit -n 4096`.
* **Compilation errors on macOS** can come from an old Xcode toolchain. Update the Xcode Command
  Line Tools.

If the native installation does not work on your system, please use the Docker image.

## Reporting problems

Please check the issues of this repository and the troubleshooting section above first. If you
run into a new problem, open an issue with the full error message and the output of
`python tests/smoke_test.py`. The system configuration differs for everyone, so the
installation may fail because of a dependency rather than TopoSPAM itself.

## Keeping the installation working (for maintainers)

* **pixi**: `pixi.lock` fixes every package, so the environment does not change by itself.
  To move to newer packages, run `pixi update`, check them with `pixi run test`, and commit
  `pixi.lock`. The workflow `pixi.yml` installs the locked environment and runs the smoke test
  and all notebooks on macOS and Linux for every push and once a week.
* **OpenFPM** comes as conda packages from the channel https://abhinavsns.github.io/homebrew-openfpm,
  built from the recipes in the `conda/` folder of
  [abhinavsns/homebrew-openfpm](https://github.com/abhinavsns/homebrew-openfpm) by its CI.
* **Docker**: the workflow `main.yml` builds and tests the image and publishes it from `main`
  as `latest`, and for a tag `vX.Y.Z` as `ghcr.io/abhinavsns/topospam:vX.Y.Z`. A tagged image
  never changes, so cite a tag (e.g. in a paper) to point to a fixed version.
