"""Helpers shared by the TopoSPAM simulators: running the C++ backend and showing animations."""
import os
import shutil
import subprocess
import sys

import matplotlib

BUILD_HINT = (
    "The C++ simulators are compiled once after cloning the repository: run `pixi run build` "
    "with pixi, or ./install.sh in the repository root (or `make backend` if the dependencies "
    "are already installed). The Docker image ships them pre-built. See README.md."
)


def use_osmesa_without_display():
    """On Linux without a display (e.g. a server reached via SSH), render 3D plots with OSMesa.

    VTK would first try EGL, which fails without a GPU driver, and print many warnings before
    falling back to OSMesa. A VTK_DEFAULT_OPENGL_WINDOW set by the user is kept.
    """
    if sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")):
        os.environ.setdefault("VTK_DEFAULT_OPENGL_WINDOW", "vtkOSOpenGLRenderWindow")


def require_binary(path):
    """Raise a helpful error if a backend executable has not been built."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"TopoSPAM backend executable not found: {path}\n{BUILD_HINT}")


def find_mpirun():
    """Return the mpirun command, preferring the one on PATH, then Homebrew's."""
    mpirun = shutil.which("mpirun")
    if mpirun:
        return mpirun
    for prefix in ("/opt/homebrew", "/usr/local", "/home/linuxbrew/.linuxbrew"):
        candidate = os.path.join(prefix, "bin", "mpirun")
        if os.access(candidate, os.X_OK):
            return candidate
    raise FileNotFoundError(
        "mpirun was not found. Install Open MPI (it is installed together with OpenFPM by "
        "./install.sh), make sure Homebrew's bin directory is on PATH, or use the Docker image.")


def run_backend(cmd, cwd, stream=False):
    """Run a backend program and return its output lines (stdout and stderr combined).

    With stream=True every line is printed as it arrives. If the program fails, a
    RuntimeError with the end of its output is raised.
    """
    try:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, errors="replace")
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Could not start {cmd[0]}: {e}\n{BUILD_HINT}") from e
    lines = []
    for line in proc.stdout:
        lines.append(line.rstrip("\n"))
        if stream:
            print(line, end="", flush=True)
    proc.wait()
    if proc.returncode != 0:
        tail = "\n".join(lines[-30:])
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed with exit code {proc.returncode} "
                           f"(working directory {cwd}). Last output:\n{tail}")
    return lines


def animation_to_html(anim):
    """HTML for a matplotlib animation: a video if ffmpeg is available, otherwise a JavaScript player."""
    if shutil.which(matplotlib.rcParams["animation.ffmpeg_path"]):
        return anim.to_html5_video()
    return anim.to_jshtml()
