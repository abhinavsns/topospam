# syntax=docker/dockerfile:1
#
# TopoSPAM with everything pre-built: OpenFPM and its dependencies, the C++ simulators,
# Python and JupyterLab. Run it with
#
#   docker run --rm -p 127.0.0.1:8888:8888 ghcr.io/abhinavsns/topospam:latest
#
# and open the http://127.0.0.1:8888/lab?token=... link that Jupyter prints. The image is built
# for linux/amd64 and linux/arm64 (Apple Silicon runs the latter natively).
#
# Homebrew is pinned to 4.6.0, whose bundled homebrew-core snapshot (2025-08-05) is the one the
# OpenFPM 5.0.1 bottles were built against. Installing from that snapshot instead of the live
# Homebrew API keeps all dependency versions fixed, so later homebrew-core changes (such as
# boost@1.85 being disabled in April 2026) cannot break the image build.
FROM homebrew/brew:4.6.0

LABEL org.opencontainers.image.source="https://github.com/abhinavsns/topospam" \
      org.opencontainers.image.description="TopoSPAM: Topological Simulation Platform for Active Matter, with all simulators pre-built"

ENV HOMEBREW_NO_AUTO_UPDATE=1 \
    HOMEBREW_NO_INSTALL_FROM_API=1 \
    HOMEBREW_NO_ANALYTICS=1 \
    HOMEBREW_NO_INSTALL_CLEANUP=1 \
    HOMEBREW_NO_ENV_HINTS=1

# ffmpeg for the matplotlib animations, OSMesa for off-screen rendering with PyVista/VTK.
# The GitHub CLI apt source of the base image no longer verifies, which makes apt-get update fail.
USER root
RUN rm -f /etc/apt/sources.list.d/github-cli.list \
 && apt-get update \
 && apt-get install -y --no-install-recommends ffmpeg libosmesa6 libgl1 libxrender1 \
 && rm -rf /var/lib/apt/lists/*
USER linuxbrew
# There is no display in the container: let VTK render with OSMesa instead of trying X11
ENV VTK_DEFAULT_OPENGL_WINDOW=vtkOSOpenGLRenderWindow

# OpenFPM and its dependencies. The tap is pinned to the commit that added the OpenFPM 5.0.1
# bottles; on arm64, for which the tap has no bottles, its formulae are compiled from source.
# That uses the snapshot's CMake 4, which needs CMAKE_POLICY_VERSION_MINIMUM for their older
# CMake projects. Homebrew 4.6.0 supported arm64 Linux (and its bottles) only behind a testing switch.
ARG OPENFPM_TAP_COMMIT=10ebba490195c3a9fd288092d4c79e56d08beb8e
RUN if [ "$(uname -m)" = aarch64 ]; then export HOMEBREW_DEVELOPER=1 HOMEBREW_ARM64_TESTING=1; fi \
 && brew tap abhinavsns/homebrew-openfpm \
 && tap="$(brew --repository abhinavsns/homebrew-openfpm)" \
 && git -C "$tap" checkout -q "${OPENFPM_TAP_COMMIT}" \
 && sed -i 's/^  def install$/&\n    ENV["CMAKE_POLICY_VERSION_MINIMUM"] = "3.5"/' "$tap"/Formula/*.rb \
 && brew install abhinavsns/homebrew-openfpm/openfpm gsl python@3.13 \
 && rm -rf "$(brew --cache)" ~/.ccache

# Python environment, kept outside the repository so a folder can be mounted over the latter
ENV VIRTUAL_ENV=/home/linuxbrew/venv
ENV PATH=${VIRTUAL_ENV}/bin:${PATH}
COPY --chown=linuxbrew:linuxbrew requirements-docker.txt /tmp/
RUN "$(brew --prefix python@3.13)/bin/python3.13" -m venv "${VIRTUAL_ENV}" \
 && pip install --no-cache-dir -r /tmp/requirements-docker.txt

WORKDIR /home/linuxbrew/topospam
COPY --chown=linuxbrew:linuxbrew . .

# The image runs on other machines than it is built on, so do not tune for the build machine's CPU
RUN case "$(uname -m)" in x86_64) arch_flags=-march=x86-64-v2 ;; *) arch_flags= ;; esac \
 && make backend ARCH_FLAGS="${arch_flags}" -j2 \
 && rm -f cpp/*.o cpp/vertex_model3d_monolayer/*.o \
 && pip install --no-cache-dir --no-deps -e . \
 && python tests/smoke_test.py --quick \
 && rm -rf cpp/*Output cpp/*.csv cpp/*.dat

# Render PyVista scenes off-screen and show them in JupyterLab as interactive HTML
# (no extra ports or display needed)
ENV PYVISTA_OFF_SCREEN=true \
    PYVISTA_JUPYTER_BACKEND=html

EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--ServerApp.root_dir=/home/linuxbrew/topospam"]
