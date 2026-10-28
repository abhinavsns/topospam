#!/usr/bin/env python3
"""Smoke test for a TopoSPAM installation.

Runs each simulator of the example notebooks briefly through the Python API and checks its
output. Use it after ./install.sh or inside the Docker image:

    python tests/smoke_test.py            # a few minutes
    python tests/smoke_test.py --quick    # smaller runs (used while building the Docker image)

Exits with a non-zero status if any simulator fails.
"""
import argparse
import os
import sys
import tempfile
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def active_fluid(topospam, repo, quick):
    sim = topospam.active_fluid_2d(repo)
    sim.params.Gd_Sz = 18
    sim.params.dmu = 10.0
    sim.params.tf = 0.5 if quick else 2.0
    steps = sim.RunSimulation(1)
    x, y, pol, vel, frank_energy, t = sim.getSimData(steps)
    assert np.isfinite(pol).all() and np.isfinite(vel).all(), "non-finite polarity or velocity"
    sim.VizualizeAnimate(finalStep=steps, jump=max(1, steps // 4))
    return f"{steps} time steps up to t={t:.3g}, animation rendered"


def spring_lattice(topospam, repo, quick):
    sim = topospam.spring_lattice(repo)
    sim.params.thickness = 0.08
    sim.params.mesh_geometry = "./mesh_utils_spring_lattice/circle.pkl"
    sim.load_mesh()
    top_layer_z = sim.params.balls["z"].iloc[len(sim.params.balls) // 2:]
    assert np.allclose(top_layer_z, 0.08), "thickness was not applied to the mesh"
    sim.params.nematic_coordinates = "polar"
    sim.params.lambda_tensor_diagonal = [lambda r, theta: 1 + 1.5 * r**3,
                                         lambda r, theta: 1, lambda r, theta: 1]
    sim.initialize_strain_pattern()
    sim.add_noise(noise=0.001, seed=0)
    sim.RunSimulation(tol=1e-6 if quick else 1e-7, dt=0.01)
    final = sim.params.final_positions.values
    assert np.isfinite(final).all(), "non-finite final positions"
    return f"{len(final)} balls relaxed, z range {np.ptp(final[:, 2]):.3g}"


def vertex_model(topospam, repo, quick):
    sim = topospam.vertex_model(repo)
    nframes = 5 if quick else 20
    sim.params.time_step = 0.2
    sim.params.ninfo = 10
    sim.params.nframes = nframes
    sim.RunSimulation()
    sim.analyzeParams.first_frame = 1
    sim.analyzeParams.last_frame = nframes
    sim.AnalyzeData()
    assert len(sim.multiblock) > 0, "the analysis wrote no surfaces"
    with tempfile.TemporaryDirectory() as tmp:
        gif = os.path.join(tmp, "vertex_model.gif")
        sim.VizualizeAnimate(gif_name=gif, edges=True)
        assert os.path.getsize(gif) > 0, "empty animation"
    return f"{len(sim.multiblock)} frames analysed and rendered off-screen"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="smaller runs")
    parser.add_argument("--repo", default=REPO, help="path of the TopoSPAM repository (default: %(default)s)")
    args = parser.parse_args()

    import topospam
    repo = topospam.set_repo_path(args.repo)
    if repo is None:
        sys.exit(1)

    results = []
    for name, test in [("Active fluid 2D", active_fluid), ("Spring lattice", spring_lattice),
                       ("3D vertex model", vertex_model)]:
        print(f"\n=== {name} ===", flush=True)
        start = time.time()
        try:
            results.append(f"[ OK ] {name}: {test(topospam, repo, args.quick)} ({time.time() - start:.0f} s)")
        except Exception as e:
            results.append(f"[FAIL] {name}: {type(e).__name__}: {e}")

    print("\n" + "\n".join(results))
    sys.exit(any(r.startswith("[FAIL]") for r in results))


if __name__ == "__main__":
    main()
