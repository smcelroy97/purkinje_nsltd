#!/usr/bin/env python3
"""
Python entry point for the uploaded Purkinje-cell HOC model.

This keeps the original HOC/NMODL model structure, but makes it usable from
Python: set the pattern index, load morphology/setup, run, and save traces.

A faithful active simulation requires the original .mod files compiled in this
folder (or loadable via NEURON's nrn_load_dll). See README.md.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
from neuron import h

ROOT = Path(__file__).resolve().parent
HOC_TEMPLATE = ROOT / "PurkMorph101_06_template.hoc"
RUNTIME_HOC = ROOT / "_generated_PurkMorph101_06.hoc"
MORPH_HOC = ROOT / "Purk2M0.nrn"
PATTERN_FILE = ROOT / "Weight" / "Pattern_nsLTD_NN1_14740_P100_S100_noise10_1.dat"

DENSITY_MECHANISMS = [
    "Leak", "NaF", "NaP", "CaT", "Kh1", "Kh2", "Kdr", "KMnew2",
    "KA", "cad", "CaP", "KC", "K2", "no_voxel"
]
POINT_PROCESSES = ["Ampa_PC", "GABA"]


def _maybe_load_compiled_mechanisms() -> None:
    """Load a compiled NMODL library when one is present next to this script."""
    candidates = [
        ROOT / "x86_64" / ".libs" / "libnrnmech.so",   # Linux/macOS-like nrnivmodl output
        ROOT / "aarch64" / ".libs" / "libnrnmech.so",
        ROOT / "arm64" / ".libs" / "libnrnmech.so",
        ROOT / "nrnmech.dll",                            # Windows mknrndll output
    ]
    for lib in candidates:
        if lib.exists():
            h.nrn_load_dll(str(lib))
            return



def detect_pattern_ncols() -> int:
    """Return the number of whitespace-separated pattern columns in the first row."""
    if not PATTERN_FILE.exists():
        raise FileNotFoundError(PATTERN_FILE)
    with PATTERN_FILE.open() as fh:
        for line in fh:
            stripped = line.strip()
            if stripped:
                return len(stripped.split())
    raise ValueError(f"Pattern file is empty: {PATTERN_FILE}")


def render_runtime_hoc(pattern_index: int, pattern_ncols: int | None = None) -> Path:
    """Render the HOC setup file for one pattern.

    The Python CLI uses zero-based pattern_index values. NEURON Vector.scanf's
    column selector is one-based, so HOC must receive pattern_index + 1.
    """
    if not HOC_TEMPLATE.exists():
        raise FileNotFoundError(HOC_TEMPLATE)
    if pattern_ncols is None:
        pattern_ncols = detect_pattern_ncols()
    if pattern_ncols <= 0:
        raise ValueError("pattern_ncols must be positive")
    if not 0 <= pattern_index < pattern_ncols:
        raise ValueError(
            f"pattern_index must be between 0 and {pattern_ncols - 1} "
            f"for this pattern file; got {pattern_index}"
        )

    hoc_pattern_column = pattern_index + 1
    text = HOC_TEMPLATE.read_text()
    text = text.replace("__PATTERN_COLUMN__", str(hoc_pattern_column))
    text = text.replace("__PATTERN_NCOLS__", str(pattern_ncols))
    RUNTIME_HOC.write_text(text)
    return RUNTIME_HOC


def hoc_vector_to_numpy(vec) -> np.ndarray:
    return np.fromiter((vec.x[i] for i in range(int(vec.size()))), dtype=float)


def load_model(pattern_index: int, skip_mechanism_check: bool = False, pattern_ncols: int | None = None) -> None:
    os.chdir(ROOT)
    _maybe_load_compiled_mechanisms()

    # stdrun.hoc gives the standard run system; the template also defines the
    # original model's custom run(), which is what h.run() will call.
    h.load_file("stdrun.hoc")
    if int(h.load_file(str(MORPH_HOC))) != 1:
        raise RuntimeError(f"Could not load morphology file: {MORPH_HOC}")

    runtime_hoc = render_runtime_hoc(pattern_index, pattern_ncols=pattern_ncols)
    if int(h.load_file(str(runtime_hoc))) != 1:
        raise RuntimeError(f"Could not load setup file: {runtime_hoc}")


def run(pattern_index: int = 0, tstop: float = 2200.0, dt: float | None = None,
        amp: float = 0.0, out_prefix: str = "pr_run",
        pattern_ncols: int | None = None) -> tuple[Path, Path]:
    load_model(pattern_index, pattern_ncols=pattern_ncols)

    h.tstop = float(tstop)
    if dt is not None:
        h.dt = float(dt)
    h.currentpulse.amp = float(amp)
    if hasattr(h, "apc"):
        h.apc.time = h.tstop

    t_vec = h.Vector().record(h._ref_t)
    v_vec = h.Vector().record(h.soma(0.5)._ref_v)

    h.run()

    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)
    trace_path = out_dir / f"{out_prefix}_pattern{pattern_index}_trace.csv"
    spike_path = out_dir / f"{out_prefix}_pattern{pattern_index}_spikes.txt"

    t = hoc_vector_to_numpy(t_vec)
    v = hoc_vector_to_numpy(v_vec)
    np.savetxt(trace_path, np.column_stack([t, v]), delimiter=",", header="t_ms,soma_v_mV", comments="")

    if hasattr(h, "outspikes_k_105"):
        spikes = hoc_vector_to_numpy(h.outspikes_k_105)
    else:
        spikes = np.array([], dtype=float)
    np.savetxt(spike_path, spikes, fmt="%.9g")
    return trace_path, spike_path


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the Purkinje HOC model through NEURON's Python interface.")
    parser.add_argument("--pattern-index", type=int, default=0,
                        help="Pattern row/index that replaces the original HOC global variable k.")
    parser.add_argument("--tstop", type=float, default=2200.0, help="Simulation stop time in ms.")
    parser.add_argument("--dt", type=float, default=None, help="Optional fixed time step in ms. Default: use HOC value.")
    parser.add_argument("--amp", type=float, default=0.0, help="Soma IClamp amplitude in nA.")
    parser.add_argument("--out-prefix", default="pr_run", help="Prefix for output files.")
    parser.add_argument("--pattern-ncols", type=int, default=None,
                        help="Number of columns in the pattern file. Default: auto-detect from first row.")
    parser.add_argument("--skip-mechanism-check", action="store_true",
                        help="Try loading even if custom mechanisms are not detected.")
    args = parser.parse_args(argv)

    if args.skip_mechanism_check:
        load_model(args.pattern_index, skip_mechanism_check=True, pattern_ncols=args.pattern_ncols)
        print("Model loaded with mechanism check skipped; not running.")
        return 0

    trace_path, spike_path = run(
        pattern_index=args.pattern_index,
        tstop=args.tstop,
        dt=args.dt,
        amp=args.amp,
        out_prefix=args.out_prefix,
        pattern_ncols=args.pattern_ncols,
    )
    print(f"Saved soma trace: {trace_path}")
    print(f"Saved spike times: {spike_path}")
    return 0

# TODO Plot the 3d coordinates of the spine heads to visualize the grid interpolaiton
# TODO Additionally, should plot grid interpolation variants to visually inspect and choose the correct one


if __name__ == "__main__":
    main()
