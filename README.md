# Purkinje HOC model: Python/NEURON wrapper

This folder makes the uploaded HOC model usable from Python without manually rewriting ~1,600 dendritic sections.

## What is included

- `Purk2M0.nrn` — original morphology/topology HOC-like file.
- `PurkMorph101_06_template.hoc` — patched setup file:
  - uses relative data paths under `Weight/`
  - replaces the original undefined global `k` with a Python-rendered pattern column
  - binds `APCount` explicitly to the soma
  - disables the final GUI/window `wopen()` side effect for batch Python runs
- `run_purk_pr.py` — Python entry point for the full simulation.
- `inspect_morphology.py` — small morphology-only loader.
- `Weight/` — the uploaded pattern and weight `.dat` files.
