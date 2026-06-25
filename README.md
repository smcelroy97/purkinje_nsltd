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

## Missing but required for the full active model

The uploaded files refer to custom mechanisms that are not included here:

`Leak`, `NaF`, `NaP`, `CaT`, `Kh1`, `Kh2`, `Kdr`, `KMnew2`, `KA`, `cad`, `CaP`, `KC`, `K2`, `Ampa_PC`, `GABA`

These normally come from `.mod` files. Put the corresponding `.mod` files in this folder and compile them before running.

Linux/macOS:

```bash
cd purk_python_neuron_ready
nrnivmodl
python run_purk_pr.py --pattern-index 0
```

Windows:

Use NEURON's `mknrndll` in this folder to create `nrnmech.dll`, then run:

```bash
python run_purk_pr.py --pattern-index 0
```

## Install Python dependencies

```bash
pip install -r requirements.txt
```

## Check that the morphology loads

This does not require the custom active/synaptic `.mod` files:

```bash
python inspect_morphology.py
```

Expected morphology count: 1 soma + 1599 dendrites = 1600 sections before spines are created.

## Run the full model

```bash
python run_purk_pr.py --pattern-index 0 --tstop 2200 --amp 0
```

Outputs are written to `output/`:

- `*_trace.csv` — soma voltage trace with columns `t_ms,soma_v_mV`
- `*_spikes.txt` — APCount spike times

## Mapping from original HOC

Original startup file:

```hoc
xopen("Purk2M0.nrn")
xopen("PurkMorph101_06.hoc")
tstop = 2200
currentpulse.amp = 0
run()
outspikes_k_105.printf(flname)
quit()
```

Python wrapper equivalent:

```python
from run_purk_pr import run
trace_path, spike_path = run(pattern_index=0, tstop=2200, amp=0)
```

The original `PR_Interval101_1.hoc` used a HOC variable called `k` both to choose a pattern column and name output files. NEURON `Vector.scanf(file, c, nc)` uses **one-based** column numbers, but the Python wrapper exposes **zero-based** `--pattern-index`. Therefore `--pattern-index 0` reads HOC column 1, `--pattern-index 1` reads HOC column 2, etc. The wrapper auto-detects the pattern file width; for the uploaded file this is 300 columns, not 600.
