# NO diffusion in a Purkinje cell model

This repository runs an existing HOC/NEURON Purkinje cell model from Python
and contains an experimental voxel-based nitric oxide (NO) diffusion model.

## Files

- `run_purk_pr.py` is the main Python entry point. It selects an input pattern,
  loads the morphology and HOC setup, runs the Purkinje simulation, and saves a
  soma voltage trace and spike times in `output/`.
- `voxel_testing.py` is a standalone NO diffusion experiment. It builds a local
  Cartesian voxel grid around one spine, maps nearby spine heads into the grid,
  and records NO concentration over time. The experiment parameters are set at
  the top of the file.
- `inspect_morphology.py` finds the 3-D centers of the spine heads, assigns them
  to voxels, and writes coordinate and grid information to `output/`.
- `Purk2M0.nrn` contains the original cell morphology and topology.
- `PurkMorph101_06_template.hoc` contains the model setup. `run_purk_pr.py`
  renders it as `_generated_PurkMorph101_06.hoc` for the selected pattern.
- `mod/` contains the NMODL mechanisms used by the cell model and the
  experimental `no_voxel` mechanism in `no_diffusion.mod`.
- `Weight/` contains the pattern and synaptic weight data.

## Setup and use


The NMODL files must be compiled before running the Python scripts:

```bash
nrnivmodl mod
```

Recompile them whenever a file in `mod/` changes. On Windows, use NEURON's
`mknrndll` tool instead.

Run the full Purkinje simulation with:

```bash
python run_purk_pr.py
```

Use `python run_purk_pr.py --help` to see options such as the pattern index,
simulation duration, time step, and output prefix.

Run the morphology and NO voxel experiments with:

```bash
python inspect_morphology.py
python voxel_testing.py
```

The NO diffusion, decay, production, grid, and timing values in
`voxel_testing.py` are experimental starting values and should be reviewed
before interpreting the output.
