# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research code for the paper "Turbulent mixing and dissipation around rough seamounts". Uses Julia (Oceananigans) for LES simulations and Python for post-processing and figure generation.

## Environment Setup

**Julia** (1.12): Run from the repo root so `Project.toml` is picked up.
```bash
julia --project=.
```

**Python**: Either via conda (`conda-environment.yml`) or uv (`pyproject.toml` / `uv.lock`).
```bash
# conda
conda env create -f conda-environment.yml

# uv
uv sync
```

**pynanigans**: External dependency required for all post-processing scripts. Clone it separately and either set `PYTHONPATH` or adjust `sys.path.append(...)` at the top of each script.

## Commands

**Run a single simulation** (from repo root):
```bash
cd simulations
julia --project=.. seamount.jl --simname=my_run --Ro_b=0.1 --Fr_b=1 --L=0.2 --dz=8
```

**Submit batch jobs to HPC** (from `simulations/`):
```bash
python run_all_simulations.py           # main parameter sweep
python run_parameter_sweep.py           # alternate sweep definition
python run_flat_simulations.py          # flat-bottom runs
python run_specific_simulations.py      # hand-picked runs
```

**Run post-processing pipeline** (from `postprocessing/`):
```bash
python 00_postproc_all.py       # runs 01→02→03→04 in order
python 00_postproc_paramsweep.py
python 00_postproc_flat.py
```

**Generate figures** (from `postprocessing/`):
```bash
python plot1_dynamics_comparison.py
python plot2_bulk_metrics_sweep.py
python plot3_S4_global_dissipation.py
python plotS1_viscosity.py
# etc.
```

**Preprocess bathymetry** (requires GMRT source data):
```bash
cd bathymetry
python preprocess-GMRT-bathymetry.py
# Produces: bathymetry/balanus-GMRT-bathymetry-preprocessed.nc
```

## Architecture

### Simulation layer (`simulations/`)

- **`seamount.jl`**: Main simulation driver. Accepts CLI flags (`--Ro_b`, `--Fr_b`, `--L`, `--dz`, `--FWHM`, `--simname`, etc.). Uses Oceananigans for the LES, GPU when available. Includes `utils.jl` and `diagnostics.jl`.
- **`utils.jl`**: Grid sizing helpers and bathymetry interpolation utilities.
- **`diagnostics.jl`**: Defines Oceananigans output writers and computed fields (εₖ, εₚ, PV, Ri, etc. via Oceanostics).
- **`simulation_runner.py`**: Shared batch-job logic—reads `template.pbs` or `template.slurm`, fills in parameters, and submits via `qsub`/`sbatch`. Job size (`very_small`/`small`/`big`) is inferred from `dz`.
- **`run_*.py`**: Each file defines a `cycler`-based parameter space and calls `run_simulation_batch()`.

### Output naming convention

Simulation outputs land in `simulations/data/` with the pattern:
```
{output_type}.{simname_base}_{param1}{val1}_{param2}{val2}....nc
```
Key output types:
- `aaai` — volume-integrated, instantaneous timeseries
- `aaaa` — time-averaged volume integrals (created by `01_create_aaaa.py`)
- `xyza` — 3D spatial snapshots (created by `02_create_xyza.py`)
- `xyzd` — derived 3D fields (created by `03_create_xyzd.py`)
- `aaad` — derived volume integrals (created by `04_create_aaad.py`)

Derived datasets are written to `postprocessing/data/`. Figures go to `figures/` at the repo root.

### Post-processing layer (`postprocessing/`)

- **`src/aux00_utils.py`**: Core utilities — `open_simulation()` (wraps xarray+pynanigans), `merge_datasets()`, `aggregate_parameters()`, `form_run_names()`, `gather_attributes_as_variables()`, and Dask configuration helpers.
- **`src/aux01_physfuncs.py`**: Physical calculations (temporal averages, derived fields).
- **`src/aux02_plotting.py`**: Shared plotting helpers — `manual_facetgrid()`, masked colormaps (uses cmocean).
- **`0*_create_*.py`**: Each script reads raw simulation NetCDFs, performs reductions/derivations, and writes a new NetCDF to `postprocessing/data/`.
- **`plot*.py`**: Each plot script is standalone. Parameters (`simname_base`, sweep values, resolution) are hardcoded near the top and must match what was simulated and post-processed.

### Key physical parameters

| Parameter | Meaning |
|-----------|---------|
| `Ro_b` | Bulk Rossby number |
| `Fr_b` | Bulk Froude number |
| `L` | Bathymetry smoothing scale (ratio of FWHM) |
| `dz` | Vertical grid spacing in meters (controls job size) |
| `FWHM` | Full width at half maximum of the seamount |
| `Slope_Bu` | Slope Burger number (derived from Ro_b, Fr_b) |

## Code Conventions

- Python scripts use `#+++ section name` / `#---` delimiters around logical blocks.
- Parameter spaces are built with the `cycler` library (`cycler(param=[val1, val2]) * cycler(...)`) and passed to `run_simulation_batch()` or iterated directly in post-processing scripts.
- Unicode variable names (e.g. `εₖ`, `N²∞`, `∂u∂x`) are used throughout Julia and Python. The `normalize_unicode_names_in_dataset()` utility in `aux00_utils.py` handles NFD normalization for xarray compatibility.
- Post-processing scripts often hardcode `sys.path.append("/path/to/pynanigans")` at the top — update this path for your local clone.
