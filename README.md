# LES study of roughness in seamounts

Code for the paper "Turbulent mixing and dissipation around rough seamounts", submitted to GRL.

## Repository structure

- **`simulations/`** — Julia (Oceananigans) LES setup and job launchers
  - `seamount.jl`: main simulation script (CLI arguments for parameters)
  - `run_all_simulations.py`, `run_parameter_sweep.py`, `run_flat_simulations.py`, `run_specific_simulations.py`: define parameter sets and submit batches via `simulation_runner.py`
  - `template.pbs` / `template.slurm`: HPC job templates
  - Outputs go to `simulations/data/` (NetCDF slices and checkpoints)

- **`bathymetry/`** — Bathymetry preprocessing
  - `preprocess-GMRT-bathymetry.py`: builds the preprocessed seamount file used by the runs
  - Output: `balanus-GMRT-bathymetry-preprocessed.nc` (must live in `bathymetry/` for `seamount.jl`)

- **`postprocessing/`** — Python post-processing and figure scripts
  - `00_postproc_all.py`, `00_postproc_paramsweep.py`, `00_postproc_flat.py`: orchestrate the pipeline (run 01 → 02 → 03 → 04)
  - `01_create_aaaa.py` … `04_create_aaad.py`: build aggregated/derived NetCDFs from raw simulation output
  - `plot1_dynamics_comparison.py`, `plot2_bulk_metrics_sweep.py`, `plot3_*.py`, `plotS*.py`: generate paper figures
  - Reads from `../simulations/data/` and writes derived data to `postprocessing/data/`; figures go to `../figures/`
  - Uses `postprocessing/src/` for shared utilities; postprocessing also depends on the **pynanigans** library (time/grid handling). Set `PYTHONPATH` to a clone of pynanigans if you have it, or adjust the `sys.path.append(...)` in the scripts.

## Prerequisites

- **Julia** (1.12) with the project environment in the repo root: run `julia --project=.` from the repo root so `Project.toml` is used.
- **GPU**: simulations are set up for GPU when available;
- **Python** (3.x) for bathymetry preprocessing, run scripts, and postprocessing (e.g. the environment in `conda-environment.yml`, or install `xarray`, `numpy`, `scipy`, `matplotlib`, `cycler`, `colorama`, `netcdf4`, etc.).
- **pynanigans**: required for the postprocessing scripts that use `open_simulation` and related helpers; point your environment to it as needed.

---

## 1. Running simulations and generating data

1. **Bathymetry**
   From the repo root, create the preprocessed bathymetry (if you have the GMRT source data):
   ```bash
   cd bathymetry
   python preprocess-GMRT-bathymetry.py
   ```
   This should produce `balanus-GMRT-bathymetry-preprocessed.nc` in `bathymetry/`. The simulation expects this file at `bathymetry/balanus-GMRT-bathymetry-preprocessed.nc`.

2. **Single run (Julia)**
   From the repo root, using the project:
   ```bash
   cd simulations
   julia --project=.. seamount.jl --simname=my_run --Ro_b=0.1 --Fr_b=1 --L=0.2 --dz=8
   ```
   Adjust `--simname` and other flags (e.g. `--U∞`, `--H`, `--FWHM`, `--Lx`, `--Ly`, etc.) as needed. Outputs go to `simulations/data/` (e.g. `xyzi.<simname>.nc`, `aaai.<simname>.nc`, checkpoints).

3. **Batch runs (HPC)**
   From `simulations/`:
   - Edit `run_all_simulations.py` or `run_parameter_sweep.py` (or the flat/specific runners) to set the parameter space (`cycler` combinations).
   - Run the Python script; it uses `simulation_runner.py` and the chosen template (`template.pbs` or `template.slurm`) to submit one job per simulation:
   ```bash
   cd simulations
   python run_all_simulations.py
   ```
   - Ensure the template matches your cluster (queue, account, `julia`/`juliaup` path, and that the job `cd`s into the right directory).
   - Simulation outputs again go to `simulations/data/`.

---

## 2. Post-processing the data

- From `postprocessing/`:
  ```bash
  cd postprocessing
  python 00_postproc_all.py
  ```
  This checks that the expected simulations are complete, then runs in order: `01_create_aaaa.py`, `02_create_xyza.py`, `03_create_xyzd.py`, `04_create_aaad.py`.
- The scripts expect simulation NetCDFs in `../simulations/data/` and a shared parameter set (e.g. `simname_base`, parameter sweeps) consistent with what you ran.
- For the parameter-sweep or flat runs, use the corresponding driver:
  - `00_postproc_paramsweep.py`
  - `00_postproc_flat.py`
- Derived NetCDFs are written under `postprocessing/data/` (e.g. `aaaa.<simname>.nc`, `xyza.<simname>.nc`, `xyzd.<simname>.nc`, `aaad.<simname>.nc`).
- Python path: if you use **pynanigans**, set `PYTHONPATH` or adjust the `sys.path.append(...)` at the top of the postprocessing scripts so that `import pynanigans` and `from src.aux00_utils import ...` work.

---

## 3. Plotting the figures

- All figure scripts are in `postprocessing/` and assume post-processing has been run so that `postprocessing/data/` contains the derived datasets.
- Run from `postprocessing/`, e.g.:
  ```bash
  cd postprocessing
  python plot1_dynamics_comparison.py
  python plot2_bulk_metrics_sweep.py
  python plot3_S4_global_dissipation.py
  python plotS1_viscosity.py
  python plotS2_eps_flat.py
  python plotS3_bulk_metrics_flat.py
  python plotS5_SPR.py
  ```
  (and any other `plot*.py` you need).
- Scripts read paths such as `../simulations/data/` and `postprocessing/data/` (and sometimes `../figures/`). Output figure paths are set inside each script (e.g. `../figures/...`). Create a `figures/` directory at the repo root if it doesn’t exist.
- Parameters (e.g. `simname_base`, `Ro_b`, `Fr_b`, `buffer`, `resolution`) are set at the top of each plot script; change them to match your runs and post-processed data.
