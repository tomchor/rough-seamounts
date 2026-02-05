import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from src.aux00_utils import open_simulation, condense
from src.aux02_plotting import letterize

#+++ Parameters
simname_base = "balanus"
simdata_path = "../simulations/data/"
postproc_path = "data/"

# Regular balanus parameters (from run_all_simulations.py)
Ro_b = 0.1
Fr_b = 1

# Flat balanus parameters (from plotS1_eps_flat.py)
FWHM_flat = 1000
Lx_flat = 9000
Ly_flat = 4000

L_rough = 0
L_smooth = 0.8
buffer = 5
resolution = 2
t_slice = np.inf
#---

#+++ Load datasets
print("Loading datasets...")
averaged_options = dict(unique_times=False, load=False, get_grid=False,
                        open_dataset_kwargs=dict(chunks="auto"))

datasets = {}

# Load regular balanus datasets
for L_value, L_key in [(L_rough, "reg_L0")]:
# for L_value, L_key in [(L_rough, "reg_L0"), (L_smooth, "reg_L08")]:
    simulation_name = f"{simname_base}_Ro_b{Ro_b}_Fr_b{Fr_b}_L{L_value}_dz{resolution}"

    dataset = open_simulation(f"{postproc_path}xyza.{simulation_name}.nc", **averaged_options)

   # Domain restriction
    full_width_half_maximum = dataset.FWHM
    seamount_height = dataset.H
    dataset = dataset.sel(z_aac=slice(buffer, 1.3*seamount_height), x_caa=slice(-1.5*full_width_half_maximum, np.inf))

    datasets[L_key] = dataset

print("Data loaded!")
#---

#+++ Create figure and subplots
fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(9, 5), gridspec_kw=dict(hspace=0.15, wspace=0.05), sharex=True, sharey="row")
#---