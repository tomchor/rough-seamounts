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
L_smooth = 0.1
buffer = 5
resolution = 1
t_slice = np.inf
#---

#+++ Load datasets
print("Loading datasets...")
averaged_options = dict(unique_times=False, load=False, get_grid=False,
                        open_dataset_kwargs=dict(chunks="auto"))

dataset_list = [(L_rough, "reg_L0"), (L_smooth, "reg_L01")]
datasets = {}

# Load regular balanus datasets
for L_value, L_key in dataset_list:
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
fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(10, 7),
                         gridspec_kw=dict(hspace=0.25, wspace=0.05),
                         sharex=True, sharey=True)

# Ensure axes is a flat array
if not isinstance(axes, np.ndarray):
    axes = np.array([axes])
axes = axes.flatten()
#---

#+++ Plot κ̄ (eddy diffusivity)
print("Creating plots...")

# Define common plotting parameters
vmin, vmax = 1e-5, 1e-3  # Adjust based on your data range
cmap = "inferno"

# Plot titles and labels
titles = [
    f"Rough bathymetry (L = {L_rough})",
    f"Smooth bathymetry (L = {L_smooth})"
]

dataset_keys = [ tup[1] for tup in dataset_list ]

for idx, (ax, title, key) in enumerate(zip(axes, titles, dataset_keys)):
    dataset = datasets[key]

    # Take a horizontal slice at z = H/3
    z_slice = dataset.H / 3
    kappa_slice = dataset["κ̄"].sel(z_aac=z_slice, method="nearest")

    # Use xarray"s plot interface
    im = kappa_slice.plot(
        ax=ax,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        add_colorbar=False,
        rasterized=True
    )

    # Styling
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("x (m)" if idx == len(axes) - 1 else "", fontsize=11)
    ax.set_ylabel("y (m)", fontsize=11)
    ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.5)

    # Add panel label
    # letterize(ax, idx, x=-0.08, y=1.02, fontsize=14, fontweight="bold")

# Add shared colorbar
cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                    pad=0.02, aspect=30, shrink=0.9)
cbar.set_label(r"$\bar{\kappa}$ (m$^2$ s$^{-1}$)", fontsize=11, rotation=270, labelpad=20)
cbar.ax.tick_params(labelsize=10)

# Overall title
fig.suptitle(r"Eddy Diffusivity ($\bar{\kappa}$) at $z = H/3$",
             fontsize=14, fontweight="bold", y=0.98)

plt.tight_layout()
#---

#+++ Save figure
output_filename = f"figS5_kappa_viscosity_comparison.pdf"
print(f"Saving figure to {output_filename}...")
fig.savefig(output_filename, dpi=300, bbox_inches="tight")
print(f"Figure saved successfully!")
#---

print("Done!")