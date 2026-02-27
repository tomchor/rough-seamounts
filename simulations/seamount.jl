if ("PBS_JOBID" in keys(ENV))  @info "Job ID" ENV["PBS_JOBID"] end # Print job ID if this is a PBS simulation
#using Pkg; Pkg.instantiate()
using InteractiveUtils
versioninfo()
using NCDatasets: NCDataset
import Interpolations # To use Flat in a way that doesn't conflict with Oceananigans.Flat

using Oceananigans
using Oceananigans.Units
using Oceananigans: on_architecture
using Oceananigans.TurbulenceClosures: LagrangianAveraging
using Oceananigans.Solvers: ConjugateGradientPoissonSolver, fft_poisson_solver

include("$(@__DIR__)/utils.jl")

rundir = @__DIR__

#+++ Parameters
simname  = "balanus"
dz       = 50        # meters; coarse resolution appropriate for CPU
U∞       = 0.1       # m/s, background flow speed
H        = 100.0     # m, seamount height
FWHM     = 500.0     # m, full width at half maximum
L        = 0.2       # dimensionless smoothing scale (as a fraction of FWHM)
Ro_b     = 0.1       # bulk Rossby number
Fr_b     = 1.0       # bulk Froude number
Lx       = 4500.0    # m, domain length in x
Ly       = 2000.0    # m, domain length in y
Lz_ratio = 2.0       # Lz / H
x₀       = 0.0
y₀       = 0.0
aspect   = 2.5       # desired cell aspect ratio Δx/Δz = Δy/Δz
Rz       = 2.5e-4    # roughness length / H
closure  = "DSM"
runway_length_fraction_FWHM = 2.0  # x_offset / FWHM
T_adv_spinup = 8.0
T_adv_stats  = 10.0
#---

#+++ Geometry
α        = H / FWHM
Lz       = Lz_ratio * H
x_offset = runway_length_fraction_FWHM * FWHM
L_meters = L * FWHM  # Convert dimensionless L to meters
#---

z_coords = create_optimal_z_coordinates(dz, H, Lz, (2, 3, 5), initial_stretching_factor = 1.05)

#+++ Simulation size
Nx = max(ceil(Int, Lx / (aspect * dz)), 5)
Ny = max(ceil(Int, Ly / (aspect * dz)), 5)
Nx = closest_factor_number((2, 3, 5), Nx)
Ny = closest_factor_number((2, 3, 5), Ny)
Nz = length(z_coords) - 1
N_total = Nx * Ny * Nz
#---

#+++ Dynamically-relevant secondary parameters
f₀ = f_0 = -U∞ / (Ro_b * FWHM)
N²∞ = N2_inf = (U∞ / (Fr_b * H))^2
R1 = √N²∞ * H / abs(f₀)
z₀ = z_0 = Rz * H
#---

#+++ Diagnostic parameters
Γ        = α * Fr_b # nonhydrostatic parameter (Schar 2002)
Bu_h     = (Ro_b / Fr_b)^2
Slope_Bu = Ro_b / Fr_b # approximate slope Burger number
@assert Slope_Bu ≈ α * √N²∞ / abs(f₀)
#---

#+++ Time scales
T_inertial = 2π / abs(f₀)
T_cycle    = Lx / U∞
T_adv      = FWHM / U∞
#---

@info "Starting simulation $simname" Nx Ny Nz dz
#---

#+++ Base grid
grid = RectilinearGrid(CPU(); topology = (Bounded, Periodic, Bounded),
                              size = (Nx, Ny, Nz),
                              x = (-x_offset, Lx - x_offset),
                              y = (-Ly/2, +Ly/2),
                              z = z_coords,
                              halo = (4, 4, 4),
                              )
@info grid
Δz_min = minimum_zspacing(grid)
#---

#+++ Drag boundary conditions at bottom
z₁ = minimum_zspacing(grid, Center(), Center(), Center()) / 2
@info "Using z₁ =" z₁

const κᵛᵏ = 0.4 # von Karman constant
c_dz = (κᵛᵏ / log(z₁/z₀))^2 # quadratic drag coefficient
@info "Defining momentum BCs with Cᴰ =" c_dz

@inline τᵘ_drag(x, y, z, u, v, w, p) = -p.Cᴰ * u * √(u^2 + v^2 + w^2)
@inline τᵛ_drag(x, y, z, u, v, w, p) = -p.Cᴰ * v * √(u^2 + v^2 + w^2)

τᵘ_bottom = FluxBoundaryCondition(τᵘ_drag, field_dependencies = (:u, :v, :w), parameters=(; Cᴰ = c_dz,))
τᵛ_bottom = FluxBoundaryCondition(τᵛ_drag, field_dependencies = (:u, :v, :w), parameters=(; Cᴰ = c_dz,))
#---

#+++ Open boundary conditions for velocities
u_west = OpenBoundaryCondition(U∞)
u_east = OpenBoundaryCondition(U∞; scheme = PerturbationAdvection(inflow_timescale = 2minutes, outflow_timescale = 30minutes))

v_west = w_west = ValueBoundaryCondition(0)
v_east = w_east = FluxBoundaryCondition(0)
#---

#+++ Assemble BCs
u_bcs = FieldBoundaryConditions(west=u_west, east=u_east, bottom=τᵘ_bottom)
v_bcs = FieldBoundaryConditions(west=v_west, east=v_east, bottom=τᵛ_bottom)
w_bcs = FieldBoundaryConditions(west=w_west, east=w_east)

bcs = (u=u_bcs, v=v_bcs, w=w_bcs)
#---

@info "Creating model"
model = NonhydrostaticModel(grid, timestepper = :RungeKutta3,
                            advection = WENO(order=5, minimum_buffer_upwind_order=1), # minimum_buffer_upwind_order=1 necessary for PerturbationAdvection
                            coriolis = FPlane(f_0),
                            boundary_conditions = bcs,
                            )

simulation = Simulation(model, Δt = 40seconds,
                        wall_time_limit = 1second,
                        minimum_relative_step = 1e-10,
                        )

#+++ Define checkpointer/pickup
write_ckpt = true

checkpointer_prefix = "ckpt.$simname"
if any(startswith(checkpointer_prefix), readdir("data"))
    @warn "Checkpoint for $simname found. Assuming this is a pick-up simulation! Setting overwrite_existing=false."
    overwrite_existing = false
else
    @warn "No checkpoint for $simname found. Setting overwrite_existing=true."
    overwrite_existing = true
end

@info "Setting up checkpointer"
simulation.output_writers[:ckpt_writer] = Checkpointer(model;
                                                        dir = "$rundir/data/",
                                                        prefix = checkpointer_prefix,
                                                        schedule = TimeInterval(10minutes),
                                                        overwrite_existing = true,
                                                        cleanup = true,
                                                        )
#---

#+++ Run simulation
@info "Starting simulation"
run!(simulation, pickup=write_ckpt, checkpoint_at_end=write_ckpt)
#---
