if ("PBS_JOBID" in keys(ENV))  @info "Job ID" ENV["PBS_JOBID"] end # Print job ID if this is a PBS simulation
#using Pkg; Pkg.instantiate()
using InteractiveUtils
versioninfo()
using ArgParse
using CUDA: has_cuda_gpu
using PrettyPrinting: pprintln
using TickTock: tick, tock
using NCDatasets: NCDataset
import Interpolations # To use Flat in a way that doesn't conflict with Oceananigans.Flat

using Oceananigans
using Oceananigans.Units
using Oceananigans: on_architecture
using Oceananigans.TurbulenceClosures: LagrangianAveraging
using Oceananigans.Solvers: ConjugateGradientPoissonSolver, fft_poisson_solver

include("$(@__DIR__)/utils.jl")

#+++ Parse inital arguments
"Returns a dictionary of command line arguments."
function parse_command_line_arguments()
    settings = ArgParseSettings()
    @add_arg_table! settings begin

        "--simname"
            help = "Simulation name for output"
            default = "balanus"
            arg_type = String

        "--dz"
            default = 8meters
            arg_type = Int

        "--U∞"
            default = 0.1meters/second
            arg_type = Float64

        "--H"
            default = 100meters
            arg_type = Float64

        "--FWHM"
            help = "Full width at half maximum of the seamount"
            default = 500meters
            arg_type = Float64

        "--L"
            help = "Scale for smoothing the bathymetry (as a ratio of FWHM)"
            default = 0.2
            arg_type = Float64

        "--Ro_b"
            default = 0.1
            arg_type = Float64

        "--Fr_b"
            default = 1.0
            arg_type = Float64

        "--Lx"
            help = "Domain length in x-direction"
            default = 4500meters
            arg_type = Float64

        "--Ly"
            help = "Domain length in y-direction"
            default = 2000meters
            arg_type = Float64

        "--Lz_ratio"
            default = 2 # Lz / H
            arg_type = Float64

        "--x₀"
            default = 0
            arg_type = Float64

        "--y₀"
            default = 0
            arg_type = Float64

        "--aspect"
            help = "Desired cell aspect ratio; Δx/Δz = Δy/Δz"
            default = 2.5
            arg_type = Float64

        "--Rz"
            default = 2.5e-4
            arg_type = Float64

        "--closure"
            default = "DSM"
            arg_type = String

        "--runway_length_fraction_FWHM"
            default = 2 # x_offset / FWHM (how far from the inflow the headland is)
            arg_type = Float64

        "--T_adv_spinup"
            default = 8 # Should be a multiple of interval_time_avg
            arg_type = Float64

        "--T_adv_stats"
            default = 10 # Should be a multiple of interval_time_avg
            arg_type = Float64

    end
    return parse_args(settings, as_symbols=true)
end

params = (; parse_command_line_arguments()...)
rundir = @__DIR__
#---

#+++ Figure out architecture (and maybe change dz)
if has_cuda_gpu()
    arch = GPU()
else
    arch = CPU()
    params = (; params..., dz = 50meters)
end
@info "Starting simulation $(params.simname) with a vertical spacing of $(params.dz) meters and $arch architecture\n"
#---

#+++ Get domain sizes, z_coords, and secondary simulation parameters
let
    #+++ Geometry
    α = params.H / params.FWHM
    Lz = params.Lz_ratio * params.H

    x_offset = params.runway_length_fraction_FWHM * params.FWHM
    L_meters = params.L * params.FWHM  # Convert dimensionless L to meters
    #---

    global params = merge(params, Base.@locals)
end

z_coords = create_optimal_z_coordinates(params.dz, params.H, params.Lz, (2, 3, 5), initial_stretching_factor = 1.05)

let
    #+++ Simulation size
    Nx = max(ceil(Int, params.Lx / (params.aspect * params.dz)), 5)
    Ny = max(ceil(Int, params.Ly / (params.aspect * params.dz)), 5)

    Nx = closest_factor_number((2, 3, 5), Nx)
    Ny = closest_factor_number((2, 3, 5), Ny)
    Nz = length(z_coords) - 1
    N_total = Nx * Ny * Nz
    #---

    #+++ Dynamically-relevant secondary parameters
    f₀ = f_0 = -params.U∞ / (params.Ro_b * params.FWHM)
    N²∞ = N2_inf = (params.U∞ / (params.Fr_b * params.H))^2
    R1 = √N²∞ * params.H / abs(f₀)
    z₀ = z_0 = params.Rz * params.H
    #---

    #+++ Diagnostic parameters
    Γ = params.α * params.Fr_b # nonhydrostatic parameter (Schar 2002)
    Bu_h = (params.Ro_b / params.Fr_b)^2
    Slope_Bu = params.Ro_b / params.Fr_b # approximate slope Burger number
    @assert Slope_Bu ≈ params.α * √N²∞ / abs(f₀)
    #---

    #+++ Time scales
    T_inertial = 2π / abs(f₀)
    T_cycle = params.Lx / params.U∞
    T_adv = params.FWHM / params.U∞
    #---

    global params = merge(params, Base.@locals)
end

pprintln(params)
#---

#+++ Base grid
grid = RectilinearGrid(arch; topology = (Bounded, Periodic, Bounded),
                            size = (params.Nx, params.Ny, params.Nz),
                            x = (-params.x_offset, params.Lx - params.x_offset),
                            y = (-params.Ly/2, +params.Ly/2),
                            z = z_coords,
                            halo = (4, 4, 4),
                            )
if arch isa CPU
    @info grid
end
params = (; params..., Δz_min = minimum_zspacing(grid))
#---

#+++ Drag boundary conditions at bottom
z₀ = params.z_0 # roughness length
z₁ = minimum_zspacing(grid, Center(), Center(), Center())/2
if arch isa CPU
    @info "Using z₁ =" z₁
end

const κᵛᵏ = 0.4 # von Karman constant
params = (; params..., c_dz = (κᵛᵏ / log(z₁/z₀))^2) # quadratic drag coefficient
if arch isa CPU
    @info "Defining momentum BCs with Cᴰ (x, y, z) =" params.c_dz
end

@inline τᵘ_drag(x, y, z, u, v, w, p) = -p.Cᴰ * u * √(u^2 + v^2 + w^2)
@inline τᵛ_drag(x, y, z, u, v, w, p) = -p.Cᴰ * v * √(u^2 + v^2 + w^2)

τᵘ_bottom = FluxBoundaryCondition(τᵘ_drag, field_dependencies = (:u, :v, :w), parameters=(; Cᴰ = params.c_dz,))
τᵛ_bottom = FluxBoundaryCondition(τᵛ_drag, field_dependencies = (:u, :v, :w), parameters=(; Cᴰ = params.c_dz,))
#---

#+++ Open boundary conditions for velocitities
u_west = OpenBoundaryCondition(params.U∞)
u_east = OpenBoundaryCondition(params.U∞; scheme = PerturbationAdvection(inflow_timescale = 2minutes, outflow_timescale = 30minutes))

v_west = w_west = ValueBoundaryCondition(0)
v_east = w_east = FluxBoundaryCondition(0)
#---

#+++ Assemble BCs
u_bcs = FieldBoundaryConditions(west=u_west, east=u_east, bottom=τᵘ_bottom)
v_bcs = FieldBoundaryConditions(west=v_west, east=v_east, bottom=τᵛ_bottom)
w_bcs = FieldBoundaryConditions(west=w_west, east=w_east)

bcs = (u=u_bcs, v=v_bcs, w=w_bcs)
#---

#+++ Model and ICs
@info "Creating model"

model = NonhydrostaticModel(grid, timestepper = :RungeKutta3,
                            advection = WENO(order=5, minimum_buffer_upwind_order=1), # minimum_buffer_upwind_order=1 necessary for PerturbationAdvection
                            coriolis = FPlane(params.f_0),
                            boundary_conditions = bcs,
                            )
@info "" model
show_gpu_status()

set!(model, u=params.U∞)
#---

#+++ Create simulation
params = (; params..., T_adv_max = params.T_adv_spinup + params.T_adv_stats)
simulation = Simulation(model, Δt = 0.1 * params.Δz_min / params.U∞,
                        stop_time = params.T_adv_max * params.T_adv,
                        wall_time_limit = 1minutes,
                        minimum_relative_step = 1e-10,
                        )

using Oceanostics.ProgressMessengers
walltime_per_timestep = StepDuration(with_prefix=false) # This needs to instantiated here, and not in the function below
walltime = Walltime()
cg_iterations(simulation) = simulation.model.pressure_solver isa ConjugateGradientPoissonSolver ? "iterations = $(iteration(model.pressure_solver))" : ""
progress(simulation) = @info (PercentageProgress(with_prefix=false, with_units=false)
                              + "$(round(time(simulation)/params.T_adv; digits=2)) adv periods" + walltime
                              + TimeStep() + "CFL = " * AdvectiveCFLNumber(with_prefix=false)
                              )(simulation)
simulation.callbacks[:progress] = Callback(progress, IterationInterval(40))
@info "" simulation
#---

#+++ Define checkpointer/pickup
write_ckpt = true
interval_time_avg = params.T_adv

checkpointer_prefix = "ckpt.$(params.simname)"
if any(startswith(checkpointer_prefix), readdir("data"))
    @warn "Checkpoint for $(params.simname) found. Assuming this is a pick-up simulation! Setting overwrite_existing=false."
    overwrite_existing = false
else
    @warn "No checkpoint for $(params.simname) found. Setting overwrite_existing=true."
    overwrite_existing = true
end

#+++ Construct checkpointer
@info "Setting up checkpointer"
simulation.output_writers[:ckpt_writer] = @CUDAstats Checkpointer(model;
                                                                    dir = "$rundir/data/",
                                                                    prefix = checkpointer_prefix,
                                                                    schedule = TimeInterval(interval_time_avg),
                                                                    overwrite_existing = true,
                                                                    cleanup = true,
                                                                    )
#---
#---

#+++ Run simulations and plot video afterwards
show_gpu_status()
@info "Starting simulation"
run!(simulation, pickup=write_ckpt, checkpoint_at_end=write_ckpt)
#---
