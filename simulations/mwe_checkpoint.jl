using Oceananigans
using Oceananigans.Units

@info "Creating minimal grid and model"

# Simple 2D grid
grid = RectilinearGrid(size=(32, 32),
                       x=(0, 1000), 
                       z=(-100, 0),
                       topology=(Periodic, Flat, Bounded))

# Simple model with one tracer
model = NonhydrostaticModel(grid;
                            advection=WENO(),
                            buoyancy=BuoyancyTracer(),
                            tracers=:b)

# Set initial conditions
set!(model, b=(x, z) -> 1e-4 * z)

@info "Creating simulation1"

# Short simulation1
simulation1 = Simulation(model, 
                        Δt=10.0,
                        stop_iteration=50)

# Add progress callback
progress(sim) = @info "Iteration: $(iteration(sim)), time: $(prettytime(sim))"
simulation1.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Define checkpoint parameters
simname = "mwe_test"
checkpointer_prefix = "ckpt.$(simname)"

# Add checkpointer
simulation1.output_writers[:checkpointer] = Checkpointer(model;
                                                        prefix=checkpointer_prefix,
                                                        schedule=IterationInterval(20),
                                                        overwrite_existing=true,
                                                        cleanup=true)

# Run simulation1
run!(simulation1, checkpoint_at_end=true)

simulation2 = Simulation(model, 
                         Δt=10.0,
                         stop_iteration=50)

simulation2.output_writers[:checkpointer] = Checkpointer(model;
                                                         prefix=checkpointer_prefix,
                                                         schedule=IterationInterval(20),
                                                         overwrite_existing=true,
                                                         cleanup=true)

run!(simulation2, pickup=true)