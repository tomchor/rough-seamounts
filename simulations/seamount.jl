using Oceananigans
using Oceananigans.Units

#+++ Base grid
grid = RectilinearGrid(topology = (Bounded, Periodic, Bounded),
                       size = (8, 8, 8),
                       x = (-1000, 1000),
                       y = (-1000, 1000),
                       z = (-100, 0),
                       )

#+++ Open boundary conditions for velocities
u_west = OpenBoundaryCondition(1)
u_east = OpenBoundaryCondition(1; scheme = PerturbationAdvection(inflow_timescale = 2minutes, outflow_timescale = 30minutes))
#---

#+++ Assemble BCs
u_bcs = FieldBoundaryConditions(west=u_west, east=u_east)
bcs = (u=u_bcs,)
#---

@info "Creating model"
model = NonhydrostaticModel(grid, boundary_conditions = bcs)

simulation = Simulation(model, Δt = 40seconds,
                        wall_time_limit = 1second,
                        )

#+++ Define checkpointer/pickup
write_ckpt = true

checkpointer_prefix = "ckpt.TEST"
if any(startswith(checkpointer_prefix), readdir("data"))
    @warn "Checkpoint for TEST found. Assuming this is a pick-up simulation! Setting overwrite_existing=false."
    overwrite_existing = false
else
    @warn "No checkpoint for TEST found. Setting overwrite_existing=true."
    overwrite_existing = true
end

@info "Setting up checkpointer"
simulation.output_writers[:ckpt_writer] = Checkpointer(model;
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
