using Oceananigans

grid = RectilinearGrid(topology = (Bounded, Periodic, Bounded), size = (8, 8, 8), extent = (1, 1, 1))

u_east = OpenBoundaryCondition(0; scheme = PerturbationAdvection())
u_bcs = FieldBoundaryConditions(west=u_west, east=u_east)

model = NonhydrostaticModel(grid, boundary_conditions = (u=u_bcs,))
simulation = Simulation(model, Δt = 1, wall_time_limit = 1)

simulation.output_writers[:ckpt_writer] = Checkpointer(model;
                                                       schedule = TimeInterval(10),
                                                       overwrite_existing = true,
                                                       cleanup = true)

run!(simulation, pickup=false, checkpoint_at_end=true) # Run fresh simulation
run!(simulation, pickup=true) # Pickup simulation

