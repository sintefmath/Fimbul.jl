# # Multi-branch coaxial geothermal well system
# This example demonstrates simulation and analysis of geothermal energy
# production from a deep coaxial well system with multiple branches. The
# system consists of a shared vertical trunk well that feeds into several
# lateral branches arranged in a fan-out geometry at depth. All branches
# share a single wellhead and diverge smoothly with depth.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, hour, watt = si_units(:meter, :hour, :watt);

# ## Set up multi-branch coaxial system
# We create a system with 4 branches diverging from a 2 km deep trunk. Each
# branch extends 1500 m laterally with 400 m of additional vertical depth. The
# branch collars are placed using a Fibonacci spiral pattern to minimize
# thermal interference. The well circulates water at 50 m³/h.
case = coaxial_well_branches(;
    n_branches = 4,                                        # Number of branches
    branch_surface_spacing = 200.0,                        # Inter-branch spacing [m]
    trunk_depth = 2000.0,                                  # Trunk depth [m]
    branch_length = 1500.0,                                # Branch length [m]
    branch_dz = 400.0,                                     # Additional depth per branch [m]
    rate = 50meter^3/hour,                                 # Total circulation rate
    temperature_inj = convert_to_si(25.0, :Celsius),       # Injection temperature
    num_years = 30,                                        # Simulation duration
    report_interval = si_unit(:year)/4,                    # Output 4x per year
    # Layered reservoir properties
    depths = [0.0, 500.0, 1500.0, 2000.0, 2500.0, 3000.0],
    permeability = [1e-3, 1e-3, 1e-3, 1e-2, 1e-3]*si_unit(:darcy),
    porosity = [0.01, 0.01, 0.01, 0.05, 0.01],
    rock_thermal_conductivity = [2.5, 2.5, 2.8, 3.5, 2.5]*watt/(meter*si_unit(:Kelvin)),
    rock_heat_capacity = [900, 900, 900, 900, 900]*si_unit(:joule)/(si_unit(:kilogram)*si_unit(:Kelvin)),
);

# ## Inspect model
# Visualize the computational mesh and well configuration. The mesh is refined
# around the wells to accurately capture thermal and hydraulic processes.
msh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(msh)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, aspect = :data, perspectiveness = 0.5,
    title = "Multi-branch coaxial geothermal system")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
wells = get_model_wells(case.model)
for (name, well) in wells
    plot_well!(ax, msh, well)
end
fig

# ### Plot reservoir properties
# Visualize the layered reservoir properties interactively.
plot_reservoir(case.model; aspect = :data)

# ## Simulate geothermal energy production
# We simulate the multi-branch system for 30 years. The system injects cooled
# water at 25°C and extracts heated water from the branches.
sim, cfg = setup_reservoir_simulator(case;
    output_substates = true,
    info_level = 0,
    initial_dt = 5.0,
    presolve_wells = true,
    relaxation = true);

# Add temperature-based timestep control
sel = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel);
sel = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :CoaxialBranch_supply)
push!(cfg[:timestep_selectors], sel);

# Run the simulation
results = simulate_reservoir(case; simulator = sim, config = cfg)

# ## Visualize results
# Examine the thermal depletion pattern around the multi-branch system
Δstates = JutulDarcy.delta_state(results.states, case.state0[:Reservoir])
plot_reservoir(case.model, Δstates;
    resolution = (600, 800), aspect = :data,
    colormap = :seaborn_icefire_gradient, key = :Temperature,
    well_arg = (markersize = 0.0,),
)

# ### Well performance
# Examine well output including flow rates, pressures, and temperatures.
plot_well_results(results.wells)
