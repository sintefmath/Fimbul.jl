# # Deep coaxial geothermal well
# This example demonstrates simulation and analysis of geothermal energy
# production from a deep coaxial closed-loop well. The well is defined by a
# general trajectory (m×3 matrix) and uses coaxial heat exchange with the
# surrounding rock formation.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, hour, watt = si_units(:meter, :hour, :watt);

# ## Set up coaxial well system
# We create a single coaxial well extending to 2500 m depth. The default
# trajectory is a vertical well. The well circulates water at 50 m³/h with
# an injection temperature of 25°C.

# Define a custom deviated trajectory (or use the default vertical one) The
# trajectory extends vertically down to 2500 m, then bends laterally with a
# radius of 150 m.

trajectory = [
    0.0    0.0    0.0;
    0.0    0.0   2500.0;
]

case = coaxial_well_branches(;
    well_trajectory = trajectory,                          # m×3 well path
    rate = 50meter^3/hour,                                 # Circulation rate
    temperature_inj = convert_to_si(25.0, :Celsius),       # Injection temperature
    num_years = 30,                                        # Simulation duration
    report_interval = si_unit(:year)/4,                    # Output 4x per year
    # Layered reservoir properties
    depths = [0.0, 500.0, 1500.0, 2000.0, 2500.0, 3000.0],
    # depths = [0.0, 50.0, 150.0, 200.0, 250.0, 1200.0],
    permeability = [1e-3, 1e-3, 1e-3, 1e-2, 1e-3]*si_unit(:darcy),
    porosity = [0.01, 0.01, 0.01, 0.05, 0.01],
    rock_thermal_conductivity = [2.5, 2.5, 2.8, 3.5, 2.5]*watt/(meter*si_unit(:Kelvin)),
    rock_heat_capacity = [900, 900, 900, 900, 900]*si_unit(:joule)/(si_unit(:kilogram)*si_unit(:Kelvin)),
    mesh_args = (offset = 25.0, offset_rel=missing), # Mesh refinement parameters
);

# ## Inspect model
# Visualize the computational mesh and well configuration. The mesh is refined
# around the well to accurately capture thermal and hydraulic processes.
msh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(msh)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, perspectiveness = 0.5,
    title = "Deep coaxial geothermal well")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
wells = get_model_wells(case.model)
for (name, well) in wells
    plot_well!(ax, msh, well)
end
fig

# ### Plot reservoir properties
# Visualize the layered reservoir properties interactively.
plot_reservoir(case.model)

# ## Simulate geothermal energy production
# We simulate the coaxial well system for 30 years. The system injects cooled
# water at 25°C and extracts heated water.
sim, cfg = setup_reservoir_simulator(case;
    tol_dp_well = 1e-2,
    output_substates = true,
    info_level = 2,
    initial_dt = 5.0,
    # presolve_wells = true,
    relaxation = true);

# Add temperature-based timestep control
sel = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel);
sel = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :CoaxialWell_supply)
push!(cfg[:timestep_selectors], sel);

# Run the simulation
results = simulate_reservoir(case; simulator = sim, config = cfg)

# ## Visualize results
# Examine the thermal depletion pattern around the coaxial well
Δstates = JutulDarcy.delta_state(results.states, case.state0[:Reservoir])
plot_reservoir(case.model, Δstates;
    resolution = (600, 800),
    colormap = :seaborn_icefire_gradient, key = :Temperature,
    well_arg = (markersize = 0.0,),
)

# ### Well performance
# Examine well output including flow rates, pressures, and temperatures.
plot_well_results(results.wells)
