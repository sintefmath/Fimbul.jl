# # Deep coaxial geothermal well
# This example demonstrates simulation and analysis of geothermal energy
# production from a deep coaxial closed-loop well. The well is defined by a
# general trajectory (m×3 matrix) and uses coaxial heat exchange with the
# surrounding rock formation. We compare two flow configurations: injection
# into the supply (inner pipe) side vs. the return (outer annulus) side.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, hour, watt = si_units(:meter, :hour, :watt);

# ## Set up coaxial well system
# We create a single coaxial well extending to 2500 m depth. The default
# trajectory is a vertical well. The well circulates water at 25 m³/h with
# an injection temperature of 25°C.

# Define common parameters
common_args = (
    well_trajectory = [
        0.0    0.0      0.0;
        0.0    0.0   2500.0;
    ],
    rate = 25meter^3/hour,
    temperature_inj = convert_to_si(25.0, :Celsius),
    num_years = 30,
    report_interval = si_unit(:year)/4,
    depths = [0.0, 500.0, 1500.0, 2000.0, 3000.0],
    permeability = [1e-1, 1e-2, 1e-2, 1e-3]*si_unit(:darcy),
    porosity = [0.1, 0.01, 0.05, 0.01],
    rock_thermal_conductivity = [2.5, 2.8, 3.5, 2.5]*watt/(meter*si_unit(:Kelvin)),
    rock_heat_capacity = [900, 900, 900, 900]*si_unit(:joule)/(si_unit(:kilogram)*si_unit(:Kelvin)),
    rock_density = [2600, 2600, 2600, 2600]*si_unit(:kilogram)/si_unit(:meter)^3,
    hxy_min = 5.0,
    mesh_args = (offset = 250.0, offset_rel=missing),
);

# ### Case 1: Inject into supply side (inner pipe)
# Cold water goes down the inner pipe, heats up at the bottom, and returns
# through the outer annulus.
case_supply = coaxial_well_branches(; inject_into = :supply, common_args...);

# ### Case 2: Inject into return side (outer annulus)
# Cold water goes down the outer annulus, heats up along the way, and returns
# through the inner pipe.
case_return = coaxial_well_branches(; inject_into = :return, common_args...);

# ## Inspect model
# Visualize the computational mesh and well configuration. The mesh is refined
# around the well to accurately capture thermal and hydraulic processes.
msh = physical_representation(reservoir_model(case_supply.model).data_domain)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, perspectiveness = 0.5, aspect=(1,1,4),
    title = "Deep coaxial geothermal well")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
wells = get_model_wells(case_supply.model)
for (name, well) in wells
    plot_well!(ax, msh, well)
end
fig

# ### Plot reservoir properties
# Visualize the layered reservoir properties interactively.
plot_reservoir(case_supply.model)

# ## Simulate both configurations
# Common simulator settings
function run_case(case)
    sim, cfg = setup_reservoir_simulator(case;
        tol_dp_well = 1e-2,
        output_substates = true,
        info_level = 2,
        initial_dt = 5.0,
        presolve_wells = true,
        relaxation = true)
    sel = VariableChangeTimestepSelector(:Temperature, 5.0;
        relative = false, model = :Reservoir)
    push!(cfg[:timestep_selectors], sel)
    sel = VariableChangeTimestepSelector(:Temperature, 5.0;
        relative = false, model = :CoaxialWell_supply)
    push!(cfg[:timestep_selectors], sel)
    return simulate_reservoir(case; simulator = sim, config = cfg)
end

# ### Simulate injection into supply side
results_supply = run_case(case_supply)

# ### Simulate injection into return side
results_return = run_case(case_return)

# ## Compare results
# ### Thermal depletion patterns
fig_cmp = Figure(size = (1200, 800))
for (i, (results, case, label)) in enumerate(zip(
    [results_supply, results_return],
    [case_supply, case_return],
    ["Inject into supply", "Inject into return"]))
    Δstates = JutulDarcy.delta_state(results.states, case.state0[:Reservoir])
    ax = Axis3(fig_cmp[1, i]; zreversed = true, perspectiveness = 0.5,
        aspect=(1,1,4), title = label)
    Jutul.plot_cell_data!(ax, msh, Δstates[end][:Temperature];
        colormap = :seaborn_icefire_gradient)
end
fig_cmp

# ### Well performance comparison
for (results, label) in zip(
    [results_supply, results_return],
    ["Supply injection", "Return injection"])
    plot_well_results(results.wells)
end
