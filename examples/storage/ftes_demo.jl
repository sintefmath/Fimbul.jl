# # Fractured Thermal Energy Storage (FTES)
# This example demonstrates how to set up and simulate a Fractured Thermal
# Energy Storage (FTES) system using Fimbul. FTES systems exploit natural or
# induced fractures in low-permeability rock to store and recover thermal energy
# by circulating water through a fracture network connected by an injector–
# producer well doublet arrangement.
#
# The system consists of one central injector surrounded by multiple producer
# wells. Horizontal fracture planes connect the injector to the producers,
# enabling thermal transport through the fracture network even in tight rock.
# During charging, hot water is injected through the central well and produced
# at the outer producers; during discharging the flow direction is reversed.

using Jutul, JutulDarcy, Fimbul
using HYPRE
using Random
using GLMakie

# ## Set up simulation case
# We create an FTES system with 8 producer wells arranged in a circle of 25 m
# radius around the central injector. The wells extend to 300 m depth and the
# fracture network consists of 25 near-horizontal fractures distributed within
# the well interval. The system is charged from April to November and
# discharged from December to March over a 3-year period.
Random.seed!(20260225)
case = Fimbul.ftes(
    (num_producers = 8, radius = 25.0, depth = 300.0),
    25;
    rate_charge = 50si"litre/second",
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = (num_years = 3,),
);

# ## Visualize the FTES system
# We first inspect the computational mesh and the embedded fracture network.
# The fracture mesh (gray surfaces) captures the codimension-one fracture
# planes that are resolved inside the 3D matrix mesh.
matrix_mesh = physical_representation(reservoir_model(case.model).data_domain)
fracture_mesh = physical_representation(case.model.models[:Fractures].data_domain)

fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1]; perspectiveness = 0.5, zreversed = true, aspect = :data,
    title = "FTES system: matrix mesh and fracture network")
Jutul.plot_mesh!(ax, fracture_mesh; color = :gray, alpha = 0.6)
Jutul.plot_mesh_edges!(ax, matrix_mesh; alpha = 0.1)
wells = get_model_wells(case.model)
colors = [:red; fill(:blue, length(wells) - 1)]
for (i, (k, w)) in enumerate(wells)
    plot_well!(ax, matrix_mesh, w; color = colors[i], linewidth = 5)
end
fig

# ## Set up reservoir simulator
# We configure solver tolerances suited to the thermal DFM system. The
# `ControlChangeTimestepSelector` is used to take very small steps when well
# controls switch between charging and discharging to maintain convergence.
simulator, config = setup_reservoir_simulator(case;
    tol_cnv = 1e-2,
    tol_mb = 1e-5,
    tol_dp_well = 1e-2,
    tol_cnv_well = Inf,
    tol_cnve_well = Inf,
    inc_tol_dT = 1e-2,
    inc_tol_dp_abs = 1e-2 * si_unit(:bar),
    initial_dt = 5.0,
    relaxation = true,
);

sel = JutulDarcy.ControlChangeTimestepSelector(case.model, 0.1, 5.0)
push!(config[:timestep_selectors], sel)
sel_T = VariableChangeTimestepSelector(:Temperature, 10.0; model = :Fractures, relative = false)
push!(config[:timestep_selectors], sel_T)
config[:timestep_max_decrease] = 1e-6;

# ## Simulate the FTES system
# Run the full multi-year simulation. Setting `info_level = 0` will show a
# progress bar rather than per-iteration solver output.
results = simulate_reservoir(case; simulator = simulator, config = config, info_level = 0);

# ## Visualize results
# ### Matrix temperature distribution
# Inspect the temperature distribution in the matrix reservoir after the final
# simulated timestep to see how thermal energy has spread around the fracture
# network.
states_m = [s[:Reservoir] for s in results.states]
plot_reservoir(case.model, states_m;
    key = :Temperature,
    aspect = :data,
    colormap = :seaborn_icefire_gradient)

# ### Fracture temperature distribution
# The fractures are the primary heat transport pathway. We visualize the
# temperature field on the fracture mesh at the end of the simulation.
states_f = [s[:Fractures] for s in results.states]
plot_reservoir(case.model.models[:Fractures], states_f;
    key = :Temperature,
    aspect = :data,
    colormap = :seaborn_icefire_gradient)

# ### Well performance over time
# Plot injection/production temperatures and flow rates throughout the
# operational schedule to assess thermal efficiency and system performance.
plot_well_results(results.wells, field = :temperature)
