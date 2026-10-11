# # High-enthalpy geothermal system above a magmatic intrusion
# <tags: Production>
# This example demonstrates simulation of a high-enthalpy geothermal system
# heated by a cooling magmatic intrusion, inspired by the production and
# reinjection scenarios of [yapparova_cold_2023](@cite). The fluid is pure water,
# modelled with the pressure-enthalpy formulation in Fimbul (`H2OSystem`), so
# the reservoir may contain liquid, vapor, two-phase and supercritical water.
#
# The example consists of two parts:
# 1. **Natural state**: A 900 °C intrusion is emplaced 2–3 km below the surface
#    and left to develop a hydrothermal convection system for 2400 years.
# 2. **Production**: Starting from the natural state, we compare 100 years of
#    production with and without reinjection of cold water.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, year = si_units(:meter, :year);

# ## Set up the natural-state case
# The domain is 10 km × 10 km × 5 km, with an ellipsoidal intrusion of
# horizontal radius 1.5 km centered horizontally. We use a relatively coarse
# mesh to keep the runtime moderate. Fluid can flow along the wellbore of a shut
# well, so we simulate the natural state without wells by setting `num_years =
# 0`. The wells are added when we set up the production scenarios.
dims = (45, 45, 40)
case_args = (
    dims = dims,
    refinement_margin = 1000.0*meter,
    injector_position = (0.0, 0.0).*meter,
)
case_ns = magmatic_intrusion(; case_args...,
    natural_state_time = 2400.0year,
    num_years = 0,
)
model = reservoir_model(case_ns.model)
mesh = physical_representation(model.data_domain);

# ### Plot reservoir properties
# We first inspect the model interactively. A low-permeability cap rock covers
# the top 500 m of the domain.
xc = model.data_domain[:cell_centroids]
x_mid = sum(extrema(xc[1, :]))/2
y_mid = sum(extrema(xc[2, :]))/2
cutaway = .!(xc[1, :] .< x_mid .&& xc[2, :] .< y_mid)
axis_args = (zreversed = true, aspect = :data, perspectiveness = 0.75, elevation = π/16)
plot_res_args = (
    resolution = (1000, 800), aspect = :data,
    # well_arg = (markersize = 0.0, ),
    axis_args = axis_args,
    cells = cutaway,
    fancy = false
)
plot_reservoir(case_ns; key = :Temperature, plot_res_args...)

# ### Initial temperature
# The initial temperature follows a linear geothermal gradient of 37.5 °C/km,
# with the intrusion at 900 °C. For the 3D plots in this example, we cut away
# the quadrant of the domain facing the viewer, so that the cut planes pass
# through the center of the intrusion.

fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1]; title = "Initial temperature", axis_args...)
T0 = convert_from_si.(case_ns.state0[:Reservoir][:Temperature], :Celsius)
plt = plot_cell_data!(ax, mesh, T0;
    cells = cutaway, colormap = :seaborn_icefire_gradient)
Colorbar(fig[1, 2], plt; label = "T (°C)")
fig

# ## Simulate the natural state
# We use the same solver settings for all simulations in this example: we remove the default maximum
# timestep of one year, since the system evolves slowly once the convection
# pattern is established, use relaxation of the Newton updates, and use
# slightly relaxed convergence tolerances. Close to the critical point of
# water, the fluid properties change rapidly, and these settings help the
# nonlinear solver through this region.
function run_case(case)
    sim, cfg = setup_reservoir_simulator(case;
        max_timestep = Inf,
        relaxation = true,
        tol_cnv = 1e-2,
        tol_mb = 1e-5,
        info_level = 2
    )
    tss = VariableChangeTimestepSelector(:Enthalpy, 1e5;
    relative = false, model = :Reservoir)
    push!(cfg[:timestep_selectors], tss)
    return simulate_reservoir(case; simulator = sim, config = cfg)
end

results_ns = run_case(case_ns);

# ### Interactive visualization of the natural state
# The interactive viewer shows the temperature at each report step of the
# natural state. Filtering out low temperatures shows how the hot plume
# develops above the intrusion.
plot_reservoir(case_ns.model, results_ns.states;
    key = :Temperature, colormap = :seaborn_icefire_gradient, cells = cutaway, plot_res_args...)

# ### Development of the convection plume
# The hot intrusion heats the surrounding water, which rises buoyantly and is
# replaced by colder water flowing in from the sides. We plot the temperature
# at selected times.
t_ns = cumsum(case_ns.dt)
nearest_step(t) = argmin(abs.(t_ns .- t))
plot_times = [300.0, 600.0, 1200.0].*year
T_range = (10.0, 900.0)
fig = Figure(size = (1000, 750))
for (i, t) in enumerate(plot_times)
    step = nearest_step(t)
    T = convert_from_si.(results_ns.states[step][:Temperature], :Celsius)
    ax = Axis3(fig[1, i];
        title = "$(round(Int, t_ns[step]/year)) years", axis_args...)
    plot_cell_data!(ax, mesh, T; cells = cutaway,
        colormap = :seaborn_icefire_gradient, colorrange = T_range)
    hidedecorations!(ax)
    S = vec(results_ns.states[step][:Saturations][1,:])
    ax = Axis3(fig[2, i]; axis_args...)
    plot_cell_data!(ax, mesh, S; cells = cutaway,
        colormap = :batlowK, colorrange = (0.0, 1.0))
    hidedecorations!(ax)
end
Colorbar(fig[3, 1:2]; colormap = :seaborn_icefire_gradient, colorrange = T_range,
    label = "T (°C)", vertical = false)
fig

# ### Phase diagram
# The evolution of the fluid can be shown in a pressure-enthalpy phase diagram,
# with temperature contours and the two-phase envelope. We follow three cells
# on the vertical line through the center of the intrusion: one at the center
# of the intrusion, one just above it, and one higher up in the convection
# plume. Each line shows the state of one cell at every timestep of the natural
# state, starting from the initial state (circle). Early in the natural state,
# the fluid passes close to the critical point of water (22.064 MPa,
# 2.085 MJ/kg), where the fluid properties change rapidly. This is the most
# challenging part of the simulation.
x_c = (sum(extrema(xc[1, :]))/2, sum(extrema(xc[2, :]))/2)
probe_depths = [2500.0, 1900.0, 1300.0, 700.0].*meter
probe_cells = [argmin((xc[1, :] .- x_c[1]).^2 .+ (xc[2, :] .- x_c[2]).^2 .+
    (xc[3, :] .- d).^2) for d in probe_depths]
probe_labels = ["Intrusion center", "Above intrusion", "Upper plume", "Near surface"]
probe_colors = [:red, :lightgreen, :cyan, :magenta]
states_ph = vcat([case_ns.state0[:Reservoir]], results_ns.states)
tables = model.system.pvt_tables
fig = Figure(size = (900, 650), fontsize = 14)
ax = Axis(fig[1, 1])
handles = plot_phase_diagram_contours!(ax, tables;
    pressure_limits = (0.0, 30e6), enthalpy_limits = (0.0, 4.5e6), lines = true)
for (i, (c, label)) in enumerate(zip(probe_cells, probe_labels))
    h = [state[:Enthalpy][c] for state in states_ph]./1e3
    p = [state[:Pressure][c] for state in states_ph]./1e6
    depth = round(Int, xc[3, c])
    lines!(ax, h, p; color = probe_colors[i], linewidth = 3, label = "$label ($depth m)")
    scatter!(ax, h, p; color = probe_colors[i], markersize = 10,
        strokecolor = :black, strokewidth = 1)
    scatter!(ax, h[1:1], p[1:1]; color = :yellow, markersize = 12,
        strokecolor = :black, strokewidth = 1, marker = :star5)
end

axislegend(ax; position = :rb)
Colorbar(fig[1, 2], handles.filled; label = "T (°C)")
fig

# ## Production scenarios
# We set up the production scenarios with wells, starting from the end of the
# natural state (`initial_state`), and simulate 100 years of production. The
# wells are multisegment wells that extend from the surface, and the well
# pressures are prescribed at the wellhead. The producer is open between 750
# and 1750 m depth, above the intrusion, and is operated at a wellhead pressure
# of 10 bar. In the second scenario, cold water (80 °C) is reinjected at a
# wellhead pressure of 80 bar (about 300 bar at the bottom of the well) in an
# injector that is open between 2050 and 2300 m depth, in the top of the
# intrusion. With
# `injection_years = 0`, the injector is not included in the model. We use
# report steps of one year.
state_ns = results_ns.states[end]
case_inj = magmatic_intrusion(; case_args...,
    natural_state_time = 0.0,
    num_years = 100,
    report_interval = 1.0year,
    initial_state = state_ns,
)
case_noinj = magmatic_intrusion(; case_args...,
    natural_state_time = 0.0,
    num_years = 100,
    injection_years = 0,
    report_interval = 1.0year,
    initial_state = state_ns,
)
wells = get_model_wells(case_inj.model)
function plot_wells!(ax)
    for well in values(wells)
        plot_well!(ax, mesh, well;
            color = :black, linewidth = 3, markersize = 0.0, fontsize = 0.0)
    end
end

# ### Simulate
results_noinj = run_case(case_noinj);
results_inj = run_case(case_inj);

# ### Well performance
# The produced fluid boils as it rises in the producer, and reaches the
# wellhead as a mixture of water and steam. The production temperature is
# therefore fixed at the saturation temperature at the wellhead pressure of
# 10 bar (about 180 °C) throughout the simulation. The specific enthalpy of the produced
# fluid, which determines the steam fraction, varies only slightly. We
# therefore show the power output, computed as the produced mass rate times
# the specific enthalpy of the produced fluid, which closely follows the
# production rate. The power output drops quickly in the first years as the
# pressure around the producer declines. Reinjection supports the reservoir
# pressure, which sustains a higher production rate and power output.
fig = Figure(size = (800, 450))
ax = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "Power output (MW)")
colors = Makie.wong_colors()
for (i, (res, label)) in enumerate(zip(
        [results_noinj, results_inj], ["Production only", "With reinjection"]))
    t = res.time./year
    q = -res.wells[:Producer][:mass_rate]
    h = [state[:Producer][:Enthalpy][1] for state in res.result.states]
    lines!(ax, t, q.*h./1e6; color = colors[i], linewidth = 2, label = label)
end
axislegend(ax; position = :rt)
fig

# ### Interactive visualization of temperature changes
# The interactive viewer shows the temperature change relative to the natural
# state for the scenario with reinjection. Filtering out values close to zero
# shows the cold zone developing around the injector.
Δstates_inj = JutulDarcy.delta_state(results_inj.states, results_ns.states[end])
plot_reservoir(case_inj.model, Δstates_inj;
    key = :Temperature, colormap = :seaborn_icefire_gradient, plot_res_args...)

# ### Reservoir temperature change
# Finally, we plot the temperature change after 100 years of production for
# both scenarios, showing only cells where the temperature change exceeds 10% of
# the largest change in that scenario. The region along the producer cools, while hot fluid drawn up from
# the plume warms the rock next to it. Reinjection creates a cold zone around
# the injector, which has not reached the producer after 100 years.
T_ns = results_ns.states[end][:Temperature]
ΔT_all = [res.states[end][:Temperature] .- T_ns for res in (results_noinj, results_inj)]
ΔT_max = maximum(maximum.(abs, ΔT_all))
ΔT_range = (-ΔT_max, ΔT_max)
lims = (extrema(xc[1, :])..., extrema(xc[2, :])..., extrema(xc[3, :])...)
fig = Figure(size = (1200, 600))
for (i, (ΔT, title)) in enumerate(zip(ΔT_all, ["Production only", "With reinjection"]))
    ax = Axis3(fig[1, i]; title = title, limits = lims, axis_args...)
    plot_cell_data!(ax, mesh, ΔT; cells = abs.(ΔT) .> 0.1*maximum(abs, ΔT) .&& cutaway,
        colormap = :seaborn_icefire_gradient) #, colorrange = ΔT_range)
    plot_wells!(ax)
end
Colorbar(fig[2, 1:2]; colormap = :seaborn_icefire_gradient, colorrange = ΔT_range,
    label = "ΔT (°C)", vertical = false)
fig
