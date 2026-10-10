# # High-enthalpy geothermal system above a magmatic intrusion
# <tags: Production>
# This example demonstrates simulation of a high-enthalpy geothermal system
# heated by a cooling magmatic intrusion, inspired by the production and
# reinjection scenarios of [yapparova_cold_2023](@cite). The fluid is pure water,
# modelled with the pressure-enthalpy formulation in Fimbul (`H2OSystem`), so
# the reservoir may contain liquid, vapor, two-phase and supercritical water.
#
# The example consists of two parts:
# 1. **Natural state**: A 500 °C intrusion is emplaced 2–3 km below the surface
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
# The domain is 10 km × 10 km × 5 km, with a cylindrical intrusion of radius
# 1.5 km centered horizontally. We use a relatively coarse mesh to keep the
# runtime moderate. Fluid can flow along the wellbore of a shut well, so we
# simulate the natural state without wells by setting `num_years = 0`. The
# wells are added when we set up the production scenarios.
dims = (51, 51, 40)
case_args = (
    dims = dims,
    injector_position = (-100.0, 250.0).*meter,
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
plot_res_args = (
    resolution = (1000, 800), aspect = :data,
    well_arg = (markersize = 0.0, ),
    axis_args = (perspectiveness = 0.5, ),
    fancy = false
)
plot_reservoir(case_ns.model; key = :permeability, plot_res_args...)

# ### Initial temperature
# The initial temperature follows a linear geothermal gradient of 37.5 °C/km,
# with the intrusion at 500 °C. For the 3D plots in this example, we cut away
# the quadrant of the domain facing the viewer, so that the cut planes pass
# through the center of the intrusion.
xc = model.data_domain[:cell_centroids]
x_mid = sum(extrema(xc[1, :]))/2
y_mid = sum(extrema(xc[2, :]))/2
cutaway = .!(xc[1, :] .< x_mid .&& xc[2, :] .< y_mid)
axis_args = (zreversed = true, aspect = :data, perspectiveness = 0.75, elevation = π/8)

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
    return simulate_reservoir(case;
        max_timestep = Inf,
        relaxation = true,
        tol_cnv = 1e-2,
        tol_mb = 1e-5,
        info_level = 2)
end

results_ns = run_case(case_ns);

# ### Interactive visualization of the natural state
# The interactive viewer shows the temperature at each report step of the
# natural state. Filtering out low temperatures shows how the hot plume
# develops above the intrusion.
plot_reservoir(case_ns.model, results_ns.states;
    key = :Temperature, colormap = :seaborn_icefire_gradient, plot_res_args...)

# ### Development of the convection plume
# The hot intrusion heats the surrounding water, which rises buoyantly and is
# replaced by colder water flowing in from the sides. We plot the temperature
# at selected times.
t_ns = cumsum(case_ns.dt)
nearest_step(t) = argmin(abs.(t_ns .- t))
plot_times = [100.0, 300.0, 600.0, 2400.0].*year
T_range = (10.0, 500.0)
fig = Figure(size = (1000, 900))
for (i, t) in enumerate(plot_times)
    step = nearest_step(t)
    T = convert_from_si.(results_ns.states[step][:Temperature], :Celsius)
    ax = Axis3(fig[(i-1)÷2 + 1, (i-1)%2 + 1];
        title = "$(round(Int, t_ns[step]/year)) years", axis_args...)
    plot_cell_data!(ax, mesh, T; cells = cutaway,
        colormap = :seaborn_icefire_gradient, colorrange = T_range)
    hidedecorations!(ax)
end
Colorbar(fig[3, 1:2]; colormap = :seaborn_icefire_gradient, colorrange = T_range,
    label = "T (°C)", vertical = false)
fig

# ### Phase diagram
# The reservoir states can be shown in a pressure-enthalpy phase diagram,
# with temperature contours and the two-phase envelope. Early in the natural
# state, the fluid in and above the intrusion passes close to the critical
# point of water (22.064 MPa, 2.085 MJ/kg), where the fluid properties change
# rapidly. This is the most challenging part of the simulation.
fig = Figure(size = (1200, 500))
handles = nothing
for (i, t) in enumerate([300.0, 2400.0].*year)
    step = nearest_step(t)
    ax = Axis(fig[1, i]; title = "t = $(round(Int, t_ns[step]/year)) years")
    global handles = plot_reservoir_state_phase_diagram!(ax, case_ns.model, results_ns.states[step];
        pressure_limits = (0.0, 50e6), enthalpy_limits = (0.0, 3.2e6),
        state_kwargs = (type = :scatter, markersize = 4, color = :black))
end
Colorbar(fig[1, 3], handles.contours.filled; label = "T (°C)")
fig

# ## Production scenarios
# We set up the production scenarios with wells, starting from the end of the
# natural state (`initial_state`), and simulate 100 years of production. The
# producer is open between 750 and 1750 m depth, above the intrusion, and is
# operated at a bottom-hole pressure of 50 bar. In the second scenario, cold
# water (80 °C) is reinjected at 300 bar in an injector that is open between
# 2050 and 2300 m depth, in the top of the intrusion. With
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
# The production rate drops quickly in the first years as the pressure around
# the producer declines. Reinjection supports the reservoir pressure, which
# sustains a higher production rate. The injector is open below the producer
# and some distance away, so the injected cold water has only a minor effect
# on the production temperature within 100 years.
fig = Figure(size = (1200, 450))
ax_q = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "Production rate (kg/s)")
ax_T = Axis(fig[1, 2]; xlabel = "Time (years)", ylabel = "Production temperature (°C)")
colors = Makie.wong_colors()
for (i, (res, label)) in enumerate(zip(
        [results_noinj, results_inj], ["Production only", "With reinjection"]))
    t = res.time./year
    q = -res.wells[:Producer][:mass_rate]
    T = convert_from_si.(res.wells[:Producer][:temperature], :Celsius)
    lines!(ax_q, t, q; color = colors[i], linewidth = 2, label = label)
    lines!(ax_T, t, T; color = colors[i], linewidth = 2, label = label)
end
axislegend(ax_q; position = :rt)
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
    plot_cell_data!(ax, mesh, ΔT; cells = abs.(ΔT) .> 0.1*maximum(abs, ΔT),
        colormap = :seaborn_icefire_gradient, colorrange = ΔT_range)
    plot_wells!(ax)
end
Colorbar(fig[2, 1:2]; colormap = :seaborn_icefire_gradient, colorrange = ΔT_range,
    label = "ΔT (°C)", vertical = false)
fig
