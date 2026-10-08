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
meter, year, bar = si_units(:meter, :year, :bar);

# ## Set up the case
# The domain is 10 km × 10 km × 5 km, with a cylindrical intrusion of radius
# 1.5 km centered horizontally. We use a coarse mesh to keep the runtime
# moderate. The producer is open between 750 and 1750 m depth, above the
# intrusion, and the injector is open between 2050 and 2300 m depth, in the top
# of the intrusion.
dims = (51, 51, 40)
case = magmatic_intrusion(;
    dims = dims,
    injector_position = (-100.0, 250.0).*meter,
    natural_state_time = 2400.0year,
    num_years = 100,
)
n_ns = case.input_data[:natural_state_steps]
model = reservoir_model(case.model)
mesh = physical_representation(model.data_domain);

# ### Inspect initial conditions
# The initial temperature follows a linear geothermal gradient of 37.5 °C/km,
# with the intrusion at 500 °C. A low-permeability cap rock covers the top
# 500 m of the domain.
nx, ny, nz = dims
Lx, Ly, Lz = 10000.0, 10000.0, 5000.0
x = range(Lx/nx/2, Lx - Lx/nx/2, length = nx)
z = range(Lz/nz/2, Lz - Lz/nz/2, length = nz)
j_mid = cld(ny, 2)
cross_section(v) = reshape(v, nx, ny, nz)[:, j_mid, :]

fig = Figure(size = (1200, 450))
T0 = convert_from_si.(case.state0[:Reservoir][:Temperature], :Celsius)
ax = Axis(fig[1, 1]; title = "Initial temperature", xlabel = "x (m)",
    ylabel = "Depth (m)", yreversed = true)
hm = heatmap!(ax, x, z, cross_section(T0); colormap = :inferno)
Colorbar(fig[1, 2], hm; label = "T (°C)")
K = model.data_domain[:permeability]
K = K isa AbstractMatrix ? K[1, :] : K
ax = Axis(fig[1, 3]; title = "Permeability", xlabel = "x (m)",
    ylabel = "Depth (m)", yreversed = true)
hm = heatmap!(ax, x, z, cross_section(log10.(K)); colormap = :viridis)
Colorbar(fig[1, 4], hm; label = "log₁₀ K (m²)")
fig

# ## Simulate the natural state
# During the natural state, both wells are shut. We remove the default
# maximum timestep of one year, since the system evolves slowly once the
# convection pattern is established.

sim, cfg = setup_reservoir_simulator(case; max_timestep = Inf, info_level = 2, relaxation = true, tol_cnv = 1e-2, tol_mb = 1e-5)

case_ns = case[1:n_ns]
results_ns = simulate_reservoir(case_ns; simulator = sim, config = cfg);

# ### Development of the convection plume
# The hot intrusion heats the surrounding water, which rises buoyantly and is
# replaced by colder water flowing in from the sides. We plot temperature in a
# vertical cross-section through the intrusion center at selected times.
t_ns = cumsum(case_ns.dt)
plot_times = [100.0, 250.0, 500.0, 2400.0].*year
fig = Figure(size = (1200, 800))
for (i, t) in enumerate(plot_times)
    step = findfirst(t_ns .>= t - 1.0)
    T = convert_from_si.(results_ns.states[step][:Temperature], :Celsius)
    ax = Axis(fig[(i-1)÷2 + 1, (i-1)%2 + 1];
        title = "t = $(round(Int, t_ns[step]/year)) years",
        xlabel = "x (m)", ylabel = "Depth (m)", yreversed = true)
    heatmap!(ax, x, z, cross_section(T); colormap = :inferno, colorrange = (10, 500))
end
Colorbar(fig[1:2, 3]; colormap = :inferno, colorrange = (10, 500), label = "T (°C)")
fig

# ### Phase diagram
# The reservoir states can be shown in a pressure-enthalpy phase diagram,
# with temperature contours and the two-phase envelope. Early in the natural
# state, the fluid in and above the intrusion passes close to the critical
# point of water (22.064 MPa, 2.085 MJ/kg), where the fluid properties change
# rapidly. This is the most challenging part of the simulation.
fig = Figure(size = (1200, 500))
handles = nothing
for (i, t) in enumerate([250.0, 2400.0].*year)
    step = findfirst(t_ns .>= t - 1.0)
    ax = Axis(fig[1, i]; title = "t = $(round(Int, t_ns[step]/year)) years")
    global handles = plot_reservoir_state_phase_diagram!(ax, case.model, results_ns.states[step];
        pressure_limits = (0.0, 50e6), enthalpy_limits = (0.0, 3.2e6),
        state_kwargs = (type = :scatter, markersize = 4, color = :black))
end
Colorbar(fig[1, 3], handles.contours.filled; label = "T (°C)")
fig

# ## Production scenarios
# We restart from the end of the natural state and simulate 100 years of
# production. The producer is operated at a bottom-hole pressure of 50 bar.
# In the second scenario, cold water (80 °C) is reinjected at 300 bar. Both
# scenarios use the same model, so we only need to replace the forces.
state_ns = results_ns.result.states[end]
case_inj = case[(n_ns+1):length(case.dt)]
case_inj = JutulCase(case_inj.model, case_inj.dt, case_inj.forces;
    state0 = state_ns, parameters = case.parameters)
case_noinj = magmatic_intrusion(;
    dims = dims,
    injector_position = (-100.0, 250.0).*meter,
    natural_state_time = 2400.0year,
    num_years = 100,
    injection_years = 0,
)
case_noinj = case_noinj[(n_ns+1):length(case_noinj.dt)]
case_noinj = JutulCase(case_inj.model, case_noinj.dt, case_noinj.forces;
    state0 = state_ns, parameters = case.parameters);

# ### Simulate
results_noinj = simulate_reservoir(case_noinj; max_timestep = Inf, info_level = -1);
results_inj = simulate_reservoir(case_inj; max_timestep = Inf, info_level = -1);

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

# ### Reservoir temperature change
# Finally, we plot the temperature change after 100 years of production in
# vertical cross-sections through the producer (top) and the injector
# (bottom). Production draws hot fluid from the plume upwards and colder water
# in from the sides, while reinjection creates a cold zone around the
# injector.
T_ns = convert_from_si.(results_ns.states[end][:Temperature], :Celsius)
ΔT_all = [convert_from_si.(res.states[end][:Temperature], :Celsius) .- T_ns
    for res in (results_noinj, results_inj)]
ΔT_max = maximum(maximum.(abs, ΔT_all))
fig = Figure(size = (1200, 800))
for (i, (dT, title)) in enumerate(zip(ΔT_all, ["Production only", "With reinjection"]))
    ΔT = reshape(dT, nx, ny, nz)
    ## Rows through the producer and injector, respectively
    yp, yi = Ly/2 - 250.0, Ly/2 + 250.0
    for (k, (yw, wname)) in enumerate(zip([yp, yi], ["producer", "injector"]))
        jw = clamp(Int(floor(yw/(Ly/ny))) + 1, 1, ny)
        ax = Axis(fig[k, i]; title = "$title – through $wname",
            xlabel = "x (m)", ylabel = "Depth (m)", yreversed = true)
        heatmap!(ax, x, z, ΔT[:, jw, :]; colormap = Reverse(:RdBu),
            colorrange = (-ΔT_max, ΔT_max))
    end
end
Colorbar(fig[1:2, 3]; colormap = Reverse(:RdBu),
    colorrange = (-ΔT_max, ΔT_max), label = "ΔT (°C)")
fig
