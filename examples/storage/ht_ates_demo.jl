# # High-Temperature Aquifer Thermal Energy Storage (HT-ATES) with H2OSystem
# <tags: Storage, ATES, HT-ATES>
# This example demonstrates setting up and simulating a High-Temperature Aquifer
# Thermal Energy Storage (HT-ATES) system using the `H2OSystem` two-phase water
# model.
#
# Unlike the standard `ates()` function (which uses a simplified single-phase
# water model with constant PVT properties), `ht_ates()` uses steam tables to
# compute accurate fluid properties at elevated temperatures via the
# pressure-enthalpy formulation. This enables physically consistent simulations
# for HT-ATES systems operating at temperatures up to and beyond 100°C.
#
# The system is charged at 200°C and discharged at 50°C with a rate of 100 l/s.
# The aquifer is located at 500 m depth.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using GLMakie
# Useful SI units
Kelvin, joule, watt = si_units(:Kelvin, :joule, :watt)
kilogram = si_unit(:kilogram)
meter = si_unit(:meter)
darcy = si_unit(:darcy);

# ## Set up the HT-ATES case
# The `ht_ates()` function is based on `ates()` but uses the `H2OSystem`
# two-phase water model (pressure-enthalpy formulation). The injection and
# boundary enthalpies are automatically computed from the steam tables to be
# consistent with the specified temperatures and hydrostatic pressures.
#
# Key parameters:
# - Charge temperature: 200°C (hot injection during summer)
# - Discharge temperature: 50°C (cold reinjection during winter)
# - Rate: 100 l/s
# - Aquifer at 500 m depth, 100 m thick
num_years = 5
case = Fimbul.ht_ates(;
    temperature_charge    = convert_to_si(200, :Celsius),  # 200°C charge temperature
    temperature_discharge = convert_to_si(50, :Celsius),   # 50°C discharge temperature
    rate_charge           = 100*si_unit(:litre)/si_unit(:second), # 100 l/s
    depths = [0.0, 450.0, 500.0, 600.0, 650.0, 900.0],    # Aquifer at 500 m depth
    aquifer_layer = 3,
    utes_schedule_args = (num_years = num_years,),
);

# ## Inspect the model
# We visualize the mesh structure and well configuration.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], zreversed = true, aspect = :data,
    title = "HT-ATES system: mesh structure and well configuration")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
wells = get_model_wells(case.model)
colors = [:red, :blue]
for (i, (k, w)) in enumerate(wells)
    plot_well!(ax, msh, w, color = colors[i], linewidth = 6)
end
plot_cell_data!(ax, msh, case.input_data[:layers],
    colormap = :rainbow, alpha = 0.3)
fig

# ### Check initial state
# Verify that the initial enthalpy and temperature are physically consistent.
# The enthalpy is derived from the steam tables at hydrostatic pressure and
# geothermal temperature.
T_init = convert_from_si.(case.state0[:Reservoir][:Temperature], :Celsius)
@show extrema(T_init)

# ## Simulate the HT-ATES system
sim, cfg = setup_reservoir_simulator(case; info_level = 0);
sel = JutulDarcy.ControlChangeTimestepSelector(case.model)
push!(cfg[:timestep_selectors], sel)
cfg[:timestep_max_decrease] = 1e-3
results = simulate_reservoir(case, simulator = sim, config = cfg)

# ## Visualize results
# Interactive visualization of the temperature field evolution.
plot_reservoir(case, results.states,
    key = :Temperature,
    aspect = :data,
    colormap = :seaborn_icefire_gradient)

# ### Inspect well output
plot_well_results(results.wells)

# ### Aquifer temperature evolution
# Track the mean aquifer temperature over time to assess thermal impact.
layers = case.input_data[:layers]
aq_layer = case.input_data[:aquifer_layer]
cells = layers .== aq_layer
times = convert_from_si.(cumsum(case.dt), :year)
T_aquifer = [state[:Temperature][cells] for state in results.states]
T_mean = [mean(convert_from_si.(T, :Celsius)) for T in T_aquifer]
T_min  = [minimum(convert_from_si.(T, :Celsius)) for T in T_aquifer]
T_max  = [maximum(convert_from_si.(T, :Celsius)) for T in T_aquifer]

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1],
    title = "Aquifer temperature evolution",
    xlabel = "Time (years)",
    ylabel = "Temperature (°C)",
    limits = ((times[1], times[end]), nothing),
)
band!(ax, times, T_min, T_max, color = (:steelblue, 0.3), label = "Min–Max range")
lines!(ax, times, T_mean, linewidth = 3, color = :steelblue, label = "Mean")
axislegend(ax)
fig
