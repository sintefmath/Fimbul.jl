# # Aquifer Thermal Energy Storage (ATES) System Simulation
# This example demonstrates how to simulate and visualize an Aquifer Thermal 
# Energy Storage (ATES) system. ATES systems store thermal energy by injecting
# hot water into an aquifer during charging periods and extracting it during
# discharging periods. The setup consists of two wells: a "hot" well for 
# injecting/producing hot water and a "cold" well for pressure support and
# cold water injection.

# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul
using HYPRE
using Statistics
using GLMakie
# Useful SI units
watt = si_unit(:watt)
joule = si_unit(:joule)
kilogram = si_unit(:kilogram)
meter = si_unit(:meter)
Kelvin = si_unit(:Kelvin)

# ## Set up ATES case
# We create a full 3D ATES case with the default five-layer system consisting of
# multiple geological layers. The aquifer layer (layer 3) is where thermal energy
# is stored through cyclic injection and production operations.
num_years = 5
case = Fimbul.ates(;
    well_distance = 400.0,                    # Distance between wells [m]
    temperature_charge = convert_to_si(85, :Celsius),   # Hot injection temperature
    temperature_discharge = convert_to_si(20, :Celsius), # Cold injection temperature
    depths = [0.0, 850.0, 900.0, 1000.0, 1050.0, 1300.0], # Depths delineating layers
    porosity = [0.01, 0.05, 0.35, 0.05, 0.01], # Layer porosities
    permeability = [1.0, 5.0, 1000.0, 5.0, 1.0].*1e-3.*si_unit(:darcy), # Layer permeabilities
    rock_thermal_conductivity = [2.5, 2.0, 1.9, 2.0, 2.5].*watt/(meter*Kelvin), # Layer thermal conductivities
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin), # Heat capacities
    aquifer_layer = 3, # Target aquifer for thermal storage (layer 3)
    utes_schedule_args = (num_years = num_years,)
);

# ## Inspect the ATES model
# First, we visualize the computational mesh and well locations. The mesh is
# refined around the wells and in the aquifer layer to accurately capture
# thermal transport processes.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 600))
ax = Axis3(fig[1, 1], zreversed = true, aspect = :data, 
    title = "ATES Mesh and Wells")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
wells = get_model_wells(case.model)
colors = [:red, :blue]  # Hot well = red, Cold well = blue
for (i, (k, w)) in enumerate(wells)
    plot_well!(ax, msh, w, color = colors[i], linewidth = 6)
end
fig

# ### Visualize reservoir properties
# Plot the reservoir properties interactively to understand the geological setup
plot_reservoir(case.model, key = :porosity, aspect = :data)

# ## ATES operational schedule
# ATES systems operate in cycles with distinct phases:
# - **Charging**: Inject hot water through hot well, produce through cold well
# - **Discharging**: Produce hot water through hot well, inject cold water through cold well  
# - **Rest**: No well activity, allows thermal equilibration
# 

# ## Simulate the ATES system
# We set up a simulator with appropriate tolerances for thermal problems and
# add a control change timestep selector to handle the transitions between
# operational phases smoothly.
sim, cfg = setup_reservoir_simulator(case; 
    info_level = 1,
    max_nonlinear_iterations = 15
);

# Add timestep selector for control changes (critical for ATES stability)
sel = JutulDarcy.ControlChangeTimestepSelector(case.model)
push!(cfg[:timestep_selectors], sel)
cfg[:timestep_max_decrease] = 1e-3

# Run the simulation (this may take several minutes for 3D case)
println("Starting ATES simulation...")
results = simulate_reservoir(case, simulator = sim, config = cfg)
println("Simulation completed successfully!")

# ## Visualize ATES results
# Plot the reservoir state evolution interactively. You can see how thermal
# plumes develop around each well during different operational phases.
plot_reservoir(case, results.states, 
    key = :Temperature, 
    aspect = :data,
    colormap = :seaborn_icefire_gradient)

# ### Plot well performance
# Examine the well responses including rates, pressures, and temperatures
plot_well_results(results.wells)

# ### Temperature evolution analysis
# Create a detailed analysis of temperature changes in the aquifer during
# different operational phases
states = results.states
times = convert_from_si.(cumsum(case.dt), :year)
geo = tpfv_geometry(msh)

# Extract aquifer temperatures over time based on layer information
domain = reservoir_model(case.model).data_domain
porosity = domain[:porosity]
aquifer_cells = porosity .> 0.3  # High porosity indicates aquifer layer

aquifer_temps = []
for state in states
    T_aquifer = state[:Temperature][aquifer_cells]
    push!(aquifer_temps, T_aquifer)
end

# Calculate thermal storage efficiency metrics
initial_temp = convert_from_si(states[1][:Temperature][aquifer_cells][1], :Celsius)
max_temp = maximum([maximum(convert_from_si.(T, :Celsius)) for T in aquifer_temps])
min_temp = minimum([minimum(convert_from_si.(T, :Celsius)) for T in aquifer_temps])

println("Aquifer temperature analysis:")
println("  Initial temperature: $(round(initial_temp, digits=1))°C")
println("  Maximum temperature: $(round(max_temp, digits=1))°C")
println("  Minimum temperature: $(round(min_temp, digits=1))°C")
println("  Temperature swing: $(round(max_temp - min_temp, digits=1))°C")

# ### Well temperature profiles over depth
# Plot temperature evolution in both wells as a function of depth
fig = Figure(size = (1400, 600))

# Temperature profiles in wells over time

# Plot temperature profiles for selected timesteps
states = results.result.states
n_steps = length(states)
ch_start = findall([true; diff([f[:Facility].control[:Hot] isa InjectorControl for f in case.forces]) .> 0])
ch_stop = findall([false; diff([f[:Facility].control[:Hot] isa InjectorControl for f in case.forces]) .< 0]).-1
dch_start = findall([false; diff([f[:Facility].control[:Hot] isa ProducerControl for f in case.forces]) .> 0])
dch_stop = findall([false; diff([f[:Facility].control[:Hot] isa ProducerControl for f in case.forces]) .< 0]).-1

key_steps = sort(vcat(ch_start, dch_stop))
colors = cgrad(:seaborn_icefire_gradient, 10, categorical = true)
color_hot = colors[8]
color_cold = colors[3]

# for (i, step) in enumerate(key_steps)
#     if step <= length(states)
#         T_hot = convert_from_si.(states[step][:Hot][:Temperature], :Celsius)
#         T_cold = convert_from_si.(states[step][:Cold][:Temperature], :Celsius)
        
#         plot_mswell_values!(ax1, case.model, :Hot, T_hot;
#             geo = geo, linewidth = 3, color = colors[i], 
#             label = "Hot ($(round(times[step], digits=1)) yr)")
#         plot_mswell_values!(ax1, case.model, :Cold, T_cold;
#             geo = geo, linewidth = 3, color = colors[i], linestyle = :dash,
#             label = "Cold ($(round(times[step], digits=1)) yr)")
#     end
# end

ax_rate = Axis(fig[1, 1], 
    xlabel = "Time (years)", 
    ylabel = "Rate (m³/h)",
    limits = ((times[1], times[end]), nothing),
)

rate = abs.(results.wells[:Hot][:rate]*si_unit(:hour))
lines!(ax_rate, times, rate, 
    color = color_hot, linewidth = 3, label = "Hot well")

ax_temp = Axis(fig[2, 1], 
    xlabel = "Time (years)", 
    ylabel = "Temperature (°C)",
    limits = ((times[1], times[end]), nothing),
)

temp = convert_from_si.(results.wells[:Hot][:temperature], :Celsius)

lines!(ax_temp, times, temp, 
    color = color_hot, linewidth = 3, label = "Hot well")

mrate = results.wells[:Hot][:mass_rate]
Cp = mean(reservoir_model(case.model).data_domain[:component_heat_capacity])
effect = mrate.*Cp.*temp
effect_mwh = abs.(effect).*1e-6  # Convert to MW

ax_effect = Axis(fig[3, 1], 
    xlabel = "Time (years)", 
    ylabel = "Effect (MW)",
    limits = ((times[1], times[end]), nothing),
)

lines!(ax_effect, times, effect_mwh,  
    color = color_hot, linewidth = 3, label = "Hot well")

ax_eta = Axis(fig[4, 1], 
    xlabel = "Time (years)", 
    ylabel = "Efficiency (–)",
    limits = ((times[1], times[end]), nothing),
)
energy = effect.*case.dt
for k in eachindex(ch_start)
    # start = ch_start[k]
    # stop = dch_end[k]
    energy_ch = sum(energy[ch_start[k]:ch_stop[k]])
    energy_dch = abs.(cumsum(energy[dch_start[k]:dch_stop[k]]))
    η = energy_dch./energy_ch
    lines!(ax_eta, times[dch_start[k]:dch_stop[k]], η, 
        color = color_hot, linewidth = 3, label = "Cumulative charge")
end

fig

# ### 3D thermal plume evolution
# Show how thermal plumes develop and interact during different phases
# We'll plot temperature changes from initial state at key time points

# Calculate temperature changes from initial state
Δstates = []
T0 = case.state0[:Reservoir][:Temperature]
for state in results.states
    ΔT = state[:Temperature] .- T0
    push!(Δstates, Dict(:Temperature => ΔT))
end

# Create 3D visualization of thermal plume evolution
fig = Figure(size = (1200, 800))
T_range = (Inf, -Inf)  # Temperature change range for consistent colorscale

geo = tpfv_geometry(msh)
cells_back = geo.cell_centroids[2,:] .> 0

for (i, step) in enumerate(key_steps[end-1:end])
    ΔT = Δstates[step][:Temperature]
    # Create a mask to show only significant temperature changes
    cells_to_show = cells_back .&& abs.(ΔT) .> 2.0  # Only show cells with >2°C change
    T_range = (
        min(T_range[1], minimum(ΔT[cells_to_show])),
        max(T_range[2], maximum(ΔT[cells_to_show]))
    )
end

for (i, step) in enumerate(key_steps[end-1:end])
    row = (i-1) ÷ 2 + 1
    col = (i-1) % 2 + 1
    
    ax = Axis3(fig[row, col], 
        title = "Year $(round(times[step], digits=1))",
        zreversed = true, 
        azimuth = -1.1*pi/2,
        aspect = :data)
    
    ΔT = Δstates[step][:Temperature]
    
    # Create a mask to show only significant temperature changes
    cells_to_show = cells_back .&& abs.(ΔT) .> 2.0  # Only show cells with >2°C change
    
    plot_cell_data!(ax, msh, ΔT, 
        cells = cells_to_show,
        colorrange = T_range,
        colormap = :seaborn_icefire_gradient)
    
    # Plot wells
    plot_well!(ax, msh, wells[:Hot], color = :black, linewidth = 4)
    plot_well!(ax, msh, wells[:Cold], color = :black, linewidth = 4)
    
    hidedecorations!(ax)


end

# Add colorbar
Colorbar(fig[2, 1:2], 
    colormap = :seaborn_icefire_gradient, 
    colorrange = T_range,
    label = "Temperature Change (°C)", 
    vertical = false)

fig