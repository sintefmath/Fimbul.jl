# # Enhanced Geothermal System (EGS) Energy Production
# This example demonstrates how to simulate and visualize energy production from
# an Enhanced Geothermal System (EGS). EGS systems create artificial geothermal
# reservoirs by hydraulically fracturing hot dry rock and circulating water
# through the created fracture network to extract thermal energy.
#
# The setup consists of horizontal wells connected by a fracture network:
# - An injector well injects cold water into the fracture system
# - A producer well extracts heated water after thermal exchange with hot rock
# - Fractures provide high-permeability pathways for fluid circulation

# Add required modules to namespace 
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# ## Set up EGS case
# We create an EGS system with two horizontal wells connected by fractures.
# The wells are positioned 200 m apart horizontally and extend laterally
# through the fractured rock formation at depth.

# Define EGS system geometry
well_spacing = 200.0  # meters between wells
fracture_radius = 300.0  # radius of fracture network
fracture_spacing = 200.0  # spacing between fractures
well_depth = 2500.0  # vertical depth to start of horizontal section
well_lateral = 2000.0  # horizontal length of wells
fracture_aperture = 1e-3  # fracture aperture in meters

# Well coordinate matrices: [x y z] for each well segment
# Each row represents a point along the well trajectory
ws, wd, wl = well_spacing, well_depth, well_lateral
well_coords = [
    [-ws/2 0.0  0.0; -ws/2 0.0 wd; -ws/2 wl wd], 
    [ ws/2 0.0  0.0;  ws/2 0.0 wd;  ws/2 wl wd]
]

# Create EGS case with fracture network
case = Fimbul.egs(
    well_coords, 
    fracture_radius, 
    fracture_spacing; # Domain size parameter
    num_years = 30, # Simulate 30 years of operation
    fracture_aperture = fracture_aperture,
    schedule_args = (report_interval = 7.0si_unit(:day),)
);

# ## Inspect EGS model
# Visualize the computational mesh, wells, and fracture network. The mesh is
# refined around wells and fractures to capture thermal and hydraulic processes
# accurately in these critical regions.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true, aspect = :data,
    title = "EGS System: Wells and Fracture Network")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.3)

# Plot wells with distinct colors
wells = get_model_wells(case.model)
colors = [:red, :blue]  # Injector = red, Producer = blue
for (i, (name, well)) in enumerate(wells)
    color = i == 1 ? :red : :blue
    label = i == 1 ? "Injector" : "Producer"
    plot_well!(ax, msh, well, color = color, linewidth = 6)
end

# Highlight fracture network (high porosity cells)
domain = reservoir_model(case.model).data_domain
is_fracture = domain[:porosity] .> 0.1  # Fractures have high porosity
fracture_cells = findall(is_fracture)
if !isempty(fracture_cells)
    plot_cell_data!(ax, msh, Float64.(is_fracture), 
        cells = fracture_cells, 
        colormap = :heat, 
        alpha = 0.6,
        colorrange = (0, 1))
end
fig

# ### Plot reservoir properties
# Examine the geological properties.
plot_reservoir(case.model, key = :porosity, aspect = :data, colormap = :hot)

# ## Simulate EGS energy production
# We simulate the EGS system for 30 years of operation. The injector is set to
# inject cold water at 20°C, while the producer extracts heated water. The 
# fracture network facilitates heat exchange between the circulating fluid and
# the surrounding hot rock matrix.

# Note: EGS simulations can be computationally intensive due to the coupled
# thermal-hydraulic processes in fractures. Setting info_level = 0 will show
# a progress bar during simulation.
results = simulate_reservoir(case; info_level = 2)

# ## Visualize EGS results
# First, plot the temperature field evolution throughout the simulation.
# The visualization shows how thermal depletion progresses through the
# fracture network over time.
plot_reservoir(case.model, results.states;
    colormap = :seaborn_icefire_gradient, 
    key = :Temperature, 
    aspect = :data)

# ### Plot well performance
# Examine the well responses including flow rates, pressures, and temperatures
# to understand EGS operational behavior and thermal decline over time.
plot_well_results(results.wells)

# ### Production temperature evolution
# For detailed analysis of thermal decline, we plot the temperature profile
# along the production well for the entire simulation period. This shows how
# thermal depletion affects different parts of the wellbore over time.
states = results.result.states
times = convert_from_si.(cumsum(case.dt), :year)
geo = tpfv_geometry(msh)

fig = Figure(size = (1200, 800))
ax = Axis(fig[1, 1], 
    xlabel = "Temperature (°C)", 
    ylabel = "Depth (m)", 
    yreversed = true,
    title = "Producer Well Temperature Evolution")

colors = cgrad(:seaborn_icefire_gradient, length(states), categorical = true)

# Plot temperature profiles for all timesteps with transparency
for (n, state) in enumerate(states)
    # Get producer well temperatures
    well_name = :Producer
    if haskey(state, well_name)
        T = convert_from_si.(state[well_name][:Temperature], :Celsius)
        plot_mswell_values!(ax, case.model, well_name, T;
            geo = geo, linewidth = 3, color = colors[n], alpha = 0.2)
    end
end

# Highlight selected timesteps with solid lines and labels
highlight_years = [1, 5, 10, 20, 30]
highlight_indices = [findfirst(t -> t >= year, times) for year in highlight_years]
highlight_indices = filter(!isnothing, highlight_indices)

for idx in highlight_indices
    if haskey(states[idx], :Producer)
        T = convert_from_si.(states[idx][:Producer][:Temperature], :Celsius)
        plot_mswell_values!(ax, case.model, :Producer, T;
            geo = geo, linewidth = 4, color = colors[idx], 
            label = "$(round(times[idx], digits=1)) years")
    end
end

axislegend(ax; position = :rt, fontsize = 14)
fig

# The temperature profiles show the progression of thermal depletion as cold
# water circulation gradually cools the rock matrix around the fractures.

# ### Temperature change analysis
# Compute and visualize the change in reservoir temperature from initial
# conditions to understand thermal depletion patterns in the EGS system.
Δstates = []
for state in results.states
    Δstate = Dict{Symbol, Any}()
    for (k, v) in state
        if haskey(case.state0[:Reservoir], k)
            v0 = case.state0[:Reservoir][k]
            Δstate[k] = v .- v0
        end
    end
    push!(Δstates, Δstate)
end

# Interactive visualization of temperature changes
plot_reservoir(case.model, Δstates;
    colormap = :seaborn_icefire_gradient, 
    key = :Temperature, 
    aspect = :data)

# ### 3D visualization of thermal depletion
# Create 3D plots showing temperature changes at key timesteps to visualize
# the spatial evolution of thermal depletion in the fracture network.
if !isempty(Δstates)
    T_min = minimum([minimum(get(s, :Temperature, [0.0])) for s in Δstates])
    T_max = maximum([maximum(get(s, :Temperature, [0.0])) for s in Δstates])
    
    # Create a mask to show interesting parts of the domain
    cells_mask = geo.cell_centroids[1, :] .> -500.0
    cells_mask = cells_mask .&& geo.cell_centroids[2, :] .> -200.0
    cells_mask = cells_mask .&& geo.cell_centroids[3, :] .> well_depth - 200.0
    
    fig = Figure(size = (1000, 800))
    
    # Plot temperature changes at selected timesteps
    plot_timesteps = min(4, length(highlight_indices))
    for (i, idx) in enumerate(highlight_indices[1:plot_timesteps])
        row = (i-1) ÷ 2 + 1
        col = (i-1) % 2 + 1
        
        ax = Axis3(fig[row, col]; 
            title = "$(round(times[idx], digits=1)) years",
            zreversed = true, 
            elevation = π/8, 
            aspect = :data)
        
        if haskey(Δstates[idx], :Temperature)
            ΔT = Δstates[idx][:Temperature]
            
            # Plot temperature changes
            plot_cell_data!(ax, msh, ΔT; 
                cells = cells_mask, 
                colormap = :seaborn_icefire_gradient, 
                colorrange = (T_min, T_max),
                alpha = 0.8)
            
            # Plot wells
            for (name, well) in wells
                color = name == :Injector ? :black : :gray
                plot_well!(ax, msh, well; 
                    color = color, 
                    linewidth = 4, 
                    markersize = 0.0)
            end
        end
        
        hidedecorations!(ax)
    end
    
    # Add colorbar
    if plot_timesteps > 2
        Colorbar(fig[3, 1:2]; 
            colormap = :seaborn_icefire_gradient, 
            colorrange = (T_min, T_max), 
            label = "Temperature Change (°C)", 
            vertical = false)
    else
        Colorbar(fig[:, 3]; 
            colormap = :seaborn_icefire_gradient, 
            colorrange = (T_min, T_max), 
            label = "Temperature Change (°C)")
    end
    
    fig
end

# ### EGS performance metrics
# Calculate and display key performance indicators for the EGS system
# including thermal power output and energy recovery efficiency.
if haskey(results.wells, :Producer)
    producer_data = results.wells[:Producer]
    
    # Calculate thermal power output over time
    mass_rate = producer_data[:mass_rate]  # kg/s
    temperature = producer_data[:temperature]  # K
    T_injection = convert_to_si(20.0, :Celsius)  # Injection temperature
    cp_water = 4200.0  # J/(kg·K) - specific heat capacity of water
    
    thermal_power = mass_rate .* cp_water .* (temperature .- T_injection) ./ 1e6  # MW
    
    fig = Figure(size = (1000, 600))
    
    # Plot thermal power evolution
    ax1 = Axis(fig[1, 1], 
        xlabel = "Time (years)", 
        ylabel = "Thermal Power (MW)",
        title = "EGS System Performance")
    lines!(ax1, times, thermal_power, color = :red, linewidth = 3)
    
    # Plot production temperature
    ax2 = Axis(fig[2, 1], 
        xlabel = "Time (years)", 
        ylabel = "Production Temperature (°C)")
    T_prod_celsius = convert_from_si.(temperature, :Celsius)
    lines!(ax2, times, T_prod_celsius, color = :blue, linewidth = 3)
    
    fig
    
    println("\nEGS Performance Summary:")
    println("Initial thermal power: $(round(thermal_power[1], digits=2)) MW")
    println("Final thermal power: $(round(thermal_power[end], digits=2)) MW")
    println("Initial production temperature: $(round(T_prod_celsius[1], digits=1)) °C")
    println("Final production temperature: $(round(T_prod_celsius[end], digits=1)) °C")
    println("Thermal decline over $(times[end]) years: $(round((1 - thermal_power[end]/thermal_power[1])*100, digits=1))%")
end
