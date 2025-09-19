# # Enhanced Geothermal System (EGS)
# This example demonstrates simulation and analysis of energy production from an
# Enhanced Geothermal System (EGS). EGS technology enables geothermal energy
# extraction from hot dry rock formations where natural permeability is
# insufficient for fluid circulation.
#
# ## EGS System Overview
# EGS systems create artificial geothermal reservoirs through hydraulic
# stimulation, creating a fracture network that connects injection and
# production wells.

# ## System Configuration
# The simulated EGS consists of:
# - **Injector well**: Injects cold water (25°C) at controlled flow rate
# - **Producer well**: Extracts heated water for surface energy conversion
# - **Fracture network**: Discrete high-permeability zones enabling fluid circulation
# - **Hot rock matrix**: Provides thermal energy source (initial temperature gradient)

# ## Required Packages and Dependencies
#
# The simulation requires several specialized Julia packages:

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul  # Core reservoir simulation framework
using HYPRE # High-performance linear algebra solvers
using GLMakie # 3D visualization and plotting capabilities

# Useful SI units
meter = si_unit(:meter)
hour = si_unit(:hour)
kilogram = si_unit(:kilogram)
Kelvin, watt, joule = si_units(:Kelvin, :watt, :joule)

# ## EGS System Geometry and Configuration
#
# Design the geometric layout of the EGS including well trajectories, fracture
# network dimensions, and operational parameters. The configuration represents a
# small, commercial-scale EGS with horizontal wells and engineered fracture
# network.

# ### Geometric Parameters
#
# Define EGS system geometry
well_spacing_x = 100.0 # Horizontal separation between wells [m]
well_spacing_z = sqrt(3/4*well_spacing_x^2) # Vertical well offset
well_depth = 2500.0 # Vertical depth to horizontal well section [m] - typical EGS depth
well_lateral = 500.0 # Length of horizontal well section [m]
fracture_radius = 200.0 # Radius of stimulated fracture network [m]
fracture_spacing = well_lateral/8 # Spacing between discrete fractures [m] - controls resolution

# ### Well Trajectory Definition
#
# Wells start vertically from surface, then transition to horizontal at depth to
# maximize contact with fracture network. The setup function supports arbitrary
# number of wells, where the first is the injector and subsequent wells are
# producers. All producers are coupled together at the surface.
ws_x, ws_z, wd, wl = well_spacing_x, well_spacing_z, well_depth, well_lateral
well_coords = [
    [0.0 0.0 0.0; 0.0 0.0 wd; 0.0 wl wd], # Injector
    [-ws_x/2 0.0  0.0; -ws_x/2 0.0 wd-ws_z; -ws_x/2 wl wd-ws_z], # Producer leg 1
    [ ws_x/2 0.0  0.0;  ws_x/2 0.0 wd-ws_z;  ws_x/2 wl wd-ws_z] # Producer leg 2
]

# ### EGS Case Generation
#
# Create the complete EGS model including reservoir domain, fracture network,
# wells, and simulation schedule using Fimbul's specialized EGS constructor.

# Create EGS case with fracture network
reports_per_year = 4 # Output frequency for results analysis
case = Fimbul.egs(well_coords, fracture_radius, fracture_spacing;
    fracture_aperture = fracture_aperture,
    rate = 100meter^3/hour, # Water injection rate
    temperature_inj = convert_to_si(25.0, :Celsius), # Injection temperature
    num_years = 10, # Years of operation
    schedule_args = (report_interval = si_unit(:year)/reports_per_year,)
);

# ## 3D Visualization of EGS System
#
# Generate comprehensive 3D visualization of the computational mesh, well trajectories,
# and fracture network. The mesh is adaptively refined around wells and fractures
# to accurately resolve thermal and hydraulic gradients in these critical regions.

# ### Extract Mesh and Geometry Information
#
# Access the computational grid and geometric properties for visualization:

# Inspect EGS model
# Visualize the computational mesh, wells, and fracture network. The mesh is
# refined around wells and fractures to capture thermal and hydraulic processes
# accurately in these critical regions.
msh = physical_representation(reservoir_model(case.model).data_domain)  # Extract mesh geometry
geo = tpfv_geometry(msh)                                                # Compute geometric properties

# ### Create 3D System Overview Plot
#
# Generate main visualization showing mesh structure and well locations:

fig = Figure(size = (1200, 800))
ax = Axis3(fig[1, 1], zreversed = true, aspect = :data,
    title = "EGS System: Wells and Fracture Network")
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.3)  # Show computational grid with transparency

# ### Well Trajectory Visualization
#
# Plot wells with distinct colors to show injection and production roles:

# Plot wells with distinct colors
wells = get_model_wells(case.model)
# colors = cgrad(:seaborn_icefire_gradient, 10, categorical = true)[[1,end]]
function plot_egs_wells(ax; colors = [:red, :blue])
    for (i, (name, well)) in enumerate(wells)
        color = colors[i]
        label = i == 1 ? "Injector" : "Producer"
        cells = well.perforations.reservoir
        xy = geo.cell_centroids[1:2, cells]
        xy = hcat(xy[:,1], xy)'
        plot_mswell_values!(ax, case.model, name, xy;
            geo = geo, linewidth = 3, color = color, label = label)
    end
end

plot_egs_wells(ax)

# ### Fracture Network Visualization
#
# Highlight high-permeability fracture zones in the computational domain:

# Highlight fracture network (high porosity cells)
domain = reservoir_model(case.model).data_domain
is_fracture = domain[:porosity] .> 0.1  # Fractures identified by enhanced porosity
fracture_cells = findall(is_fracture)
if !isempty(fracture_cells)
    plot_cell_data!(ax, msh, Float64.(is_fracture),
        cells = fracture_cells,
        colorrange = (0, 1))
end
fig

# ### Reservoir Property Visualization
#
# Examine the spatial distribution of geological properties, particularly
# porosity which distinguishes fracture zones from rock matrix:

# Plot reservoir properties
# Examine the geological properties.
plot_reservoir(case.model, key = :porosity, aspect = :data, colormap = :hot)

# ## Numerical Simulation Configuration
#
# Configure the reservoir simulator with appropriate solver settings, convergence criteria,
# and timestep control for robust simulation of coupled thermal-hydraulic processes in EGS.

# ### Primary Solver Setup
#
# Initialize the reservoir simulator with tolerance settings optimized for EGS:

# Simulate EGS energy production
# We simulate the EGS system for 30 years of operation. The injector is set to
# inject cold water at 20°C, while the producer extracts heated water. The
# fracture network facilitates heat exchange between the circulating fluid and
# the surrounding hot rock matrix.
sim, cfg = setup_reservoir_simulator(case;
    info_level = 2,          # Detailed convergence information (0=progress bar, 1=basic, 2=detailed)
    tol_cnv = Inf,          # Disable CNV tolerance (material balance) - use other criteria
    inc_tol_dT = 1e-2,      # Temperature increment tolerance [K] - ensures thermal accuracy
    inc_tol_dp_rel = 1e-3,  # Relative pressure increment tolerance [-] - controls pressure changes
    initial_dt = 5.0,       # Initial timestep [days] - conservative start for stability
    relaxation = true);     # Enable solution relaxation for numerical stability

# ### Adaptive Timestep Control
#
# Add specialized timestep selectors to control solution quality during thermal transients.
# These selectors monitor temperature changes and adjust timesteps to maintain accuracy:

sel = VariableChangeTimestepSelector(:Temperature, 5.0; relative = false, model = :Injector)
push!(cfg[:timestep_selectors], sel)
sel = VariableChangeTimestepSelector(:Temperature, 5.0; relative = false, model = :Producer)
push!(cfg[:timestep_selectors], sel)

##
# Note: EGS simulations can be computationally intensive due to the coupled
# thermal-hydraulic processes in fractures. Setting info_level = 0 will show
# a progress bar during simulation.
results = simulate_reservoir(case; simulator = sim, config = cfg)
# ## Results Analysis and Visualization
#
# Analyze simulation results to understand EGS performance, thermal depletion patterns,
# and energy production characteristics throughout the operational period.

# ### Temperature Field Evolution
#
# Visualize how thermal depletion progresses through the fracture network over time.
# This shows the cooling front advancement and thermal interaction between fractures:

# Visualize EGS results

# First, plot the temperature field evolution throughout the simulation.
# The visualization shows how thermal depletion progresses through the
# fracture network over time.
plot_reservoir(case.model, results.states;
    colormap = :seaborn_icefire_gradient,
    key = :Temperature,
    aspect = :data)

# ### Well Performance Analysis
#
# Examine operational metrics including flow rates, pressures, and temperatures
# to understand EGS system behavior and identify thermal breakthrough timing:

# Plot well performance
# Examine the well responses including flow rates, pressures, and temperatures
# to understand EGS operational behavior and thermal decline over time.
plot_well_results(results.wells)

# ### Temperature Change Analysis
#
# Calculate and visualize temperature changes relative to initial conditions
# to quantify thermal depletion magnitude and spatial distribution:

# Calculate temperature changes from initial conditions
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

# Plot temperature changes to highlight thermal depletion zones
plot_reservoir(case.model, Δstates;
    colormap = :seaborn_icefire_gradient,
    key = :Temperature,
    aspect = :data)

# ### Fracture-Level Thermal Evolution
#
# Create detailed time-series visualization of temperature changes within 
# the fracture network to understand thermal depletion progression:

# Extract fracture zones and prepare time-series data
is_frac = isapprox.(domain[:porosity], maximum(domain[:porosity]))  # Identify fracture cells
time = convert_from_si.(cumsum(dt), :year)                          # Convert timesteps to years
ΔT = [Δstate[:Temperature][is_frac] for Δstate in Δstates]         # Temperature changes in fractures
colorrange = extrema(vcat(ΔT...))                                   # Color scale for all timesteps

# Set up 3D visualization parameters for fracture temperature evolution
x = geo.cell_centroids[:, is_frac]    # Fracture cell coordinates
xlim = extrema(x, dims=2)             # Determine spatial bounds
xlim = [xl .+ (xl[2]-xl[1]).*(-0.1, 0.1) for xl in xlim]  # Add 10% margin
limits = Tuple(xlim)

# Create multi-panel time-series visualization
fig = Figure(size = (650, 800))
steps = Int.(round.(range(1, length(ΔT), 3)))  # Select 3 representative timesteps
for (n, ΔT_n) in enumerate(ΔT[steps])
    ax_n = Axis3(fig[n, 1],
    zreversed = true, aspect = (1,6,1),
    azimuth = 1.2pi,
    elevation = pi/20,
    limits = limits,
    title = "$(round(time[steps[n]], digits=1)) years",
    titlegap = -10
    )
    ## Plot temperature changes in fracture network
    plot_cell_data!(ax_n, msh, ΔT_n;
        cells = is_frac,
        colorrange = colorrange,
        colormap = :seaborn_icefire_gradient)
    plot_egs_wells(ax_n; colors = [:black, :black])
    hidedecorations!(ax_n)
end

Colorbar(fig[length(steps)+1, 1];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange,
    label = "ΔT (°C)", vertical = false, flipaxis = false)

fig

# ### Fracture Data Extraction and Analysis
#
# Extract detailed flow and thermal data from individual fractures to analyze
# their relative contributions to overall EGS energy production performance.

## Fracture data extraction
function get_fracture_data(states, model, well)
    """
    Extract flow rates, temperatures, and energy fluxes from individual fractures
    connected to the specified well. This function processes simulation results
    to compute fracture-level performance metrics.
    
    Args:
        states: Simulation state results for all timesteps
        model: EGS model containing well and fracture geometry
        well: Well symbol (:Injector or :Producer) to analyze
    
    Returns:
        Dictionary containing fracture-level data arrays:
        - :y: Y-coordinates of fracture locations
        - :Temperature: Temperature evolution in each fracture
        - :MassFlux: Mass flow rate through each fracture
        - :EnergyFlux: Thermal energy flux through each fracture
    """

    ## Extract well perforation data and connectivity information
    perf = model.models[well].data_domain.representation.perforations  # Well-fracture connections
    N = model.models[well].domain.representation.neighborship           # Cell connectivity matrix
    is_frac = findall(isapprox.(perf.WI, maximum(perf.WI)))            # High-productivity fracture connections
    jj = [cell_ijk(msh, c)[2] for c in perf.reservoir[is_frac]]        # Y-indices of fracture cells

    ## Initialize data arrays for fracture analysis
    JJ = unique(jj)              # Unique fracture Y-locations
    nstep = length(states)       # Number of simulation timesteps
    nfrac = length(JJ)           # Number of discrete fractures
    Q, Qh, T = zeros(nstep, nfrac), zeros(nstep, nfrac), zeros(nstep, nfrac)  # Pre-allocate arrays
    y = geo.cell_centroids[2, perf.reservoir[is_frac]]  # Y-coordinates of fracture centers

    ## Process each simulation timestep
    for (sno, state) in enumerate(states)
        ## Extract mass fluxes and enthalpies from well model
        Qn = state[well][:TotalMassFlux]    ## Mass flux per well segment [kg/s]
        h = state[well][:FluidEnthalpy]     ## Fluid enthalpy [J/kg]
        
        ## Apply upwind scheme for enthalpy (use upstream cell value)
        h = [q > 0 ? h[N[1,s]] : h[N[2,s]] for (s,q) in enumerate(Qn)]
        
        ## Compute energy fluxes combining mass and enthalpy
        Qhn = Qn.*h                         ## Energy flux [W]

        ## Calculate net fluxes per fracture segment (flow into fracture)
        Qn = .-diff(Qn)[is_frac]           ## Net mass flux into fracture [kg/s]
        Qhn = .-diff(Qhn)[is_frac]         ## Net energy flux into fracture [W]
        Tn = state[well][:Temperature][is_frac]  ## Temperature in fracture [K]
        h = state[well][:FluidEnthalpy][is_frac] ## Enthalpy in fracture [J/kg]
        
        ## Aggregate data by fracture Y-location
        for (fno, j) in enumerate(JJ)
            ix = jj .== j                   ## Cells belonging to this fracture
            Q[sno, fno] = sum(Qn[ix])       ## Total mass flow for this fracture
            Qh[sno, fno] = sum(Qhn[ix])     ## Total energy flow for this fracture  
            T[sno, fno] = mean(Tn[ix])      ## Average temperature for this fracture
        end
    end

    ## Return organized fracture data
    data = Dict()
    data[:y] = y                    ## Spatial coordinates
    data[:Temperature] = T          ## Temperature evolution [K]
    data[:MassFlux] = Q            ## Mass flux evolution [kg/s]
    data[:EnergyFlux] = Qh         ## Energy flux evolution [W]

    return data
end

# ### Apply Fracture Data Extraction
#
# Extract detailed performance data from injection and production wells
# to analyze individual fracture contributions and thermal depletion patterns:

##
states = results.result.states;                              # Extract simulation states
fdata_inj = get_fracture_data(states, case.model, :Injector);   # Injector fracture data
fdata_prod = get_fracture_data(states, case.model, :Producer);  # Producer fracture data

# ### Energy Production Analysis Setup
#
# Prepare arrays and visualization settings for comprehensive energy analysis:

nsteps = length(dt)                                         # Number of simulation timesteps
energy, cat, dodge = Float64[], Int[], Int[]               # Arrays for energy analysis plots

# Define color scheme for fracture identification
colors = reverse(cgrad(:seaborn_icefire_gradient, size(fdata_inj[:Temperature], 2), categorical = true))
function plot_fracture_data(ax, data, stacked = false)
    df_prev = zeros(nsteps)
    for (fno, df) in enumerate(eachcol(data))
        if stacked
            ## plot values on top of each other to illustrate realtive contribution
            x = vcat(time, reverse(time))
            df .+= df_prev
            y = vcat(df_prev, reverse(df))
            poly!(ax, x, y; color = colors[fno],
            strokecolor = :black, strokewidth = 1, label = "Fracture $fno")
            df_prev = df
        else
            lines!(ax, time, df; color = colors[fno],
            linewidth = 2, label = "Fracture $fno")
        end
    end
end

# ### Comprehensive Energy Performance Visualization
#
# Create multi-panel plot showing temperature evolution, thermal power production,
# and energy analysis across all fractures to understand EGS system performance:

##

fig = Figure(size = (1000, 800))
xmax = round(maximum(time))
limits = ((0, xmax).+(-0.1, 0.1).*xmax, nothing)
xticks = 0:xmax
function make_axis(title, ylabel, rno; kwargs...)
    ax = Axis(fig[rno,1];
    title = title, xlabel = "Time (years)", ylabel = ylabel, limits = limits,
    xticks = xticks, kwargs...)
    return ax
end

# ### Panel 1: Fracture Temperature Evolution
#
# Show how temperature changes in each fracture over the operational period:

# Plot temperature
ax = make_axis("Temperature", "T (°C)", 1)
temperature = fdata_prod[:Temperature]                         # Producer temperatures [K]
plot_fracture_data(ax, convert_from_si.(temperature, :Celsius), false)  # Convert to Celsius
hidexdecorations!(ax, grid = false)                           # Clean appearance for multi-panel

# ### Panel 2: Thermal Power Production  
#
# Display net thermal power extracted from the EGS system with stacked contributions:

ax_pwr = make_axis("Thermal power", "Power (MW)", 2)
effect = .-(fdata_prod[:EnergyFlux] .+ fdata_inj[:EnergyFlux])  # Net energy extraction [W]
plot_fracture_data(ax_pwr, effect./1e6, true)                # Convert to MW, use stacked plot
hidexdecorations!(ax_pwr, grid = false)

# ### Panel 3: Annual Energy Production per Fracture
#
# Calculate and display cumulative energy production for each fracture annually:

height, cat, dodge = Float64[], Int[], Int[]                  # Arrays for grouped bar chart
n = reports_per_year                                          # Reports per year
energy_per_year = []
for (fno, effect_f) in enumerate(eachcol(effect))
    ## Integrate power over annual periods to get energy [GWh]
    energy_f = [sum(effect_f[k:k+n-1].*dt[k:k+n-1])./(si_unit(:giga)*si_unit(:hour)) for k = 1:n:length(dt)-n+1]
    push!(energy_per_year, energy_f)
    push!(cat, 1:length(energy_f)...)
    push!(dodge, fill(fno, length(energy_f))...)
end
ax = make_axis("Energy per year per fracture", "Energy (GWh)", 3)
barplot!(ax, cat, vcat(energy_per_year...);
dodge = dodge, color = colors[dodge], strokecolor = :black, strokewidth = 1)
hidexdecorations!(ax, grid = false)

# ### Panel 4: Relative Fracture Contribution Analysis
#
# Show the relative contribution of each fracture to total energy production:

η = reduce(hcat, energy_per_year)                             # Combine energy arrays
η = η ./ sum(η, dims=2)                                       # Normalize to fractions
ax = make_axis("Energy fraction per year per fracture", "Fraction (-)", 4)
barplot!(ax, cat, η[:];
dodge = dodge, color = colors[dodge], strokecolor = :black, strokewidth = 1)

# Add legend showing fracture identification
Legend(fig[2:3,2], ax_pwr)
fig