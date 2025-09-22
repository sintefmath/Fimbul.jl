# # Enhanced Geothermal System (EGS)
# This example demonstrates simulation and analysis of energy production from an
# Enhanced Geothermal System (EGS). EGS technology enables geothermal energy
# extraction from hot dry rock formations where natural permeability is
# insufficient for fluid circulation.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul # Core reservoir simulation framework
using HYPRE # High-performance linear algebra solvers
using GLMakie # 3D visualization and plotting capabilities
using Statistics # For computing statistics

# Useful SI units
meter = si_unit(:meter)
day = si_unit(:day)
watt = si_unit(:watt);

# ## EGS setup
# We consider an EGS system with one injection well and two production wells.
# The wells extend 2500 m vertically before they continue horizontally. The
# horizontal sections are arranged in a triangular pattern, connected by a
# stimulated fracture network comprising eight fractures intersecting the wells
# at a right angle. Thermal energy is produced by circulating cold water through
# this fracture network, which extracts heat from the surrounding hot rock
# matrix by conduction. To leverage buoyancy effects, the injection well is
# placed at a lower elevation than the production wells, forcing the colder (and
# therefore denser) water to sweep a larger volume of the fracture network,
# thereby enhancing heat extraction.

# ### Define EGS geometry
# We will use the setup function `egs`, which takes as input the well
# coordinates, fracture radius, and fracture spacing to create the complete EGS.
well_depth = 2500.0meter # Vertical depth to horizontal well section [m]
well_spacing_x = 100.0meter # Horizontal separation between production wells [m]
well_spacing_z = sqrt(3/4*well_spacing_x^2) # Vertical offset between injector and producers [m]
well_lateral = 500.0meter # Length of horizontal well section [m]

# Well coordinates are defined as a vector of nx3 arrays, each containing the
# (x,y,z) coordinates of a well trajectory. The first well is the injector,
# while preceding wells are interpreted as producers. All producers are coupled
# together in a single well at the top to allow enable controlling/monitoring
# the total production rate.
ws_x, ws_z, wd, wl = well_spacing_x, well_spacing_z, well_depth, well_lateral
well_coords = [
    [0.0 0.0 0.0; 0.0 0.0 wd; 0.0 wl wd], # Injector
    [-ws_x/2 0.0  0.0; -ws_x/2 0.0 wd-ws_z; -ws_x/2 wl wd-ws_z], # Producer leg 1
    [ ws_x/2 0.0  0.0;  ws_x/2 0.0 wd-ws_z;  ws_x/2 wl wd-ws_z] # Producer leg 2
]
# The fracture network is defined by the fracture radius and spacing along the
# wellbore, starting at the horizontal section of the well.
fracture_radius = 200.0meter # Radius of stimulated fracture network [m]
fracture_spacing = well_lateral/8 # Spacing between discrete fractures [m]

# ### Create EGS VariableChangeTimestepSelector
# We set up a scenario describing 10 years of operation with a water injection
# rate of 9250 m³/day (approximately 107 liters/second) at a temperature of
# 25°C. The simulation will output results four times per year for analysis.
reports_per_year = 4 # Output frequency for results analysis
case = Fimbul.egs(well_coords, fracture_radius, fracture_spacing;
    rate = 9250meter^3/day, # Water injection rate
    temperature_inj = convert_to_si(25.0, :Celsius), # Injection temperature
    num_years = 10, # Years of operation
    schedule_args = (report_interval = si_unit(:year)/reports_per_year,)
);

# ### Inspect EGS model
# Visualize the computational mesh, wells, and fracture network. The mesh is
# refined around wells and fractures to capture thermal and hydraulic processes
# accurately in these critical regions.
msh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(msh)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, aspect = :data,
    title = "EGS System: Wells and Fracture Network")
# Show computational mesh with transparency
Jutul.plot_mesh_edges!(ax, msh)
# Plot wells
wells = get_model_wells(case.model)
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
# Highlight fracture network (high porosity cells)
domain = reservoir_model(case.model).data_domain
is_fracture = isapprox.(domain[:porosity], maximum(domain[:porosity]))
fracture_cells = findall(is_fracture)
plot_mesh!(ax, msh; cells = fracture_cells)
fig

# ## Simulate system
# To simulate the system, we configure the reservoir simulator with settings
# optimized for EGS applications.
sim, cfg = setup_reservoir_simulator(case;
    info_level = 2, # 0=progress bar, 1=basic, 2=detailed
    tol_cnv = Inf, # Disable CNV tolerance - use other criteria
    inc_tol_dT = 1e-2, # Temperature increment tolerance [K]
    inc_tol_dp_rel = 1e-3, # Relative pressure increment tolerance [-]
    initial_dt = 5.0, # Initial timestep [s]
    relaxation = true); # Enable relaxation in Newton solver

# We add a specialized timestep selector to control solution quality during
# thermal transients. These selectors monitor temperature changes and adjust
# timesteps aiming at a maximum change of 5 K per timestep.
sel = VariableChangeTimestepSelector(:Temperature, 5.0; 
relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel)
# sel = VariableChangeTimestepSelector(:Temperature, 5.0; relative = false, model = :Producer)
# push!(cfg[:timestep_selectors], sel)

# Note: EGS simulations can be computationally intensive. Depending on your
# system, the simulation may take several minutes to complete.
results = simulate_reservoir(case; simulator = sim, config = cfg)

# ## Analysis and Visualization
# Next, we analyze and visualize the simulation results to understand the EGS
# performance, thermal depletion patterns, and energy production characteristics
# throughout the operational period.


# ### Reservoir state evolution
# First, plot the reservoir state throughout the simulation.
plot_reservoir(case.model, results.states;
    colormap = :seaborn_icefire_gradient,
    key = :Temperature,
    aspect = :data)

# ### Deviation from initial conditions
# It is often more informative to visualize the deviation from the inital
# conditions to highlight thermal depletion zones. We therefore compute the
# change in reservoir variables to the initial state for all timesteps.
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
plot_reservoir(case.model, Δstates;
    colormap = :seaborn_icefire_gradient,
    key = :Temperature,
    aspect = :data)

# ### Fracture-level thermal evolution
# Before we inspect well output, we visualize the temperature within the
# fracture network. This helps us understand how effectively the cold water
# propagates through the fractures.

# Extract fracture zones and prepare time-series data
dt = case.dt
time = convert_from_si.(cumsum(dt), :year)
ΔT = [Δstate[:Temperature][is_fracture] for Δstate in Δstates]
colorrange = extrema(vcat(ΔT...))

# Define plot limits
x = geo.cell_centroids[:, is_fracture]
xlim = extrema(x, dims=2)
limits = Tuple([xl .+ (xl[2]-xl[1]).*(-0.1, 0.1) for xl in xlim])

# Visualize fracture temperature at three representative timesteps
fig = Figure(size = (650, 800))
steps = Int.(round.(range(1, length(ΔT), 3)))
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
# Add colorbar
Colorbar(fig[length(steps)+1, 1];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange,
    label = "ΔT (°C)", vertical = false, flipaxis = false)
fig

# ### Well Performance Analysis
# Next, we examine the well responses including flow rates, pressures, and
# temperatures. Notice how the flow rate declines over time due to cooling of
# the fracture network and consequent increase in fluid viscosity.
plot_well_results(results.wells)
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
GWh = si_unit(:giga)*si_unit(:watt)*si_unit(:hour)  # GWh unit conversion
for (fno, eff_f) in enumerate(eachcol(effect))
    ## Integrate power over annual periods to get energy [GWh]
    energy_f = [sum(eff_f[k:k+n-1].*dt[k:k+n-1])./GWh for k = 1:n:length(dt)-n+1]
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