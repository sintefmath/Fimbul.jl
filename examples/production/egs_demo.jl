# # Enhanced Geothermal System (EGS)
# This example demonstrates simulation and analysis of energy production from an
# Enhanced Geothermal System (EGS). EGS technology enables geothermal energy
# extraction from hot dry rock formations where natural permeability is
# insufficient for fluid circulation.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul # Core reservoir simulation framework
using HYPRE # High-performance linear solvers
using GLMakie # 3D visualization and plotting capabilities
using Statistics # For computing statistics

# Useful SI units
meter, day, watt = si_units(:meter, :day, :watt)

# ## EGS setup
# We consider an EGS system with one injection well and two production wells.
# The wells extend 2500 m vertically before they continue 500 m horizontally.
# The horizontal sections are arranged in a triangular pattern, connected by a
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
# together into a single well at the top to enable controlling/monitoring the
# total production rate.
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

# ### Create simulation case
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

# ### Inspect model
# Visualize the computational mesh, wells, and fracture network. The mesh is
# refined around wells and fractures to capture thermal and hydraulic processes
# accurately in these critical regions.
msh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(msh)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, aspect = :data,
    title = "EGS System: Wells and Fracture Network")
# Show computational mesh with transparency
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.3)
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
# We configure the simulator to improve nonlinear convergence. Standard CNV
# criteria are disabled because fractures create numerical challenges: the
# fluid/energy flux through thin, high-permeability fracture cells can be orders
# of magnitude larger than the accumulation within those cells, making material
# balance errors appear artificially large despite physical accuracy. Instead,
# we use increment-based criteria that directly monitor changes in primary
# variables (temperature and pressure) between Newton iterations.
sim, cfg = setup_reservoir_simulator(case;
    info_level = 0, # 0=progress bar, 1=basic, 2=detailed
    tol_cnv = Inf, # Disable CNV tolerance - use other criteria
    inc_tol_dT = 1e-2, # Temperature increment tolerance [K]
    inc_tol_dp_rel = 1e-3, # Relative pressure increment tolerance [-]
    initial_dt = 5.0, # Initial timestep [s]
    relaxation = true); # Enable relaxation in Newton solver

# We add a specialized timestep selector to control solution quality during
# thermal transients. These selectors monitor temperature changes and adjust
# timesteps aiming at a maximum change of 5°C per timestep.
sel = VariableChangeTimestepSelector(:Temperature, 5.0; 
relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel)

# Note: EGS simulations can be computationally intensive. Depending on your
# system, the simulation may take several minutes to complete.
results = simulate_reservoir(case; simulator = sim, config = cfg)

# ## Interactive Visualization
# Next, we analyze and visualize the simulation results interactively to
# understand the EGS performance, thermal depletion patterns, and energy
# production characteristics throughout the operational period.

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

# ### Well Performance Analysis
# Next, we examine the well responses including flow rates, pressures, and
# temperatures. The injector is configured to mirror the production rate to
# avoid pressure buildup. Notice how the flow rate declines over time due to
# cooling of the fracture network and consequent increase in fluid viscosity.
plot_well_results(results.wells)

# ## Fracture-level performance
# Understanding the performance of each individual fracture is crucial for
# optimizing EGS systems. We first visualize the temperature within the fracture
# network at selected timesteps to understand thermal depletion patterns.

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
        cells = is_fracture,
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

# ### Fracture metrics
# Finally, we analyze key performance metrics for each individual fracture over
# the entire simulation period, including temperature evolution, thermal power
# production, annual energy production. Notice how the first three fractures
# contribute more than 50 % of the total energy production throughout the
# simulation period. The last three fractures contribute only 21 % in the first
# year, but their relative contribution increases to 27 % in the final year.

# Extract fracture data from simulation results
states = results.result.states;
fdata_inj = Fimbul.get_egs_fracture_data(states, case.model, :Injector; geo=geo);
fdata_prod = Fimbul.get_egs_fracture_data(states, case.model, :Producer; geo=geo);
nsteps = length(dt)
energy, cat, dodge = Float64[], Int[], Int[]

# Utility for plotting fracture data
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
# Utility for making axes
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

# Panel 1: fracture temperature
ax = make_axis("Temperature", "T (°C)", 1)
temperature = fdata_prod[:Temperature]
plot_fracture_data(ax, convert_from_si.(temperature, :Celsius), false)
hidexdecorations!(ax, grid = false)

# Panel 2: thermal power production
ax_pwr = make_axis("Thermal power", "Power (MW)", 2)
effect = .-(fdata_prod[:EnergyFlux] .+ fdata_inj[:EnergyFlux])
plot_fracture_data(ax_pwr, effect./1e6, true)
hidexdecorations!(ax_pwr, grid = false)

# Panel 3: annual energy production per fracture
energy_per_year, cat, dodge = [], Int[], Int[]
n = reports_per_year
GWh = si_unit(:giga)*si_unit(:watt)*si_unit(:hour)
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

# Panel 4: relative fracture contribution
η = reduce(hcat, energy_per_year)
η = η ./ sum(η, dims=2)
ax = make_axis("Energy fraction per year per fracture", "Fraction (-)", 4)
barplot!(ax, cat, η[:];
dodge = dodge, color = colors[dodge], strokecolor = :black, strokewidth = 1)

# Add legend
Legend(fig[2:3,2], ax_pwr)
fig