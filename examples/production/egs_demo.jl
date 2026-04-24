# # Enhanced Geothermal System (EGS)
# This example demonstrates simulation and analysis of energy production from an
# Enhanced Geothermal System (EGS). EGS technology enables geothermal energy
# extraction from hot dry rock formations where natural permeability is
# insufficient for fluid circulation.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul # Core reservoir simulation framework
using HYPRE # High-performance linear solvers
using GLMakie # 3D visualization and plotting capabilities

# Useful SI units
meter, day, watt = si_units(:meter, :day, :watt);

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
# We will use `egs_well_coordinates` to generate smooth deviated well
# trajectories with a curved bend from vertical to horizontal, and `egs` to
# assemble the complete simulation case using a Discrete Fracture Model (DFM).
fracture_radius  = 200.0meter  # Radius of stimulated fracture disks [m]
fracture_spacing = 125.0meter  # Spacing between discrete fractures [m]

# `egs_well_coordinates` returns `(injector_coords, producer_coords)` as vectors
# of n×3 trajectory matrices. The injector runs vertically to 2500 m depth then
# horizontally; the two producer legs are offset 100 m in x and raised slightly
# to leverage buoyancy-driven flow.
inj, prod = Fimbul.egs_well_coordinates(
    well_depth      = 2500.0meter,
    well_spacing_x  = 100.0meter,
    well_lateral    = 1000.0meter,
    bend_radius     = 200.0meter,
);

# ### Create simulation case
# We set up a scenario describing 10 years of operation with a water injection
# rate of 9250 m³/day (approximately 107 liters/second) at a temperature of
# 25°C. The simulation will output results four times per year for analysis.
num_years = 10 # Total simulation period [years]
case = Fimbul.egs(inj, prod, fracture_radius, fracture_spacing;
    rate            = 9250meter^3/day,                   # Water injection rate
    temperature_inj = convert_to_si(25.0, :Celsius),     # Injection temperature
    num_years       = num_years,
    schedule_args   = (report_interval = si_unit(:year)/4,)
);

# ### Inspect model
# Visualize the computational mesh, DFM fracture network, and wells.
msh = physical_representation(reservoir_model(case.model).data_domain)
geo = tpfv_geometry(msh)

# Get DFM fracture geometry
frac_domain = case.model.models[:Fractures].data_domain
frac_mesh   = physical_representation(frac_domain)
frac_geo    = tpfv_geometry(frac_mesh)

fig = Figure(size = (800, 800))
# Get aspect ration from mesh size
x_range = extrema(geo.cell_centroids[1, :])
y_range = extrema(geo.cell_centroids[2, :])
xy_aspect = (x_range[2]-x_range[1]) / (y_range[2]-y_range[1])

ax = Axis3(fig[1, 1]; zreversed = true, aspect = (xy_aspect,1,3), perspectiveness = 0.5,
    title = "EGS System: Wells and DFM Fracture Network")
Jutul.plot_mesh_edges!(ax, msh; alpha = 0.2)  # Background matrix mesh
Jutul.plot_mesh!(ax, frac_mesh; color = :lightgray)  # DFM fractures

# Plot wells
wells_dict = get_model_wells(case.model; data_domain=true)
function plot_egs_wells(ax; colors = [:red, :blue])
    for (i, (name, well)) in enumerate(wells_dict)
        xy = permutedims(well[:cell_centroids][1:2, :])
        plot_mswell_values!(ax, case.model, name, xy;
            geo = geo, linewidth = 3, color = colors[i])
    end
end
plot_egs_wells(ax)
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
    info_level = 2, # 0=progress bar, 1=basic, 2=detailed
    output_substates = true, # Output results at each Newton iteration
    initial_dt = 5.0, # Initial timestep [s]
    relaxation = true # Enable relaxation in Newton solver
);

# We add a specialized timestep selector to control solution quality during
# thermal transients. These selectors monitor temperature changes and adjust
# timesteps aiming at a maximum change of 5°C per timestep.
sel = VariableChangeTimestepSelector(:Temperature, 5.0; 
relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel);

# Note: EGS simulations can be computationally intensive. Depending on your
# system, the simulation may take several minutes to complete.
results = simulate_reservoir(case; simulator = sim, config = cfg)

# ## Interactive Visualization
# Next, we analyze and visualize the simulation results interactively to
# understand the EGS performance, thermal depletion patterns, and energy
# production characteristics throughout the operational period.

# ### Reservoir state evolution
# First, plot the reservoir state throughout the simulation.
plot_res_args = (
    resolution = (600, 800), aspect = :data, 
    colormap = :seaborn_icefire_gradient, key = :Temperature,
    well_arg = (markersize = 0.0, ),
)
plot_reservoir(case.model, results.states; plot_res_args...)

# ### Deviation from initial conditions
# It is often more informative to visualize the deviation from the inital
# conditions to highlight thermal depletion zones. We therefore compute the
# change in reservoir variables to the initial state for all timesteps.
Δstates = JutulDarcy.delta_state(results.states, case.state0[:Reservoir])
plot_reservoir(case.model, Δstates; plot_res_args...)

# ### Well Performance Analysis
# Next, we examine the well responses including flow rates, pressures, and
# temperatures. The injector is configured to mirror the production rate to
# avoid pressure buildup. Notice how the flow rate declines over time due to
# cooling of the fracture network and consequent increase in fluid viscosity.
plot_well_results(results.wells)

# ## Fracture-level performance
# Understanding the performance of each individual fracture is crucial for
# optimizing EGS systems. We first visualize the DFM fracture temperature at
# selected timesteps.

# Expand to all internal sub-steps to obtain the full time axis with dt
states_full, dt_full, _ = Jutul.expand_to_ministeps(results.result)
time_full = cumsum(dt_full) ./ si_unit(:year)

# Compute fracture temperature change (ΔT) vs. initial state
T0_frac = case.state0[:Fractures][:Temperature]
ΔT_frac = [state[:Fractures][:Temperature] .- T0_frac for state in states_full]
colorrange = extrema(vcat(ΔT_frac...))
fxf, fyf, fzf = frac_geo.cell_centroids[1,:], frac_geo.cell_centroids[2,:], frac_geo.cell_centroids[3,:]
xlim_f = [(extrema(fxf) .+ diff(collect(extrema(fxf))).*[-0.3, 0.3])...]
ylim_f = [(extrema(fyf) .+ diff(collect(extrema(fyf))).*[-0.1, 0.1])...]
zlim_f = [(extrema(fzf) .+ diff(collect(extrema(fzf))).*[-0.3, 0.3])...]
limits_f = (xlim_f, ylim_f, zlim_f)

n_steps_f = length(ΔT_frac)
steps = Int.(round.([0.125, 0.25, 1.0] .* n_steps_f))

fig = Figure(size = (650, 800))
for (n, ΔT_n) in enumerate(ΔT_frac[steps])
    ax_n = Axis3(fig[n, 1];
        perspectiveness = 0.5, zreversed = true, aspect = (1, 6, 1),
        azimuth = 1.2π, elevation = π/20, limits = limits_f,
        title = "$(round(time_full[steps[n]], digits=1)) years", titlegap = -10)
    scatter!(ax_n, fxf, fyf, fzf; color = ΔT_n, colorrange = colorrange,
        colormap = :seaborn_icefire_gradient, markersize = 6)
    plot_egs_wells(ax_n; colors = [:black, :black])
    hidedecorations!(ax_n)
end
Colorbar(fig[length(steps)+1, 1];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange,
    label = "ΔT (°C)", vertical = false, flipaxis = false)
fig

# ### Fracture metrics
# Finally, we analyze key performance metrics for each individual fracture over
# the entire simulation period: temperature evolution, thermal power production,
# and annual energy production.

# Extract fracture data from simulation results using the DFM-aware helper.
# `get_egs_fracture_data` groups fracture cells by y-position and returns mean
# temperature and total thermal energy per fracture for every timestep.
states, dt, _ = Jutul.expand_to_ministeps(results.result)
time = cumsum(dt) ./ si_unit(:year)

fdata = Fimbul.get_egs_fracture_data(states, case.model)
n_frac = length(fdata[:y])
colors = cgrad(:BrBg, n_frac, categorical = true)

function plot_fracture_data(ax, time, data; stacked = false)
    df_prev = zeros(length(time))
    for (fno, df) in enumerate(eachcol(data))
        df = copy(df)
        if stacked
            df .+= df_prev
            poly!(ax, vcat(time, reverse(time)), vcat(df_prev, reverse(df));
                color = colors[fno], strokecolor = :black, strokewidth = 1,
                label = "Fracture $fno")
            df_prev = df
        else
            lines!(ax, time, df; color = colors[fno], linewidth = 2,
                label = "Fracture $fno")
        end
    end
end

fig = Figure(size = (1000, 800))
xmax   = round(maximum(time))
limits = ((0, xmax) .+ (-0.1, 0.1) .* xmax, nothing)
xticks = 0:xmax
function make_axis(title, ylabel, rno; kwargs...)
    Axis(fig[rno, 1]; title = title, xlabel = "Time (years)", ylabel = ylabel,
        limits = limits, xticks = xticks, kwargs...)
end

first_step = findfirst(time .> 1/104)

ax = make_axis("Temperature", "T (°C)", 1)
temperature = fdata[:Temperature][first_step:end, :]
plot_fracture_data(ax, time[first_step:end],
    convert_from_si.(temperature, :Celsius))
hidexdecorations!(ax, grid = false)

# Thermal power approximated from rate of change of stored thermal energy
Er     = fdata[:TotalThermalEnergy]
power  = .-diff(vcat(Er[1,:]', Er), dims = 1) ./ dt  # negative: energy leaving

ax_pwr = make_axis("Thermal power", "Power (MW)", 2)
plot_fracture_data(ax_pwr, time[first_step:end],
    power[first_step:end, :] ./ 1e6; stacked = true)
hidexdecorations!(ax_pwr, grid = false)

ax = make_axis("Annual energy production per fracture", "Energy (GWh)", 3)
energy_per_year, cat, dodge = [], Int[], Int[]
GWh = si_unit(:giga) * si_unit(:watt) * si_unit(:hour)
ix = vcat(0, [findfirst(isapprox.(time, y; atol = 1e-2)) for y in 1:num_years]) .+ 1
for (fno, pwr_f) in enumerate(eachcol(power))
    energy_f = [sum(pwr_f[ix[k]:ix[k+1]-1] .* dt[ix[k]:ix[k+1]-1]) / GWh
                for k in 1:length(ix)-1]
    push!(energy_per_year, energy_f)
    push!(cat, 1:length(energy_f)...)
    push!(dodge, fill(fno, length(energy_f))...)
end
barplot!(ax, cat, vcat(energy_per_year...);
    dodge = dodge, color = colors[dodge], strokecolor = :black, strokewidth = 1)
hidexdecorations!(ax, grid = false)

ax = make_axis("Annual energy fraction per fracture", "Fraction (-)", 4)
η = reduce(hcat, energy_per_year)
η = η ./ sum(η, dims = 2)
barplot!(ax, cat, η[:];
    dodge = dodge, color = colors[dodge], strokecolor = :black, strokewidth = 1)

Legend(fig[2:3, 2], ax_pwr)
fig
