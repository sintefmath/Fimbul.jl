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

# Shared helper: visualize mesh edges, DFM fractures and wells for any EGS case.
function plot_egs_mesh(c, title)
    m_  = physical_representation(reservoir_model(c.model).data_domain)
    g_  = tpfv_geometry(m_)
    fm_ = physical_representation(c.model.models[:Fractures].data_domain)
    xr  = diff(vcat(extrema(g_.cell_centroids[1,:])...))[1]
    yr  = diff(vcat(extrema(g_.cell_centroids[2,:])...))[1]
    zr  = diff(vcat(extrema(g_.cell_centroids[3,:])...))[1]
    asp = (xr, yr, zr) ./ max.(xr, yr, zr)
    fig_ = Figure(size = (800, 800))
    ax_  = Axis3(fig_[1, 1]; zreversed = true, aspect = asp,
        perspectiveness = 0.0, title = title)
    Jutul.plot_mesh_edges!(ax_, m_; alpha = 0.2)
    Jutul.plot_mesh!(ax_, fm_; color = :lightgray)
    for (i, (name, well)) in enumerate(get_model_wells(c.model; data_domain = true))
        xy = permutedims(well[:cell_centroids][1:2, :])
        plot_mswell_values!(ax_, c.model, name, xy;
            geo = g_, linewidth = 3, color = [:red, :blue][i])
    end
    return fig_
end

# Plot wells closure (used in fracture ΔT visualizations below)
wells_dict = get_model_wells(case.model; data_domain=true)
function plot_egs_wells(ax; colors = [:red, :blue])
    for (i, (name, well)) in enumerate(wells_dict)
        xy = permutedims(well[:cell_centroids][1:2, :])
        plot_mswell_values!(ax, case.model, name, xy;
            geo = geo, linewidth = 3, color = colors[i])
    end
end
plot_egs_mesh(case, "Run 1 – Uniform fractures, well spacing 100 m")

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
cfg[:tolerances][:Fractures][:energy_conservation] = (CNV = Inf, EB = 1e-5, increment_dT = 1e-2)
cfg[:tolerances][:Fractures][:mass_conservation] = (CNV = Inf, MB = 1e-5, increment_dp_abs = 1e-2*si_unit(:bar))

# We add a specialized timestep selector to control solution quality during
# thermal transients. These selectors monitor temperature changes and adjust
# timesteps aiming at a maximum change of 5°C per timestep.
sel = VariableChangeTimestepSelector(:Temperature, 10.0; 
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
    resolution = (600, 800), aspect = aspect, 
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
# Well results for all runs are shown together after the third simulation.

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
steps = Int.(round.([0.125, 0.5, 1.0] .* n_steps_f))

fig = Figure(size = (650, 800))
for (n, ΔT_n) in enumerate(ΔT_frac[steps])
    ax_n = Axis3(fig[n, 1];
        perspectiveness = 0.5, zreversed = true, aspect = (1, 6, 1),
        azimuth = 1.2π, elevation = π/20, limits = limits_f,
        title = "$(round(time_full[steps[n]], digits=1)) years", titlegap = -10)
    plot_cell_data!(ax_n, frac_mesh, ΔT_n; colorrange = colorrange,
        colormap = :seaborn_icefire_gradient)
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
# Passing the full `case` lets `get_egs_fracture_data` use the `:cut_no` field
# from `cut_mesh` to directly assign each fracture cell to its originating
# fracture plane, avoiding the heuristic y-coordinate grouping.
states, dt, _ = Jutul.expand_to_ministeps(results.result)
time = cumsum(dt) ./ si_unit(:year)

fdata = Fimbul.get_egs_fracture_data(states, case)
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

# Thermal power from the fracture energy balance:
#   power = dE/dt + q_out - q_in = heat extracted from rock matrix
# where q_in is the energy carried into the fracture by the injected cold fluid
# and q_out is the energy carried out by the produced warm fluid. This is more
# accurate than the bare -dE/dt approximation because it isolates the
# matrix-to-fracture conductive heat extraction from changes in stored fluid
# energy caused by the circulating fluid.
Er    = fdata[:TotalThermalEnergy]
dEdt  = diff(vcat(Er[1,:]', Er), dims = 1) ./ dt  # rate of change of fracture stored energy
q_in  = fdata[:q_in]   # injector → fracture energy flux [W]
q_out = fdata[:q_out]  # fracture → producer energy flux [W]
power = dEdt .+ q_out .- q_in  # heat extracted from rock matrix [W]

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

# ## Second run: random fracture angles and apertures
# To understand how fracture orientation and aperture variability affect EGS
# performance, we run a second scenario where both properties are drawn from
# normal distributions. Fracture tilt angles are sampled from N(0°, 5°) and
# apertures from N(0.5 mm, 0.1 mm). All other parameters are kept identical.
case2 = Fimbul.egs(inj, prod, fracture_radius, fracture_spacing;
    rate              = 9250meter^3/day,
    temperature_inj   = convert_to_si(25.0, :Celsius),
    num_years         = num_years,
    fracture_theta    = (0.0, deg2rad(5.0)),    # N(0°, 5°) tilt
    fracture_aperture = (0.5e-3, 1e-4),           # N(0.5 mm, 0.1 mm)
    schedule_args     = (report_interval = si_unit(:year)/4,)
);

sim2, cfg2 = setup_reservoir_simulator(case2;
    info_level       = 2,
    output_substates = true,
    initial_dt       = 5.0,
    relaxation       = true
);
cfg2[:tolerances][:Fractures][:energy_conservation] = (CNV = Inf, EB = 1e-5, increment_dT = 1e-2)
cfg2[:tolerances][:Fractures][:mass_conservation]   = (CNV = Inf, MB = 1e-5, increment_dp_abs = 1e-2*si_unit(:bar))
sel2 = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :Reservoir)
push!(cfg2[:timestep_selectors], sel2);

results2 = simulate_reservoir(case2; simulator = sim2, config = cfg2)

states2, dt2, _ = Jutul.expand_to_ministeps(results2.result)
time2  = cumsum(dt2) ./ si_unit(:year)
fdata2 = Fimbul.get_egs_fracture_data(states2, case2)

frac_domain2 = case2.model.models[:Fractures].data_domain
frac_mesh2   = physical_representation(frac_domain2)
frac_geo2    = tpfv_geometry(frac_mesh2)
msh2 = physical_representation(reservoir_model(case2.model).data_domain)
geo2 = tpfv_geometry(msh2)
plot_egs_mesh(case2, "Run 2 – Random fractures, well spacing 100 m")

# ### Fracture ΔT visualization – random case
states_full2, dt_full2, _ = Jutul.expand_to_ministeps(results2.result)
time_full2 = cumsum(dt_full2) ./ si_unit(:year)

T0_frac2 = case2.state0[:Fractures][:Temperature]
ΔT_frac2 = [state[:Fractures][:Temperature] .- T0_frac2 for state in states_full2]
colorrange2 = extrema(vcat(ΔT_frac2...))

fxf2, fyf2, fzf2 = frac_geo2.cell_centroids[1,:], frac_geo2.cell_centroids[2,:], frac_geo2.cell_centroids[3,:]
xlim_f2 = [(extrema(fxf2) .+ diff(collect(extrema(fxf2))).*[-0.3, 0.3])...]
ylim_f2 = [(extrema(fyf2) .+ diff(collect(extrema(fyf2))).*[-0.1, 0.1])...]
zlim_f2 = [(extrema(fzf2) .+ diff(collect(extrema(fzf2))).*[-0.3, 0.3])...]
limits_f2 = (xlim_f2, ylim_f2, zlim_f2)

n_steps_f2 = length(ΔT_frac2)
steps2 = Int.(round.([0.125, 0.5, 1.0] .* n_steps_f2))

wells_dict2 = get_model_wells(case2.model; data_domain = true)
function plot_egs_wells2(ax; colors = [:red, :blue])
    for (i, (name, well)) in enumerate(wells_dict2)
        xy = permutedims(well[:cell_centroids][1:2, :])
        plot_mswell_values!(ax, case2.model, name, xy;
            geo = geo2, linewidth = 3, color = colors[i])
    end
end

fig = Figure(size = (650, 800))
for (n, ΔT_n) in enumerate(ΔT_frac2[steps2])
    ax_n = Axis3(fig[n, 1];
        perspectiveness = 0.5, zreversed = true, aspect = (1, 6, 1),
        azimuth = 1.2π, elevation = π/20, limits = limits_f2,
        title = "$(round(time_full2[steps2[n]], digits=1)) years", titlegap = -10)
    plot_cell_data!(ax_n, frac_mesh2, ΔT_n; colorrange = colorrange2,
        colormap = :seaborn_icefire_gradient)
    plot_egs_wells2(ax_n; colors = [:black, :black])
    hidedecorations!(ax_n)
end
Colorbar(fig[length(steps2)+1, 1];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange2,
    label = "ΔT (°C)", vertical = false, flipaxis = false)
fig

# ### Fracture metrics – random case
n_frac2  = length(fdata2[:y])
colors2  = cgrad(:BrBg, n_frac2, categorical = true)

function plot_fracture_data2(ax, time, data; stacked = false)
    df_prev = zeros(length(time))
    for (fno, df) in enumerate(eachcol(data))
        df = copy(df)
        if stacked
            df .+= df_prev
            poly!(ax, vcat(time, reverse(time)), vcat(df_prev, reverse(df));
                color = colors2[fno], strokecolor = :black, strokewidth = 1,
                label = "Fracture $fno")
            df_prev = df
        else
            lines!(ax, time, df; color = colors2[fno], linewidth = 2,
                label = "Fracture $fno")
        end
    end
end

xmax2   = round(maximum(time2))
limits2 = ((0, xmax2) .+ (-0.1, 0.1) .* xmax2, nothing)
xticks2 = 0:xmax2

first_step2 = findfirst(time2 .> 1/104)

fig = Figure(size = (1000, 800))
function make_axis2(title, ylabel, rno; kwargs...)
    Axis(fig[rno, 1]; title = title, xlabel = "Time (years)", ylabel = ylabel,
        limits = limits2, xticks = xticks2, kwargs...)
end

ax2 = make_axis2("Temperature", "T (°C)", 1)
temperature2 = fdata2[:Temperature][first_step2:end, :]
plot_fracture_data2(ax2, time2[first_step2:end],
    convert_from_si.(temperature2, :Celsius))
hidexdecorations!(ax2, grid = false)

Er2    = fdata2[:TotalThermalEnergy]
dEdt2  = diff(vcat(Er2[1,:]', Er2), dims = 1) ./ dt2
q_in2  = fdata2[:q_in]
q_out2 = fdata2[:q_out]
power2 = dEdt2 .+ q_out2 .- q_in2

ax2_pwr = make_axis2("Thermal power", "Power (MW)", 2)
plot_fracture_data2(ax2_pwr, time2[first_step2:end],
    power2[first_step2:end, :] ./ 1e6; stacked = true)
hidexdecorations!(ax2_pwr, grid = false)

energy_per_year2, cat2, dodge2 = [], Int[], Int[]
ix2 = vcat(0, [findfirst(isapprox.(time2, y; atol = 1e-2)) for y in 1:num_years]) .+ 1
ax2_e = make_axis2("Annual energy production per fracture", "Energy (GWh)", 3)
for (fno, pwr_f) in enumerate(eachcol(power2))
    energy_f = [sum(pwr_f[ix2[k]:ix2[k+1]-1] .* dt2[ix2[k]:ix2[k+1]-1]) / GWh
                for k in 1:length(ix2)-1]
    push!(energy_per_year2, energy_f)
    push!(cat2, 1:length(energy_f)...)
    push!(dodge2, fill(fno, length(energy_f))...)
end
barplot!(ax2_e, cat2, vcat(energy_per_year2...);
    dodge = dodge2, color = colors2[dodge2], strokecolor = :black, strokewidth = 1)
hidexdecorations!(ax2_e, grid = false)

ax2_frac = make_axis2("Annual energy fraction per fracture", "Fraction (-)", 4)
η2 = reduce(hcat, energy_per_year2)
η2 = η2 ./ sum(η2, dims = 2)
barplot!(ax2_frac, cat2, η2[:];
    dodge = dodge2, color = colors2[dodge2], strokecolor = :black, strokewidth = 1)

Legend(fig[2:3, 2], ax2_pwr)
fig

# ## Third run: wider well spacing (200 m)
# We repeat the uniform fracture scenario with a 200 m injector–producer
# spacing to study how increased well-to-well distance affects heat extraction.
inj3, prod3 = Fimbul.egs_well_coordinates(
    well_depth      = 2500.0meter,
    well_spacing_x  = 200.0meter,
    well_lateral    = 1000.0meter,
    bend_radius     = 200.0meter,
)

case3 = Fimbul.egs(inj3, prod3, fracture_radius, fracture_spacing;
    rate              = 9250meter^3/day,
    temperature_inj   = convert_to_si(25.0, :Celsius),
    num_years         = num_years,
    schedule_args     = (report_interval = si_unit(:year)/4,)
);
plot_egs_mesh(case3, "Run 3 – Uniform fractures, well spacing 200 m")

sim3, cfg3 = setup_reservoir_simulator(case3;
    info_level       = 2,
    output_substates = true,
    initial_dt       = 5.0,
    relaxation       = true
);
cfg3[:tolerances][:Fractures][:energy_conservation] = (CNV = Inf, EB = 1e-5, increment_dT = 1e-2)
cfg3[:tolerances][:Fractures][:mass_conservation]   = (CNV = Inf, MB = 1e-5, increment_dp_abs = 1e-2*si_unit(:bar))
sel3 = VariableChangeTimestepSelector(:Temperature, 5.0;
    relative = false, model = :Reservoir)
push!(cfg3[:timestep_selectors], sel3);

results3 = simulate_reservoir(case3; simulator = sim3, config = cfg3)

states3, dt3, _ = Jutul.expand_to_ministeps(results3.result)
time3  = cumsum(dt3) ./ si_unit(:year)
fdata3 = Fimbul.get_egs_fracture_data(states3, case3)

# ### Well performance – all three runs
plot_well_results([results.wells, results2.wells, results3.wells];
    names = ["Uniform 100 m", "Random 100 m", "Uniform 200 m"])

frac_domain3 = case3.model.models[:Fractures].data_domain
frac_mesh3   = physical_representation(frac_domain3)
frac_geo3    = tpfv_geometry(frac_mesh3)
msh3 = physical_representation(reservoir_model(case3.model).data_domain)
geo3 = tpfv_geometry(msh3)

# ### Fracture ΔT visualization – wide spacing
states_full3, dt_full3, _ = Jutul.expand_to_ministeps(results3.result)
time_full3 = cumsum(dt_full3) ./ si_unit(:year)

T0_frac3 = case3.state0[:Fractures][:Temperature]
ΔT_frac3 = [state[:Fractures][:Temperature] .- T0_frac3 for state in states_full3]
colorrange3 = extrema(vcat(ΔT_frac3...))

fxf3, fyf3, fzf3 = frac_geo3.cell_centroids[1,:], frac_geo3.cell_centroids[2,:], frac_geo3.cell_centroids[3,:]
xlim_f3 = [(extrema(fxf3) .+ diff(collect(extrema(fxf3))).*[-0.3, 0.3])...]
ylim_f3 = [(extrema(fyf3) .+ diff(collect(extrema(fyf3))).*[-0.1, 0.1])...]
zlim_f3 = [(extrema(fzf3) .+ diff(collect(extrema(fzf3))).*[-0.3, 0.3])...]
limits_f3 = (xlim_f3, ylim_f3, zlim_f3)

n_steps_f3 = length(ΔT_frac3)
steps3 = Int.(round.([0.125, 0.5, 1.0] .* n_steps_f3))

wells_dict3 = get_model_wells(case3.model; data_domain = true)
function plot_egs_wells3(ax; colors = [:red, :blue])
    for (i, (name, well)) in enumerate(wells_dict3)
        xy = permutedims(well[:cell_centroids][1:2, :])
        plot_mswell_values!(ax, case3.model, name, xy;
            geo = geo3, linewidth = 3, color = colors[i])
    end
end

fig = Figure(size = (650, 800))
for (n, ΔT_n) in enumerate(ΔT_frac3[steps3])
    ax_n = Axis3(fig[n, 1];
        perspectiveness = 0.5, zreversed = true, aspect = (1, 6, 1),
        azimuth = 1.2π, elevation = π/20, limits = limits_f3,
        title = "$(round(time_full3[steps3[n]], digits=1)) years", titlegap = -10)
    plot_cell_data!(ax_n, frac_mesh3, ΔT_n; colorrange = colorrange3,
        colormap = :seaborn_icefire_gradient)
    plot_egs_wells3(ax_n; colors = [:black, :black])
    hidedecorations!(ax_n)
end
Colorbar(fig[length(steps3)+1, 1];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange3,
    label = "ΔT (°C)", vertical = false, flipaxis = false)
fig

# ### Fracture metrics – wide spacing
n_frac3  = length(fdata3[:y])
colors3  = cgrad(:BrBg, n_frac3, categorical = true)

function plot_fracture_data3(ax, time, data; stacked = false)
    df_prev = zeros(length(time))
    for (fno, df) in enumerate(eachcol(data))
        df = copy(df)
        if stacked
            df .+= df_prev
            poly!(ax, vcat(time, reverse(time)), vcat(df_prev, reverse(df));
                color = colors3[fno], strokecolor = :black, strokewidth = 1,
                label = "Fracture $fno")
            df_prev = df
        else
            lines!(ax, time, df; color = colors3[fno], linewidth = 2,
                label = "Fracture $fno")
        end
    end
end

xmax3   = round(maximum(time3))
limits3 = ((0, xmax3) .+ (-0.1, 0.1) .* xmax3, nothing)
xticks3 = 0:xmax3
first_step3 = findfirst(time3 .> 1/104)

fig = Figure(size = (1000, 800))
function make_axis3(title, ylabel, rno; kwargs...)
    Axis(fig[rno, 1]; title = title, xlabel = "Time (years)", ylabel = ylabel,
        limits = limits3, xticks = xticks3, kwargs...)
end

ax3 = make_axis3("Temperature", "T (°C)", 1)
temperature3 = fdata3[:Temperature][first_step3:end, :]
plot_fracture_data3(ax3, time3[first_step3:end],
    convert_from_si.(temperature3, :Celsius))
hidexdecorations!(ax3, grid = false)

Er3    = fdata3[:TotalThermalEnergy]
dEdt3  = diff(vcat(Er3[1,:]', Er3), dims = 1) ./ dt3
q_in3  = fdata3[:q_in]
q_out3 = fdata3[:q_out]
power3 = dEdt3 .+ q_out3 .- q_in3

ax3_pwr = make_axis3("Thermal power", "Power (MW)", 2)
plot_fracture_data3(ax3_pwr, time3[first_step3:end],
    power3[first_step3:end, :] ./ 1e6; stacked = true)
hidexdecorations!(ax3_pwr, grid = false)

energy_per_year3, cat3, dodge3 = [], Int[], Int[]
ix3 = vcat(0, [findfirst(isapprox.(time3, y; atol = 1e-2)) for y in 1:num_years]) .+ 1
ax3_e = make_axis3("Annual energy production per fracture", "Energy (GWh)", 3)
for (fno, pwr_f) in enumerate(eachcol(power3))
    energy_f = [sum(pwr_f[ix3[k]:ix3[k+1]-1] .* dt3[ix3[k]:ix3[k+1]-1]) / GWh
                for k in 1:length(ix3)-1]
    push!(energy_per_year3, energy_f)
    push!(cat3, 1:length(energy_f)...)
    push!(dodge3, fill(fno, length(energy_f))...)
end
barplot!(ax3_e, cat3, vcat(energy_per_year3...);
    dodge = dodge3, color = colors3[dodge3], strokecolor = :black, strokewidth = 1)
hidexdecorations!(ax3_e, grid = false)

ax3_frac = make_axis3("Annual energy fraction per fracture", "Fraction (-)", 4)
η3 = reduce(hcat, energy_per_year3)
η3 = η3 ./ sum(η3, dims = 2)
barplot!(ax3_frac, cat3, η3[:];
    dodge = dodge3, color = colors3[dodge3], strokecolor = :black, strokewidth = 1)

Legend(fig[2:3, 2], ax3_pwr)
fig

# ## Comparison: uniform (100 m) vs. random (100 m) vs. uniform (200 m)
# We compare all three runs by looking at the aggregate power and cumulative
# energy, as well as the spread across individual fractures.

total_power1 = sum(power,  dims = 2)[:]
total_power2 = sum(power2, dims = 2)[:]
total_power3 = sum(power3, dims = 2)[:]

fig = Figure(size = (1000, 700))
xmax_c   = round(max(maximum(time), maximum(time2), maximum(time3)))
lim_c    = ((0, xmax_c) .+ (-0.1, 0.1) .* xmax_c, nothing)
xtick_c  = 0:xmax_c

function make_cax(title, ylabel, rno)
    Axis(fig[rno, 1]; title = title, xlabel = "Time (years)", ylabel = ylabel,
        limits = lim_c, xticks = xtick_c)
end

# ── Total power ──────────────────────────────────────────────────────────────
ax_c1 = make_cax("Total thermal power", "Power (MW)", 1)
lines!(ax_c1, time[first_step:end],  total_power1[first_step:end]  ./ 1e6;
    color = :steelblue,   linewidth = 2, label = "Uniform 100 m")
lines!(ax_c1, time2[first_step:end], total_power2[first_step:end]  ./ 1e6;
    color = :orangered,   linewidth = 2, linestyle = :dash, label = "Random 100 m")
lines!(ax_c1, time3[first_step:end], total_power3[first_step:end]  ./ 1e6;
    color = :forestgreen, linewidth = 2, linestyle = :dot,  label = "Uniform 200 m")
axislegend(ax_c1; position = :rb)
hidexdecorations!(ax_c1, grid = false)

# ── Per-fracture power spread ────────────────────────────────────────────────
ax_c2 = make_cax("Per-fracture power (shaded band = min–max)", "Power (MW)", 2)
function band_fractures!(ax, t, pwr; color = :steelblue, label = "")
    n_f        = size(pwr, 2)
    total_pwr  = sum(pwr, dims = 2)[:]
    pmin = [minimum(r) for r in eachrow(pwr)]
    pmax = [maximum(r) for r in eachrow(pwr)]
    band!(ax, t, pmin ./ 1e6, pmax ./ 1e6; color = (color, 0.3))
    lines!(ax, t, total_pwr ./ n_f ./ 1e6; color = color, linewidth = 2,
        label = label)
end
fs = first_step
band_fractures!(ax_c2, time[fs:end],  power[fs:end,  :];  color = :steelblue,   label = "Uniform 100 m (mean)")
band_fractures!(ax_c2, time2[fs:end], power2[fs:end, :];  color = :orangered,   label = "Random 100 m (mean)")
band_fractures!(ax_c2, time3[fs:end], power3[fs:end, :];  color = :forestgreen, label = "Uniform 200 m (mean)")
axislegend(ax_c2; position = :rb)
hidexdecorations!(ax_c2, grid = false)

# ── Cumulative energy per year ───────────────────────────────────────────────
ax_c3 = make_cax("Annual energy production", "Energy (GWh)", 3)
function annual_energy(pwr, t, dt_v)
    ix = vcat(0, [findfirst(isapprox.(t, y; atol = 1e-2)) for y in 1:num_years]) .+ 1
    [sum(pwr[ix[k]:ix[k+1]-1] .* dt_v[ix[k]:ix[k+1]-1]) / GWh for k in 1:length(ix)-1]
end
ann1 = annual_energy(total_power1, time,  dt)
ann2 = annual_energy(total_power2, time2, dt2)
ann3 = annual_energy(total_power3, time3, dt3)
years = 1:num_years
barplot!(ax_c3, years .- 0.28, ann1; width = 0.25, color = :steelblue,   label = "Uniform 100 m")
barplot!(ax_c3, years,         ann2; width = 0.25, color = :orangered,   label = "Random 100 m")
barplot!(ax_c3, years .+ 0.28, ann3; width = 0.25, color = :forestgreen, label = "Uniform 200 m")
axislegend(ax_c3; position = :rt)

fig
