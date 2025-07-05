# # Borehole Thermal Energy Storage (BTES)
# This example demonstrates how to set up and simulate a Borehole Thermal Energy
# Storage (BTES) system using Fimbul. The BTES system consists of multiple
# boreholes coupled in series and parallel.
using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

# ## Set up simulation case
# We consider a domain with 48 BTES wells that all reach 100 m depth. The BTES
# wells are charged during the summer at 90 °C and discharged during the winter
# months at 10 °C. This cycle is simulated for 4 years.
case, sections = btes(num_wells = 48, depths = [0.0, 0.5, 100, 125],
    charge_months = ["April", "May", "June", "July", "August", "September"],
    discharge_months = ["October", "November", "December", "January", "February", "March"],
    num_years = 4,
);

# ### Visualize BTES system
# The 48 BTES wells are placed in a pattern so that all wells are approximately
# the same distance apart, and divided into 6 sections of 8 wells each. During
# charging, water is injected into the inner-most well of each section at 0.5
# l/s, with water flowing from the inner-most well to the outer-most well.
# During discharging, water flows in the opposite direction, from the outer-most
# well to the inner-most well.
msh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1]; zreversed = true, aspect = :data,
azimuth = 0, elevation = π/2)
Jutul.plot_mesh_edges!(ax, msh, alpha = 0.2)
colors = Makie.wong_colors()
lns, labels = [], String[]
for (sno, (section, wells)) in enumerate(sections)
    for (wno, wname) in enumerate(wells)
        well = case.model.models[wname].domain.representation
        l = plot_well!(ax, msh, well; color = colors[sno], fontsize = 0)
        if wno == 1
            push!(lns, l)
            push!(labels, "Section $sno")
        end
    end
end
Legend(fig[1, 2], lns, labels);
fig

# ## Set up reservoir simulator
simulator, config = setup_reservoir_simulator(case;
    tol_cnv = 1e-2,
    tol_mb = 1e-6,
    timesteps = :auto,
    initial_dt = 5.0,
    relaxation = true);

# Changing from charging to discharging results in a thermal shock that is
# challenging to resolve for the nonlinear solver. We therefore use a timestep
# selector that reduces the timestep to 5 seconds when the control changes.
sel = JutulDarcy.ControlChangeTimestepSelector(
    case.model, 0.0, convert_to_si(5.0, :second))
push!(config[:timestep_selectors], sel)
config[:timestep_max_decrease] = 1e-6
for ws in well_symbols(case.model)
    config[:tolerances][ws][:default] = 1e-2
end

# ## Simulate the case
results = simulate_reservoir(case;
simulator=simulator, config=config, info_level=0);

# ### Interactive visualization
# We plot the reservoir state and the well output interactively.
plot_reservoir(case.model, results.states;
well_fontsize = 0, key = :Temperature, step = length(case.dt),
colormap = :seaborn_icefire_gradient)

plot_well_results(results.wells, field = :temperature)

# ## Temperature evolution in the reservoir
# We plot the temperature evolution in the reservoir after each of the four
# charging stages. To enhance the visualization, we only plot the cells below 50
# m depth that have a temperature above 15 °C.
dstep = Int64.(length(case.dt)/(4*2))
steps = dstep:2*dstep:length(case.dt)
fig = Figure(size = (1000, 750))
geo = tpfv_geometry(msh)
bottom = geo.cell_centroids[3,:] .>= 50.0
T_min, T_max = Inf, -Inf
for (sno, step) in enumerate(steps)
    ax = Axis3(fig[(sno-1)÷2+1, (sno-1)%2+1];
    limits = (-50, 50, -50, 50, 40, 125),
    title = "Charge $sno", zreversed = true, aspect = :data)
    T = convert_from_si.(results.states[step][:Temperature], :Celsius)
    cells = bottom .& (T .>= 15.0)
    T_min = min(T_min, minimum(T[cells]))
    T_max = max(T_max, maximum(T[cells]))
    plot_cell_data!(ax, msh, T; cells = cells, colormap = :seaborn_icefire_gradient)
    hidedecorations!(ax)
end
Colorbar(fig[1:2, 3]; 
colormap = :seaborn_icefire_gradient, colorrange = (T_min, T_max), 
label = "T ≥ 15.0 °C", vertical = true)
fig

# ## Temperature evolution in the wells
# Designing a BTES system requires in-depth understanding of how the temperature
# evolves as water runs through the pipes. We plot the temperature in the first
# section for each report step during the first charge and discharge stages.
# This kind of visualization can be used when designing a BTES system to
# determine e.g., how many wells are needed for a given depth, and how wells
# should be coupled in series/parallel.

# We set up a function that plots the temperature evolution in a given section
# of the BTES system for a given set of timesteps. The function plots the
# temperature as a function of well depth for all wells in the section.
function plot_btes_temperature(ax, section, timesteps)
    colors = cgrad(:ice, length(timesteps), categorical = true)
    msh = physical_representation(reservoir_model(case.model).data_domain)
    geo = tpfv_geometry(msh)
    time = convert_from_si.(cumsum(case.dt), :day)
    out = []
    for (sno, step) in enumerate(timesteps)
        label = String("$(time[step]) days")
        for (wno, well) in enumerate(sections[Symbol("S$section")])
            lbl = wno == 1 ? label : nothing
            T = results.result.states[step][well][:Temperature]
            T = convert_from_si.(T, :Celsius)
            wm = case.model.models[well]
            num_nodes = length(wm.domain.representation.perforations.reservoir)
            plot_mswell_values!(ax, case.model, well, T;
            nodes = 1:num_nodes, geo = geo, color = colors[sno], linewidth = 6, label = lbl)
        end
    end
end

fig = Figure(size = (1000, 750))
section = 1
time = convert_from_si.(cumsum(case.dt), :year)

# Temperature evolution during the first charging
well = sections[Symbol("S$section")][1]
T_ch = convert_to_si(90.0, :Celsius)
is_charge = [f[:Facility].control[well].temperature == T_ch for f in case.forces]
stop = findfirst(diff(is_charge) .< 0)
ax = Axis(fig[1, 1]; title = "Temperature evolution, section $section (charging)",
xlabel = "T [°C]", ylabel = "Depth [m]", yreversed = true)
plot_btes_temperature(ax, section, 1:stop)
Legend(fig[1,2], ax; fontsize = 20)

# Temperature evolution during the first discharging
well = sections[Symbol("S$section")][end-1]
T_dch = convert_to_si(10.0, :Celsius)
is_discharge = [f[:Facility].control[well].temperature == T_dch for f in case.forces]
start = findfirst(diff(is_discharge) .> 0)+1
stop = findfirst(diff(is_discharge) .< 0)
ax = Axis(fig[2, 1]; title = "Temperature evolution, section $section (discharging)",
xlabel = "T [°C]", ylabel = "Depth [m]", yreversed = true)
plot_btes_temperature(ax, section, start:stop)
Legend(fig[2,2], ax; fontsize = 20)

fig