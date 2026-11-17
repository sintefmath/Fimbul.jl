# # Borehole Thermal Energy Storage (BTES)
# <tags: Storage, BTES>
# This example demonstrates how to set up and simulate seasonal Borehole Thermal
# Energy Storage (BTES) using Fimbul, and compares three patterns for placing
# the boreholes. The BTES system consists of multiple boreholes coupled in
# series and parallel.
using Jutul, JutulDarcy
using Fimbul
using HYPRE
using GLMakie

# ## Set up simulation cases
# We create BTES systems with 30 wells that extend to 100 m depth, divided into
# 5 sectors of 6 wells each. The wells are placed in a sunflower, a circular
# and a rectangular (6 x 5) pattern, all with 5 m between neighboring wells.
# The systems operate on a seasonal cycle: wells are charged during summer
# months (April to September) with hot water at 90°C, and discharged during
# winter months (October to March) with cold water at 10°C. This seasonal
# operation is simulated over a 4-year period to study the thermal energy
# storage performance.
#
# The wells of the sunflower pattern are grouped into sectors along its spiral
# spines (`sector_grouping = :spiral`, the default). For the circular pattern,
# the wells are grouped into sectors from the center outwards
# (`sector_grouping = :sequential`), instead of into circle sectors, which is
# the default. For the rectangular pattern, the sectors are blocks of 3 x 2
# wells, with the last sector in the top row (`sector_grouping = :cartesian`).
#
# !!! tip "Adjusting the well fields"
#     The size of the fields is set by `num_wells` and `num_sectors`. Instead
#     of the number of wells, the spacing between them can be given
#     (`well_spacing`), and the size of the field can be set with `radius`
#     (sunflower and circular) or `extent = (Lx, Ly)` (rectangular), see
#     [`btes`](@ref). For the rectangular pattern, the blocks that make up the
#     sectors are as square as possible, or set by `sector_dims = (nx, ny)`;
#     for example, 45 wells in 5 sectors gives four blocks of 3 x 3 wells and a
#     last sector of 9 wells in the two top rows. The sunflower sectors are
#     made of whole spines when the number of sectors is 2, 3, 5 or 8;
#     otherwise, a few outer wells are moved between neighboring sectors to
#     even them out. Larger fields take longer to simulate.
patterns = [:sunflower, :circular, :rectangular]
pattern_names = Dict(:sunflower => "Sunflower", :circular => "Circular",
    :rectangular => "Rectangular")
num_wells, num_sectors = 30, 5
pattern_args = Dict(:sunflower => (sector_grouping = :spiral,),
    :circular => (sector_grouping = :sequential,),
    :rectangular => (sector_grouping = :cartesian,))
case_args = (
    num_wells = num_wells,
    num_sectors = num_sectors,
    depths = [0.0, 0.5, 100, 125],
    charge_period = ["April", "September"],
    discharge_period = ["October", "March"],
    num_years = 4,
    reversed_discharge = true,
)
cases = Dict(p => btes(p; pattern_args[p]..., case_args...) for p in patterns);

# ### Visualize the BTES layouts
# During charging, water is injected at 0.5 l/s into the first well of each
# sector and flows sequentially through all wells in series to the last well.
# During discharging, the flow is reversed (`reversed_discharge = true`), with
# water flowing from the last well back to the first. Each sector is modelled
# as a single closed-loop well holding its boreholes in series, which is much
# faster to simulate than one well per borehole (`well_grouping = :borehole`).
#
# For the sunflower pattern, the sectors are bundles of adjacent spiral spines
# (see below), with the flow running from the center outwards. With 30 wells in
# 5 sectors, each sector is a single spine of the 5-spine family. For the
# circular pattern, the wells are numbered from the center outwards,
# counter-clockwise on each ring, and the sectors are consecutive groups in
# this order. For the rectangular pattern, the sectors are blocks of 3 x 2
# wells, tiling the field from the bottom up, and the wells left over in the
# top row form the last sector. The flow runs through each sector column by
# column.
#
# The structure underlying each pattern is shown in light gray:
# - Sunflower: well ``n`` sits at radius ``R\sqrt{(n - 1/2)/(N - 1/2)}`` and
#   angle ``nΔθ``, where ``Δθ = 2π/φ^2`` and ``φ`` is the golden ratio. Wells
#   ``F`` apart in this order, for consecutive Fibonacci numbers ``F``, line up
#   along spirals, the "spines" of the sunflower.
# - Circular: a center well surrounded by concentric rings, one well spacing
#   apart.
# - Rectangular: the Cartesian grid on which the wells are placed.
structure_color = (:gray, 0.4)

function plot_sunflower_spines!(ax, N, R; num_points = 200)
    Δθ = 2π/Base.MathConstants.golden^2
    ## The two largest spine families with at least two wells per spine, which
    ## include the family that the sectors are made of
    fib = [1, 2]
    while fib[end] + fib[end-1] <= N/2
        push!(fib, fib[end] + fib[end-1])
    end
    for F in fib[end-1:end]
        ## Turning angle between wells F apart, taken as the smallest equivalent
        ## angle, so that the spine through them is a smooth spiral
        δ = mod(F*Δθ + π, 2π) - π
        for n in 1:F
            t = range(0.0, (N - n)/F, length = num_points)
            r = R.*sqrt.((n .+ F.*t .- 0.5)./(N - 0.5))
            θ = n*Δθ .+ δ.*t
            lines!(ax, r.*cos.(θ), r.*sin.(θ); color = structure_color, linewidth = 1)
        end
    end
end

function plot_circular_rings!(ax, R, spacing)
    for k in 1:round(Int, R/spacing)
        arc!(ax, Point2f(0, 0), k*spacing, 0, 2π; color = structure_color, linewidth = 1)
    end
end

function plot_cartesian_grid!(ax, xy)
    xs, ys = sort(unique(xy[1, :])), sort(unique(xy[2, :]))
    for x in xs
        lines!(ax, [x, x], [ys[1], ys[end]]; color = structure_color, linewidth = 1)
    end
    for y in ys
        lines!(ax, [xs[1], xs[end]], [y, y]; color = structure_color, linewidth = 1)
    end
end

# We plot the well positions of each pattern, colored by sector, with the wells
# of each sector connected in the order the fluid passes them during charging.
# The first well of each sector is marked with a larger symbol.
colors = Makie.wong_colors()
sector_colors = Makie.to_colormap(:tab10)
fig = Figure(size = (1200, 450))
for (pno, p) in enumerate(patterns)
    field, spacing = Fimbul.btes_well_field(p;
        num_wells = num_wells, num_sectors = num_sectors, pattern_args[p]...,
        depths = case_args.depths)
    xy_all = hcat([wc[1:2, 1] for sector in field for wc in sector]...)
    R = maximum(sqrt.(sum(xy_all.^2, dims = 1)))
    ax = Axis(fig[1, pno]; title = pattern_names[p], aspect = DataAspect(),
        xlabel = "x [m]", ylabel = "y [m]")
    if p == :sunflower
        plot_sunflower_spines!(ax, num_wells, R)
    elseif p == :circular
        plot_circular_rings!(ax, R, spacing)
    else
        plot_cartesian_grid!(ax, xy_all)
    end
    for (sno, sector) in enumerate(field)
        xy = hcat([wc[1:2, 1] for wc in sector]...)
        lines!(ax, xy[1, :], xy[2, :]; color = sector_colors[sno], linewidth = 2)
        scatter!(ax, xy[1, :], xy[2, :]; color = sector_colors[sno], markersize = 10,
            label = "Sector $sno")
        scatter!(ax, xy[1, 1:1], xy[2, 1:1]; color = sector_colors[sno], markersize = 20)
    end
    pno == 1 && Legend(fig[2, 1:3], ax; orientation = :horizontal)
end
fig

# ## Set up reservoir simulators
# We configure the simulator with specific tolerances and timestep controls to
# handle the challenging thermal transients in the BTES system. We also store
# the states of all ministeps (`output_substates = true`), which we use to
# compute the energy exchanged with the ground accurately.
#
# The transition from charging to discharging creates a thermal shock that is
# numerically challenging for the nonlinear solver. We use a specialized
# timestep selector that reduces the timestep to 5 seconds during control
# changes to maintain numerical stability and convergence.
function setup_btes_simulator(case)
    simulator, config = setup_reservoir_simulator(case;
        tol_cnv = 1e-2,
        tol_mb = 1e-6,
        tol_cnve_well = 1e-2,
        tol_cnv_well = 1e-2,
        timesteps = :auto,
        initial_dt = 5.0,
        target_its = 10,
        relaxation = true,
        output_substates = true,
        info_level = -1,
    )
    sel = JutulDarcy.ControlChangeTimestepSelector(
        case.model, 0.0, convert_to_si(5.0, :second))
    push!(config[:timestep_selectors], sel)
    config[:timestep_max_decrease] = 1e-6
    return simulator, config
end

# ## Simulate the BTES systems
# Run the full 4-year simulation for each pattern. Note that these simulations
# can take a few minutes to run. Setting `info_level = 0` will show a progress
# bar while the simulations run.
results = Dict(map(patterns) do p
    case = cases[p]
    simulator, config = setup_btes_simulator(case)
    p => simulate_reservoir(case;
        simulator = simulator, config = config, info_level = 0)
end);

# Since we stored the ministep states, `result.states` holds the reservoir
# state after every ministep, while `result.result.states` holds the full state
# after each report step. Below, we use the latter whenever we need the state
# at a given report step.
#
# Charging and discharging steps are found from the injection temperature. The
# sector wells are injected into at their supply part during charging, and at
# their return part during discharging, since `reversed_discharge = true`. All
# patterns share the same schedule.
T_ch = convert_to_si(90.0, :Celsius)
T_dch = convert_to_si(10.0, :Celsius)
function is_injecting(forces, T_inj)
    return any(c -> c isa InjectorControl && c.temperature == T_inj,
        values(forces[:Facility].control))
end
schedule = cases[patterns[1]].forces
is_charge = [is_injecting(f, T_ch) for f in schedule]
is_discharge = [is_injecting(f, T_dch) for f in schedule]

# Contiguous ranges of steps where `mask` holds, one per cycle. The state after
# a step is the state at its end, so the last step of each range gives the
# state at the end of the period.
function step_ranges(mask)
    starts = findall(i -> mask[i] && (i == 1 || !mask[i-1]), eachindex(mask))
    stops = findall(i -> mask[i] && (i == length(mask) || !mask[i+1]), eachindex(mask))
    return [a:b for (a, b) in zip(starts, stops)]
end
charge_ranges = step_ranges(is_charge)
discharge_ranges = step_ranges(is_discharge)

# ## Results for the sunflower pattern
# ### Interactive visualization of results
# We first plot the final temperature distribution in the reservoir and well
# performance over time for the sunflower pattern. The reservoir plot shows the
# thermal plumes developed around each well sector, while the well results show
# the temperature and flow rate evolution throughout the simulation.
case, result = cases[:sunflower], results[:sunflower]
plot_reservoir(case.model, result.states;
    well_fontsize = 0,
    key = :Temperature,
    step = length(result.states),
    colormap = :seaborn_icefire_gradient)

#
plot_well_results(result.wells, field = :temperature)

# ### Temperature evolution within the boreholes
# Understanding the temperature profiles within individual boreholes is crucial
# for BTES system design. We analyze how temperature varies with depth along the
# wellbore during both charging and discharging phases. This information can
# help determine optimal well depth, number of wells needed, and the most
# effective series/parallel coupling configuration for a given thermal storage
# capacity requirement.
#
# To this end, we define a utility function to plot temperature evolution in a
# BTES sector. This function visualizes temperature as a function of well depth
# for all wells in a given sector at specified timesteps. Each sector is a
# single well, and the `borehole` field gives the borehole of each well node.
function plot_btes_temperature(ax, case, result, sector, timesteps)
    colors = cgrad(:ice, length(timesteps), categorical = true)
    msh = physical_representation(reservoir_model(case.model).data_domain)
    geo = tpfv_geometry(msh)
    time = convert_from_si.(cumsum(case.dt), :day)
    well = only(case.input_data[:sectors][Symbol("S$sector")])
    d = case.model.models[well].data_domain
    borehole, tag = d[:borehole, Cells()], d[:tag, Cells()]
    for (sno, step) in enumerate(timesteps)
        label = String("$(round((time[step]),digits=1)) days")
        T = result.result.states[step][well][:Temperature]
        T = convert_from_si.(T, :Celsius)
        for b in 1:maximum(borehole)
            lbl = b == 1 ? label : nothing
            nodes = findall(i -> borehole[i] == b &&
                tag[i] in (:pipe_left, :pipe_right), eachindex(borehole))
            plot_mswell_values!(ax, case.model, well, T;
            nodes = nodes, geo = geo, color = colors[sno], linewidth = 6, label = lbl)
        end
    end
end

# #### Temperature evolution during the first charging cycle
# We first visualize how thermal energy propagates through the wellbore network
# as hot water flows from the innermost to the outermost well in the sector.
fig = Figure(size = (1000, 750))
sector = 1
ax = Axis(fig[1, 1]; title = "Temperature evolution, sector $sector (charging)",
xlabel = "T [°C]", ylabel = "Depth [m]", yreversed = true)
plot_btes_temperature(ax, case, result, sector, charge_ranges[1])
Legend(fig[1,2], ax; fontsize = 20)

# #### Temperature evolution during the first discharging cycle
# Next, we show how stored thermal energy is extracted as cold water flows
# through the sector in the opposite direction, from the outermost to the
# innermost well.
ax = Axis(fig[2, 1]; title = "Temperature evolution, sector $sector (discharging)",
xlabel = "T [°C]", ylabel = "Depth [m]", yreversed = true)
plot_btes_temperature(ax, case, result, sector, discharge_ranges[1])
Legend(fig[2,2], ax; fontsize = 20)

fig

# ## Comparison of the patterns
# ### Temperature evolution in the reservoir
# Next, we compare how thermal energy storage develops in the subsurface for
# the three patterns. We examine the temperature distribution at the end of
# each of the four charging periods to see how the thermal plume grows. For
# better visualization, we only display cells below 50 m depth with
# temperatures above 15°C (warmer than the natural ground temperature), seen
# from above with x to the right and y upwards.
axis_args = (
    zreversed = true,
    aspect = :data,
    elevation = π/2,
    azimuth = -π/2,
    perspectiveness = 0.0,
)
T_range = (15.0, convert_from_si(T_ch*0.95, :Celsius))
steps = last.(charge_ranges)
fig = Figure(size = (1100, 850))
for (pno, p) in enumerate(patterns)
    case, result = cases[p], results[p]
    msh = physical_representation(reservoir_model(case.model).data_domain)
    geo = tpfv_geometry(msh)
    bottom = geo.cell_centroids[3,:] .>= 50.0
    for (cno, step) in enumerate(steps)
        ax = Axis3(fig[pno, cno]; axis_args...,
            limits = (-40, 40, -40, 40, 40, 125))
        T = result.result.states[step][:Reservoir][:Temperature]
        T = convert_from_si.(T, :Celsius)
        plot_cell_data!(ax, msh, T; cells = bottom,
            colormap = :seaborn_icefire_gradient, colorrange = T_range)
        hidedecorations!(ax)
        hidespines!(ax)
    end
    Label(fig[pno, 0], pattern_names[p]; rotation = π/2, fontsize = 18,
        font = :bold, tellheight = false)
end
for cno in eachindex(steps)
    Label(fig[0, cno], "Charge $cno"; fontsize = 18, font = :bold,
        tellwidth = false)
end
Colorbar(fig[1:length(patterns), length(steps)+1];
    colormap = :seaborn_icefire_gradient, colorrange = T_range,
    label = "T [°C]", vertical = true)
fig

# ### Charged and discharged energy
# Finally, we compare how much energy each pattern stores and recovers. The net
# energy exchanged with the ground over a ministep of length Δt is
# ``Δt \sum_w q_w h_w``, where ``q_w`` is the surface mass rate of well part
# ``w`` (positive for injection) and ``h_w`` its specific enthalpy at the
# surface. Summing over the charging steps of a cycle gives the charged energy,
# and summing over the discharging steps gives the discharged energy (with
# opposite sign). The ratio of the two is the recovery factor of the cycle.
#
# The substates hold the ministeps of each report step but the last, whose
# state is the report state itself.
function exchanged_energy(result, steps)
    E = 0.0
    for step in steps
        report_state = result.result.states[step]
        ministep_states = [report_state[:substates]..., report_state]
        dt = [m[:dt] for m in result.result.reports[step][:ministeps] if m[:success]]
        @assert length(ministep_states) == length(dt)
        for (state, Δt) in zip(ministep_states, dt)
            f = state[:Facility]
            E += Δt*sum(f[:TotalSurfaceMassRate].*f[:SurfaceEnthalpy])
        end
    end
    return E
end

MWh = si_unit(:mega)*si_unit(:watt)*si_unit(:hour)
E_charged = Dict(p =>
    [exchanged_energy(results[p], r) for r in charge_ranges]./MWh
    for p in patterns)
E_discharged = Dict(p =>
    [-exchanged_energy(results[p], r) for r in discharge_ranges]./MWh
    for p in patterns)

# We plot the charged energy as light bars and the discharged energy as dark
# bars on top, with the recovery factor of each cycle above.
num_cycles = min(length(charge_ranges), length(discharge_ranges))
num_patterns = length(patterns)
width = 0.8/num_patterns
E_max = maximum(maximum.(values(E_charged)))
fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1]; xlabel = "Cycle", ylabel = "Energy [MWh]",
    xticks = 1:num_cycles,
    title = "Charged (light) and discharged (dark) energy per cycle")
for (pno, p) in enumerate(patterns)
    x = (1:num_cycles) .+ (pno - (num_patterns + 1)/2)*width
    E_ch = E_charged[p][1:num_cycles]
    E_dch = E_discharged[p][1:num_cycles]
    barplot!(ax, x, E_ch; width = width, color = (colors[pno], 0.35),
        strokecolor = :black, strokewidth = 1)
    barplot!(ax, x, E_dch; width = width, color = colors[pno],
        strokecolor = :black, strokewidth = 1, label = pattern_names[p])
    text!(ax, x, E_ch .+ 0.01*E_max;
        text = ["$(round(Int, 100*r))%" for r in E_dch./E_ch],
        align = (:center, :bottom), fontsize = 14)
end
ylims!(ax, 0, 1.1*E_max)
Legend(fig[2, 1], ax; orientation = :horizontal)
fig
