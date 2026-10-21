# # Deep coaxial borehole heat exchanger (BHE) demo
# <tags: Production, Workflow>
# This example demonstrates simulation and calibration of a deep coaxial
# closed-loop well, also known as a coaxial borehole heat exchanger (BHE). The
# well is defined by a general trajectory (m×3 matrix) and consists of two
# concentric pipes: an inner pipe and an outer annulus, separated from the
# rock by a layer of grout.
#
# We explore three scenarios:
# 1. **Layered reservoir**: Three geological layers, where the middle layer
#    has a distinctly lower thermal conductivity than the other two. This
#    model represents the real system in the rest of the example.
# 2. **Inner pipe conductivity study**: The same reservoir with varying inner
#    pipe wall thermal conductivity to demonstrate the effect of heat loss
#    between the supply and return pipes.
# 3. **Calibration of a digital twin**: A model of the same system with
#    homogeneous geology – a digital twin (DT) – is calibrated against
#    observations from the layered model. We calibrate the thermal
#    conductivity to two different types of data: the production temperature
#    at the wellhead over time, and the temperature along the well at the end
#    of the simulation.

# Add required modules to namespace
using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie

# Useful SI units
meter, hour, watt, Kelvin, joule, kilogram, year = si_units(
    :meter, :hour, :watt, :Kelvin, :joule, :kilogram, :year);

# ## Common parameters
# All scenarios share the same well trajectory, operational parameters,
# schedule and mesh settings. The well is vertical and 2500 m deep. Water is
# injected into the outer annulus at 10 °C and produced from the inner pipe,
# which gives a higher production temperature than the opposite flow
# direction: The cold water flowing down the annulus is in direct thermal
# contact with the rock, where the temperature difference is largest. We
# simulate five years of operation with monthly reports.
#
# The mesh is refined horizontally towards the well, with 25 m cells
# vertically along the well. The mesh has an additional interval below the
# bottom of the well with coarser cells. All rock properties except the thermal
# conductivity are the same everywhere, so that the thermal conductivity is the
# property that determines the heat extraction.
base_args = (
    well_trajectory = [0.0 0.0 0.0; 0.0 0.0 2500.0],
    inject_into = :outer,
    rate = 25meter^3/hour,
    temperature_inj = convert_to_si(10.0, :Celsius),
    num_years = 5,
    report_interval = year/12,
    hz = [25.0, 25.0, 25.0, 100.0],
    hxy_min = 5.0,
    mesh_args = (offset = 250.0, offset_rel = missing),
    permeability = 1e-3*si_unit(:darcy),
    porosity = 0.01,
    rock_heat_capacity = 900*joule/(kilogram*Kelvin),
    rock_density = 2600*kilogram/meter^3,
);

# ## Simulator helper
# A reusable function to set up the simulator with consistent settings, which
# is also used for the simulations performed during the calibration. We limit
# the temperature change per time step in the reservoir and the well. For the
# forward simulations, we only store the states at the report steps. The
# gradients computed during the calibration require the states of all time
# steps, so these are stored during the calibration.
function setup_simulator(case; output_substates = false)
    sim, cfg = setup_reservoir_simulator(case;
        output_substates = output_substates,
        info_level = -1,
        initial_dt = 120.0,
        presolve_wells = true,
        relaxation = true)
    for model in (:Reservoir, :CoaxialWell_supply)
        sel = VariableChangeTimestepSelector(:Temperature, 5.0;
            relative = false, model = model)
        push!(cfg[:timestep_selectors], sel)
    end
    return sim, cfg
end

function run_case(case)
    sim, cfg = setup_simulator(case)
    return simulate_reservoir(case; simulator = sim, config = cfg)
end

# ## Scenario 1 – Layered reservoir
# The reservoir consists of three geological layers. The middle layer has a
# distinctly lower thermal conductivity than the other two:
# - **Upper layer** (0–800 m): 3.0 W/(m·K)
# - **Middle layer** (800–1700 m): 1.5 W/(m·K)
# - **Lower layer** (1700–3000 m): 3.5 W/(m·K)
#
# The lower layer is split in two at the bottom of the well, to allow for a
# coarser vertical resolution below the well.
λ_unit = watt/(meter*Kelvin)
depths_layers = [0.0, 800.0, 1700.0, 3000.0]
λ_layers = [3.0, 1.5, 3.5]*λ_unit
depths = [0.0, 800.0, 1700.0, 2500.0 + 1e-3, 3000.0]
case_ref = coaxial_bhe(; base_args...,
    depths = depths, rock_thermal_conductivity = [λ_layers; λ_layers[end]]);

# ### Thermal conductivity and well
# The layers are shown with a quarter of the domain cut away, so that the well
# in the center is visible. The well is plotted using its trajectory.
function plot_trajectory!(ax, trajectory; color = :black, linewidth = 4, kwarg...)
    return lines!(ax, trajectory[:, 1], trajectory[:, 2], trajectory[:, 3];
        color = color, linewidth = linewidth, kwarg...)
end

rdomain_ref = reservoir_model(case_ref.model).data_domain
msh = physical_representation(rdomain_ref)
x = rdomain_ref[:cell_centroids]
fig = Figure(size = (700, 800))
ax = Axis3(fig[1, 1]; zreversed = true, perspectiveness = 0.75,
    aspect = (1, 1, 3))
colormap = cgrad(:BrBG_5, 3, categorical = true)
cutout = .!(x[1, :] .< 10.0 .&& x[2, :] .< 10.0)
for i in 1:3
    cells = findall(cutout .&& (rdomain_ref[:rock_thermal_conductivity] .== λ_layers[i]))
    Jutul.plot_mesh!(ax, msh; cells = cells, color = colormap[i])
end
plot_trajectory!(ax, base_args.well_trajectory)
elements = [PolyElement(color = colormap[i], strokecolor = :black, strokewidth = 1)
    for i in 1:3]
labels = ["$name: $(λ_layers[i]/λ_unit) W/(m·K)"
    for (i, name) in enumerate(("Upper layer", "Middle layer", "Lower layer"))]
Legend(fig[1, 2], elements, labels, "Thermal conductivity", valign = :center)
fig

# ### Simulate
results_ref = run_case(case_ref);

# ### Temperature change
# We plot the temperature after five years relative to the initial
# temperature, with half of the domain cut away. The rock around the well is
# cooled, most strongly in the lower layer, where the rock temperature is
# highest and heat is efficiently supplied from the surrounding rock due to
# the high conductivity. The color scale is saturated at ±10 K to show the
# extent of the region affected by the well, since the largest changes are
# limited to the cells closest to the well.
cutout = .!(x[1, :] .< 0.0 .&& x[2, :] .< 0.0)
ΔT = results_ref.states[end][:Temperature] .- case_ref.state0[:Reservoir][:Temperature]
fig = Figure(size = (700, 800))
ax = Axis3(fig[1, 1]; zreversed = true, perspectiveness = 0.75,
    aspect = (1, 1, 3), title = "Temperature change after 5 years")
plt = Jutul.plot_cell_data!(ax, msh, ΔT; cells = cutout,
    colormap = :vik, colorrange = (-10.0, 10.0), shading = NoShading)
Colorbar(fig[1, 2], plt; label = "ΔT (K)")
fig

# ### Well temperature profiles
# The temperature with depth in the different parts of the well after five
# years: the inner pipe, the outer annulus and the grout. The black line shows
# the reservoir temperature in the cells along the well for reference. The
# grout temperature follows the temperature of the rock next to the well,
# which is drawn down less in layers with high conductivity, where heat is
# more efficiently supplied from the surrounding rock. This gives distinct
# steps in the temperature profile at the layer interfaces.
colors_well = cgrad(:BrBG_4, 4, categorical = true)[[1, 2, 4]]
function plot_well_profiles!(ax, case, results)
    well = case.model.models[:CoaxialWell_supply].data_domain
    T_well = convert_from_si.(
        results.result.states[end][:CoaxialWell_supply][:Temperature], :Celsius)
    for (j, tag) in enumerate(unique(well[:tag]))
        name = titlecase(replace(string(tag), "_" => " "))
        cells = well[:tag] .== tag
        lines!(ax, T_well[cells], well[:cell_centroids][3, cells];
            color = colors_well[j], linewidth = 4, label = name)
    end
    res_cells = well.representation.perforations.reservoir
    T_res = convert_from_si.(
        results.result.states[end][:Reservoir][:Temperature], :Celsius)
    z_res = reservoir_model(case.model).data_domain[:cell_centroids][3, res_cells]
    lines!(ax, T_res[res_cells], z_res; color = :black, linewidth = 2, label = "Reservoir")
end

fig = Figure(size = (600, 600))
ax = Axis(fig[1, 1]; title = "Well temperature after 5 years",
    xlabel = "Temperature (°C)", ylabel = "Depth (m)", yreversed = true)
plot_well_profiles!(ax, case_ref, results_ref)
axislegend(ax; position = :rt)
fig

# ## Scenario 2 – Effect of inner pipe wall thermal conductivity
# Using the layered reservoir, we vary the inner pipe wall thermal
# conductivity from 0.0 (perfectly insulating inner pipe) to 4× the default
# value of 0.38 W/(m·K). A higher conductivity increases heat transfer between
# the supply and return pipes, causing more thermal short-circuiting and a
# lower production temperature.
λ_pipe_default = 0.38  # default inner pipe wall thermal conductivity [W/(m·K)]
λ_pipe = [0.0, 1.0, 4.0]*λ_pipe_default
labels_pipe = ["λ = 0.00", "λ = 0.38 (default)", "λ = 1.52"]

cases_pipe = []
results_pipe = []
for λ in λ_pipe
    case_i = coaxial_bhe(; base_args...,
        depths = depths, rock_thermal_conductivity = [λ_layers; λ_layers[end]],
        well_args = (
            inner_pipe_thermal_conductivity = λ,
            ## Keep the outer pipe conductivity constant
            outer_pipe_thermal_conductivity = λ_pipe_default,
        ))
    push!(cases_pipe, case_i)
    push!(results_pipe, run_case(case_i))
end

# ### Well temperature profiles for different inner pipe conductivities
# Each panel shows the temperature–depth profile for a given inner pipe
# conductivity value. Higher conductivity leads to more heat leakage from
# the warm produced water in the inner pipe to the cold injected water in the
# annulus, reducing the net temperature gain at the surface.
fig = Figure(size = (600, 800))
for (i, (case, results, label)) in enumerate(zip(cases_pipe, results_pipe, labels_pipe))
    ax = Axis(fig[i, 1]; title = label,
        xlabel = "Temperature (°C)", ylabel = "Depth (m)", yreversed = true)
    plot_well_profiles!(ax, case, results)
    i < length(λ_pipe) && hidexdecorations!(ax, grid = false)
    i == 1 && axislegend(ax; position = :rt)
end
linkaxes!(filter(c -> c isa Axis, fig.content)...)
fig

# ### Power output over time
# Compare the thermal power output for each inner pipe conductivity. Power is
# computed as mass flow rate × heat capacity × (production temperature −
# injection temperature). The power output drops considerably with increasing
# inner pipe conductivity.
function thermal_power(results; c_p = 4184.0joule/(kilogram*Kelvin))
    w = results.wells[:CoaxialWell_supply]
    return abs.(w[:mass_rate]) .* c_p .* (w[:temperature] .- base_args.temperature_inj)
end

fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1]; title = "Effect of inner pipe thermal conductivity",
    xlabel = "Time (years)", ylabel = "Power (MW)")
colors_pipe = Makie.wong_colors()
for (i, (results, label)) in enumerate(zip(results_pipe, labels_pipe))
    lines!(ax, results.time./year, thermal_power(results)./1e6;
        color = colors_pipe[i], linewidth = 4, label = label)
end
axislegend(ax; position = :rt)
fig

# ## Scenario 3 – Calibration of a digital twin
# In practice, the geology around the well is only partially known. We now use
# the layered model from Scenario 1 as a stand-in for the real system (the
# *reference*), and calibrate a digital twin (DT) of the system against
# observations from the reference. The DT has the same mesh as the reference,
# but homogeneous geology. As a prior guess for the thermal conductivity, we
# use 3.5 W/(m·K), which is a typical value for sedimentary rocks.
λ_prior = 3.5λ_unit
case_prior = coaxial_bhe(; base_args...,
    depths = depths, rock_thermal_conductivity = λ_prior);
results_prior = run_case(case_prior);

# ### Utilities for comparison
# The production temperature is the temperature at the top node of the
# production side of the well. The downhole temperatures are the grout
# temperatures along the well, as measured by a fibre-optic cable for
# distributed temperature sensing (DTS) installed behind the casing. We sort
# them by depth.
production_temperature(results) = convert_from_si.(
    results.wells[:CoaxialWell_supply][:temperature], :Celsius)

function downhole_temperature(case, results)
    well = case.model.models[:CoaxialWell_supply].data_domain
    cells = findall(well[:tag] .== :grout)
    z = well[:cell_centroids][3, cells]
    order = sortperm(z)
    T = results.result.states[end][:CoaxialWell_supply][:Temperature][cells]
    return (z = z[order], T = convert_from_si.(T[order], :Celsius))
end

# In the plots, the reference is shown as a thick, dashed black line, drawn
# on top of the DT results, which are shown as solid lines with a fixed color
# for each DT.
reference_style = (color = :black, linewidth = 6, linestyle = :dash)
dt_colors = Dict(zip(
    ["DT (prior)", "DT (producer data)", "DT (downhole temperatures)"],
    Makie.wong_colors()))
line_style(label) = label == "Reference" ? reference_style :
    (color = dt_colors[label], linewidth = 3, linestyle = :solid)

function plot_comparison(cases, results, labels; title = "")
    fig = Figure(size = (1200, 500))
    ax_t = Axis(fig[1, 1]; title = "Production temperature",
        xlabel = "Time (years)", ylabel = "Temperature (°C)")
    ax_z = Axis(fig[1, 2]; title = "Downhole temperature after 5 years",
        xlabel = "Temperature (°C)", ylabel = "Depth (m)", yreversed = true)
    for (c, r, l) in zip(cases, results, labels)
        z_order = l == "Reference" ? 1.0 : 0.0
        plt = lines!(ax_t, r.time./year, production_temperature(r);
            line_style(l)..., label = l)
        translate!(plt, 0, 0, z_order)
        prof = downhole_temperature(c, r)
        plt = lines!(ax_z, prof.T, prof.z; line_style(l)..., label = l)
        translate!(plt, 0, 0, z_order)
    end
    axislegend(ax_t; position = :rt)
    Label(fig[0, :], title; fontsize = 20, font = :bold)
    return fig
end

# ### Compare reference and DT before calibration
# The DT with the prior conductivity overestimates the production temperature
# of the reference. The downhole temperatures also show that the DT cannot
# represent the effect of the layers, which gives distinct steps in the
# reference temperature profile at the layer interfaces.
plot_comparison([case_ref, case_prior], [results_ref, results_prior],
    ["Reference", "DT (prior)"]; title = "Reference vs. DT before calibration")

# ### Calibration 1: Producer data
# The first calibration uses only the production temperature. This is the
# data that is most easily available in practice. With only a single
# measurement location, we cannot expect to resolve the layering, and we
# therefore calibrate a single (effective) thermal conductivity for the entire
# reservoir by lumping all cells into one group.
#
# We use the dict-based optimization framework in JutulDarcy to calibrate the
# DT. The function `setup_reservoir_dict_optimization` extracts the parameters
# of the DT into a dictionary, together with a function that sets up a new case
# from the (modified) dictionary. The objective function is the squared
# mismatch in production temperature, integrated over time, which is computed
# by `well_mismatch_thermal`.
states_ref = results_ref.result.states
time_scale = sum(case_prior.dt)
objective_producer = (model, state, dt, step_info, forces) ->
    well_mismatch_thermal(model, [:CoaxialWell_supply],
        states_ref, state, dt, step_info, forces;
        w_bhp = 0.0, w_energy = 0.0, w_temp = 1.0/Kelvin, scale = time_scale)

num_cells = number_of_cells(rdomain_ref)
dprm_producer = setup_reservoir_dict_optimization(case_prior)
free_optimization_parameter!(dprm_producer, [:model, :rock_thermal_conductivity];
    abs_min = 0.5λ_unit, abs_max = 6.0λ_unit, lumping = ones(Int, num_cells))
sim, cfg = setup_simulator(case_prior; output_substates = true);
prm_producer = optimize_reservoir(dprm_producer, objective_producer;
    simulator = sim, config = cfg, max_it = 25)
case_producer = dprm_producer.setup_function(prm_producer)
λ_producer = first(prm_producer[:model][:rock_thermal_conductivity])
println("Calibrated effective conductivity: $(round(λ_producer/λ_unit, digits = 2)) W/(m·K)")
results_producer = run_case(case_producer);

# #### Compare reference and DT after calibration to producer data
# The production temperature of the calibrated DT matches the reference well.
# The calibrated conductivity is close to, but slightly higher than, the
# thickness-weighted mean of the reference along the well (about 2.6
# W/(m·K)): Most of the heat is extracted in the deeper part of the well,
# where the temperature difference between the rock and the fluid is largest,
# and this part includes the layer with the highest conductivity. The downhole
# temperatures still cannot represent the variation due to the layering.
plot_comparison(
    [case_ref, case_prior, case_producer],
    [results_ref, results_prior, results_producer],
    ["Reference", "DT (prior)", "DT (producer data)"];
    title = "Calibration to producer data")

# ### Calibration 2: Downhole temperatures
# Distributed temperature sensing gives temperature measurements along the
# entire well, which can be used to resolve the vertical variation of the
# thermal conductivity. We assume that the depths of the geological layers are
# known, e.g., from geological logs, but not their thermal conductivities. We
# therefore calibrate one thermal conductivity per layer, by lumping all cells
# in the same layer into one group. As for the calibration to producer data,
# the DT with the prior conductivity is used as initial guess.
z = x[3, :]
num_layers = length(depths_layers) - 1
layers = [clamp(searchsortedlast(depths_layers, zi), 1, num_layers) for zi in z]

# The DTS measurements are the grout temperatures of the reference at the end
# of the simulation. Since the DT has the same mesh as the reference, the well
# cells coincide and the measurements can be compared directly. The objective
# function is the mean squared mismatch between the measurements and the
# grout temperatures of the DT, evaluated at the final time step only. Note
# that the mismatch is computed at all time steps and multiplied by zero for
# all but the final step. This ensures that the dependence of the objective
# function on the state is detected when the sparsity of the objective
# function is determined, which is done using the first time step.
well = case_ref.model.models[:CoaxialWell_supply].data_domain
grout = findall(well[:tag] .== :grout)
T_dts = states_ref[end][:CoaxialWell_supply][:Temperature][grout]

function objective_dts(model, state, dt, step_info, forces)
    is_final = step_info[:substep_global] == step_info[:Nsubstep_global]
    T = state[:CoaxialWell_supply][:Temperature]
    mismatch = 0.0
    for (i, c) in enumerate(grout)
        mismatch += ((T[c] - T_dts[i])/Kelvin)^2
    end
    return (is_final ? 1.0 : 0.0)*mismatch/length(grout)
end

dprm_dts = setup_reservoir_dict_optimization(case_prior)
free_optimization_parameter!(dprm_dts, [:model, :rock_thermal_conductivity];
    abs_min = 0.5λ_unit, abs_max = 6.0λ_unit, lumping = layers)
sim, cfg = setup_simulator(case_prior; output_substates = true);
prm_dts = optimize_reservoir(dprm_dts, objective_dts;
    simulator = sim, config = cfg, max_it = 50)
case_dts = dprm_dts.setup_function(prm_dts)
results_dts = run_case(case_dts);

# #### Comparison reference vs DT after calibration to downhole temperatures
# The DT calibrated to the downhole temperatures captures the effect of the
# layering on the temperature profile. The production temperature, which was
# not part of the objective function, is also close to the reference.
plot_comparison(
    [case_ref, case_prior, case_dts],
    [results_ref, results_prior, results_dts],
    ["Reference", "DT (prior)", "DT (downhole temperatures)"];
    title = "Calibration to downhole temperatures")

# ### Calibrated thermal conductivities
# Finally, we compare the calibrated thermal conductivities with the layered
# conductivity of the reference. Calibration to the production temperature
# gives a single effective conductivity, whereas calibration to the downhole
# temperatures recovers the conductivity of each layer. The dotted line marks
# the bottom of the well.
λ_dts = [prm_dts[:model][:rock_thermal_conductivity][findfirst(==(layer), layers)]
    for layer in 1:num_layers]
println("Calibrated layer conductivities: $(round.(λ_dts./λ_unit, digits = 2)) W/(m·K)")

function stairs_depth(depths, values)
    z = vcat([[depths[i], depths[i+1]] for i in eachindex(values)]...)
    v = vcat([[v, v] for v in values]...)
    return v, z
end

fig = Figure(size = (600, 700))
ax = Axis(fig[1, 1]; title = "Thermal conductivity",
    xlabel = "Thermal conductivity (W/(m·K))", ylabel = "Depth (m)",
    yreversed = true)
v, zv = stairs_depth(depths_layers, λ_layers./λ_unit)
plt = lines!(ax, v, zv; reference_style..., label = "Reference")
translate!(plt, 0, 0, 1)
for (label, λ) in (("DT (prior)", λ_prior), ("DT (producer data)", λ_producer))
    lines!(ax, [λ, λ]./λ_unit, [0.0, 3000.0]; line_style(label)..., label = label)
end
v, zv = stairs_depth(depths_layers, λ_dts./λ_unit)
lines!(ax, v, zv; line_style("DT (downhole temperatures)")...,
    label = "DT (downhole temperatures)")
hlines!(ax, [2500.0]; color = :black, linestyle = :dot)
axislegend(ax; position = :lt)
fig
