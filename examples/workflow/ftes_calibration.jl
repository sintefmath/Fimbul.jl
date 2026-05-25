# # Calibrating an idealized FTES model
# This example demonstrates a simple history-matching workflow for Fractured
# Thermal Energy Storage (FTES) with Fimbul. We first generate synthetic data
# from a heterogeneous reference model and then calibrate a simpler idealized
# model against the producer temperature response.
#
# The reference model contains gently dipping fractures with varying strike and
# aperture, while the calibration model uses a perfectly horizontal fracture
# set. We calibrate the matrix rock thermal conductivity together with the
# aperture of each idealized fracture.

using Dates
using Random

using GLMakie
using HYPRE
using Jutul, JutulDarcy, Fimbul

# ## Helper functions
# We collect a few helper functions up front to keep the workflow sections
# focused on the calibration procedure itself.
function setup_ftes_simulator(case; kwargs...)
    simulator, config = setup_reservoir_simulator(case;
        relaxation = true,
        initial_dt = 5.0,
        info_level = 0,
        kwargs...,
    )

    control_selector = JutulDarcy.ControlChangeTimestepSelector(
        case.model, 0.0, 60.0)
    push!(config[:timestep_selectors], control_selector)

    reservoir_temperature_selector = VariableChangeTimestepSelector(
        :Temperature, 15.0; model = :Reservoir, relative = false)
    push!(config[:timestep_selectors], reservoir_temperature_selector)

    fracture_temperature_selector = VariableChangeTimestepSelector(
        :Temperature, 15.0; model = :Fractures, relative = false,
    )
    push!(config[:timestep_selectors], fracture_temperature_selector)

    config[:timestep_max_decrease] = 1e-6
    return simulator, config
end

function get_well_observables(well_results)
    t_days = convert_from_si.(well_results.time, :day)
    producer_temp = convert_from_si.(well_results[:Producer, :temperature], :Celsius)
    return (time = t_days, producer_temp = producer_temp)
end

function plot_match!(ax, reference, simulated::AbstractVector)
    colors = Makie.wong_colors(6)[[1, 2]]

    lines!(ax, reference.time, reference.producer_temp;
        label = "Reference", linewidth = 5, linestyle = :dash, color = :black)
    lines!(ax, simulated[1].time, simulated[1].producer_temp;
        label = "Initial", linewidth = 2, color = colors[1])
    if length(simulated) > 2
        for observable in simulated[2:end-1]
            lines!(ax, observable.time, observable.producer_temp;
                label = "Optimization iterate", linewidth = 2, color = :gray)
        end
    end
    if length(simulated) > 1
        lines!(ax, simulated[end].time, simulated[end].producer_temp;
            label = "Final", linewidth = 2, color = colors[2])
    end

    axislegend(ax, position = :rt, merge = true)
    return ax
end

function plot_ftes_wells!(ax, case)
    colors = Makie.wong_colors(6)[[2, 6]]
    for (i, xw) in enumerate(case.input_data[:well_coordinates])
        color = ifelse(i == 1, colors[1], colors[2])
        lines!(ax, xw[1, :], xw[2, :], xw[3, :]; color = color, linewidth = 3)
    end
    return ax
end

# ## Shared FTES configuration
# The reference and idealized cases share the same well pattern, computational
# domain, and operating schedule. This makes it easier to isolate the impact of
# structural simplifications in the fracture description.
wells = (num_producers = 8, radius = 25.0, depth = 220.0)
depth_window = (z_min = 50.0, z_max = 190.0)
mesh_args = (
    hxy_min = 12.0,
    hxy_max = 60.0,
    offset = 140.0,
    offset_rel = missing,
)

T_charge = convert_to_si(90.0, :Celsius)
T_discharge = convert_to_si(20.0, :Celsius)
controls = Fimbul.ftes_controls(
    rate_charge = 12.0liter / second,
    temperature_charge = T_charge,
    temperature_discharge = T_discharge,
    producer_bhp_fraction = 0.2,
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = (num_years = 1,),
)

# ## Build a synthetic reference case
# The reference model uses a heterogeneous fracture description with varying
# strike, dip, and aperture. We simulate this model to create synthetic
# production-temperature observations for calibration.
Random.seed!(2026)
fractures_reference = Fimbul.setup_ftes_fractures(
    10,
    depth_window.z_min,
    depth_window.z_max;
    uniform_spacing = true,
    strike = (0.0, 20.0),
    dip = (0.0, 2.5),
    radius = 55.0,
    aperture = (0.5e-3, 1.5e-4),
    porosity = 0.5,
)

discretization_reference = Fimbul.ftes_discretization(
    wells,
    fractures_reference;
    info_level = 1,
    mesh_args...,
)

parameters_reference = Fimbul.ftes_parameters(
    fracture_properties = (
        aperture = fractures_reference[:aperture],
        porosity = fractures_reference[:porosity],
    ),
)

case_reference = Fimbul.ftes(
    discretization_reference,
    parameters_reference,
    controls;
    info_level = 1,
)

# ## Inspect the reference geometry
# We begin by visualizing the matrix mesh, embedded fracture network, and well
# layout used to generate the synthetic observations.
matrix_mesh = physical_representation(reservoir_model(case_reference.model).data_domain)
fracture_mesh = physical_representation(case_reference.model.models[:Fractures].data_domain)

axis_args = (
    perspectiveness = 0.75,
    zreversed = true,
    aspect = :data,
    elevation = 0.025π,
    azimuth = 1.35π,
)
fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1]; axis_args..., title = "FTES reference system")
Jutul.plot_mesh!(ax, fracture_mesh; color = :gray)
Jutul.plot_mesh_edges!(ax, matrix_mesh; alpha = 0.1)
plot_ftes_wells!(ax, case_reference)
fig

# ## Simulate the reference case
# We use a conservative timestep strategy around control changes and enable
# `output_substates` so that reservoir states can be visualized later.
sim, cfg = setup_ftes_simulator(case_reference; output_substates = true)
result_reference = simulate_reservoir(
    case_reference; simulator = sim, config = cfg)

# ## Build the idealized calibration model
# The calibration model keeps the same well layout and fracture count, but uses
# perfectly horizontal fractures with a uniform initial aperture. This creates a
# deliberately simpler proxy that must absorb the mismatch through calibration.
fractures_idealized = Fimbul.setup_ftes_fractures(
    10,
    depth_window.z_min,
    depth_window.z_max;
    uniform_spacing = true,
    strike = 0.0,
    dip = 0.0,
    radius = 55.0,
    aperture = 0.5e-3,
    porosity = 0.5,
)

discretization_idealized = Fimbul.ftes_discretization(
    wells,
    fractures_idealized;
    info_level = 1,
    mesh_args...,
)

parameters_initial = Fimbul.ftes_parameters(
    matrix_properties = (
        rock_thermal_conductivity = 4.0 * si"watt/(meter*Kelvin)",
    ),
    fracture_properties = (
        aperture = fractures_idealized[:aperture],
        porosity = fractures_idealized[:porosity],
    ),
)

case_initial = Fimbul.ftes(
    discretization_idealized, parameters_initial, controls; info_level = 1)
sim, cfg = setup_ftes_simulator(case_initial; output_substates = true)
result_initial = simulate_reservoir(
    case_initial; simulator = sim, config = cfg)

# ## Compare the initial mismatch
# Before calibration, the simplified model does not reproduce the producer
# temperature curve from the heterogeneous reference case.
reference_observables = get_well_observables(result_reference.wells)
initial_observables = get_well_observables(result_initial.wells)

fig = Figure(size = (900, 500))
ax = Axis(fig[1, 1];
    xlabel = "Time (days)",
    ylabel = "Producer temperature (°C)",
    title = "Initial mismatch",
)
plot_match!(ax, reference_observables, [initial_observables])
fig

# ## Define the calibration objective
# We match the producer temperature over the full schedule using a least-squares
# misfit integrated over time.
reference_well_data = result_reference.wells
producer_temperature_reference = Jutul.get_1d_interpolator(
    reference_well_data.time,
    reference_well_data[:Producer, :temperature],
)
total_time = reference_well_data.time[end]

function setup_idealized_case(parameters, step_info = missing)
    return Fimbul.ftes(discretization_idealized, parameters, controls)
end

import JutulDarcy: compute_well_qoi

function mismatch_objective(model, state, dt, step_info, forces)
    time = step_info[:time]
    producer_temperature = compute_well_qoi(model, state, forces, :Producer, :temperature)
    ΔT = (producer_temperature_reference(time) - producer_temperature) /
        convert_to_si(1.0, :Kelvin)
    return dt * ΔT^2 / total_time
end

# ## Run the calibration
# We use dict-based optimization and free one matrix parameter together with the
# aperture assigned to each idealized fracture. The optimization is intentionally
# short to keep the example lightweight; increase `max_it` for a tighter match.
opt = JutulDarcy.setup_reservoir_dict_optimization(
    parameters_initial,
    setup_idealized_case,
)

free_optimization_parameter!(
    opt,
    [:reservoir, :matrix, :rock_thermal_conductivity],
    abs_min = 1.0 * si"watt/(meter*Kelvin)",
    abs_max = 4.0 * si"watt/(meter*Kelvin)",
)
free_optimization_parameter!(
    opt,
    [:reservoir, :fractures, :aperture],
    abs_min = 1.0e-4 * si"meter",
    abs_max = 1.0e-3 * si"meter",
)

simulator, config = setup_ftes_simulator(case_initial)
parameters_optimized = JutulDarcy.optimize_reservoir(
    opt,
    mismatch_objective;
    deps = :case,
    max_it = 2,
    optimizer = :lbfgsb_qp,
    simulator = simulator,
    config = config,
    solution_history = true,
)

# ## Inspect optimization history
# The objective history provides a compact view of how much the mismatch is
# reduced during the calibration.


# ## Compare the calibrated response
# We now simulate the optimized idealized model and compare the resulting
# producer temperature curve against the initial proxy and the synthetic
# reference data.
simulated_observables = []
push!(simulated_observables, initial_observables)
# For illustraional purposes, we simulate the intermediate optimization iterates
# to show how the response evolves during calibration. Adjust the range of
# `simulate_int` to control how many intermediate iterates to recompute.
simulate_int = 2:length(opt.history.solutions)-1
simulate_int = setdiff(simulate_int, [1, length(opt.history.solutions)])
for h in opt.history.solutions[simulate_int]
	case = setup_idealized_case(h.parameters)
    sim, cfg = setup_ftes_simulator(case; output_substates = true)
	result = simulate_reservoir(case; simulator = sim, config = cfg)
	observables = get_well_observables(result.wells)
	push!(simulated_observables, observables)
end
case = setup_idealized_case(opt.history.solutions[end].parameters)
sim, cfg = setup_ftes_simulator(case; output_substates = true)
result = simulate_reservoir(case; simulator = sim, config = cfg)
push!(simulated_observables, get_well_observables(result.wells))

# ### 
fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1];
    xlabel = "Time (days)",
    title = "Producer temperature (°C)",
)
plot_match!(ax, reference_observables, simulated_observables)

ax = Axis(fig[1, 2];
    xlabel = "Optimization iteration",
    title = "Temperature mismatch objective (-)",
)

colors = Makie.wong_colors(6)[[1, 2]]
lines!(ax, opt.history.objectives; linewidth = 2, color = :black)
scatter!(ax, opt.history.objectives, markersize = 8, color = :gray)
scatter!(ax, 1, opt.history.objectives[1];
    markersize = 10, color = colors[1], label = "Initial")
scatter!(ax, length(opt.history.objectives), opt.history.objectives[end];
    markersize = 10, color = colors[2], label = "Final")
fig

# ## Inspect the calibrated parameters
# The optimized parameter dictionary can be printed directly for further
# inspection or reused in subsequent FTES studies.
println("Initial parameters:")
println(parameters_initial)
println()
println("Optimized parameters:")
println(parameters_optimized)
