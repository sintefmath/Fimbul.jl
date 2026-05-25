using Jutul, JutulDarcy, Fimbul
using HYPRE
using GLMakie
using Random

meter, day, year, darcy, bar, liter, second = si_units(
	:meter,
	:day,
	:year,
	:darcy,
	:bar,
	:liter,
	:second,
)

function merge_fracture_sets(fracture_sets...)
	merged = Dict{Symbol, Any}()
	all_keys = reduce(union, (collect(keys(fr)) for fr in fracture_sets))
	for key in all_keys
		values = Any[]
		for fractures in fracture_sets
			haskey(fractures, key) || continue
			append!(values, fractures[key])
		end
		merged[key] = values
	end
	return merged
end

function setup_ftes_simulator(case)
	simulator, config = setup_reservoir_simulator(case;
		output_substates = true,
		relaxation = true,
		initial_dt = 5.0,
		info_level = 0,
	)
	sel_cc = JutulDarcy.ControlChangeTimestepSelector(case.model, 0.0, 60.0*si_unit(:day))
	push!(config[:timestep_selectors], sel_cc)
	sel_vc = VariableChangeTimestepSelector(:Temperature, 15.0;
		model = :Reservoir,
		relative = false,
	)
	push!(config[:timestep_selectors], sel_vc)
	sel_vc_f = VariableChangeTimestepSelector(:Temperature, 15.0;
		model = :Fractures,
		relative = false,
	)
	push!(config[:timestep_selectors], sel_vc_f)
	config[:timestep_max_decrease] = 1e-6
	return simulator, config
end

function run_case(case; info_level = 1)
	simulator, config = setup_ftes_simulator(case)
	return simulate_reservoir(case;
		simulator = simulator,
		config = config,
		info_level = info_level,
	)
end

function get_well_observables(well_results)
	t_days = convert_from_si.(well_results.time, :day)
	producer_temp = convert_from_si.(well_results[:Producer, :temperature], :Celsius)
	return (time = t_days, producer_temp = producer_temp)
end

function plot_match!(ax, reference, simulated::Vector)

	colors = Makie.wong_colors(6)[[1,2]]

	lines!( # Plot reference
		ax, reference.time, reference.producer_temp;
		label = "Reference", linewidth = 6, linestyle=:dash, color = :black)
	lines!( # Plot initial guess
		ax, simulated[1].time, simulated[1].producer_temp;
		label = "Initial", linewidth = 2, color=colors[1])
	for k in 2:length(simulated)-1
		lines!( # Plot optimization iterations
			ax, simulated[k].time, simulated[k].producer_temp;
			label = "Opt. Iteration", linewidth = 2, color=:gray)
	end
	if length(simulated) > 1
		lines!( # Plot final optimized result
			ax, simulated[end].time, simulated[end].producer_temp;
			label = "Final", linewidth = 2, color=colors[2])
	end
	axislegend(ax, position = :rb, merge=true)

	return fig
end

# ## Shared FTES setup
wells = (num_producers = 8, radius = 25.0, depth = 220.0)
depth_window = (z_min = 50.0, z_max = 190.0)
mesh_args = (
	hxy_min = 12.0,
	hxy_max = 60.0,
	offset = 140.0,
	offset_rel = missing,
)

T_charge = convert_to_si(90, :Celsius)
T_discharge = convert_to_si(20, :Celsius)
controls = Fimbul.ftes_controls(
	rate_charge = 12.0liter/second,
	temperature_charge = T_charge,
	temperature_discharge = T_discharge,
	producer_bhp_fraction = 0.2,
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
	utes_schedule_args = (num_years = 1,),
)

# ## Step 1: "Real" FTES case
Random.seed!(2026)
fractures_horizontal = Fimbul.setup_ftes_fractures(
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

Δz = depth_window.z_max - depth_window.z_min
fractures_conductive = Fimbul.setup_ftes_fractures(
	1,
	depth_window.z_min + 0.3Δz,
	depth_window.z_max - 0.3Δz;
	strike = (0.0, 15.0),
	dip = (0.0, 2.5),
	radius = Inf,
	aperture = 0.5e-3,
	porosity = 0.4,
)

fractures_real = merge_fracture_sets(fractures_horizontal, fractures_conductive)
fractures_real = fractures_horizontal
real_fracture_permeability = vcat(fill(2.0e-4*darcy, 10), fill(1.0e-3*darcy, 5))

discretization_real = Fimbul.ftes_discretization(
	wells,
	fractures_real;
	info_level = 1,
    mesh_args...
)

parameters_real = Fimbul.ftes_parameters(
	# matrix_properties = (
	# 	permeability = 5.0e-5*darcy,
	# 	porosity = 0.08,
	# ),
	fracture_properties = (
		aperture = fractures_real[:aperture],
		porosity = fractures_real[:porosity],
	),
)

case_real = Fimbul.ftes(discretization_real, parameters_real, controls; info_level = 1)

##
matrix_mesh = physical_representation(reservoir_model(case_real.model).data_domain)
fracture_mesh = physical_representation(case_real.model.models[:Fractures].data_domain)

axis_args = (perspectiveness = 0.75, zreversed = true, aspect = :data,
    elevation = 0.025π, azimuth = 1.35π)
fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1]; axis_args...,
    title = "FTES system: matrix mesh and fracture network")
Jutul.plot_mesh!(ax, fracture_mesh; color = :gray)
Jutul.plot_mesh_edges!(ax, matrix_mesh; alpha = 0.1)
colors = Makie.wong_colors(6)[[2,6]]

function plot_ftes_wells(ax)
    for (i, xw) in enumerate(case_real.input_data[:well_coordinates])
        color = ifelse(i == 1, colors[1], colors[2])
        lines!(ax, xw[1,:], xw[2,:],  xw[3,:],
            color = color, linewidth = 3)
    end
end
plot_ftes_wells(ax)
fig

##
result_real = run_case(case_real; info_level = 2)

##
msh = physical_representation(reservoir_model(case_real.model).data_domain)
geo = tpfv_geometry(msh)
x_range = diff(vcat(extrema(geo.cell_centroids[1, :])...))[1]
y_range = diff(vcat(extrema(geo.cell_centroids[2, :])...))[1]
z_range = diff(vcat(extrema(geo.cell_centroids[3, :])...))[1]
aspect  = (x_range, y_range, z_range) ./ max.(x_range, y_range, z_range)

##
using Dates
timestamps = case_real.input_data[:timestamps][2:end]

steps = findall([Dates.monthname(t) ∈ ["December", "April"] .&&
    Dates.day(t) == 1 for t in timestamps])
steps = steps[[1, 2, end-1, end]] # Select first two and last two cycles for better visualization
cells = .!(geo.cell_centroids[1,:] .< 0.0 .&& geo.cell_centroids[2,:] .< 0.0)
colorrange = convert_from_si.((T_discharge, T_charge), :Celsius)
fig = Figure(size = (800, 900))
for (k, step) in enumerate(steps)
    row = (k-1)÷2 + 1
    col = (k-1)%2 + 1
    month = Dates.monthname(case_real.input_data[:timestamps][step])
    year = Dates.year(case_real.input_data[:timestamps][step])
    ax = Axis3(fig[row, col];
        title = "$month $year",
        zreversed = true, aspect = aspect, axis_args...,
        azimuth = 1.25π, titlegap = -50)
    T = convert_from_si.(result_real.result.states[step][:Reservoir][:Temperature], :Celsius)
    plot_cell_data!(ax, msh, T;
        cells = cells, colormap = :seaborn_icefire_gradient, colorrange = colorrange)
    hidedecorations!(ax)
end
Colorbar(fig[3, 1:2];
    colormap = :seaborn_icefire_gradient, colorrange = colorrange,
    label = "Temperature (°C)", vertical = false, flipaxis = false)
colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)
fig


# ## Step 2: Idealized FTES case
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
    mesh_args...
)

parameters_initial = Fimbul.ftes_parameters(
	matrix_properties = (
		rock_thermal_conductivity = 4.0*si"watt/(meter*Kelvin)",
	),
	fracture_properties = (
		aperture = fractures_idealized[:aperture],
		porosity = fractures_idealized[:porosity],
	),
)

case_initial = Fimbul.ftes(discretization_idealized, parameters_initial, controls; info_level = 1)
result_initial = run_case(case_initial; info_level = 2)

##
reference_ws = get_well_observables(result_real.wells)
initial_ws = get_well_observables(result_initial.wells)
fig = Figure(size = (900, 700))
ax = Axis(fig[1, 1])
plot_match!(ax, reference_ws, [initial_ws])
fig

# ## Step 3: Calibrate the idealized case
reference_well_data = result_real.wells
prod_temp_ref = get_1d_interpolator(reference_well_data.time, reference_well_data[:Producer, :temperature])
total_time = reference_well_data.time[end]

function setup_idealized_case(prm, step_info = missing)
	return Fimbul.ftes(discretization_idealized, prm, controls)
end

import JutulDarcy: compute_well_qoi
function mismatch_objective(model, state, dt, step_info, forces)
	time = step_info[:time]

	producer_temp = compute_well_qoi(model, state, forces, :Producer, :temperature)
	dtemp = (prod_temp_ref(time) - producer_temp)/convert_to_si(1.0, :Kelvin)

	return dt*dtemp^2/total_time
end

opt = JutulDarcy.setup_reservoir_dict_optimization(parameters_initial, setup_idealized_case)

free_optimization_parameter!(opt, [:reservoir, :matrix, :rock_thermal_conductivity],
	abs_min = 1.0*si"watt/(meter*Kelvin)",
	abs_max = 4.0*si"watt/(meter*Kelvin)",
)
free_optimization_parameter!(opt, [:reservoir, :fractures, :aperture],
	abs_min = 1.0e-4*si"meter",
	abs_max = 1.0e-3*si"meter",
)

sim, cfg = setup_ftes_simulator(case_initial)
output_path = jutul_output_path("ftes-calibration", subfolder = joinpath("papers", "get-2026"))
parameters_optimized = JutulDarcy.optimize_reservoir(
	opt,
	mismatch_objective;
	deps = :case,
	max_it = 2,
	solution_history = true,
	optimizer = :lbfgsb_qp,
	simulator = sim,
	config = cfg,
	output_path = output_path,
)

##
history = opt.history[4]
simulated_opt = []
for h in history
	case = setup_idealized_case(h.parameters)
	result = run_case(case; info_level = 0)
	observables = get_well_observables(result.wells)
	push!(simulated_opt, observables)
end

##
case_optimized = Fimbul.ftes(discretization_idealized, parameters_optimized, controls; info_level = 1)
result_optimized = run_case(case_optimized)

##
fig = Figure(size = (900, 700))
ax = Axis(fig[1, 1])
plot_match!(ax, observed_real, simulated_opt)
ax = Axis(fig[1, 2], yscale = log10)
lines!(ax, opt.history.objectives; label = "Objective value", linewidth = 2)
fig

##


observed_real = get_well_observables(result_real.wells)
observed_initial = get_well_observables(result_initial.wells)
observed_optimized = get_well_observables(result_optimized.wells)


match_figure = plot_match(observed_real, simulated_opt)

lines(opt.history.objectives)

display(match_figure)

println("Initial parameters:")
println(parameters_initial)
println()
println("Optimized parameters:")
println(parameters_optimized)
