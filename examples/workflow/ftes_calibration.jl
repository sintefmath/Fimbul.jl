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
	# linsolve = GenericKrylov(
	# 	preconditioner = ILUZeroPreconditioner(),
	# 	rtol = 1e-3,
	# )
	simulator, config = setup_reservoir_simulator(case;
		output_substates = true,
		relaxation = true,
		initial_dt = 5.0,
		info_level = 0,
	)
	for well in well_symbols(case.model)
		config[:tolerances][well][:mass_conservation] = (
			CNV = Inf, MB = 1e-5,
			increment_dp_abs = 1e-3*si"atm", increment_dp_rel = Inf,
			increment_dz = Inf, increment_saturation = Inf)
		config[:tolerances][well][:energy_conservation] = (
			CNV = Inf, EB = 1e-5, increment_dT = 1e-2)
	end
	# config[:tolerances][:Fractures][:mass_conservation] = (
	# 	CNV = Inf, MB = 1.0e-5,
	# 	increment_dp_abs = 1e-3*si"atm", increment_dp_rel = Inf,
	# 	increment_dz = Inf, increment_saturation = Inf)
	# config[:tolerances][:Fractures][:energy_conservation] = (
	# 	CNV = Inf, EB = 1e-5, increment_dT = 1e-2)
	push!(config[:timestep_selectors], JutulDarcy.ControlChangeTimestepSelector(case.model, 0.0, 5.0))
	push!(config[:timestep_selectors], VariableChangeTimestepSelector(:Temperature, 15.0;
		model = :Reservoir,
		relative = false,
	))
	push!(config[:timestep_selectors], VariableChangeTimestepSelector(:Temperature, 15.0;
		model = :Fractures,
		relative = false,
	))
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
	producer_rate = -well_results[:Producer, :wrat]/(liter/second)
	injector_bhp = convert_from_si.(well_results[:Injector, :bhp], :bar)
	return (time = t_days, producer_temp = producer_temp, producer_rate = producer_rate, injector_bhp = injector_bhp)
end

function plot_match(reference, initial, optimized)
	fig = Figure(size = (1100, 800))

	ax1 = Axis(fig[1, 1],
		title = "Producer temperature",
		xlabel = "Time (days)",
		ylabel = "Temperature (C)",
	)
	scatter!(ax1, reference.time, reference.producer_temp, label = "Reference")
	lines!(ax1, initial.time, initial.producer_temp, label = "Initial")
	lines!(ax1, optimized.time, optimized.producer_temp, label = "Optimized")
	axislegend(ax1, position = :rb)

	ax2 = Axis(fig[2, 1],
		title = "Producer rate",
		xlabel = "Time (days)",
		ylabel = "Rate (L/s)",
	)
	scatter!(ax2, reference.time, reference.producer_rate, label = "Reference")
	lines!(ax2, initial.time, initial.producer_rate, label = "Initial")
	lines!(ax2, optimized.time, optimized.producer_rate, label = "Optimized")
	axislegend(ax2, position = :rb)

	ax3 = Axis(fig[3, 1],
		title = "Injector bottom-hole pressure",
		xlabel = "Time (days)",
		ylabel = "Pressure (bar)",
	)
	scatter!(ax3, reference.time, reference.injector_bhp, label = "Reference")
	lines!(ax3, initial.time, initial.injector_bhp, label = "Initial")
	lines!(ax3, optimized.time, optimized.injector_bhp, label = "Optimized")
	axislegend(ax3, position = :rb)

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

controls = Fimbul.ftes_controls(
	rate_charge = 12.0liter/second,
	temperature_charge = convert_to_si(90.0, :Celsius),
	temperature_discharge = convert_to_si(20.0, :Celsius),
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
	porosity = (0.5, 0.2),
	# porosity = 0.5,
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
using Statistics

# names = [:Transmissibilities, :RockThermalConductivities, :FluidThermalConductivities]
# for name in names
# 	values = case_real.parameters[:Fractures][name]
# 	fix = values .== 0.0
# 	if any(fix)
# 		mean_value = mean(values[.!fix])
# 		values[fix] .= mean_value
# 		@info "Replaced $(sum(fix)) zero values in $name with mean value $mean_value"
# 	end
# end


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
	porosity = 0.3,
)

discretization_idealized = Fimbul.ftes_discretization(
	wells,
	fractures_idealized;
	info_level = 1,
    mesh_args...
)

parameters_initial = Fimbul.ftes_parameters(
	# matrix_properties = (
	# 	permeability = 1.5e-4*darcy,
	# 	porosity = 0.12,
	# ),
	fracture_properties = (
		aperture = fractures_idealized[:aperture],
		porosity = fractures_idealized[:porosity],
		# permeability = fill(3.0e-4*darcy, 10),
		# porosity = fill(0.2, 10),
	),
)

case_initial = Fimbul.ftes(discretization_idealized, parameters_initial, controls; info_level = 1)
result_initial = run_case(case_initial; info_level = 2)

# ## Step 3: Calibrate the idealized case
truth_well_data = result_real.wells
prod_temp_ref = get_1d_interpolator(truth_well_data.time, truth_well_data[:Producer, :temperature])
prod_rate_ref = get_1d_interpolator(truth_well_data.time, truth_well_data[:Producer, :wrat])
inj_bhp_ref = get_1d_interpolator(truth_well_data.time, truth_well_data[:Injector, :bhp])
total_time = truth_well_data.time[end]

function setup_idealized_case(prm, step_info = missing)
	return Fimbul.ftes(discretization_idealized, prm, controls)
end

import JutulDarcy: compute_well_qoi
function mismatch_objective(model, state, dt, step_info, forces)
	time = step_info[:time]

	producer_temp = compute_well_qoi(model, state, forces, :Producer, :temperature)
	producer_rate = compute_well_qoi(model, state, forces, :Producer, :wrat)
	injector_bhp = compute_well_qoi(model, state, forces, :Injector, :bhp)

	dtemp = (prod_temp_ref(time) - producer_temp)/convert_to_si(1.0, :Kelvin)
	# drate = (prod_rate_ref(time) - producer_rate)/(liter/second)
	# dbhp = (inj_bhp_ref(time) - injector_bhp)/bar

	return dt*(dtemp^2 + 0.05*drate^2 + 0.1*dbhp^2)/total_time
end

opt = JutulDarcy.setup_reservoir_dict_optimization(parameters_initial, setup_idealized_case)
# free_optimization_parameter!(opt, [:reservoir, :matrix, :properties, :permeability],
# 	abs_min = 1.0e-5*darcy,
# 	abs_max = 5.0e-4*darcy,
# )
# free_optimization_parameter!(opt, [:reservoir, :matrix, :properties, :porosity],
# 	abs_min = 0.03,
# 	abs_max = 0.2,
# )
free_optimization_parameter!(opt, [:reservoir, :fractures, :properties, :aperture],
	abs_min = 1.0e-4,
	abs_max = 1.0e-3,
)
# free_optimization_parameter!(opt, [:reservoir, :fractures, :properties, :permeability],
# 	abs_min = 5.0e-5*darcy,
# 	abs_max = 2.0e-3*darcy,
# )
free_optimization_parameter!(opt, [:reservoir, :fractures, :properties, :porosity],
	abs_min = 0.05,
	abs_max = 0.6,
)

parameters_optimized = JutulDarcy.optimize_reservoir(
	opt,
	mismatch_objective;
	deps = :case,
	max_it = 12,
	optimizer = :lbfgsb_qp,
)

case_optimized = Fimbul.ftes(discretization_idealized, parameters_optimized, controls; info_level = 1)
result_optimized = run_case(case_optimized)

observed_real = get_well_observables(result_real.well_results)
observed_initial = get_well_observables(result_initial.well_results)
observed_optimized = get_well_observables(result_optimized.well_results)

match_figure = plot_match(observed_real, observed_initial, observed_optimized)
display(match_figure)

println("Initial parameters:")
println(parameters_initial)
println()
println("Optimized parameters:")
println(parameters_optimized)

##
Tf = case_real.parameters[:Fractures][:Transmissibilities]
is_zero = Tf .== 0.0
fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1])
plot_mesh!(ax, fracture_mesh, outer=false)
plot_mesh_edges!(ax, fracture_mesh, outer=false, alpha = 0.1)
fig
geo = tpfv_geometry(fracture_mesh)
xf = geo.face_centroids[:, is_zero]
scatter!(ax, xf[1,:], xf[2,:], xf[3,:], markersize = 10)
fig