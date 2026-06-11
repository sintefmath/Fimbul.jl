using Fimbul, JutulDarcy, Jutul, CoolProp, HYPRE, GLMakie

##
const BTES_CONFIG = (
	num_years = 1,
    num_wells = 16,
	num_sectors = 2,
	geothermal_gradient = 0.03,
	temperature_charge = convert_to_si(200.0, :Celsius),
	temperature_discharge = convert_to_si(50.0, :Celsius),
	rate_charge = 0.5 * si_unit(:litre) / si_unit(:second),
	rate_discharge = 0.5 * si_unit(:litre) / si_unit(:second),
	temperature_surface = convert_to_si(10.0, :Celsius),
	depths = [0.0, 10.0, 250.0, 300.0],
	density = [2580, 2580, 2580]*si"kilogram/meter^3",
    thermal_conductivity = [3.7, 3.7, 3.7]*si"watt/meter/Kelvin",
    heat_capacity = [900, 900, 900]*si"joule/kilogram/Kelvin",
	n_z = [5, 25, 3],
	well_layers = [2],
    charge_period = ["April", "September"],
    discharge_period = ["October", "March"],
	report_interval = 14 * si"day",
	utes_schedule_args = NamedTuple(),
)

bhp = 20.0 * si_unit(:atm)
function add_enthalpy_to_injectors(controls, enthalpy_of_pt)
	converted = Dict{Symbol, Any}()
	if controls[:B1_supply].target isa TotalRateTarget #|| true
		temperature = BTES_CONFIG.temperature_charge
	elseif controls[:B16_supply].target isa TotalRateTarget
		temperature = BTES_CONFIG.temperature_discharge
	else
		error("Unexpected control targets")
	end
	for (name, ctrl) in pairs(controls)
		if ctrl isa InjectorControl
			target = TotalRateTarget(BTES_CONFIG.rate_charge)
			# temperature = BTES_CONFIG.temperature_charge
			# target = ctrl.target
			converted[name] = InjectorControl(
				target,
				ctrl.injection_mixture;
				density = ctrl.mixture_density,
				phases = ctrl.phases,
				temperature = temperature,
				enthalpy = enthalpy_of_pt(name, temperature),
				tracers = ctrl.tracers,
				factor = ctrl.factor,
				check = false,
			)
		elseif ctrl isa ProducerControl
			# converted[name] = ctrl
			# ctrl.target
			converted[name] = ProducerControl(BottomHolePressureTarget(bhp))
				# ctrl.target
				# ctrl = ProducerControl(BTES_CONFIG.rate_discharge)
		else
			converted[name] = ctrl
		end
	end
	return converted
end

function rebuild_btes_with_h2o(case0, tables; config = BTES_CONFIG)
	domain = reservoir_model(case0.model).data_domain
	wells = collect(values(get_model_wells(case0; data_domain = true)))
	# for well in wells
	# 	if contains(string(well.representation.name), "_supply")
	# 		tag = well[:tag]
	# 		N = well.representation.neighborship
	# 		ix = vec(tag[N][1,:] .== :pipe_left .&& tag[N][2,:] .== :grout_left)
	# 		well[:material_thermal_conductivity][findall(ix)[1:4]] .= 0.0
	# 		ix = vec(tag[N][1,:] .== :pipe_right .&& tag[N][2,:] .== :grout_right)
	# 		well[:material_thermal_conductivity][findall(ix)[end-3:end]] .= 0.0
	# 	end
	# 	# well[:thermal_well_index][1:min(5,n)] .= 0
	# end
	sys = Fimbul.H2OSystem(tables)
	model, parameters = setup_reservoir_model(domain, sys;
		wells = wells,
		thermal = true,
		extra_out = true,
        # dp_max_abs_well = 1e4,
		# block_backend = false,
	)

	rmodel = reservoir_model(model)
	append!(rmodel.output_variables, [
		:Saturations,
		:PhaseMassDensities,
		:PhaseViscosities,
		:FluidEnthalpy,
		:Temperature,
		:Enthalpy,
	])
	unique!(rmodel.output_variables)

	p0 = case0.state0[:Reservoir][:Pressure]
	T0 = case0.state0[:Reservoir][:Temperature]
	h0 = tables[:enthalpy].(p0, T0)
	wellhead_pressure = Dict(
		well => p0[first(model.models[well].data_domain.representation.perforations.reservoir)]
		for well in well_symbols(model)
	)

	state0 = setup_reservoir_state(model,
		Pressure = p0,
		Enthalpy = h0,
	)

	# if haskey(case0.state0, :Facility) && haskey(state0, :Facility)
	# 	state0[:Facility][:SurfaceTemperature] .= case0.state0[:Facility][:SurfaceTemperature]
	# end

    # p0 = 10.0si"atm"
	# for well in well_symbols(model)
	# 	if haskey(case0.state0, well)
	# 		well_state0 = case0.state0[well]
	# 		state0[well][:Pressure] .= well_state0[:Pressure] .+ p0
	# 			# state0[well][:Pressure] .= wellhead_pressure[well] .+ p0
	# 		if haskey(state0[well], :Temperature) && haskey(well_state0, :Temperature)
	# 			state0[well][:Temperature] .= well_state0[:Temperature]
	# 		end
	# 		if haskey(state0[well], :Enthalpy) && haskey(well_state0, :Pressure) && haskey(well_state0, :Temperature)
	# 			state0[well][:Enthalpy] .= tables[:enthalpy].(well_state0[:Pressure], well_state0[:Temperature])
	# 		end
	# 	end
	# end
	# state0[:Facility][:BottomHolePressure] = [state0[well][:Pressure][1] for well in well_symbols(model)]
	# state0[:Facility][:SurfaceEnthalpy] = [state0[well][:Enthalpy][1] for well in well_symbols(model)]

	geo = tpfv_geometry(physical_representation(domain))
	z_bc = geo.boundary_centroids[3, :]
	bottom = isapprox.(z_bc, maximum(z_bc))
	z_hat = z_bc[.!bottom] .- minimum(z_bc[.!bottom])
	bc_cells = geo.boundary_neighbors[.!bottom]
	rho_ref = reservoir_model(case0.model).system.rho_ref[1]
	# p_bc = minimum()
	# 5.0 * si_unit(:atm) .+ rho_ref * gravity_constant .* z_hat
	# T_bc = config.temperature_surface .+ config.geothermal_gradient .* z_hat
	bc_cells = Int[]
	T_bc = Float64[]
	p_bc = Float64[]
	for bc in case0.forces[1][:Reservoir].bc
		push!(bc_cells, bc.cell)
		push!(T_bc, bc.temperature)
		push!(p_bc, bc.pressure)
	end
	h_bc = tables[:enthalpy].(p_bc, T_bc)
	bc = flow_boundary_condition(bc_cells, domain, p_bc, T_bc; enthalpy = h_bc)

	control_charge, control_discharge, sectors = Fimbul.setup_controls(
		model,
		config.num_sectors,
		config.rate_charge,
		config.rate_discharge,
		config.temperature_charge,
		config.temperature_discharge,
	)
	enthalpy_of_pt = function (well, temperature)
		if !isfinite(temperature)
			# return NaN
			temperature = config.temperature_charge
		end
		pressure = wellhead_pressure[well]
		p = 55si"bar"
		h = tables[:enthalpy](p, temperature)
		h = (p, T) -> tables[:enthalpy](p, T)
		return h
	end
	control_charge_h2o = add_enthalpy_to_injectors(control_charge, enthalpy_of_pt)
	control_discharge_h2o = add_enthalpy_to_injectors(control_discharge, enthalpy_of_pt)

	forces_charge = setup_reservoir_forces(model; control = control_charge_h2o, bc = bc)
	forces_discharge = setup_reservoir_forces(model; control = control_discharge_h2o, bc = bc)
	forces_rest = setup_reservoir_forces(model; bc = bc)
	dt, forces, timestamps = make_utes_schedule(
		forces_charge,
		forces_discharge,
		forces_rest;
		charge_period = config.charge_period,
		discharge_period = config.discharge_period,
		num_years = config.num_years,
		report_interval = config.report_interval,
		config.utes_schedule_args...,
	)

	info = Dict{Symbol, Any}(
		:description => "BTES case rebuilt with Fimbul.H2OSystem",
		:sectors => sectors,
		:timestamps => timestamps,
	)

	return JutulCase(model, dt, forces;
		state0 = state0,
		parameters = parameters,
		input_data = info,
	)
end

##
case0 = Fimbul.btes(;BTES_CONFIG...)

##
using JLD2
fn = joinpath(dirname(pathof(Fimbul)), "..", "coolprop_tables.jld2")
tables = load(fn)["data"]
case = rebuild_btes_with_h2o(case0, tables; config = BTES_CONFIG)

##
sim, cfg = setup_reservoir_simulator(case;
	tol_cnv = 1e-2,
	tol_mb = 1e-6,
	info_level = 2,
	# tol_cnve_well = Inf,
	tol_dp_well = 1e-3,
	# tolerances = (inc_tol_dh_rel = 1e-1,),
	output_substates = true,
    max_timestep_cuts=50,
    relaxation = true,
);
sel = JutulDarcy.ControlChangeTimestepSelector(
    case.model, 0.0, convert_to_si(5.0, :second))
push!(cfg[:timestep_selectors], sel)
cfg[:timestep_max_decrease] = 1e-6;
sel = VariableChangeTimestepSelector(:Temperature, 10.0; 
    relative = false, model = :Reservoir)
push!(cfg[:timestep_selectors], sel);
sel = VariableChangeTimestepSelector(:Temperature, 10.0; 
    relative = false, model = :B1_supply)
push!(cfg[:timestep_selectors], sel);
cfg
# for (w, cfg_w) in cfg[:tolerances]
# 	w ∈ well_symbols(case.model) || continue
# 	cfg_w[:energy_conservation] = (CNV = Inf, EB = 1.0e-2,
# 	increment_dT = Inf, increment_dh_abs = Inf, increment_dh_rel = 0.1)
# 	# cfg[:tolerances][w] = cfg_w
# end

##
results = simulate_reservoir(case; simulator = sim, config = cfg)

##
well = :B1_supply
pipe_left = case.model.models[well].data_domain[:tag] .== :pipe_left
pipe_right = case.model.models[well].data_domain[:tag] .== :pipe_right
grout_left = case.model.models[well].data_domain[:tag] .== :grout_left
grout_right = case.model.models[well].data_domain[:tag] .== :grout_right
colors = cgrad(:Paired_10, 10, categorical=true)[3:end]
time = cumsum(case.dt)./si_unit(:day)

function plot_btes_temperature_profiles(results, steps; legendposition = :cb, legend_axno = 1)
	fig = Figure(size = (1200, 500))

	cells = case.model.models[:B1_supply].data_domain.representation.perforations.reservoir
	z = case.model.models[:B1_supply].data_domain[:cell_centroids][3, :]

	for (sno, step) in enumerate(steps)
		
		ax = Axis(fig[1, sno], title = "Time = $(Int(round(time[step], digits = 0))) days", xlabel = "Temperature (°C)", ylabel = "Depth (m)", yreversed = true, aspect = AxisAspect(1))
		state_w = results.result.states[step][:B1_supply]
		for (sno, (section, label)) in enumerate(zip((pipe_left, pipe_right, grout_left, grout_right), ["Pipe left", "Pipe right", "Grout left", "Grout right"]))
				ix = section
				T = convert_from_si.(state_w[:Temperature][ix], :Celsius)
				lines!(ax, T, z[ix]; color=colors[sno], linewidth = 4, label = label)
		end

		state_r = results.result.states[step][:Reservoir]
		T_reservoir = convert_from_si.(state_r[:Temperature][cells], :Celsius)
		lines!(ax, T_reservoir, z[grout_right .|| grout_left]; color=:black, linewidth=4, label="Reservoir")
		sno > 1 && hideydecorations!(ax, ticks=false, grid=false)
		sno == legend_axno && axislegend(ax; position=legendposition)
	end
	linkaxes!(filter(c -> c isa Axis, fig.content)...)
	fig

end

##
using CairoMakie
fig_path = joinpath(@__DIR__, "figures")
GLMakie.closeall()

fig = plot_btes_temperature_profiles(results, [1, 7, 14]; legendposition = :lb, legend_axno = 3)
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_temperature_charge.png"), fig; backend = CairoMakie)

fig = plot_btes_temperature_profiles(results, [15, 21, 27]; legendposition = :rt, legend_axno = 3)
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_temperature_discharge.png"), fig; backend = CairoMakie)

##

function plot_btes_phase_diagram(results, steps)
	pressure_limits = [1e4, 7e6]
	enthalpy_limits = [1e3, 1e6]

	tables = case.model.models[:Reservoir].system.pvt_tables

	fig, handles = plot_phase_diagram_contours(tables;
		contour_kwargs = (labels=true,),
		variable=:temperature,
		filled=true, lines=true,
		pressure_limits = pressure_limits, enthalpy_limits = enthalpy_limits,
		n_pressure = 1000, n_enthalpy = 1000,
		axis_kwargs = (aspect = AxisAspect(1),),
	)

	cells = case.model.models[:B1_supply].data_domain.representation.perforations.reservoir
	for step in steps
        state_w = results.result.states[step][:B1_supply]
		state_w = results.result.states[step][:B1_supply]
		for (sno, (section, label)) in enumerate(zip((pipe_left, pipe_right, grout_left, grout_right), ["Pipe left", "Pipe right", "Grout left", "Grout right"]))
			ix = section
			p = state_w[:Pressure][ix]./1e6
			h = state_w[:Enthalpy][ix]./1e3
			lines!(handles.axis, h, p; color=colors[sno], linewidth = 4, label = label)
		end

		state_r = results.result.states[step][:Reservoir]
		p_reservoir = state_r[:Pressure][cells]./1e6
		h_reservoir = state_r[:Enthalpy][cells]./1e3
		lines!(handles.axis, h_reservoir, p_reservoir; color=:black, linewidth=4, label="Reservoir")

	end

	Colorbar(fig[1,2], handles.filled, label = "Temperature (°C)", width = 15)

	return fig

end

##
fig = plot_btes_phase_diagram(results, [1, 7, 14])
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_phase_diagram_charge.png"), fig; backend = CairoMakie)

fig = plot_btes_phase_diagram(results, [15, 21, 27])
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_phase_diagram_discharge.png"), fig; backend = CairoMakie)

##
mesh = physical_representation(reservoir_model(case.model).data_domain)
fig = Figure(size = (1200, 1000))
ax = Axis3(fig[1,1], zreversed = true, perspectiveness = 0.75)
plot_mesh_edges!(ax, mesh)
for well in well_symbols(case.model)
	plot_well!(ax, mesh, case.model.models[well].data_domain, markersize = 0.0, fontsize = 0.0)
end
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_setup.png"), fig; backend = CairoMakie)

##
GLMakie.closeall()

states = [state[:Reservoir] for state in results.result.states]
x = tpfv_geometry(mesh).cell_centroids
remove = x[1, :] .< 0.0 .&& x[2, :] .< 0.0
keep = .!remove
function plot_reservoir_states(results, steps)

	crange = (25, 115)

	fig = Figure(size = (1200, 1000))
	for (sno, step) in enumerate(steps)
		ax = Axis3(fig[1,sno], zreversed = true, perspectiveness = 0.75, title = "Time = $(Int(round(time[step], digits = 0))) days", titlegap=-10)
		state_r = results.result.states[step][:Reservoir]
		T = convert_from_si.(state_r[:Temperature][keep], :Celsius)
		plot_cell_data!(ax, mesh, T;
		cells = keep, colormap = :seaborn_icefire_gradient, colorrange = crange)
		# for well in well_symbols(case.model)
		# 	plot_well!(ax, mesh, case.model.models[well].data_domain, markersize = 0.0, fontsize = 0.0)
		# end
		# hidedecorations!(ax)
		# hidespines!(ax)
	end
	Colorbar(fig[1,length(steps)+1], colormap = :seaborn_icefire_gradient, label = "Temperature (°C)", width = 15, limits = crange)

	return fig
end

##
GLMakie.closeall()
fig = plot_reservoir_states(results, [14])
display(GLMakie.Screen(), fig)
save(joinpath(fig_path, "ht_btes_reservoir_temperature.png"), fig; backend = CairoMakie)

fig = plot_reservoir_states(results, [15, 21, 27])
display(GLMakie.Screen(), fig)
# plot_reservoir(case, states; 
# axis_args = (perspectiveness = 0.75,), cells = keep, aspect = (1,1,1), well_fontsize = 0.0, well_arg = (markersize = 0.0,))

##




















##
dd = case.model.models[:B1_supply].data_domain
section =  pipe_left .+ 2.0*pipe_right .+ 3.0*grout_left .+ 4.0*grout_right
dd[:section, Cells()] = Int64.(section)

states, _ = Jutul.expand_to_ministeps(results.result)

JutulDarcy.plot_well_states_interactive(:B1_supply, case.model, states)

##
##
sim, cfg = setup_reservoir_simulator(case0;
	tol_cnv = 1e-2,
	tol_mb = 1e-7,
	info_level = 2,
)
results0 = simulate_reservoir(case0; simulator = sim, config = cfg)

##
st = deepcopy(case.state0)
st = Jutul.evaluate_all_secondary_variables(case.model, st)

##
reservoir_states = [state[:Reservoir] for state in results.states]
min_temperature = [minimum(state[:Temperature]) for state in reservoir_states]
max_temperature = [maximum(state[:Temperature]) for state in reservoir_states]
min_enthalpy = [minimum(state[:Enthalpy]) for state in reservoir_states]
max_enthalpy = [maximum(state[:Enthalpy]) for state in reservoir_states]

##
st = deepcopy(case.state0)
st = Jutul.evaluate_all_secondary_variables(case.model, st)

##
pressure_limits = [1e4, 6e6]
enthalpy_limits = [1e3, 3e6]

tables = case.model.models[:Reservoir].system.pvt_tables
fig, handles = plot_phase_diagram_contours(tables;
    contour_kwargs = (labels=true,),
    variable=:temperature,
    filled=true, lines=true,
    pressure_limits = pressure_limits, enthalpy_limits = enthalpy_limits,
    n_pressure = 1000, n_enthalpy = 1000,
)

# timesteps = 1:length(results.result.states)
# p = [st[:Cold][:Pressure][1] for st in results.result.states]
# h = [st[:Cold][:Enthalpy][1] for st in results.result.states]

well = :B1_supply
pipe_left = case.model.models[well].data_domain[:tag] .== :pipe_left
pipe_right = case.model.models[well].data_domain[:tag] .== :pipe_right
pipe = pipe_left .| pipe_right
if false
    n = 15
    colors = cgrad(:viridis, n, categorical=true)
    lines!(handles.axis, h[1:n]./1e3, p[1:n]./1e6; color=:black, linewidth=2)
    for (k, (pk, hk)) in enumerate(zip(p[1:n], h[1:n]))
        scatter!(handles.axis, hk/1e3, pk/1e6; color=colors[k])
    end
else
    n = length(results.result.states)
    colors = cgrad(:viridis, n, categorical=true)
    for timestep = 1:5:n
        st = results.result.states[timestep][:B1_supply]
		p_pipe = st[:Pressure][pipe]
		h_pipe = st[:Enthalpy][pipe]
		lines!(handles.axis, h_pipe./1e3, p_pipe./1e6; color=:blue, linewidth=2)
		scatter!(handles.axis, h_pipe./1e3, p_pipe./1e6; color=:blue, markersize=10)
		p_grout = st[:Pressure][.!pipe]
		h_grout = st[:Enthalpy][.!pipe]
		lines!(handles.axis, h_grout./1e3, p_grout./1e6; color=:red, linewidth=2)
		# scatter!(handles.axis, h_grout./1e3, p_grout./1e6; color=:red, markersize=10)
        # plot_reservoir_state_ph!(handles.axis, case.model, st, color=colors[timestep])
    end
end
fig