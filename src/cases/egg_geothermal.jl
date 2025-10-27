"""
    egg_geothermal(; <keyword arguments>)

# Keyword arguments

- `include_wells = true`: Include wells in the model.
- `wells_distance_ij = 20`: Distance (in cell numbers) between the wells in i
  and j directions.
- `simple_well = false`: Use simple well model instead of the full well model.
- `geothermal_gradient = 0.03si_unit(:Kelvin)/si_unit(:meter)`: Geothermal
  gradient in the reservoir.
- `pressure_top = 50si_unit(:bar)`: Pressure at the top of the reservoir.
- `temperature_top = convert_to_si(50, :Celsius)`: Temperature at the top of the
  reservoir.
- `use_bc = true`: Use fixed Dirichlet pressure and temperature boundary
  conditions at the reservoir sides.

"""
function egg_geothermal(;
    include_wells = true,
    well_distance_ij = 20,
    simple_well = false,
    geothermal_gradient = 0.03si_unit(:Kelvin)/si_unit(:meter),
    pressure_top = 50si_unit(:bar),
    temperature_top = convert_to_si(50, :Celsius),
    use_bc = true,
    )

    # ## Load EGG data from file
    egg_dir = JutulDarcy.GeoEnergyIO.test_input_file_path("EGG")
    data_pth = joinpath(egg_dir, "EGG.DATA")
    case0 = setup_case_from_data_file(data_pth)
    domain = reservoir_model(case0.model).data_domain

    # ## Set up wells
    if include_wells
        d_ij = ceil(well_distance_ij/2)
        d_ij = clamp(d_ij,1,29)
        inj_well = setup_vertical_well(domain, 30 + d_ij, 30 - d_ij;
            name = :WellA, simple_well = simple_well)
        prod_well = setup_vertical_well(domain, 30 - d_ij, 30 + d_ij;
            name = :WellB, simple_well = simple_well)
        obs_well = setup_vertical_well(domain, 30, 30;
            name = :WellObs, simple_well = simple_well)
        wells = [inj_well, prod_well, obs_well]
    else
        wells = []
    end

    # ## Set up reservoir model
    model = setup_reservoir_model(
        domain, :geothermal,
        thermal = true,
        wells = wells
    );
    rmodel = reservoir_model(model)
    push!(rmodel.output_variables, :PhaseMassDensities, :PhaseViscosities)

    # ## Get geometry and gradients
    mesh = physical_representation(rmodel.data_domain)
    geo = tpfv_geometry(mesh)
    z = geo.cell_centroids[3,:]
    depth = z .- minimum(z)
    rho_ref = rmodel.system.rho_ref[1]
    dp_dz = rho_ref.*9.81
    dT_dz = geothermal_gradient

    # ## Set up initial state
    state0 = setup_reservoir_state(model,
        Pressure = pressure_top .+ dp_dz.*depth,
        Temperature = temperature_top .+ dT_dz.*depth
    )

    # ## Set up forces
    if use_bc
        sides = abs.(geo.boundary_normals[3,:]) .== 0
        cells = mesh.boundary_faces.neighbors[sides]
        bc_cells, bc_pressure, bc_temperature = [], [], []
        for c in cells
            push!(bc_cells, c)
            push!(bc_pressure, pressure_top + dp_dz*depth[c])
            push!(bc_temperature, temperature_top + dT_dz*depth[c])
        end
        bc = flow_boundary_condition(bc_cells, domain, bc_pressure, bc_temperature)
    else
        bc = nothing
    end

    dt = 100si_unit(:year)
    forces = setup_reservoir_forces(model, bc = bc)
    
    # ## Return case
    return JutulCase(model, [dt], [forces], state0 = state0)

end

"""
    egg_geothermal_doublet(; <keyword arguments>)

# Keyword arguments
- `rate_injection = 10si_unit(:litre)/si_unit(:second)`: Injection rate.
- `rate_production = rate_injection`: Production rate.
- `temperature_injection = convert_to_si(10.0, :Celsius)`: Injection temperature.
- `rate_observation = missing`: Observation well rate. Set to
  min(rate_injection, rate_production)/1000 if missing.
- `bhp = 45*si_unit(:bar)`: Bottom hole pressure of supporting well if BCs are
  not given (default - set use_bc = true to force BCs).
- `report_interval = si_unit(:year)/4`: Report interval for simulation results.
- `production_time = 25.0si_unit(:year)`: Production time.
- `well_distance_ij = 30`: Distance between the wells in i and j directions.
- All other keyword arguments are passed to `egg_geothermal`.

""" 
function egg_geothermal_doublet(;
    rate_injection = 10si_unit(:litre)/si_unit(:second),
    rate_production = rate_injection,
    temperature_injection = convert_to_si(10.0, :Celsius),
    rate_observation = missing,
    bhp = 45*si_unit(:bar),
    report_interval = si_unit(:year)/4,
    production_time = 25.0si_unit(:year),
    kwargs...
    )

    # ## Set up base case
    case0 = egg_geothermal(;
        well_distance_ij = 30,
        temperature_top = convert_to_si(90.0, :Celsius),
        use_bc = false,
        kwargs...
        )

    # ## Set up forces
    rmodel = reservoir_model(case0.model)
    rho_ref = rmodel.system.rho_ref[1]
    bc = case0.forces[1][:Reservoir].bc
    # Set observation well control
    if ismissing(rate_observation)
        rate_observation = 1e-3*min(rate_injection, rate_production)
    end
    target = TotalRateTarget(-rate_observation)
    ctrl_obs = ProducerControl(target)
    # Set controls
    if isnothing(bc)
        target = BottomHolePressureTarget(bhp)
    else 
        target = TotalRateTarget(-rate_production)
    end
    ctrl_prod = ProducerControl(target)
    target = TotalRateTarget(rate_injection)
    ctrl_inj = InjectorControl(target, [1.0];
        density = rho_ref, temperature = temperature_injection)
    control = Dict(:WellA => ctrl_prod, :WellB => ctrl_inj, :WellObs => ctrl_obs)
    # Define reservoir forces
    forces = setup_reservoir_forces(case0.model, control = control, bc = bc)
    # Set report steps
    n_step = Int(floor(production_time/report_interval))
    dt = fill(report_interval, n_step)
    time_residual = production_time - sum(dt)
    if time_residual < 0.1*report_interval
        dt[end] += time_residual
    else
        dt = vcat(dt, time_residual)
    end
    forces = fill(forces, n_step)

    # ## Return case
    case = JutulCase(case0.model, dt, forces, 
        state0 = case0.state0, parameters = case0.parameters)
    return case

end


"""
    egg_ates(;
    temperature_charge = convert_to_si(90.0, :Celsius),
    temperature_discharge = convert_to_si(10.0, :Celsius),
    rate_charge = 25si_unit(:litre)/si_unit(:second),
    rate_discharge = rate_charge,
    rate_observation = missing,
    bhp_charge = 25.0si_unit(:bar),
    bhp_discharge = 45.0si_unit(:bar),
    charge_months = ("June", "July", "August", "September"),
    discharge_months = ("December", "January", "February", "March"),
    report_interval = si_unit(:year)/12,
    num_years = 5,
    kwargs...
    )

# Keyword arguments
- `temperature_charge = convert_to_si(90.0, :Celsius)`: Charge temperature.
- `temperature_discharge = convert_to_si(10.0, :Celsius)`: Discharge temperature.
- `rate_charge = 25si_unit(:litre)/si_unit(:second)`: Charge rate.
- `rate_discharge = rate_charge`: Discharge rate.
- `rate_observation = missing`: Observation well rate. Set to
  min(rate_charge, rate_discharge)/1000 if missing.
- `bhp_charge = 25.0si_unit(:bar)`: Bottom hole pressure of supporting well
  during charge if BCs are not given (not default - set use_bc = false to force
  no BCs).
- `bhp_discharge = 45.0si_unit(:bar)`: Bottom hole pressure of supporting well
  during discharge if BCs are not given (not default - set use_bc = false to
  force no BCs).
- `charge_period = ["June", "September"]`: Charge period.
- `discharge_period = ["December", "March"]`: Discharge period.
- `num_years = 5`: Number of years to simulate.
- `report_interval = si_unit(:year)/12`: Report interval for simulation results.
- `utes_schedule_args = NamedTuple()`: Additional arguments passed to
  `make_utes_schedule()`.
- All other keyword arguments are passed to `egg_geothermal`.

"""
function egg_ates(;
    temperature_chare = convert_to_si(90.0, :Celsius),
    temperature_dischare = convert_to_si(10.0, :Celsius),
    rate_charge = 25si_unit(:litre)/si_unit(:second),
    rate_discharge = rate_charge,
    rate_observation = missing,
    bhp_charge = 25.0si_unit(:bar),
    bhp_discharge = 45.0si_unit(:bar),
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    num_years = 5,
    report_interval = si_unit(:year)/12,
    utes_schedule_args = NamedTuple(),
    kwargs...
    )

    # ## Get base case
    case0 = egg_geothermal(; kwargs...)

    # ## Set up forces
    rmodel = reservoir_model(case0.model)
    rho_ref = rmodel.system.rho_ref[1]
    bc = case0.forces[1][:Reservoir].bc
    # Set observation well control
    if ismissing(rate_observation)
        rate_observation = 1e-3*min(rate_charge, rate_discharge)
    end
    target = TotalRateTarget(-rate_observation)
    ctrl_obs = ProducerControl(target)
    # Set charge controls
    target = TotalRateTarget(rate_charge)
    ctrl_hot = InjectorControl(target, [1.0];
        density = rho_ref, temperature = temperature_chare)
    if isnothing(bc)
        target = BottomHolePressureTarget(bhp_charge)
    else 
        target = TotalRateTarget(-rate_charge)
    end
    ctrl_cold = ProducerControl(target)
    control = Dict(:WellA => ctrl_hot, :WellB => ctrl_cold, :WellObs => ctrl_obs)
    forces_charge = setup_reservoir_forces(case0.model, control = control, bc = bc)
    # Set discharge controls
    target = TotalRateTarget(-rate_discharge)
    ctrl_hot = ProducerControl(target)
    if isnothing(bc)
        target = BottomHolePressureTarget(bhp_discharge)
    else
        target = TotalRateTarget(rate_discharge)
    end
    ctrl_cold = InjectorControl(target, [1.0];
        density = rho_ref, temperature = temperature_dischare)
    control = Dict(:WellA => ctrl_hot, :WellB => ctrl_cold, :WellObs => ctrl_obs)
    forces_discharge = setup_reservoir_forces(case0.model, control = control, bc = bc)
    # Set rest forces
    forces_rest = setup_reservoir_forces(case0.model, bc = bc)
    dt, forces = make_utes_schedule(forces_charge, forces_discharge, forces_rest;
        charge_period = charge_period,
        discharge_period = discharge_period,
        num_years = num_years,
        report_interval = report_interval
    )

    # ## Return case
    case = JutulCase(case0.model, dt, forces, 
        state0 = case0.state0, parameters = case0.parameters)
    return case

end