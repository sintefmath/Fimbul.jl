to_kelvin = T -> convert_to_si(T, :Celsius)
second, year, day = si_units(:second, :year, :day)
meter = si_unit(:meter)
litre = si_unit(:litre)
Kelvin = si_unit(:Kelvin)
darcy = si_unit(:darcy)

function btes(;
    num_wells = 50,
    well_spacing = 5.0,
    depths = [0.0, 0.5, 50, 65],
    well_layers = [1, 2],
    density = [30, 2580, 2580]*kilogram/meter^3,
    thermal_conductivity = [0.034, 3.7, 3.7]*watt/meter/Kelvin,
    heat_capacity = [1500, 900, 900]*joule/kilogram/Kelvin,
    geothermal_gradient = 0.03Kelvin/meter,
    temperature_charge = to_kelvin(90.0),
    temperature_discharge = to_kelvin(10.0),
    rate_charge = 0.5litre/second,
    rate_discharge = rate_charge,
    temperature_surface = to_kelvin(10.0),
    num_years = 5,
    years = missing,
    charge_months = ["June", "July", "August", "September"],
    discharge_months = ["December", "January", "February", "March"],
    report_interval = 14day,
    n_z = [4, 10, 3],
    n_xy = 5,
    mesh_args = NamedTuple(),
    )

    # ## Create mesh
    x = fibonacci_pattern_2d(num_wells; spacing = well_spacing)
    well_coordinates = map(x -> [x], x)
    hz = diff(depths)./n_z
    println(hz)
    hxy = well_spacing/n_xy
    mesh, layers, metrics = extruded_mesh(well_coordinates, depths;
        hxy_min = hxy, hz = hz, mesh_args...)

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite,
    # with a styrofoam layer on top
    density = density[layers]
    thermal_conductivity = thermal_conductivity[layers]
    heat_capacity = heat_capacity[layers]
    domain = reservoir_domain(mesh,
        permeability = 1e-6darcy,
        porosity = 0.01,
        rock_density = density,
        rock_heat_capacity = heat_capacity,
        rock_thermal_conductivity = thermal_conductivity,
        component_heat_capacity = 4.278e3joule/kilogram/Kelvin,
    )
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    nl = length(layers)
    geo = tpfv_geometry(mesh)
    for (wno, xw) in enumerate(well_coordinates)
        name = Symbol("B$wno")
        println("Adding well $name ($wno/$num_wells)")
        xw = xw[1]
        d = max(norm(xw, 2))
        v = (d > 0) ? xw./d : (1.0, 0.0)
        # Shift coordiates a bit to avoid being exactly on the node
        xw = xw .+ (hxy_min/2) .* v
        trajectory = [xw[1] xw[2] 0.0; xw[1] xw[2] depths[end]]
        cells = Jutul.find_enclosing_cells(mesh, trajectory, n = 100)
        filter!(c -> layers[c] ∈ well_layers, cells)
        w_sup, w_ret = setup_btes_well(domain, cells, name=name, btes_type=:u1)
        push!(well_models, w_sup, w_ret)
    end
    # Make the model
    model, parameters = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
    );

    # ## Set up initial state and boundary conditions
    geo = tpfv_geometry(mesh)
    z_bc = geo.boundary_centroids[3, :]
    bottom = map(v -> isapprox(v, maximum(z_bc)), z_bc)
    # Define pressure and temperature profiles
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = geothermal_gradient
    T = z -> temperature_surface .+ dTdz*z
    p = z -> 5atm .+ dpdz.*z
    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- minimum(z_cells)
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );
    # Set boundary conditions
    z_bc = z_bc[.!bottom]
    z_hat = z_bc .- minimum(z_bc)
    bc_cells = geo.boundary_neighbors[.!bottom]
    bc = flow_boundary_condition(bc_cells, domain, p(z_hat), T(z_hat));

    # ## Set up controls
    # Rate control for supply side
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0], 
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(p(0.0))
    ctrl_prod = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()
    for well in well_models
        if contains(String(well.name), "_supply")
            control_charge[well.name] = ctrl_charge
            control_discharge[well.name] = ctrl_discharge
        else
            control_charge[well.name] = ctrl_prod
            control_discharge[well.name] = ctrl_prod
        end
    end
    forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
    forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc);
    forces_rest = setup_reservoir_forces(model, bc=bc)
    # Make schedule
    forces, dt = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest;
        charge_months = charge_months,
        discharge_months = discharge_months,
        num_years = num_years,
        years = years,
        report_interval = report_interval
    )

    # ## Assemble and return model
    case = JutulCase(model, dt, forces, state0 = state0)
    return case

end