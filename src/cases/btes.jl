to_kelvin = T -> convert_to_si(T, :Celsius)
meter, = si_units(:meter)
year, day = si_units(:year, :day)
litre, = si_units(:litre)
Kelvin, = si_unit(:Kelvin)
second, = si_units(:second)
darcy, = si_units(:darcy)

function btes(;
    num_wells = 50,
    depth = 50,
    temperature_charge = to_kelvin(90.0),
    temperature_discharge = to_kelvin(10.0),
    rate_charge = 0.5litre/second,
    rate_discharge = 0.5litre/second,
    temperature_surface = to_kelvin(10.0),
    geothermal_gradient = 0.03Kelvin/meter,
    num_years = 5,
    years = missing,
    charge_months = ["June", "July", "August", "September"],
    discharge_months = ["December", "January", "February", "March"],
    report_interval = 14day,
    mehsing_args = NamedTuple(),
    )

    # ## Create mesh
    well_coordinates = fibonachi_pattern_2d(num_wells; spacing = 5.0)
    mesh, metrics = extruded_mesh(well_coordinates, depth)

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite
    domain = reservoir_domain(mesh,
        permeability = 1e-6darcy,
        porosity = 0.01,
        rock_density = 2650kilogram/meter^3,
        rock_heat_capacity = 790joule/kilogram/Kelvin,
        rock_thermal_conductivity = 3.0watt/meter/Kelvin,
        component_heat_capacity = 4.278e3joule/kilogram/Kelvin
    );
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    for (wno, xw) in enumerate(well_coordinates)
        d = max(norm(xw, 2))
        v = (d > 0) ? xw ./d : (1.0, 0.0)
        # Shift coordiates a bit to avoid not be exactly on the node
        xw = xw .+ (hxy_min/2) .* v
        trajectory = [xw[1] xw[2] 0.0; xw[1] xw[2] depth]
        cells = Jutul.find_enclosing_cells(mesh, trajectory)
        name = Symbol("B$wno")
        w_sup, w_ret = setup_btes_well(domain, cells, name=name, btes_type=:u1)
        push!(well_models, w_sup, w_ret)
    end
    # Make the model
    model, parameters = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models
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
    p = z -> 1atm .+ rho.*dpdz.*z
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