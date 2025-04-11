darcy = si_unit(:darcy)

function geothermal_doublet(;
    temperature_inj = convert_to_si(10.0, :Celsius),
    rate = 50si_unit(:litre)/si_unit(:second),
    temperature_surface = convert_to_si(10.0, :Celsius),
    high_months =  ["January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"],
    low_months = nothing,
    num_years = 20,
    report_interval = si_unit(:year)/12,
    kwargs...
    )

    spacing_top = 100.0
    spacing_bottom = 1500.0
    depth_1 = 800.0
    depth_2 = 2000.0

    trajectory_inj = [
        -spacing_top/2 0.0 0.0;
        -spacing_top/2 0.0 depth_1;
        -spacing_bottom 0.0 depth_2;
    ]

    trajectory_prod = [
        spacing_top/2 0.0 0.0;
        spacing_top/2 0.0 depth_1;
        spacing_bottom 0.0 depth_2;
    ]

    xw_inj = [Tuple(x) for x in eachrow(unique(trajectory_inj[:, 1:2], dims=1))]
    xw_prod = [Tuple(x) for x in eachrow(unique(trajectory_prod[:, 1:2], dims=1))]
    xw = vcat([xw_inj], [xw_prod])

    depths = [0.0, 500.0, 1900.0, 2000.0, 2500.0]

    # Create the mesh
    mesh, layers, metrics = extruded_mesh(xw, depths)

    permeability = [1e-3, 1e-1, 1e-1, 1e-3]*darcy
    porosity = [0.01, 0.2, 0.35, 0.01]
    density = [2000, 2580, 2600, 2400]*kilogram/meter^3
    thermal_conductivity = [2.0, 2.8, 3.5, 1.9]*watt/meter/Kelvin
    heat_capacity = [1500, 900, 900, 1500]*joule/kilogram/Kelvin

    permeability = permeability[layers]
    porosity = porosity[layers]
    density = density[layers]
    thermal_conductivity = thermal_conductivity[layers]
    heat_capacity = heat_capacity[layers]
    domain = reservoir_domain(mesh,
        permeability = permeability,
        porosity = porosity,
        rock_density = density,
        rock_heat_capacity = heat_capacity,
        rock_thermal_conductivity = thermal_conductivity,
        component_heat_capacity = 4.278e3joule/kilogram/Kelvin,
    )

    cells_inj = Jutul.find_enclosing_cells(mesh, trajectory_inj, n = 100)
    WI = map(c -> compute_peaceman_index(mesh, permeability[c], 0.1, c), cells_inj)
    WI[layers[cells_inj] .!== 3] .= 0.0
    well_inj = setup_well(domain, cells_inj; 
        name = :Injector, WI = WI, simple_well = false)

    cells_prod = Jutul.find_enclosing_cells(mesh, trajectory_prod, n = 100)
    WI = map(c -> compute_peaceman_index(mesh, permeability[c], 0.1, c), cells_prod)
    WI[layers[cells_prod] .!== 3] .= 0.0
    well_prod = setup_well(domain, cells_prod;
        name = :Producer, WI = WI, simple_well = false)

    model, parameters = setup_reservoir_model(
        domain, :geothermal,
        wells = [well_inj, well_prod]
    );

    bc, state0, = set_dirichlet_bcs(model;
        pressure_surface = 1atm,
        temperature_surface = temperature_surface
    )

    rho = reservoir_model(model).system.rho_ref[1]

    inj_target = TotalRateTarget(rate)
    ctrl_inj = InjectorControl(inj_target, [1.0], 
        density=rho, temperature=temperature_inj)
    # BHP control for producer
    prod_target = TotalRateTarget(-rate)
    ctrl_prod = ProducerControl(prod_target);
    # Set up forces
    control = Dict()
    control[:Injector] = ctrl_inj
    control[:Producer] = ctrl_prod
    forces = setup_reservoir_forces(model, control=control, bc=bc)
    
    # Set up the schedule
    forces, dt = make_production_schedule(forces, missing, missing;
        high_months = high_months,
        low_months = low_months,
        report_interval = report_interval,
        num_years = num_years
    )

    case = JutulCase(model, dt, forces; state0 = state0)

    return case
end