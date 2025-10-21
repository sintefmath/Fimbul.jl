meter = si_unit(:meter)
hour, day = si_units(:hour, :day)
darcy = si_unit(:darcy)

"""
    geothermal_doublet(;  <keyword arguments>)

Generic setup function for geothermal doublet case

# Keyword arguments
- `depths = [0.0, 500.0, 2400.0, 2500.0, 3000.0]`: Depths delineating geological layers.
- `permeability = [1e-3, 5e-2, 1.0, 1e-3]*darcy`: Permeability of the layers.
- `porosity = [0.01, 0.2, 0.35, 0.01]`: Porosity of the layers.
- `density = [2000, 2580, 2600, 2400]*kilogram/meter^3`: Rock density in the layers.
- `thermal_conductivity = [2.0, 2.8, 3.5, 1.9]*watt/meter/Kelvin`: Thermal conductivity in the layers.
- `heat_capacity = [1500, 900, 900, 1500]*joule/kilogram/Kelvin`: Heat capacity in the layers.
- `aquifer_layer = 3`: Index of the aquifer layer.
- `spacing_top = 100.0`: Horizontal well spacing at the surface
- `spacing_bottom = 1000.0`: Horizontal well spacing in the aquifer
- `depth_1 = 800.0`: Depth at which well starts to deviate.
- `depth_2 = 2500.0`: Depth of wells
- `temperature_inj = convert_to_si(20.0, :Celsius)`: Injection temperature.
- `rate = 300meter^3/hour`: Injection and production rate.
- `temperature_surface = convert_to_si(10.0, :Celsius)`: Temperature at the
  surface.
- `num_years = 200`: Number of years to run the simulation.
- `report_interval = si_unit(:year)`: Reporting interval for the simulation.

"""
function geothermal_doublet(;
    depths = [0.0, 500.0, 2400.0, 2500.0, 3000.0],
    permeability = [1e-3, 5e-2, 1.0, 1e-3]*darcy,
    porosity = [0.01, 0.2, 0.35, 0.01],
    density = [2000, 2580, 2600, 2400]*kilogram/meter^3,
    thermal_conductivity = [2.0, 2.8, 3.5, 1.9]*watt/meter/Kelvin,
    heat_capacity = [1500, 900, 900, 1500]*joule/kilogram/Kelvin,
    aquifer_layer = 3,
    spacing_top = 100.0,
    spacing_bottom = 1000.0,
    depth_1 = 800.0,
    depth_2 = depths[aquifer_layer+1],
    temperature_inj = convert_to_si(20.0, :Celsius),
    rate = 300meter^3/hour,
    temperature_surface = convert_to_si(10.0, :Celsius),
    num_years = 100,
    report_interval = si_unit(:year),
    )

    trajectory_inj = [
        -spacing_top/2 0.0 0.0;
        -spacing_top/2 0.0 depth_1;
        -spacing_bottom/2 0.0 depth_2;
    ]

    trajectory_prod = [
        spacing_top/2 0.0 0.0;
        spacing_top/2 0.0 depth_1;
        spacing_bottom/2 0.0 depth_2;
    ]

    well_coords = [trajectory_inj, trajectory_prod]

    xw_inj = [Tuple(x) for x in eachrow(unique(trajectory_inj[:, 1:2], dims=1))]
    xw_prod = [Tuple(x) for x in eachrow(unique(trajectory_prod[:, 1:2], dims=1))]
    xw = vcat([xw_inj], [xw_prod])

    # Create the mesh
    thickness = diff(depths)
    nz = [4,10,20,4]
    mesh, layers, metrics = extruded_mesh(xw, depths;
    hxy_min = 100/3, hxy_max = 750, hz = thickness./nz, dist_min_factor = 3.0, offset_rel = 2.5)

    permeability = permeability[layers]
    # permeability = repeat(permeability', 3, 1)
    # permeability[3,:] .*= 0.25 # Reduce vertical permeability
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

    num_layers = Int64.(number_of_cells(mesh)/metrics.nc_2d)
    mesh_layer = repeat(1:num_layers, metrics.nc_2d)

    trajectory_inj[:,2] .+= 0.5*metrics.hxy_min
    cells_inj = Jutul.find_enclosing_cells(mesh, trajectory_inj, n = 1000)
    ix = unique(i -> mesh_layer[cells_inj[i]], eachindex(cells_inj))
    cells_inj = cells_inj[ix]
    WI = [compute_peaceman_index(mesh, permeability[c], 0.1, c) for c in cells_inj]
    WI[layers[cells_inj] .!== 3] .= 0.0
    well_inj = setup_well(domain, cells_inj; 
        name = :Injector, WI = WI, simple_well = false)

    trajectory_prod[:,2] .+= 0.5*metrics.hxy_min
    cells_prod = Jutul.find_enclosing_cells(mesh, trajectory_prod, n = 1000)
    ix = unique(i -> mesh_layer[cells_prod[i]], eachindex(cells_prod))
    cells_prod = cells_prod[ix]
    WI = [compute_peaceman_index(mesh, permeability[c], 0.1, c) for c in cells_prod]
    WI[layers[cells_prod] .!== 3] .= 0.0
    well_prod = setup_well(domain, cells_prod;
        name = :Producer, WI = WI, simple_well = false)

    model = setup_reservoir_model(
        domain, :geothermal,
        wells = [well_inj, well_prod],
    );
    rmodel = reservoir_model(model)
    push!(rmodel.output_variables, :PhaseViscosities)

    bc, state0, = set_dirichlet_bcs(model, [:top, :bottom];
        pressure_surface = 5atm,
        temperature_surface = temperature_surface
    )

    rho = reservoir_model(model).system.rho_ref[1]

    inj_target = JutulDarcy.ReinjectionTarget(NaN, [:Producer])
    # inj_target = TotalRateTarget(rate)
    ctrl_inj = InjectorControl(inj_target, [1.0], 
        density=rho, temperature=temperature_inj, tracers = [1.0])
    # BHP control for producer
    prod_target = TotalRateTarget(-rate)
    ctrl_prod = ProducerControl(prod_target);
    # Set up forces
    control = Dict()
    control[:Injector] = ctrl_inj
    control[:Producer] = ctrl_prod
    forces = setup_reservoir_forces(model, control=control, bc=bc)
    
    # Set up the schedule
    @assert report_interval <= year
    n = Int64(ceil(year/report_interval))
    dt = year/n;
    dt = fill(dt, n*num_years)

    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end