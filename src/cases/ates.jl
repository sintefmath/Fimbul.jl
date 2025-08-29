function ates(;
    well_distance = 500.0,
    temperature_charge = convert_to_si(95, :Celsius),
    temperature_discharge = convert_to_si(25, :Celsius),
    rate_charge = 300meter^3/hour,
    rate_discharge = rate_charge,
    depths = [0.0, 850.0, 900.0, 1000.0, 1050.0, 1300.0],
    aquifer_layer = 3,
    porosity = [0.01, 0.05, 0.2, 0.05, 0.01],
    permeability = [1.0, 5.0, 1000.0, 5.0, 1.0].*1e-3.*si_unit(:darcy),
    thermal_conductivity = [2.5, 2.0, 1.5, 2.0, 2.5].*watt/(meter*Kelvin),
    rock_heat_capacity = fill(900.0*joule/(kilogram*Kelvin), length(depths)-1),
    utes_schedule_args = NamedTuple()
)

    Cf = 4184.0
    Cr = rock_heat_capacity[aquifer_layer]
    ϕ = porosity[aquifer_layer]
    charge_duration = 0.5si_unit(:year)
    Vin = rate_charge*charge_duration
    layer_thickness = diff(depths)
    Haq = layer_thickness[aquifer_layer]

    thermal_radius = thermal_radius_aquifer(Vin, Haq, ϕ, Cf, Cr)

    msh, layers = make_ates_cart_mesh(well_distance, depths, aquifer_layer;
        thermal_radius = thermal_radius,
    )

    porosity = porosity[layers]
    permeability = permeability[layers]
    thermal_conductivity = thermal_conductivity[layers]
    rock_heat_capacity = rock_heat_capacity[layers]
    domain = reservoir_domain(msh;
        porosity = porosity,
        permeability = permeability,
        thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity
    )

    k = cell_ijk(msh, findlast(layers .== aquifer_layer))[3]

    xw_hot  = [-well_distance/2 0.0 0.0; -well_distance/2 0.0 depths[end]]
    cell = Jutul.find_enclosing_cells(msh, xw_hot)[1]
    ij = cell_ijk(msh, cell)[1:2]
    hot_well = setup_vertical_well(domain, ij[1], ij[2]; toe=k, simple_well=false, name = :Hot)
    open = layers[hot_well.perforations.reservoir] .== aquifer_layer
    hot_well.perforations.WI[.!open] .= 0.0

    xw_cold  = [well_distance/2 0.0 0.0; well_distance/2 0.0 depths[end]]
    cell = Jutul.find_enclosing_cells(msh, xw_cold)[1]
    ij = cell_ijk(msh, cell)[1:2]
    cold_well = setup_vertical_well(domain, ij[1], ij[2]; toe=k, simple_well=false, name = :Cold)
    open = layers[cold_well.perforations.reservoir] .== aquifer_layer
    cold_well.perforations.WI[.!open] .= 0.0

    model, _ = setup_reservoir_model(
        domain, :geothermal;
        wells = [hot_well, cold_well]
    )

    bc, state0 = set_dirichlet_bcs(
        model;
        pressure_surface = 1atm,
        temperature_surface = convert_to_si(10.0, :Celsius),
        geothermal_gradient = 0.03Kelvin/meter,
    )
    
    rho = reservoir_model(model).system.rho_ref[1]

    rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Cold])
    ctrl_hot = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_charge)

    rate_target = TotalRateTarget(-rate_charge)
    ctrl_cold = ProducerControl(rate_target)

    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)

    forces_charge = setup_reservoir_forces(model; bc = bc, control = control)

    rate_target = TotalRateTarget(-rate_discharge)
    ctrl_hot  = ProducerControl(rate_target)

    rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Hot])
    ctrl_cold = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_discharge)

    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)

    forces_discharge = setup_reservoir_forces(model; bc = bc, control = control)

    forces_rest = setup_reservoir_forces(model; bc = bc)

    forces, dt = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest;
        utes_schedule_args...
    )

    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end

function make_ates_cart_mesh(well_distance, depths, aquifer_layer;
    hxy_min = missing,
    hxy_max = missing,
    hz_min = missing,
    hz_max = missing,
    thermal_radius = missing,
    offset_rel = 1.0
)

    thermal_radius = ismissing(thermal_radius) ? well_distance/4 : thermal_radius
    offset = offset_rel*maximum([well_distance, 2*thermal_radius])

    xy_hot = (-well_distance/2, 0.0)
    xy_cold = (well_distance/2, 0.0)

    hxy_min = ismissing(hxy_min) ? well_distance/35 : hxy_min
    hxy_max = ismissing(hxy_max) ? hxy_min*10 : hxy_max
    @assert hxy_max ≥ hxy_min "hxy_max must be greater than or equal to hxy_min"

    xy = []
    hxy = [hxy_max, hxy_min, hxy_min, hxy_min, hxy_min, hxy_min, hxy_max];
    for dim = 1:2
        xy_d = [
            xy_hot[dim]-offset,
            xy_hot[dim]-thermal_radius,
            xy_hot[dim]-hxy_min/2,
            xy_hot[dim]+thermal_radius,
            xy_cold[dim]-thermal_radius,
            xy_cold[dim]-hxy_min/2,
            xy_cold[dim]+thermal_radius,
            xy_cold[dim]+offset
        ]
        xy_d, _ = Fimbul.interpolate_z(xy_d, hxy)

        push!(xy, xy_d)
    end

    layer_thickness = diff(depths)
    hz_min = ismissing(hz_min) ? layer_thickness[aquifer_layer]/10 : hz_min
    hz_max = ismissing(hz_max) ? hz_min*10 : hz_max
    @assert hz_max ≥ hz_min "hz_max must be greater than or equal to hz_min"

    hz = fill(hz_max, length(layer_thickness))
    hz[[-1, 0, 1] .+ aquifer_layer] .= hz_min
    z, layers = Fimbul.interpolate_z(depths, hz)

    x = (xy[1], xy[2], z)
    sizes = map(x->diff(x), (xy[1], xy[2], z))
    dims = Tuple([length(s) for s in sizes])

    msh = CartesianMesh(dims, sizes, origin = minimum.(x))

    layers = repeat(layers, inner = dims[1]*dims[2])

    return msh, layers

end