function ates(;
    well_distance = missing,
    depths = [0.0, 850.0, 900.0, 1000.0, 1050.0, 1300.0],
    porosity = [0.01, 0.05, 0.35, 0.05, 0.01],
    permeability = [1.0, 5.0, 1000.0, 5.0, 1.0].*1e-3.*si_unit(:darcy),
    rock_thermal_conductivity = [2.5, 2.0, 1.5, 2.0, 2.5].*watt/(meter*Kelvin),
    rock_heat_capacity = fill(900.0*joule/(kilogram*Kelvin), length(depths)-1),
    aquifer_layer = 3,
    temperature_charge = convert_to_si(95, :Celsius),
    temperature_discharge = convert_to_si(25, :Celsius),
    rate_charge = missing,
    rate_discharge = rate_charge,
    balanced_injection = true,
    temperature_surface = convert_to_si(10.0, :Celsius),
    thermal_gradient = 0.03Kelvin/meter,
    utes_schedule_args = NamedTuple(),
    use_2d = false,
    mesh_args = NamedTuple(),
)

    Cf = 4184.0
    Cr = rock_heat_capacity[aquifer_layer]
    ϕ = porosity[aquifer_layer]
    Caq = Cf*ϕ + Cr*(1 - ϕ)
    layer_thickness = diff(depths)
    Haq = layer_thickness[aquifer_layer]
    charge_duration = 0.5si_unit(:year)
    thermal_radius = missing
    if ismissing(rate_charge)
        thermal_radius = 250.0
        rate_charge = (thermal_radius^2*Caq*π*Haq/Cf)/charge_duration
        if use_2d
            rate_charge *= 2/(π*thermal_radius)
        end
    end
    if ismissing(rate_discharge)
        rate_discharge = rate_charge
    end
    Vin = rate_charge*charge_duration

    if ismissing(thermal_radius)
        thermal_radius = thermal_radius_aquifer(Vin, Haq, ϕ, Cf, Cr)
    end
    if ismissing(well_distance)
        well_distance = 2*thermal_radius
    end

    msh, layers = make_ates_cart_mesh(well_distance, depths, aquifer_layer;
        thermal_radius = thermal_radius,
        use_2d = use_2d,
        mesh_args...
    )

    porosity = porosity[layers]
    permeability = permeability[layers]
    rock_thermal_conductivity = rock_thermal_conductivity[layers]
    rock_heat_capacity = rock_heat_capacity[layers]
    domain = reservoir_domain(msh;
        porosity = porosity,
        permeability = permeability,
        rock_thermal_conductivity = rock_thermal_conductivity,
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

    subset = use_2d ? [:top, :bottom] : :all
    # bc, state0 = set_dirichlet_bcs(
    #     model, subset;
    #     pressure_surface = 1atm,
    #     temperature_surface = convert_to_si(10.0, :Celsius),
    #     geothermal_gradient = 0.03Kelvin/meter,
    # )
    
    rho = reservoir_model(model).system.rho_ref[1]

    rmodel = reservoir_model(model)
    geo = tpfv_geometry(msh)

    rho = rmodel.system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = thermal_gradient
    p = z -> 1atm .+ dpdz.*z
    T = z -> temperature_surface .+ dTdz*z

    # Set boundary conditions
    z_bdr = geo.boundary_centroids[3, :]
    xy_bdr = geo.boundary_centroids[1:2, :]
    cells_bdr = geo.boundary_neighbors
    z0 = minimum(z_bdr)

    top = isapprox.(z_bdr, minimum(z_bdr))
    bottom = isapprox.(z_bdr, maximum(z_bdr))
    west = isapprox.(xy_bdr[1, :], minimum(xy_bdr[1, :]))
    east = isapprox.(xy_bdr[1, :], maximum(xy_bdr[1, :]))
    south = isapprox.(xy_bdr[2, :], minimum(xy_bdr[2, :]))
    north = isapprox.(xy_bdr[2, :], maximum(xy_bdr[2, :]))

    cells_bc = Int64[]
    z_bc = Float64[]
    if use_2d
        subset = top .|| bottom .|| west .|| east
    else
        subset = top .|| bottom .|| south .|| north .|| west .|| east
    end
    # for s in subset
    #     if s == :top
    #         ix = top
    #     elseif s == :bottom
    #         ix = bottom
    #     elseif s == :sides
    #         ix = sides
    #     else
    #         @error "Unknown boundary condition subset: $s"
    #     end
    #     push!(cells_bc, cells_bdr[ix]...)
    #     push!(z_bc, z_bdr[ix]...)
    # end
    cells_bc = cells_bdr[subset]
    z_bc = z_bdr[subset]
    z_hat = z_bc .- z0
    bc = flow_boundary_condition(cells_bc, rmodel.data_domain, p(z_hat), T(z_hat));

    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- z0
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );

    if balanced_injection
        rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Cold])
    else
        rate_target = TotalRateTarget(rate_charge)
    end
    ctrl_hot = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_charge)

    rate_target = TotalRateTarget(-rate_charge)
    ctrl_cold = ProducerControl(rate_target)

    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)

    forces_charge = setup_reservoir_forces(model; bc = bc, control = control)

    rate_target = TotalRateTarget(-rate_discharge)
    ctrl_hot  = ProducerControl(rate_target)
    if balanced_injection
        rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Hot])
    else
        rate_target = TotalRateTarget(rate_discharge)
    end
    ctrl_cold = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_discharge)

    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)

    forces_discharge = setup_reservoir_forces(model; bc = bc, control = control)

    forces_rest = setup_reservoir_forces(model; bc = bc)

    dt, forces = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest,
        utes_schedule_args...
    )

    case = JutulCase(model, dt, forces; state0 = state0)

    return case

end

function ates_simple(;
    well_distance = 500.0,
    aquifer_thickness = 100.0,
    depth = 1000.0,
    porosity = [0.2, 0.01],
    permeability = [1000.0, 1.0].*1e-3.*si_unit(:darcy),
    thermal_conductivity = [2.0, 2.0].*watt/(meter*Kelvin),
    rock_heat_capacity = [900.0, 900.0]*joule/(kilogram*Kelvin),
    kwargs...
)

    aquifer_layer = 2
    # TODO: add back-of-the-envelope calculation in how far a the plume will
    # spread downwards via conduction
    depths = [0.0, depth, depth + aquifer_thickness, (depth + aquifer_thickness)*1.1]

    make_prop = prop -> [prop[1], prop[2], prop[1]]
    porosity = make_prop(porosity)
    permeability = make_prop(permeability)
    thermal_conductivity = make_prop(thermal_conductivity)
    rock_heat_capacity = make_prop(rock_heat_capacity)

    return ates(;
        well_distance = well_distance,
        depths = depths,
        porosity = porosity,
        permeability = permeability,
        thermal_conductivity = thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity,
        aquifer_layer = aquifer_layer,
        kwargs...
    )

end

function make_ates_cart_mesh(well_distance, depths, aquifer_layer;
    hxy_min = missing,
    hxy_max = missing,
    hz_min = missing,
    hz_max = missing,
    thermal_radius = missing,
    offset_rel = 1.0,
    use_2d = false,
)

    thermal_radius = ismissing(thermal_radius) ? well_distance/4 : thermal_radius
    offset = offset_rel*maximum([well_distance, 2*thermal_radius])

    xy_hot = (-well_distance/2, 0.0)
    xy_cold = (well_distance/2, 0.0)

    hxy_min = ismissing(hxy_min) ? well_distance/25 : hxy_min
    hxy_max = ismissing(hxy_max) ? hxy_min*10 : hxy_max
    @assert hxy_max ≥ hxy_min "hxy_max must be greater than or equal to hxy_min"

    xy = []
    hxy = [hxy_max, hxy_min, hxy_min, hxy_min, hxy_min, hxy_min, hxy_max];
    dims = use_2d ? 1 : 1:2
    for dim in dims
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

    if use_2d
        push!(xy, [-0.5, 0.5])
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