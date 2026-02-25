import Jutul.CutCellMeshes: PlaneCut, PolygonalSurface, cut_mesh

function ftes(well_coordinates::Vector{Matrix{Float64}}, num_fractures::Int;
    fracture_dip = (0.0, 5.0),
    fracture_strike = (0.0, 5.0),
    fracture_radius = Inf,
    fracture_aperture = 0.5e-3,
    fracture_porosity = 0.5,
    marix_permeability = 1e-4si_unit(:darcy),
    matrix_porosity = 0.01,
    density = 2500kilogram/meter^3,
    thermal_conductivity = 2.5watt/meter/Kelvin,
    heat_capacity = 1000joule/kilogram/Kelvin,
    rate = 200si_unit(:litre)/si_unit(:second),
    temperature_charge = convert_to_si(95.0, :Celsius),
    temperature_discharge = convert_to_si(10.0, :Celsius),
    charge_period = ["April", "November"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    mesh_args = NamedTuple(),
    )

    # Make constraints from well coordinates
    collars = hcat([x[1:2, 1] for x in well_coordinates]...)
    Δx_min, Δx_max = Fimbul.min_max_distance(collars)
    hxy_min = Δx_min/3
    well_outline = Fimbul.offset_boundary(collars, Δx_max; n=24)
    well_outline = hcat(well_outline, well_outline[:, 1]) # Close the loop
    collars = [permutedims([x[1] x[2]]).+hxy_min/2 for x in eachcol(collars)]
    cell_constraints = [x for x in collars]
    push!(cell_constraints, well_outline)
    # Determine layers form the well depth
    well_depth = maximum(maximum(x[3, :]) for x in well_coordinates)
    depths = [0.0, well_depth+1e-2, well_depth*1.25]
    # Generate mesh
    matrix_mesh, layers, _ = extruded_mesh(cell_constraints, depths; hxy_min=hxy_min, mesh_args...)
    # Add fractures
    fracture_faces = Int[]
    for _ in 1:num_fractures
        is_frac = falses(number_of_faces(matrix_mesh))
        is_frac[fracture_faces] .= true
        α = fracture_strike[1] + randn()*fracture_strike[2]
        δ = fracture_dip[1] + randn()*fracture_dip[2]
        normal = strike_dip_to_normal(α, δ)
        center = [0.0, 0.0, well_depth*rand()]
        plane = PlaneCut(center, normal)
        matrix_mesh, info = cut_mesh(matrix_mesh, plane; min_cut_fraction = 0.0, extra_out=true)
        face_index = filter(x->x>0, info["face_index"])
        fracture_faces = findall(is_frac[face_index])
        new_faces = findall(info["face_index"] .== 0)
        append!(fracture_faces, new_faces)
    end
    # Generate embedded mesh for fractures
    fracture_mesh = Jutul.EmbeddedMeshes.EmbeddedMesh(matrix_mesh, fracture_faces)

    # Generate matrix and fracture domains
    matrix = reservoir_domain(matrix_mesh;
        permeability=marix_permeability, porosity=matrix_porosity)
    fractures = JutulDarcy.fracture_domain(fracture_mesh;
        aperture=fracture_aperture, porosity=fracture_porosity, matrix_faces=fracture_faces)

    frac_args = (WI = missing,)
    wells = []
    for (wno, x) in enumerate(well_coordinates)
        if wno == 1
            name = :Injector
        else
            name = Symbol("Producer_$(wno-1)")
        end
        cells = Jutul.find_enclosing_cells(matrix_mesh, permutedims(x), n=1_000_000)
        well = setup_well(matrix, cells, fractures; name = name, frac_args=frac_args)
        push!(wells, well)
    end
    
    model = JutulDarcy.setup_reservoir_model(matrix, fractures, :geothermal;
        wells=wells, block_backend=true)

    geo = tpfv_geometry(matrix_mesh)
    ρ = reservoir_model(model).system.rho_ref[1]
    dpdz = gravity_constant * ρ
    p0(z) = 10si_unit(:atm)# .+ dpdz.*z
    dTdz = 0.03si_unit(:Kelvin)/si_unit(:meter)
    T0(z) = convert_to_si(10.0, :Celsius)# .+ dTdz.*z
    z = geo.cell_centroids[3, :]

    state0 = setup_reservoir_state(model; Pressure=p0(z), Temperature=T0(z))
    
    bc_cells = geo.boundary_neighbors
    bc_temperature = state0[:Reservoir][:Temperature][bc_cells]
    bc_pressure = state0[:Reservoir][:Pressure][bc_cells]
    bc = flow_boundary_condition(bc_cells, matrix, bc_pressure, bc_temperature)

    ctrl_prod = ProducerControl(BottomHolePressureTarget(p0(0.0)*0.1))

    ctrl_inj = InjectorControl(TotalRateTarget(rate), [1.0];
        density=ρ, temperature=temperature_charge)
    control_charge = Dict()
    for w in wells
        name = w.representation.name
        if name == :Injector
            control_charge[name] = ctrl_inj
        else
            control_charge[name] = ctrl_prod
        end
    end
    forces_charge = setup_reservoir_forces(model, bc=bc, control=control_charge)
    
    ctrl_inj = InjectorControl(TotalRateTarget(rate), [1.0];
        density=ρ, temperature=temperature_discharge)
    control_discharge = Dict()
    for w in wells
        name = w.representation.name
        if name == :Injector
            control_discharge[name] = ctrl_inj
        else
            control_discharge[name] = ctrl_prod
        end
    end
    forces_discharge = setup_reservoir_forces(model, bc=bc, control=control_discharge)

    forces_rest = setup_reservoir_forces(model, bc=bc)

    dt, forces = make_utes_schedule(forces_charge, forces_discharge, forces_rest;
        charge_period=charge_period, discharge_period=discharge_period,
        utes_schedule_args...)

    case = JutulCase(model, dt, forces; state0=state0)

    return case

end

function strike_dip_to_normal(strike::Float64, dip::Float64)
    strike_rad = deg2rad(strike)
    dip_rad = deg2rad(dip)
    normal_x = -sin(dip_rad) * sin(strike_rad)
    normal_y = -sin(dip_rad) * cos(strike_rad)
    normal_z = cos(dip_rad)
    return [normal_x, normal_y, normal_z]
end

function ftes(num_producers::Int, radius::Float64, depth::Float64, num_fractures::Int; kwargs...)

    # Place producers in a circle around the injector
    Δθ = 2π/(num_producers)
    r = radius
    producer_coordinates = [[r*cos(i*Δθ), r*sin(i*Δθ), 0.0] for i in 0:num_producers-1]
    wc = pushfirst!(producer_coordinates, [0.0, 0.0, 0.0]) # Add injector at the center

    well_coordinates = Vector{Matrix{Float64}}(undef, length(wc))
    for (wno, x) in enumerate(wc)
        x = repeat(x, 1, 2)
        x[3, 2] = depth
        well_coordinates[wno] = x
    end
    well_coordinates = [wc for wc in well_coordinates]

    return ftes(well_coordinates, num_fractures; kwargs...)

end