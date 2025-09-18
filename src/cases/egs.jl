meter = si_unit(:meter)
day = si_unit(:day)

function egs(well_coords, fracture_radius, fracture_spacing;
    well_names = missing,
    fracture_aperture=1e-3,
    porosity = 0.01,
    permeability = 0.5e-3 .* darcy,
    rock_thermal_conductivity = 2.5 .* watt/(meter*Kelvin),
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    temperature_inj = convert_to_si(25, :Celsius),
    rate_prod = 10kilogram/second/(1000kilogram/meter^3),
    num_years = 20,
    mesh_args = NamedTuple(),
    schedule_args = NamedTuple(),
)

    msh = Fimbul.make_egs_cart_mesh(
        well_coords, fracture_aperture, fracture_radius, fracture_spacing;
        mesh_args...
    )
    msh = UnstructuredMesh(msh)

    Δy = [cell_dims(msh, c)[2] for c in 1:number_of_cells(msh)]

    is_frac = isapprox.(Δy, fracture_aperture)
    geo = tpfv_geometry(msh)
    xz = geo.cell_centroids[[1,3],:]
    xc = 0.0
    for wc in well_coords
        xc = xc .+ wc[2,[1,3]]
    end
    xc = xc ./ (length(well_coords))
    println(xc)
    r = vec(sum((xc .- xz).^2, dims=1).^0.5)
    is_frac = is_frac .&& (r .<= fracture_radius)

    nc = number_of_cells(msh)
    function process_prop(v, v_frac)
        v = fill(v, nc)
        v[is_frac] .= v_frac
        return v
    end

    permeability_frac = fracture_aperture^2/12
    porosity = process_prop(porosity, 0.5)
    permeability = process_prop(permeability, permeability_frac)
    domain = reservoir_domain(msh;
        porosity = porosity,
        permeability = permeability,
        rock_thermal_conductivity = rock_thermal_conductivity,
        rock_heat_capacity = rock_heat_capacity
    )

    if ismissing(well_names)
        well_names = [:Inj]
        for i in 2:length(well_coords)
            push!(well_names, Symbol("Prod$(i-1)"))
        end
    end

    # TODO: support multiple producers
    wells = []
    N = get_neighborship(UnstructuredMesh(msh))
    z = geo.cell_centroids[3, :]
    WI_max = 1e-9

    cells = find_well_cells(msh, well_coords[1], z, N, is_frac)
    well_inj = setup_well(domain, cells;
            name = :Injector,
            simple_well = false
        )
    
    open = isapprox.(z[cells], maximum(z[cells]), atol=1.0)
    well_inj.perforations.WI[.!open] .= 0.0
    well_inj.perforations.WI[open] = min.(well_inj.perforations.WI[open], WI_max)
    push!(wells, well_inj)

    println("Adding well Prod")
    cells = [cells[1]]
    neighbors = [0; 0]
    cell_no = 1;
    for (wno, wc) in enumerate(well_coords)
        wno == 1 ? continue : nothing
        println("Adding producer branch $(wno-1)/$(length(well_coords)-1)")
        cells_k = Jutul.find_enclosing_cells(msh, wc, n = 1000)
        num_cells = length(cells_k)
        push!(cells, cells_k...)
        joint_neighbors = [1; cell_no + 1]
        nodes = collect(1:num_cells) .+ cell_no
        branch_neighbors = [nodes[1:end-1]'; nodes[2:end]']
        neighbors = hcat(neighbors, joint_neighbors, branch_neighbors)
        cell_no += num_cells
    end

    neighbors = Int64.(neighbors[:, 2:end])
    cells = Int64.(cells)

    well_prod = setup_well(domain, cells;
    N = neighbors, name = :Producer, simple_well = false)

    wells = [well_inj, well_prod]

    model, _ = setup_reservoir_model(domain, :geothermal; wells = wells, split_wells = true)
    bc, state0 = set_dirichlet_bcs(model, pressure_surface = 10atm)

    rho = reservoir_model(model).system.rho_ref[1]
    # Injector: inject cooled water at a rate equal to production
    # TODO: assert that well exists at construction
    rate_target = TotalRateTarget(rate_prod)
    # rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Prod1])
    ctrl_inj = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_inj)
    # Producer
    # TODO: support multiple periods with different rates
    rate_target = TotalRateTarget(-rate_prod)
    ctrl_prod = ProducerControl(rate_target)
    control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)

    limits = Dict()
    limits[:Producer] = (bhp = 10atm,)
    forces = setup_reservoir_forces(model, control = control, bc = bc, limits = limits)

    # dt_target = si_unit(:year)/12

    # dt_0 = 1si_unit(:hour)
    # n_rampup = 20
    # # dt_0*α^(n_rampup-1) = dt_target
    # α = (dt_target/dt_0)^(1/(n_rampup-1))
    # dt = [dt_0*α^(i-1) for i in 1:n_rampup]

    # # n_step = Int(ceil(si_unit(:year)/dt_target))
    # dt = vcat(dt, fill(dt_target, num_years*12))

    dt, forces = make_schedule(
        [forces], [(1,1), (12,31)];
        num_years = num_years,
        schedule_args...)

    # ## Create and return the complete simulation case
    case = JutulCase(model, dt, forces; state0 = state0)
    return case

end

function find_well_cells(msh, wc, z, N, is_frac)
    cells = Jutul.find_enclosing_cells(msh, wc, n = 1000)
    frac_cells = Int[]
    for cell in cells
        faces = msh.faces.cells_to_faces[cell]
        c = vcat(N[:, faces]...)
        keep = is_frac[c]
        push!(frac_cells, c[keep]...)
    end
    push!(cells, frac_cells...)
    unique!(cells)
    cells = topo_sort_well(cells, msh, N, z)
end

function make_egs_cart_mesh(well_coords, fracture_aperture, fracture_radius, fracture_spacing;
    hx_min = missing,
    hx_max = missing,
    hy_min = missing,
    hy_max = missing,
    hz_min = missing,
    hz_max = missing,
    offset_rel = 2.0,
)
    # Calculate domain offset to minimize boundary effects
    well_distance = abs(well_coords[1][1,1] - well_coords[2][1,1])
    y_min = minimum(well_coords[1][:, 2])
    y_max = maximum(well_coords[1][:, 2])

    depth = well_coords[1][3,3]

    offset = offset_rel*(well_distance + fracture_radius)
    # Well coordinates
    hx_min = ismissing(hx_min) ? well_distance/10 : hx_min
    hx_max = ismissing(hx_max) ? offset/5 : hx_max
    hz_min = ismissing(hz_min) ? hx_min : hz_min
    hz_max = ismissing(hz_max) ? (depth - fracture_radius)/5 : hz_max

    function make_coords_xz(x0, xw, h_min, h_max)

        x_min, x_max = extrema(xw)
        x_mid = (x_min + x_max)/2

        xw .-= h_min/2
        xw = vcat(x_mid - fracture_radius, xw, x_mid + fracture_radius)
        xd = sort(vcat(x0, xw))
        nx = length(xd)-1
        hxz = fill(h_max, nx)
        x_mid = diff(xd)./2 .+ xd[1:end-1]
        is_frac = x_min - fracture_radius .<= x_mid .<= x_max + fracture_radius
        hxz[is_frac] .= h_min

        xd = first(Fimbul.interpolate_z(xd, hxz))

        return xd

    end

    xw = []
    for wc in well_coords
        x = unique!(wc[:, 1])
        push!(xw, x...)
    end
    x_min, x_max = extrema(xw)
    x_mid = (x_min + x_max)/2
    x0 = [
        x_mid - offset,
        x_mid + offset
    ]
    x = make_coords_xz(x0, xw, hx_min, hx_max)
    println(x)

    zw = []
    for wc in well_coords
        z = wc[2,3]
        push!(zw, z...)
    end
    z0 = [
        0.0,
        4000.0
    ]
    z = make_coords_xz(z0, zw, hz_min, hz_max)

    yw = []
    for wc in well_coords
        y = unique!(wc[:, 2])
        push!(yw, y...)
    end
    unique!(yw)

    y_min, y_max = extrema(yw)

    hy_min = ismissing(hy_min) ? fracture_spacing/3 : hy_min
    hy_max = ismissing(hy_max) ? max(offset, y_max - fracture_radius)/5 : hy_max

    fracture_start = y_min + fracture_spacing/2
    fracture_end = y_max - fracture_spacing/2

    fracture_y = fracture_start:fracture_spacing:fracture_end
    fracture_y = sort(vcat(fracture_y, fracture_y .+ fracture_aperture))

    y = [
        y_min - offset,
        y_min - hy_min/2,
        fracture_y...,
        y_max - hy_min/2,
        y_max + offset
    ]

    println(y)

    hy = fill(hy_max, length(y)-1)
    y_mid = diff(y)./2 .+ y[1:end-1]
    is_frac = y_min .< y_mid .< y_max
    hy[is_frac] .= hy_min

    y = first(Fimbul.interpolate_z(y, hy))

    xyz = (x, y, z)
    sizes = map(x->diff(x), xyz)
    dims = Tuple([length(s) for s in sizes])
    msh = CartesianMesh(dims, sizes, origin = minimum.(xyz))

    return msh

end