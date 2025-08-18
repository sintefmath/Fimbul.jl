Kelvin = si_unit(:Kelvin)
meter = si_unit(:meter)
year = si_unit(:year)

function fimbul_logo(; 
    L = 1000meter,
    blob_temperatures = convert_to_si.([10.0, 95.0, 95.0], :Celsius),
    time = 500year,
    num_steps = 100,
    mesh_args = NamedTuple()
    )
    
    msh, blob_coords = make_fimbul_mesh(mesh_args...)
    domain = reservoir_domain(msh)
    model, parameters = setup_reservoir_model(domain, :geothermal; kgrad = :avgmpfa);

    geo = tpfv_geometry(msh)

    xy = geo.cell_centroids
    xy = [(x[1], x[3]) for x in eachcol(xy)]

    surface_temperature = convert_to_si(10.0, :Celsius)
    surface_pressure = convert_to_si(5.0, :atm)

    dTdz = 0.03Kelvin/meter
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant

    x = geo.boundary_centroids[[1,3],:]
    x_ext = extrema(x[1,:])
    z_ext = extrema(x[2,:])
    T = z -> surface_temperature .+ dTdz.*z
    p = z -> surface_pressure .+ dpdz.*z
    cells = geo.boundary_neighbors
    bc_cells, bc_pressure, bc_temperature = Int64[], Float64[], Float64[]

    on_boundary = x -> 
    isapprox(x[1], x_ext[1]) || isapprox(x[1], x_ext[2]) ||
    isapprox(x[2], z_ext[1]) || isapprox(x[2], z_ext[2])

    for (k, x_k) in enumerate(eachcol(x))
        !on_boundary(x_k) ? continue : nothing
        println(x_k)
        push!(bc_cells, cells[k])
        push!(bc_pressure, p(x_k[2]))
        push!(bc_temperature, T(x_k[2]))
    end
    bc = flow_boundary_condition(bc_cells, domain, bc_pressure, bc_temperature)

    ##
    z = geo.cell_centroids[3,:]
    z = z .- minimum(z)
    T = surface_temperature .+ dTdz.*z
    p = surface_pressure .+ dpdz.*z

    for (bno, xb) in enumerate(blob_coords)

        xb = [(x[1], -x[2] + L) for x in xb]
        inside = map(p -> Fimbul.point_in_polygon(p, xb), xy)
        cells = findall(inside)
        T[cells] .= blob_temperatures[bno]

    end

    ##
    state0 = setup_reservoir_state(model; Pressure = p, Temperature = T)

    forces = setup_reservoir_forces(model; bc = bc)
    dt = fill(time/num_steps, num_steps)

    case = JutulCase(model, dt, forces; state0 = state0)

    plot_args = (aspect = :data, azimuth = π/2, elevation = 0)

    return case, plot_args

end

function make_fimbul_mesh(; L = 1000, hxy_min = L/50, hxy_max = L/10)

    circle = (x0, r, n) -> 
    [r.*(cos(θ), sin(θ)) .+ x0 for θ in range(π/2, 2π+π/2; length = n+1)]

    x0 = (0.0, 0.0)
    xb = circle((0.0, 0.0), L/3, 3)[1:end-1]
    cc = []
    for x0 in xb
        cc_k = circle(x0, L/5, 20)
        push!(cc, cc_k)
    end
    boundary = [(-L, -L), (L, -L), (L, L), (-L, L)]
    msh0, _ = Fimbul.extruded_mesh(cc, [0.0, 1.0]; 
    boundary = boundary, hxy_min = hxy_min, hxy_max = hxy_max, hz = 1.0, recombine_to_quads = true)

    pts = [(x[1], x[3], -x[2]+L) for x in msh0.node_points]
    pts, dim = Jutul.convert_coord_points(pts)
    msh = UnstructuredMesh(
        msh0.faces.cells_to_faces,
        msh0.boundary_faces.cells_to_faces,
        msh0.faces.faces_to_nodes,
        msh0.boundary_faces.faces_to_nodes,
        pts,
        msh0.faces.neighbors,
        msh0.boundary_faces.neighbors;
        z_is_depth = true
        )

    return msh, cc

end