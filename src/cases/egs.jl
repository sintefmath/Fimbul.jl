meter = si_unit(:meter)
day = si_unit(:day)

function egs(well_coords, fracture_radius, fracture_spacing;
    well_names = missing,
    fracture_aperture=0.5e-3,
    porosity = 0.01,
    permeability = 1e-3 .* darcy,
    rock_thermal_conductivity = 2.5 .* watt/(meter*Kelvin),
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    temperature_inj = convert_to_si(25, :Celsius),
    rate = 100kilogram/second/(1000kilogram/meter^3),
    balanced_injection = true,
    num_years = 20,
    model_args = NamedTuple(),
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
    xw, zw = [], []
    for wc in well_coords
        push!(xw, mean(wc[:,1]))
        push!(zw, maximum(wc[:,3]))
    end
    xmin, xmax = extrema(xw)
    zmin, zmax = extrema(zw)
    xc = ((xmin + xmax)/2, (zmin + zmax)/2)
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

    cells = [cells[1]]
    neighbors = [0; 0]
    cell_no = 1;
    for (wno, wc) in enumerate(well_coords)
        wno == 1 ? continue : nothing
        println("Adding producer branch $(wno-1)/$(length(well_coords)-1)")
        cells_k = find_well_cells(msh, wc, z, N, is_frac)
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
    open = isapprox.(z[cells], maximum(z[cells]), atol=1.0)
    well_prod.perforations.WI[.!open] .= 0.0
    well_prod.perforations.WI[open] = min.(well_prod.perforations.WI[open], WI_max)

    wells = [well_inj, well_prod]

    model, _ = setup_reservoir_model(
        domain, :geothermal; wells = wells, model_args...)
    bc, state0 = set_dirichlet_bcs(model, pressure_surface = 10atm)

    rho = reservoir_model(model).system.rho_ref[1]
    # Injector: inject cooled water at a rate equal to production
    # TODO: assert that well exists at construction
    if balanced_injection
        rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Producer])
    else
        rate_target = TotalRateTarget(rate)
    end
    ctrl_inj = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_inj)
    # Producer
    # TODO: support multiple periods with different rates
    rate_target = TotalRateTarget(-rate)
    ctrl_prod = ProducerControl(rate_target)
    control = Dict(:Injector => ctrl_inj, :Producer => ctrl_prod)

    forces = setup_reservoir_forces(model, control = control, bc = bc)#, limits = limits)

    dt, forces = make_schedule(
        [forces], [(1,1), (1,1)];
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
    hx_min = ismissing(hx_min) ? well_distance/5 : hx_min
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

    hy_min = ismissing(hy_min) ? fracture_spacing/5 : hy_min
    hy_max = ismissing(hy_max) ? max(offset, y_max - fracture_radius)/5 : hy_max

    fracture_start = y_min + fracture_spacing/2
    fracture_end = y_max - fracture_spacing/2
    fractured_length = fracture_end - fracture_start
    @assert fractured_length > 0 "Fracture spacing too large for well length"
    nfrac = Int(round(fractured_length/fracture_spacing)) + 1
    
    fracture_y = range(fracture_start, stop=fracture_end, length=nfrac)
    fracture_y = sort(vcat(fracture_y, fracture_y .+ fracture_aperture))

    y = [
        y_min - offset,
        y_min - hy_min/2,
        fracture_y...,
        y_max + hy_min/2,
        y_max + offset
    ]

    hy = fill(hy_max, length(y)-1)
    y_mid = diff(y)./2 .+ y[1:end-1]
    is_frac = y_min .< y_mid .< y_max
    hy[is_frac] .= hy_min
    dy = diff(y)
    ix = [1, length(dy)]
    n = ceil.(dy[ix]./hy[ix])
    fix = n .< 5 .&& .! is_frac[ix]
    hy[ix[fix]] .= dy[ix[fix]]./5

    y = first(Fimbul.interpolate_z(y, hy))

    xyz = (x, y, z)
    sizes = map(x->diff(x), xyz)
    dims = Tuple([length(s) for s in sizes])
    msh = CartesianMesh(dims, sizes, origin = minimum.(xyz))

    return msh

end

 """
    get_egs_fracture_data(states, model, well)

Extract flow rates, temperatures, and energy fluxes from individual fractures
connected to the specified well. This function processes simulation results
to compute fracture-level performance metrics.

Args:
    states: Simulation state results for all timesteps
    model: EGS model containing well and fracture geometry
    well: Well symbol to analyze

Returns:
    Dictionary containing fracture-level data arrays:
    - :y: Y-coordinates of fracture locations
    - :Temperature: Temperature evolution in each fracture
    - :MassFlux: Mass flow rate to/from well into each fracture
    - :EnergyFlux: Thermal energy flux to/from well into each fracture

    NOTE: Positive flux indicates flow into the fracture from the well, negative
    indicates flow out of the fracture into the well.
    
"""
function get_egs_fracture_data(states, model, well; geo = missing)

    # Mesh and geometry
    msh = physical_representation(reservoir_model(model).data_domain)
    geo  = ismissing(geo) ? tpfv_geometry(msh) : geo

    ## Extract well perforation data and connectivity information
    perf = model.models[well].data_domain.representation.perforations # Well-fracture connections
    N = model.models[well].domain.representation.neighborship # Cell connectivity matrix
    is_frac = findall(isapprox.(perf.WI, maximum(perf.WI))) # High-productivity fracture connections
    jj = [cell_ijk(msh, c)[2] for c in perf.reservoir[is_frac]] # Y-indices of fracture cells

    ## Initialize data arrays for fracture analysis
    JJ = unique(jj) # Unique fracture Y-locations
    nstep = length(states) # Number of simulation timesteps
    nfrac = length(JJ) # Number of discrete fractures
    Q, Qh, T = zeros(nstep, nfrac), zeros(nstep, nfrac), zeros(nstep, nfrac) # Pre-allocate arrays
    y = geo.cell_centroids[2, perf.reservoir[is_frac]] # Y-coordinates of fracture centers

    ## Process each simulation timestep
    for (sno, state) in enumerate(states)
        ## Extract mass fluxes and enthalpies from well model
        Qn = state[well][:TotalMassFlux] # Mass flux per well segment [kg/s]
        h = state[well][:FluidEnthalpy] # Fluid enthalpy [J/kg]

        ## Apply upwind scheme for enthalpy (use upstream cell value)
        h = [q > 0 ? h[N[1,s]] : h[N[2,s]] for (s,q) in enumerate(Qn)]
        
        ## Compute energy fluxes combining mass and enthalpy
        Qhn = Qn.*h # Energy flux [W]

        ## Calculate net fluxes per fracture segment (flow into fracture)
        Qn = .-diff(Qn)[is_frac] # Net mass flux into fracture [kg/s]
        Qhn = .-diff(Qhn)[is_frac] # Net energy flux into fracture [W]
        Tn = state[well][:Temperature][is_frac] # Temperature in fracture [K]
        h = state[well][:FluidEnthalpy][is_frac] # Enthalpy in fracture [J/kg]

        ## Aggregate data by fracture Y-location
        for (fno, j) in enumerate(JJ)
            ix = jj .== j # Cells belonging to this fracture
            Q[sno, fno] = sum(Qn[ix]) # Total mass flow for this fracture
            Qh[sno, fno] = sum(Qhn[ix]) # Total energy flow for this fracture
            T[sno, fno] = mean(Tn[ix]) # Average temperature for this fracture
        end
    end

    ## Return organized fracture data
    data = Dict()
    data[:y] = y # Spatial coordinates
    data[:Temperature] = T # Temperature evolution [K]
    data[:MassFlux] = Q # Mass flux evolution [kg/s]
    data[:EnergyFlux] = Qh # Energy flux evolution [W]

    return data
end