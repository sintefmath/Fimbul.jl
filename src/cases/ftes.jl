function ftes(;
    depths = [0.0, 100.0, 300.0, 400.0],
    permeability_matrix = convert_to_si(1e-4, :darcy),
    porosity_matrix = 0.01,
    density_matrix = 2600.0*kilogram/meter^3,
    thermal_conductivity_matrix = 2.0*watt/(meter*Kelvin),
    rock_heat_capacity_matrix = 900.0*joule/(kilogram*Kelvin),
    wellpark_radius = 17.5meter,
    num_producers = 8,
    fracture_spacing = missing,
    num_fractures = missing,
    fracture_radius = wellpark_radius*2.0,
    aperture_fracture = convert_to_si(1e-4, :meter),
    porosity_fracture = 0.5,
    pressure_surface = 10.0atm,
    temperature_surface = convert_to_si(4.0, :Celsius),
    rate_charge = 15litre/second,
    rate_discharge = rate_charge,
    temperature_charge = convert_to_si(120.0, :Celsius),
    temperature_discharge = convert_to_si(45.0, :Celsius),
    use_bc = true,
    nc_xy = 4,
    nc_z = missing,
    mesh_args = NamedTuple(),
    utes_schedule_args...
    )

    if fracture_radius < wellpark_radius
        @warn "Fracture radius is smaller than wellpark radius " *
        "- no producers will intersect the fracture network."
    end

    r = wellpark_radius
    # define well coordinates
    xw = [[r.*(sin(θ), cos(θ))] for θ in 
    range(0, stop=2π, length=num_producers+1)][1:end-1]
    xw = vcat([[(0.0, 0.0)]], xw)
    d_xy = Fimbul.min_distance(vcat(xw...))

    # define domain outline, at a distance r from wells
    dr = fracture_radius - wellpark_radius
    outline = Fimbul.offset_boundary(vcat(xw...), dr)
    push!(outline, outline[1])
    cc = vcat(xw, [outline])

    # Define number of fractures for each layer
    layer_thickness = diff(depths)
    num_layers = length(layer_thickness)

    has_num_fractures = !ismissing(num_fractures)
    has_fracture_spacing = !ismissing(fracture_spacing)
    if has_num_fractures && has_fracture_spacing
        error("Please provide only num_fractures or fracture_spacing, not both")
    elseif !has_num_fractures && !has_fracture_spacing
        num_fractures = fill(10, num_layers)
        num_fractures[[1, end]] .= 0
    end
    if ismissing(num_fractures)
        if fracture_spacing isa Real
            fracture_spacing = fill(fracture_spacing, num_layers)
            fracture_spacing[[1,end]] = Inf
        end
        num_fractures = Int(round(layer_thickness./fracture_spacing)) - 1
    end
        
    @assert length(num_fractures) == num_layers
    hxy_min = d_xy/nc_xy
    if ismissing(nc_z)
        nc_z = num_fractures.*5
        nc_z[nc_z .== 0] .= 5
    end
    hz = layer_thickness./nc_z # cell heights per 
    interpolation = []
    for (k, nf) in enumerate(num_fractures)
        top, bottom = false, false
        if num_fractures[k] > 0
            push!(interpolation, :nothing)
        else
            top = k > 1 && num_fractures[k-1] > 0
            bottom = k < num_layers && num_fractures[k+1] > 0
        end
        if top && bottom
            push!(interpolation, :both)
        elseif top
            push!(interpolation, :top)
        else
            push!(interpolation, :bottom)
        end
    end
    interpolation = fill(:nothing, num_layers)
    interpolation = [:bottom, :nothing, :top]
    # return num_fractures, depths, interpolation

    # define mesh with fractures (with a fairly large offset)
    mesh, layers, fractures, metrics = horizontal_fractured_mesh(
        cc, depths, num_fractures; 
        hxy_min = hxy_min, hz = hz, interpolation = interpolation, mesh_args...)

    # restrict fractures to domain defined by outline
    geo = tpfv_geometry(mesh)
    x = geo.cell_centroids[1:2, :]
    x = [(x[1],x[2]) for x in eachcol(x)]
    for (i, xi) = enumerate(x)
        in_domain = Fimbul.point_in_polygon(xi, outline)
        fractures[i] = fractures[i] && in_domain
    end

    # define rock and fracture properties
    num_layers = maximum(layers)
    process_property = (name, value) -> begin
        if value isa Real
            value = fill(value, num_layers)
        end
        if length(value) != num_layers
            error("Wrong size of matrix property $name. Please define as a " *
            "scalar, or a vector with one value per layer ($num_layers))")
        end
        value = value[layers]
    end
    perm = process_property(:permeability_matrix, permeability_matrix)
    poro = process_property(:porosity_matrix, porosity_matrix)
    rho_rock = process_property(:density_matrix, density_matrix)
    th_cond = process_property(:thermal_conductivity_matrix,
        thermal_conductivity_matrix)
    Cpr = process_property(:rock_heat_capacity_matrix,
        rock_heat_capacity_matrix)
    perm[fractures] .= aperture_fracture^2/12
    poro = fill(porosity_matrix, num_layers)[layers]
    poro[fractures] .= porosity_fracture

    # construct the reservoir domain
    domain = reservoir_domain(mesh;
        permeability = perm,
        porosity = poro,
        rock_density = rho_rock,
        thermal_conductivity = th_cond,
        rock_heat_capacity = Cpr,
    )

    nl = findfirst(layers .== maximum(layers))
    println(nl)
    well_models = []
    z = geo.cell_centroids[3, :]

    # add vertical injector well at center location (0,0)
    println("Adding injector")
    rcells = find_well_cells(mesh, xw[1], z, metrics.hxy_min, layers)
    w = setup_well(domain, rcells;
        name = :Injector,
        simple_well=false
    )
    push!(well_models, w)

    # add producer well with nontrivial shape covering the other 9 locations
    println("Adding producer")
    rcells = [rcells[1]]
    wcells = [1]
    neighbors = [NaN; NaN]
    cell_no = 1;
    for (i, x) in enumerate(xw)
        (i == 1) ? continue : 
        println("Adding producer branch $(i-1)/$(length(xw)-1)")
        # Find enclosing cells of well trajectory
        branch_rcells = find_well_cells(mesh, x, z, metrics.hxy_min, layers)
        num_cells = length(branch_rcells)
        push!(rcells, branch_rcells...)
        joint_neighbors = [1; cell_no + 1]
        branch_wcells = collect(1:num_cells) .+ cell_no
        push!(wcells, branch_wcells...)
        branch_neighbors = [wcells[1:end-1]'; wcells[2:end]']
        neighbors = hcat(neighbors, joint_neighbors, branch_neighbors)
        cell_no += num_cells

    end
    neighbors = Int64.(neighbors[:, 2:end])
    rcells = Int64.(rcells)

    w = setup_well(domain, rcells;
        N = neighbors,
        perforation_cells_well = wcells,
        well_cell_centers = geo.cell_centroids[:, rcells],
        name = :Producer,
        simple_well=false
    )
    push!(well_models, w)

    # Make reservoir (multi) model
    model = setup_reservoir_model(
        domain, :geothermal,
        thermal = true,
        wells = well_models,
    )

    # define initial state with pressure and temerature gradients, and zero-flow boundary
    # conditions
    bc, state0, p, T = set_dirichlet_bcs(model;
        pressure_surface = pressure_surface,
        temperature_surface = temperature_surface
    )
    bc = use_bc ? bc : nothing

    # ## Set up controls
    # Rate control for supply side
    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0], 
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    z_min = minimum(z[rcells])
    bhp_target = BottomHolePressureTarget(p(z_min))
    ctrl_prod = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()
    # control_charge[:Injector] = ctrl_charge
    # control_discharge[:Inj] = ctrl_discharge
    # control_charge[:Prod] = ctrl_prod
    # control_discharge[:Prod] = ctrl_prod

    for (name, wm) in get_model_wells(model)
        if name == :Injector
            control_charge[name] = ctrl_charge
            control_discharge[name] = ctrl_discharge
        else
            control_charge[name] = ctrl_prod
            control_discharge[name] = ctrl_prod
        end
    end
    forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
    forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc);
    forces_rest = setup_reservoir_forces(model, bc=bc)

    dt, forces, timestamps = Fimbul.make_utes_schedule(
        forces_charge, forces_discharge, forces_rest; utes_schedule_args...)

    info = Dict()
    info[:timestamps] = timestamps

    case = JutulCase(model, dt, forces, state0 = state0, input_data = info)

    return case

end

function find_well_cells(mesh, x, z, hxy_min, layers)

    xy = [x[1][1], x[1][2]]
    x = [xy[1] xy[2] 0.0; xy[1] xy[2] 25.0]
    x[:,1:2] .+= hxy_min/2
    cells = Jutul.find_enclosing_cells(mesh, x, n = 100)
    z_cells = z[cells]
    cell = cells[findall(z_cells .== minimum(z_cells))[1]]

    num_layers = maximum(layers)
    nl = findfirst(layers .== num_layers)
    cells = collect((1:nl) .+ (cell-1))
    keep = layers[cells] .< num_layers
    cells = cells[keep]

    return cells

end

# function setup_reservoir_simulator_skattoera(case, name = nothing;
#     subfolder = missing,
#     use_convergence_monitor_relaxation = true,
#     kwargs...)

#     if ismissing(subfolder)
#         subfolder = joinpath("projects", "ghost-digit")
#     end
#     if !isnothing(name)
#         output_path = jutul_output_path(name, subfolder=subfolder)
#     else
#         output_path = nothing
#     end

#     simulator, config = setup_reservoir_simulator(case;
#         method=:newton,
#         presolve_wells = true,
#         nldd_partition=missing,
#         tol_cnv=1e-2,
#         tol_mb=1e-6,
#         output_path=output_path,
#         timesteps=:auto,
#         max_nonlinear_iterations=8,
#         target_its=8,
#         initial_dt=convert_to_si(5.0, :minute),
#         relaxation=true,
#         info_level=2,
#         kwargs...
#     );
#     # Timestep selector adjusting timestep based on control changes
#     sel = JutulDarcy.ControlChangeTimestepSelector(case.model,
#         0.1, convert_to_si(1.0, :minute))
#     push!(config[:timestep_selectors], sel)
#     config[:timestep_max_decrease] = 1e-10
#     if use_convergence_monitor_relaxation
#         Jutul.ConvergenceMonitors.set_convergence_monitor_relaxation!(config)
#     end

#     for ws in well_symbols(case.model)
#         config[:tolerances][ws][:default] = 1e-2
#     end
#     config[:tolerances][:Reservoir][:default] = 1e-2

#     return simulator, config

# en