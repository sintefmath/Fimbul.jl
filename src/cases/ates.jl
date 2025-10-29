"""
    ates(; <keyword arguments>)

Create an Aquifer Thermal Energy Storage (ATES) simulation case.

# Keyword Arguments

## Geometry Parameters
- `well_distance = missing`: Distance between hot and cold wells [m]. If missing, 
  calculated as 2×thermal_radius
- `depths = [0.0, 850.0, 900.0, 1000.0, 1050.0, 1300.0]`: Layer interface depths [m]
- `aquifer_layer = 3`: Index of the aquifer layer for thermal storage
- `use_2d = false`: If true, creates a 2D model

## Rock Properties (per layer)
- `porosity = [0.01, 0.05, 0.35, 0.05, 0.01]`: Porosity for each layer [-]
- `permeability`: Permeability for each layer [m²]
- `rock_thermal_conductivity`: Thermal conductivity [W/(m·K)]
- `rock_heat_capacity`: Rock heat capacity [J/(kg·K)]

## Thermal Parameters
- `temperature_surface = 10°C`: Surface temperature [K]
- `thermal_gradient = 0.03`: Geothermal gradient [K/m]

## Operational Parameters
- `temperature_charge = 95°C`: Temperature of injected water during charging [K]
- `temperature_discharge = 25°C`: Temperature of injected water during discharging [K]
- `rate_charge = missing`: Injection/production rate during charging [m³/s]. If
  missing, calculated based on a target thermal_radius = well_distance/2 or, if
  well_distance missing, 250 m
- `rate_discharge = rate_charge`: Rate during discharging [m³/s]
- `balanced_injection = true`: If true, the injection well target will be set
  equal to the production rate

## Advanced Options
- `utes_schedule_args = NamedTuple()`: Additional arguments passed to
  make_utes_schedule
- `mesh_args = NamedTuple()`: Additional arguments passed to make_ates_cart_mesh

# Returns
A `JutulCase` for ATES
"""
function ates(;
    well_distance = missing,
    depths = [0.0, 850.0, 900.0, 1000.0, 1050.0, 1300.0],
    porosity = [0.01, 0.05, 0.35, 0.05, 0.01],
    permeability = [1.0, 5.0, 1000.0, 5.0, 1.0].*1e-3.*si_unit(:darcy),
    rock_thermal_conductivity = [2.5, 2.0, 1.5, 2.0, 2.5].*watt/(meter*Kelvin),
    rock_heat_capacity = 900.0*joule/(kilogram*Kelvin),
    aquifer_layer = 3,
    temperature_charge = convert_to_si(95, :Celsius),
    temperature_discharge = convert_to_si(25, :Celsius),
    rate_charge = missing,
    rate_discharge = rate_charge,
    balanced_injection = true,
    temperature_surface = convert_to_si(10.0, :Celsius),
    thermal_gradient = 0.03Kelvin/meter,
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    utes_schedule_args = NamedTuple(),
    use_2d = false,
    mesh_args = NamedTuple(),
)

    # ## Process input and set defaults
    # Fill layer properties given as scalars
    layer_thickness = diff(depths)
    process_value = v -> v isa Number ? fill(v, length(depths)-1) : v
    porosity = process_value(porosity)
    permeability = process_value(permeability)
    rock_thermal_conductivity = process_value(rock_thermal_conductivity)
    rock_heat_capacity = process_value(rock_heat_capacity)
    # Calculate thermal properties for the aquifer layer
    Cf, Cr = 4184.0, rock_heat_capacity[aquifer_layer]
    ϕ = porosity[aquifer_layer]
    Caq = Cf*ϕ + Cr*(1 - ϕ)
    Haq = layer_thickness[aquifer_layer]
    # Default charging duration (6 months)
    # TODO: Use charge duration from input
    charge_duration = 0.5si_unit(:year)
    ch_start = Fimbul.process_time(2025, charge_period[1])
    ch_end = Fimbul.process_time(2025, charge_period[2], true)
    charge_duration = (ch_end - ch_start).value*1e-3
    thermal_radius = missing
    # Calculate rate based on a target 250 m thermal radius if not provided
    if ismissing(rate_charge)
        thermal_radius = ismissing(well_distance) ? 250.0 : well_distance/2
        rate_charge = (thermal_radius^2*Caq*π*Haq/Cf)/charge_duration
        # Adjust for 2D simulation (per unit width)
        if use_2d
            rate_charge *= 2/(π*thermal_radius)
        end
    end
    # Set discharge rate equal to charge rate if not specified
    if ismissing(rate_discharge)
        rate_discharge = rate_charge
    end
    Vin = rate_charge*charge_duration
    # Calculate thermal radius from injected volume if not already set above
    if ismissing(thermal_radius)
        thermal_radius = thermal_radius_aquifer(Vin, Haq, ϕ, Cf, Cr)
    end
    # Set well distance based on thermal radius if not specified
    if ismissing(well_distance)
        well_distance = 2*thermal_radius
    end

    # ## Create reservoir domain
    # Create mesh
    nearwell_radius = 0.5*min(thermal_radius, well_distance/2)
    msh, layers = make_ates_cart_mesh(well_distance, depths, aquifer_layer;
        nearwell_radius = nearwell_radius,
        use_2d = use_2d,
        mesh_args...
    )
    # Set properties and create domain
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

    # ## Setup wells
    k = cell_ijk(msh, findlast(layers .== aquifer_layer))[3]
    # Setup hot well
    xw_hot  = [-well_distance/2 0.0 0.0; -well_distance/2 0.0 depths[end]]
    cell = Jutul.find_enclosing_cells(msh, xw_hot)[1]
    ij = cell_ijk(msh, cell)[1:2]
    hot_well = setup_vertical_well(domain, ij[1], ij[2]; toe=k, simple_well=false, name = :Hot)
    # Only perforate in the aquifer layer
    hot_rcells = hot_well.representation.perforations.reservoir
    WI = [compute_peaceman_index(msh, permeability[c], 0.1, c) for c in hot_rcells]
    WI[layers[hot_rcells] .!== aquifer_layer] .= 0.0
    hot_well[:well_indices] = WI
    # Setup cold well
    xw_cold  = [well_distance/2 0.0 0.0; well_distance/2 0.0 depths[end]]
    cell = Jutul.find_enclosing_cells(msh, xw_cold)[1]
    ij = cell_ijk(msh, cell)[1:2]
    cold_well = setup_vertical_well(domain, ij[1], ij[2]; toe=k, simple_well=false, name = :Cold)
    # Only perforate in the aquifer layer
    cold_rcells = cold_well.representation.perforations.reservoir
    WI = [compute_peaceman_index(msh, permeability[c], 0.1, c) for c in cold_rcells]
    WI[layers[cold_rcells] .!== aquifer_layer] .= 0.0
    cold_well[:well_indices] = WI

    # ## Setup reservoir model
    model = setup_reservoir_model(
        domain, :geothermal;
        wells = [hot_well, cold_well]
    )

    # ## Set boundary and initial conditions
    # Get fluid density for well controls
    rho = reservoir_model(model).system.rho_ref[1]
    # Setup geometry and calculate hydrostatic/geothermal profiles
    rmodel = reservoir_model(model)
    geo = tpfv_geometry(msh)
    # Hydrostatic pressure gradient
    dpdz = rho*gravity_constant
    # Geothermal temperature gradient
    dTdz = thermal_gradient
    # Pressure and temperature as functions of depth
    p = z -> 10atm .+ dpdz.*z
    T = z -> temperature_surface .+ dTdz*z
    # Set boundary conditions
    z_bdr = geo.boundary_centroids[3, :]
    xy_bdr = geo.boundary_centroids[1:2, :]
    cells_bdr = geo.boundary_neighbors
    z0 = minimum(z_bdr)
    top = isapprox.(z_bdr, z0)
    bottom = isapprox.(z_bdr, maximum(z_bdr))
    west = isapprox.(xy_bdr[1, :], minimum(xy_bdr[1, :]))
    east = isapprox.(xy_bdr[1, :], maximum(xy_bdr[1, :]))
    south = isapprox.(xy_bdr[2, :], minimum(xy_bdr[2, :]))
    north = isapprox.(xy_bdr[2, :], maximum(xy_bdr[2, :]))
    # Select boundary faces based on setup
    cells_bc, z_bc = Int64[], Float64[]
    if use_2d
        # For 2D: apply BCs on top, bottom, west (left), east (right)
        subset = top .|| bottom .|| west .|| east
    else
        # For 3D: apply BCs on all external boundaries
        subset = top .|| bottom .|| south .|| north .|| west .|| east
    end
    cells_bc = cells_bdr[subset]
    z_bc = z_bdr[subset]
    z_hat = z_bc .- z0
    # Create flow boundary conditions
    bc = flow_boundary_condition(cells_bc, rmodel.data_domain, p(z_hat), T(z_hat));
    # Setup initial state with hydrostatic pressure and geothermal temperature
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- z0
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );

    # ## Set up ATES schedule
    # Set up well controls for charging phase
    if balanced_injection
        rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Cold])
    else
        rate_target = TotalRateTarget(rate_charge)
    end
    # Hot well: inject hot water
    ctrl_hot = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_charge)
    # Cold well: produce water
    rate_target = TotalRateTarget(-rate_charge)
    ctrl_cold = ProducerControl(rate_target)
    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)
    # Set up forces
    forces_charge = setup_reservoir_forces(model; bc = bc, control = control)
    # Set up well controls for discharging phase (roles reversed)
    # Hot well: produce hot water
    rate_target = TotalRateTarget(-rate_discharge)
    ctrl_hot  = ProducerControl(rate_target)
    # Cold well: inject cold water
    if balanced_injection
        rate_target = JutulDarcy.ReinjectionTarget(NaN, [:Hot])
    else
        rate_target = TotalRateTarget(rate_discharge)
    end
    ctrl_cold = InjectorControl(
        rate_target, [1.0], density = rho, temperature = temperature_discharge)
    # Set up forces
    control = Dict(:Hot => ctrl_hot, :Cold => ctrl_cold)
    forces_discharge = setup_reservoir_forces(model; bc = bc, control = control)
    # Set up forces for rest periods (no active wells)
    forces_rest = setup_reservoir_forces(model; bc = bc)
    # Create UTES operational schedule
    dt, forces, timestamps = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest,
        charge_period = charge_period,
        discharge_period = discharge_period;
        utes_schedule_args...
    )

    # ## Additional info
    info = Dict()
    info[:description] = "Aquifer Thermal Energy Storage (ATES) case "*
    "generated using Fimbul.ates()"
    info[:well_distance] = well_distance
    info[:thermal_radius] = thermal_radius
    info[:layers] = layers
    info[:aquifer_layer] = aquifer_layer
    info[:timestamps] = timestamps

    # ## Create and return the complete simulation case
    case = JutulCase(model, dt, forces; state0 = state0, input_data = info)

    return case

end

"""
    ates_simple(; kwargs...)

Create a simplified ATES case with only aquifer and cap rock/basement layers.

This is a convenience function for quick ATES simulations with minimal setup.
It creates a three-layer system: cap rock, aquifer, and basement.

# Keyword Arguments
- `well_distance = 500.0`: Distance between wells [m]
- `aquifer_thickness = 100.0`: Thickness of aquifer layer [m]
- `depth = 1000.0`: Depth to top of aquifer [m]
- `porosity = [0.2, 0.01]`: [aquifer, cap rock/basement] porosity [-]
- `permeability`: [aquifer, cap rock/basement] permeability [m²]
- `thermal_conductivity`: [aquifer, cap rock/basement] thermal conductivity [W/(m·K)]
- `rock_heat_capacity`: [aquifer, cap rock/basement] heat capacity [J/(kg·K)]
- `kwargs...`: Additional arguments passed to `ates()`

# Returns
A `JutulCase` object for the simplified ATES system.

# Example
```julia
case = ates_simple(
    well_distance = 400.0,
    aquifer_thickness = 50.0,
    depth = 800.0
)
```
"""
function ates_simple(;
    well_distance = missing,
    aquifer_thickness = 100.0,
    depth = 1000.0,
    porosity = [0.35, 0.05],
    permeability = [1000.0, 5.0].*1e-3.*si_unit(:darcy),
    rock_thermal_conductivity = [2.0, 2.0].*watt/(meter*Kelvin),
    rock_heat_capacity = [900.0, 900.0]*joule/(kilogram*Kelvin),
    kwargs...
)

    # Middle layer is the aquifer
    aquifer_layer = 2  
    # Create simple three-layer system: cap rock - aquifer - basement
    # TODO: Add back-of-the-envelope calculation for thermal plume spread
    # downwards to set basement thickness
    depths = [0.0, depth, depth + aquifer_thickness, (depth + aquifer_thickness)*1.25]
    # Map properties to three layers
    make_prop = prop -> [prop[2], prop[1], prop[2]]  # [cap, aquifer, basement]
    porosity = make_prop(porosity)
    permeability = make_prop(permeability)
    rock_thermal_conductivity = make_prop(rock_thermal_conductivity)
    rock_heat_capacity = make_prop(rock_heat_capacity)
    # Call main ates function
    return Fimbul.ates(;
        well_distance=well_distance,
        depths=depths,
        porosity=porosity,
        permeability=permeability,
        rock_thermal_conductivity=rock_thermal_conductivity,
        rock_heat_capacity=rock_heat_capacity,
        aquifer_layer=aquifer_layer,
        kwargs...
    )

end

"""
    make_ates_cart_mesh(well_distance, depths, aquifer_layer; kwargs...)

Generate a Cartesian mesh optimized for ATES simulations.

Creates a structured mesh with refined gridding around wells and in the aquifer
layer to accurately capture thermal transport processes.

# Arguments
- `well_distance`: Distance between hot and cold wells [m]
- `depths`: Array of layer interface depths [m]
- `aquifer_layer`: Index of the aquifer layer for refined gridding

# Keyword Arguments
- `hxy_min = missing`: Minimum horizontal grid spacing [m]
- `hxy_max = missing`: Maximum horizontal grid spacing [m]  
- `hz_min = missing`: Minimum vertical grid spacing [m]
- `hz_max = missing`: Maximum vertical grid spacing [m]
- `thermal_radius = missing`: Radius of thermal influence for grid refinement [m]
- `offset_rel = 1.0`: Relative offset for domain boundaries
- `use_2d = false`: Create 2D mesh if true

# Returns
- `msh`: CartesianMesh object
- `layers`: Array mapping cells to layer indices

# Grid Design
The mesh uses graded refinement:
- Fine grid around wells and within thermal radius
- Coarse grid at domain boundaries
- Refined vertical resolution in and around aquifer layer
"""
function make_ates_cart_mesh(well_distance, depths, aquifer_layer;
    hxy_min = missing,
    hxy_max = missing,
    hz_min = missing,
    hz_max = missing,
    nearwell_radius = missing,
    offset_rel = 1.0,
    use_2d = false,
)

    # ## Process input
    # Set default thermal radius if not provided
    nearwell_radius = ismissing(nearwell_radius) ? well_distance/2 : nearwell_radius
    # Calculate domain offset to minimize boundary effects
    offset = offset_rel*maximum([well_distance, 2*nearwell_radius])
    # Well coordinates
    xy_hot = (-well_distance/2, 0.0)
    xy_cold = (well_distance/2, 0.0)
    # Set default horizontal grid spacing
    hxy_min = ismissing(hxy_min) ? nearwell_radius/6 : hxy_min
    hxy_max = ismissing(hxy_max) ? (offset - nearwell_radius)/5 : hxy_max
    @assert hxy_max ≥ hxy_min "hxy_max must be greater than or equal to hxy_min"
    # Set vertical grid spacing
    layer_thickness = diff(depths)
    hz_min = ismissing(hz_min) ? layer_thickness[aquifer_layer]/10 : hz_min
    hz_max = ismissing(hz_max) ? maximum(layer_thickness)/5 : hz_max
    @assert hz_max ≥ hz_min "hz_max must be greater than or equal to hz_min"

    # ## Create coordinate arrays with refinement
    # Create coordinates for x dimension
    x = [
        xy_hot[1]-offset,          # Left boundary
        xy_hot[1]-nearwell_radius, # Left thermal boundary
        xy_hot[1]-hxy_min/2,       # Near hot well
        xy_hot[1]+nearwell_radius, # Between wells
        xy_cold[1]-nearwell_radius,# Between wells  
        xy_cold[1]-hxy_min/2,      # Near cold well
        xy_cold[1]+nearwell_radius,# Right thermal boundary
        xy_cold[1]+offset          # Right boundary
    ]
    # Grid spacing array: coarse-fine-fine-coarse-fine-fine-coarse
    hx = [hxy_max, hxy_min, hxy_min, hxy_max, hxy_min, hxy_min, hxy_max];
    dx = diff(x)
    hx[dx .< hxy_max] .= hxy_min
    nx = round.(dx./hx)
    hx[nx .< 5] .= dx[nx .< 5]./5
    x, _ = Fimbul.interpolate_z(x, hx)
    # Create coordinates for y dimension
    if !use_2d
        y = [
            xy_hot[2]-offset,
            xy_hot[2]-nearwell_radius,
            xy_hot[2]-hxy_min/2,
            xy_hot[2]+nearwell_radius,
            xy_hot[2]+offset
        ]
        hy = [hxy_max, hxy_min, hxy_min, hxy_max]
        hy[diff(y) .< hxy_max] .= hxy_min
        y, _ = Fimbul.interpolate_z(y, hy)
    else
        y = [-0.5, 0.5]
    end
    # Create coordinates in vertical dimension
    hz = fill(hz_max, length(layer_thickness))
    hz[aquifer_layer] = hz_min
    hz[layer_thickness .< hz_max] .= hz_min
    nz = round.(layer_thickness./hz)
    hz[nz .< 5] .= layer_thickness[nz .< 5]./5
    z, layers = Fimbul.interpolate_z(depths, hz)
    # ## Create Cartesian mesh
    xyz = (x, y, z)
    sizes = map(x->diff(x), xyz)
    dims = Tuple([length(s) for s in sizes])
    msh = CartesianMesh(dims, sizes, origin = minimum.(xyz))
    # Map layers to all cells (repeat for each horizontal cell)
    layers = repeat(layers, inner = dims[1]*dims[2])

    return msh, layers
end