to_kelvin = T -> convert_to_si(T, :Celsius)
second, year, day = si_units(:second, :year, :day)
meter = si_unit(:meter)
litre = si_unit(:litre)
Kelvin = si_unit(:Kelvin)
darcy = si_unit(:darcy)

"""
    btes(; <keyword arguments>)

Setup function for borehole thermal energy storage (BTES) system.

# Keyword arguments
- `num_wells = 48`: Number of wells in the BTES system.
- `num_sections = 6`: Number of sections in the BTES system. The system is
  divided into equal circle sectors, and all wells in each sector are coupled in
  series.
- `well_spacing = 5.0`: Horizontal spacing between wells in meters.
- `depths = [0.0, 0.5, 50, 65]`: Depths delineating geological layers in meters.
- `well_layers = [1, 2]`: Layers in which the wells are placed
- `density = [30, 2580, 2580]*kilogram/meter^3`: Rock density in the layers.
- `thermal_conductivity = [0.034, 3.7, 3.7]*watt/meter/Kelvin`: Thermal
  conductivity in the layers.
- `heat_capacity = [1500, 900, 900]*joule/kilogram/Kelvin`: Heat capacity in the
  layers.
- `geothermal_gradient = 0.03Kelvin/meter`: Geothermal gradient.
- `temperature_charge = to_kelvin(90.0)`: Injection temperature during charging.
- `temperature_discharge = to_kelvin(10.0)`: Injection temperature during
  discharging.
- `rate_charge = 0.5litre/second`: Injection rate during charging.
- `rate_discharge = rate_charge`: Injection rate during discharging.
- `temperature_surface = to_kelvin(10.0)`: Temperature at the surface.
- `num_years = 5`: Number of years to run the simulation.
- `charge_period = ["June", "September"]`: Period during which
  the system is charged.
- `discharge_period = ["December", "March"]`: Period during which
  the system is discharged.
- `report_interval = 14day`: Reporting interval for the simulation.
- `utes_schedule_args = NamedTuple()`: Additional arguments for the UTES schedule.
- `n_z = [3, 8, 3]`: Number of layers in the vertical direction for each layer.
- `n_xy = 3`: Number of layers in the horizontal direction for each layer.
- `mesh_args = NamedTuple()`: Additional arguments for the mesh generation.
"""
function btes(;
    num_wells = 48,
    num_sectors = 6,
    well_spacing = 5.0,
    depths = [0.0, 0.5, 50, 65],
    well_layers = [1, 2],
    density = [30, 2580, 2580]*kilogram/meter^3,
    thermal_conductivity = [0.034, 3.7, 3.7]*watt/meter/Kelvin,
    heat_capacity = [1500, 900, 900]*joule/kilogram/Kelvin,
    geothermal_gradient = 0.03Kelvin/meter,
    temperature_charge = to_kelvin(90.0),
    temperature_discharge = to_kelvin(10.0),
    rate_charge = 0.5litre/second,
    rate_discharge = rate_charge,
    temperature_surface = to_kelvin(10.0),
    num_years = 4,
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    report_interval = 14day,
    utes_schedule_args = NamedTuple(),
    n_z = [3, 8, 3],
    n_xy = 3,
    mesh_args = NamedTuple(),
    )

    # ## Create mesh
    x = fibonacci_pattern_2d(num_wells; spacing = well_spacing)
    well_coordinates = map(x -> [x], x)
    hz = diff(depths)./n_z
    hxy = well_spacing/n_xy
    mesh, layers, metrics = extruded_mesh(well_coordinates, depths;
        hxy_min = hxy, hz = hz, mesh_args...)

    # ## Set up model
    # Set up reservoir domain with rock properties similar to that of granite,
    # with a styrofoam layer on top
    density = density[layers]
    thermal_conductivity = thermal_conductivity[layers]
    heat_capacity = heat_capacity[layers]
    domain = reservoir_domain(mesh,
        permeability = 1e-6darcy,
        porosity = 0.01,
        rock_density = density,
        rock_heat_capacity = heat_capacity,
        rock_thermal_conductivity = thermal_conductivity,
        component_heat_capacity = 4.278e3joule/kilogram/Kelvin,
    )
    # Set up BTES wells
    hxy_min = metrics.hxy_min
    well_models = []
    nl = length(layers)
    geo = tpfv_geometry(mesh)
    for (wno, xw) in enumerate(well_coordinates)
        name = Symbol("B$wno")
        println("Adding well $name ($wno/$num_wells)")
        xw = xw[1]
        d = max(norm(xw, 2))
        v = (d > 0) ? xw./d : (1.0, 0.0)
        # Shift coordiates a bit to avoid being exactly on the node
        xw = xw .+ (hxy_min/2) .* v
        trajectory = [xw[1] xw[2] 0.0; xw[1] xw[2] depths[end]]
        cells = Jutul.find_enclosing_cells(mesh, trajectory, n = 100)
        filter!(c -> layers[c] ∈ well_layers, cells)
        w_sup, w_ret = setup_btes_well(domain, cells, name=name, btes_type=:u1)
        push!(well_models, w_sup, w_ret)
    end
    # Make the model
    model, parameters = setup_reservoir_model(
        domain, :geothermal,
        wells = well_models,
    );

    # ## Set up initial state and boundary conditions
    geo = tpfv_geometry(mesh)
    z_bc = geo.boundary_centroids[3, :]
    bottom = map(v -> isapprox(v, maximum(z_bc)), z_bc)
    # Define pressure and temperature profiles
    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = geothermal_gradient
    T = z -> temperature_surface .+ dTdz*z
    p = z -> 5atm .+ dpdz.*z
    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- minimum(z_cells)
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );
    # Set boundary conditions
    z_bc = z_bc[.!bottom]
    z_hat = z_bc .- minimum(z_bc)
    bc_cells = geo.boundary_neighbors[.!bottom]
    bc = flow_boundary_condition(bc_cells, domain, p(z_hat), T(z_hat));

    control_charge, control_discharge, sections = setup_controls(model, num_sectors,
        rate_charge, rate_discharge, temperature_charge, temperature_discharge);
    forces_charge = setup_reservoir_forces(model, control=control_charge, bc=bc)
    forces_discharge = setup_reservoir_forces(model, control=control_discharge, bc=bc);
    forces_rest = setup_reservoir_forces(model, bc=bc)
    # Make schedule
    dt, forces = make_utes_schedule(
        forces_charge, forces_discharge, forces_rest;
        charge_period = charge_period,
        discharge_period = discharge_period,
        num_years = num_years,
        report_interval = report_interval,
        utes_schedule_args...,
    )

    # ## Assemble and return model
    case = JutulCase(model, dt, forces, state0 = state0)
    return case, sections

end

function setup_controls(model, number_of_sectors, 
    rate_charge, rate_discharge, temperature_charge, temperature_discharge)

    rho = reservoir_model(model).system.rho_ref[1]
    rate_target = TotalRateTarget(rate_charge)
    ctrl_charge = InjectorControl(rate_target, [1.0], 
        density=rho, temperature=temperature_charge)
    rate_target = TotalRateTarget(rate_discharge)
    ctrl_discharge = InjectorControl(rate_target, [1.0],
        density=rho, temperature=temperature_discharge);
    # BHP control for return side
    bhp_target = BottomHolePressureTarget(1.0si_unit(:atm))
    ctrl_ret = ProducerControl(bhp_target);
    # Set up forces
    control_charge = Dict()
    control_discharge = Dict()

    msh = physical_representation(reservoir_model(model).data_domain)
    geo = tpfv_geometry(msh)
    xy = geo.cell_centroids[1:2,:]
    # map from (x,y) to polar coordiates
    r = sqrt.(xy[1,:].^2 .+ xy[2,:].^2)
    θ = atan.(xy[2,:], xy[1,:]) .+ π
    assigned = []
    wells = well_symbols(model)
    filter!(well -> contains(String(well), "_supply"), wells)
    get_return = (well) -> Symbol(replace(String(well), "_supply" => "_return"))
    sectors = Dict()
    wells_per_sector = div(length(wells), number_of_sectors)
    rem = length(wells) - wells_per_sector*number_of_sectors
    wells_per_sector = fill(wells_per_sector, number_of_sectors)
    wells_per_sector[1:rem] .+= 1
    wmodels = [model.models[well] for well in wells]
    wc = [wm.domain.representation.perforations.reservoir[1] for wm in wmodels]
    θ = θ[wc]
    order_θ = sortperm(θ)
    wtot = 0
    for sno in 1:number_of_sectors
        wno = (1:wells_per_sector[sno]) .+ wtot
        wtot += wells_per_sector[sno]
        sw = wells[order_θ[wno]]
        well_radii = r[wc[order_θ[wno]]]
        order = sortperm(well_radii)
        sec_wells = Symbol[]
        for (k, wno) in enumerate(order)
            well_sup = sw[wno]
            well_ret = get_return(well_sup)
            @assert well_sup ∉ assigned
            @assert well_ret ∉ assigned
            if k == 1
                # Water is injected into innermost well during charging
                control_charge[well_sup] = ctrl_charge
                # Discharging runs from outer to inner
                well_prev = get_return(sw[order[k+1]])
                target = JutulDarcy.ReinjectionTarget(NaN, [well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            elseif k == length(sw)
                # Water is injected into outermost well during discharging
                control_discharge[well_sup] = ctrl_discharge
                # Charging runs from inner to outer
                well_prev = get_return(sw[order[k-1]])
                target = JutulDarcy.ReinjectionTarget(NaN, [well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
            else
                # Charging runs from inner to outer
                well_prev = get_return(sw[order[k-1]])
                target = JutulDarcy.ReinjectionTarget(NaN, [well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_charge[well_sup] = ctrl
                # Discharging runs from outer to inner
                well_prev = get_return(sw[order[k+1]])
                target = JutulDarcy.ReinjectionTarget(NaN, [well_prev])
                ctrl = InjectorControl(target, [1.0],
                    density=rho, temperature=NaN; check=false)
                control_discharge[well_sup] = ctrl
            end
            control_charge[well_ret] = ctrl_ret
            control_discharge[well_ret] = ctrl_ret
            push!(assigned, well_sup, well_ret)
            push!(sec_wells, well_sup, well_ret)
        end
        sectors[Symbol("S$sno")] = sec_wells
    end

    @assert sort(assigned) == sort(well_symbols(model))

    return control_charge, control_discharge, sectors

end