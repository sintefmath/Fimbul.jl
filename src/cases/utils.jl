atm = si_unit(:atm)
year, day = si_units(:year, :day)
"""
    make_utes_schedule(forces_charge, forces_discharge, forces_rest; <keyword arguments>...)

Construct a schedule for a UTES system with a cycle of charge -- rest --
discharge -- rest.

# Keyword arguments

- `charge_months::Vector{String} = ["June", "July", "August", "September"]`:
  Months in which the system is charged.
- `discharge_months::Vector{String} = ["December", "January", "February",
  "March"]`: Months in which the system is discharged.
- `start_month::Union{Missing, String}`: Month in which the schedule starts.
  Defaults to the first month of charging.
- `num_years::Int`: Number of years the schedule is repeated (starting from
  2025). If provided, keyword argument `years` must be missing.
- `years::Vector{Int}`: Years in which the schedule is repeated. Defaults to
  `2025:num_years`. If provided, keyword argument `num_years` must be missing.
- `report_interval = 14si_unit(:day)`: Interval at which the simulation output
  is reported.
"""
function make_utes_schedule(forces_charge, forces_discharge, forces_rest;
    charge_months::Union{Nothing, Vector{String}} = ["June", "July", "August", "September"],
    discharge_months::Union{Nothing, Vector{String}} = ["December", "January", "February", "March"],
    start_month::Union{Missing, String} = missing,
    num_years = 5,
    years = missing,
    report_interval = 14day
    )

    # ## Process input
    # Validate months
    charge_months = isnothing(charge_months) ? [] : charge_months
    discharge_months = isnothing(discharge_months) ? [] : discharge_months
    @assert intersect(charge_months, discharge_months) == []
        "Charge and discharge months must be disjoint"
    # TODO add more month checks
    # Porcess years/number of cycles
    if ismissing(years)
        years = 2025:2025+num_years-1
    else
        @assert ismissing(num_years) "Please provide either num_years or years"
    end
    @assert all(diff(years) .== 1) "Years must be consecutive"

    # ## Construct schedule
    # Set month order
    dt_vec, forces = Float64[], []
    if ismissing(start_month)
        if !isempty(charge_months)
            start_month = charge_months[1]
        elseif !isempty(discharge_months)
            start_month = discharge_months[1]
        else
            start_month = "January"
        end
    end
    start_monthno = findall(monthname.(1:12) .== start_month)[1]
    month_ix = ((0:11).+start_monthno.-1).%12 .+ 1
    # Set up schedule for each year
    for year in years
        for mno in month_ix
            mname = monthname(mno)
            # Determine report step length
            num_days = daysinmonth(year, mno)
            time = num_days*day
            if  report_interval > time
                @warn "Report intervall $report_interval is larger than "*
                "the length of $mname. Adjusting to $num_days days"
                n_steps = 1
            else
                n_steps = max(Int(round(time/report_interval)), 1)
            end
            dt = fill(time/n_steps, n_steps)
            # Set forces
            if mname in charge_months
                push!(dt_vec, dt...)
                push!(forces, fill(forces_charge, n_steps)...)
            elseif mname in discharge_months
                push!(dt_vec, dt...)
                push!(forces, fill(forces_discharge, n_steps)...)
            else
                push!(dt_vec, dt...)
                push!(forces, fill(forces_rest, n_steps)...)
            end
        end
    end

    return forces, dt_vec

end

function make_production_schedule(forces_high, forces_low, forces_rest;
        high_months::Union{Nothing, Vector{String}} = [
            "January", "February", "March", "April", "May", "June",
            "July", "August", "September", "October", "November", "December"],
        low_months::Union{Nothing, Vector{String}} = nothing,
        kwargs...
    )

    return make_utes_schedule(forces_high, forces_low, forces_rest;
        charge_months = high_months,
        discharge_months = low_months,
        kwargs...
    )

end

function set_dirichlet_bcs(model, subset = :all;
        pressure_surface = 1atm, 
        temperature_surface = convert_to_si(10.0, :Celsius),
        geothermal_gradient = 0.03Kelvin/meter,
    )

    if subset == :all
        subset = [:top, :sides, :bottom]
    elseif subset isa Symbol
        subset = [subset]
    end

    rmodel = reservoir_model(model)
    mesh = physical_representation(rmodel.data_domain)
    geo = tpfv_geometry(mesh)

    rho = reservoir_model(model).system.rho_ref[1]
    dpdz = rho*gravity_constant
    dTdz = geothermal_gradient
    p = z -> pressure_surface .+ dpdz.*z
    T = z -> temperature_surface .+ dTdz*z

    # Set boundary conditions
    z_bdr = geo.boundary_centroids[3, :]
    cells_bdr = geo.boundary_neighbors
    z0 = minimum(z_bdr)

    top = isapprox.(z_bdr, minimum(z_bdr))
    bottom = isapprox.(z_bdr, maximum(z_bdr))
    sides = .!top .&& .!bottom

    cells_bc = Int64[]
    z_bc = Float64[]
    for s in subset
        if s == :top
            ix = top
        elseif s == :bottom
            ix = bottom
        elseif s == :sides
            ix = sides
        else
            @error "Unknown boundary condition subset: $s"
        end
        push!(cells_bc, cells_bdr[ix]...)
        push!(z_bc, z_bdr[ix]...)
    end
    z_hat = z_bc .- z0
    bc = flow_boundary_condition(cells_bc, rmodel.data_domain, p(z_hat), T(z_hat));

    # Set initial conditions
    z_cells = geo.cell_centroids[3, :]
    z_hat = z_cells .- z0
    state0 = setup_reservoir_state(model,
        Pressure = p(z_hat),
        Temperature = T(z_hat)
    );

    return bc, state0, p, T

end