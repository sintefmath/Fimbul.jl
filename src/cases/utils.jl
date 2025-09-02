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
# function make_utes_schedule(forces_charge, forces_discharge, forces_rest;
#     charge_period::Union{Nothing, Vector{String}} = [6, 9],
#     discharge_period::Union{Nothing, Vector{String}} = [12, 3],
#     # charge_months::Union{Nothing, Vector{String}} = ["June", "July", "August", "September"],
#     # discharge_months::Union{Nothing, Vector{String}} = ["December", "January", "February", "March"],
#     start_month::Union{Missing, String} = missing,
#     num_years = 5,
#     years = missing,
#     report_interval = 14day
#     )

#     # Process years/number of cycles
#     if ismissing(years)
#         years = 2025:2025+num_years-1
#     else
#         @assert ismissing(num_years) "Please provide either num_years or years"
#     end
#     @assert all(diff(years) .== 1) "Years must be consecutive"

#     start_year = years[1]
    
#     start_month = monthname(start_monthno)

#     # ## Process input
#     # Validate months
#     charge_periods = isnothing(charge_periods) ? [] : charge_periods
#     discharge_periods = isnothing(discharge_periods) ? [] : discharge_periods

#     @assert intersect(charge_periods, discharge_periods) == []
#         "Charge and discharge months must be disjoint"
#     # TODO add more month checks
  

#     # ## Construct schedule
#     # Set month order
#     dt_vec, forces = Float64[], []
#     if ismissing(start_month)
#         if !isempty(charge_months)
#             start_month = charge_months[1]
#         elseif !isempty(discharge_months)
#             start_month = discharge_months[1]
#         else
#             start_month = "January"
#         end
#     end
#     start_monthno = findall(monthname.(1:12) .== start_month)[1]
#     month_ix = ((0:11).+start_monthno.-1).%12 .+ 1
#     # Set up schedule for each year
#     for year in years
#         for mno in month_ix
#             mname = monthname(mno)
#             # Determine report step length
#             num_days = daysinmonth(year, mno)
#             time = num_days*day
#             if  report_interval > time
#                 @warn "Report intervall $report_interval is larger than "*
#                 "the length of $mname. Adjusting to $num_days days"
#                 n_steps = 1
#             else
#                 n_steps = max(Int(round(time/report_interval)), 1)
#             end
#             dt = fill(time/n_steps, n_steps)
#             # Set forces
#             if mname in charge_months
#                 push!(dt_vec, dt...)
#                 push!(forces, fill(forces_charge, n_steps)...)
#             elseif mname in discharge_months
#                 push!(dt_vec, dt...)
#                 push!(forces, fill(forces_discharge, n_steps)...)
#             else
#                 push!(dt_vec, dt...)
#                 push!(forces, fill(forces_rest, n_steps)...)
#             end
#         end
#     end

#     return forces, dt_vec

# end

function make_schedule(forces, periods;
    start_year = missing,
    num_years = 1,
    report_interval = 14day
    )

    start_year = ismissing(start_year) ? Dates.year(now()) : start_year
    years = (0:num_years-1) .+ start_year

    num_periods = length(periods)-1
    @assert num_periods > 0 "At least one period is required"
    @assert length(forces) == num_periods "Number of forces must match number of periods"

    if length(report_interval) == 1
        report_interval = fill(report_interval, num_periods)
    end

    dt, force_vec, timestamps = Float64[], [], DateTime[]
    for year in years
        periods_year = Fimbul.process_periods(year, periods)
        for k in eachindex(periods_year)[1:end-1]
            start_time, end_time = periods_year[k], periods_year[k+1]
            Δt = Dates.value(end_time - start_time)*1e-3
            dt_k = report_interval[k]
            n_step = max(1, Int(round(Δt/dt_k)))
            dt_k = Δt/n_step
            dt_k, forces_k = fill(dt_k, n_step), fill(forces[k], n_step)
            dts = Dates.Millisecond(Int(round(Δt/n_step*1e3)))
            ts = start_time:dts:(end_time - dts)
            push!(dt, dt_k...)
            push!(force_vec, forces_k...)
            push!(timestamps, ts...)
        end
        if year == years[end]
            push!(timestamps, periods_year[end])
        end
    end

    return dt, force_vec, timestamps

end

function make_utes_schedule(forces_charge, forces_discharge, forces_rest;
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"],
    start_year = missing,
    num_years = 5,
    kwargs...
    )

    start_year = ismissing(start_year) ? Dates.year(now()) : start_year
    
    ch_start = process_time(start_year, charge_period[1])
    ch_end = process_time(start_year, charge_period[2], true)
    dch_start = process_time(start_year, discharge_period[1])
    dch_end = process_time(start_year, discharge_period[2], true)

    periods = [ch_start, ch_end, dch_start, dch_end, ch_start + Dates.Year(1)]
    periods = process_periods(start_year, periods)
    keep = diff(periods) .> Dates.Second(0)
    println("Periods: $periods")
    println("Keep: $keep")

    periods = periods[[keep; true]]

    periods = [(Dates.month(p), Dates.day(p), Dates.hour(p)) for p in periods]
    forces = [forces_charge; forces_rest; forces_discharge; forces_rest][keep]

    dt, forces, timestamps = make_schedule(forces, periods;
        start_year = start_year, num_years = num_years, kwargs...)

    return dt, forces, timestamps

end

function process_periods(year, periods)

    println("Periods: $periods")
    periods_proc = Vector{DateTime}(undef, length(periods))
    for (pno, period) in enumerate(periods)
        time = process_time(year, period)
        if pno > 1
            time < periods_proc[pno-1] ? time += Dates.Year(1) : nothing
        end
        periods_proc[pno] = time
    end

    return periods_proc

end

function process_time(year, time, is_end=false)

    if time isa DateTime
        return time
    elseif time isa String
        time = Dates.monthname_to_value(time, Dates.ENGLISH)
    end
    n = length(time)

    @assert time isa Int || time isa Tuple{Int,Int} || time isa Tuple{Int,Int,Int} "Time must be an Int (month), Tuple{Int,Int} (month, day) or Tuple{Int,Int,Int} (month, day, hour)"

    time = DateTime(year, time...)

    !is_end ? (return time) : nothing;

    if n == 1
        time += Dates.Month(1)
    elseif n == 2
        time += Dates.Day(1)
    else
        error("Invalid time specification")
    end

    return time

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