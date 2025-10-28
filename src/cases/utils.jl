function make_schedule(forces, durations::Vector{Float64};
    num_reports = missing,
    report_interval = missing,
    num_cycles = 1
    )

    num_periods = length(durations)
    @assert num_periods > 0 "At one period is required"
    @assert length(forces) == num_periods "Number of forces must match number of periods"
    @assert all(durations .> 0.0) "All durations must be positive"

    has_num_reports = !ismissing(num_reports)
    has_report_interval = !ismissing(report_interval)
    if has_num_reports && has_report_interval
        error("Specify either report_interval or num_reports, not both")
    elseif !has_num_reports && !has_report_interval
        num_reports = fill(5, num_periods)
        report_interval = fill(missing, num_periods)
    elseif has_num_reports
        if num_reports isa Real
            num_reports = fill(Int(num_reports), num_periods)
        end
        report_interval = fill(missing, num_periods)
        @assert length(num_reports) == num_periods "Length of num_reports must match number of periods"
    elseif has_report_interval
        if report_interval isa Real
            report_interval = fill(report_interval, num_periods)
        end
        num_reports = fill(missing, num_periods)
        @assert length(report_interval) == num_periods "Length of report_interval must match number of periods"
    end

    forces0 = copy(forces)
    dt, forces = Float64[], Any[]
    for (f, d, n, r) in zip(forces0, durations, num_reports, report_interval)
        if !ismissing(n)
            # Number of reports specified
        elseif !ismissing(r)
            # Report interval specified
            n = max(Int(ceil(d/r)), 1)
        else
            error("Either num_reports or report_interval must be specified")
        end
        push!(dt, fill(d/n, n)...)
        push!(forces, fill(f, n)...)
    end

    dt = repeat(dt, num_cycles)
    forces = repeat(forces, num_cycles)
    timestamps = cumsum([0.0; dt])

    return dt, forces, timestamps

end

"""
    make_schedule(forces, periods; start_year=missing, num_years=1, report_interval=14day)

Create a simulation schedule with time steps and forces for multi-period
reservoir simulations using actual dates to define periods.

This function generates a detailed schedule based on user-defined periods and
associated forcing conditions. It automatically handles multi-year simulations
with periodic patterns and provides fleible time step control within each
period.

# Arguments
- `forces`: Vector of forcing conditions, one for each period. Each element
  should contain well controls, boundary conditions, or other simulation drivers
  for the corresponding period.
- `periods`: Vector of time period definitions. Length should be one more than
  `forces` to define period boundaries. Periods can be specified with month
  numbers, (month, day) tuples, (month, day, hour) tuples, or month name
  strings.

# Keyword Arguments
- `start_year::Union{Int,Missing}=missing`: Starting year for the simulation. If
  `missing`, uses the current year.
- `num_years::Int=1`: Number of years to simulate with the periodic schedule.
- `report_interval::Union{Real,Vector}=14day`: Time step size within each
  period. Can be a single value applied to all periods or a vector with specific
  intervals for each period.

# Returns
- `dt`: Vector of time step sizes in seconds for the entire simulation
- `force_vec`: Vector of forcing conditions corresponding to each time step
- `timestamps`: Vector of DateTime objects marking the time points

# Example
```julia
# Define seasonal operations: charge in summer, rest in winter
forces = [charge_controls, rest_controls]
periods = ["June", "September", "December"]  # June-Sep charge, Sep-Dec rest
dt, forces, times = make_schedule(forces, periods;
    num_years=3, report_interval=7si_unit(:day))
```
"""
function make_schedule(forces, periods;
    start_year = missing,
    num_years = 1,
    report_interval = missing,
    num_reports = missing
    )

    start_year = ismissing(start_year) ? Dates.year(now()) : start_year
    years = (0:num_years-1) .+ start_year

    num_periods = length(periods)-1
    @assert num_periods > 0 "At least one period is required"
    @assert length(forces) == num_periods "Number of forces must match number of periods"

    if num_reports isa Real || ismissing(num_reports)
        num_reports = fill(num_reports, num_periods)
    end
    @assert length(num_reports) == num_periods "Length of num_reports must match number of periods"
    if report_interval isa Real || ismissing(report_interval)
        report_interval = fill(report_interval, num_periods)
    end
    @assert length(report_interval) == num_periods "Length of report_interval must match number of periods"

    dt, force_vec, timestamps = Float64[], [], DateTime[]
    for year in years
        periods_year = Fimbul.process_periods(year, periods)
        if num_years > 1
            @assert periods_year[end] == periods_year[1] + Year(1) "The periods must make up a full year if num_years > 1"
        end
        for k in eachindex(periods_year)[1:end-1]
            start_time, end_time = periods_year[k], periods_year[k+1]
            Δt = Dates.value(end_time - start_time)*1e-3
            dt_k, forces_k, ts = make_schedule([forces[k]], [Δt]; 
                num_reports = num_reports[k], report_interval = report_interval[k], num_cycles = 1)
            ts = start_time .+ Millisecond.(round.(ts[1:end-1]*1e3))
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

"""
    make_utes_schedule(forces_charge, forces_discharge, forces_rest; kwargs...)

Create a specialized schedule for Underground Thermal Energy Storage (UTES)

This function generates a three-phase operational schedule typical for thermal energy storage:
1. **Charging phase**: Inject hot water to store thermal energy
2. **Rest phase**: No well activity 
3. **Discharging phase**: Extract stored thermal energy for heating applications

The schedule automatically handles seasonal operations with user-defined charge
and discharge periods, typically aligned with energy availability (summer
charging) and demand (winter discharging).

NOTE: The function does not validate that the provided forcing conditions are
appropriate for UTES operations. It is up to the user to ensure that the
`forces_charge`, `forces_discharge`, and `forces_rest` inputs contain suitable
well controls.

# Arguments
- `forces_charge`: Forcing conditions during charging phase (hot water injection)
- `forces_discharge`: Forcing conditions during discharging phase (thermal energy extraction)  
- `forces_rest`: Forcing conditions during rest periods (no well activity)

# Keyword Arguments
- `charge_period=["June", "September"]`: Start and end of charging phase (see
  `make_schedule` for valid formats)
- `discharge_period=["December", "March"]`: Start and end of discharging phase (see
  `make_schedule` for valid formats)
- `start_year::Union{Int,Missing}=missing`: Starting year for simulation.
  Defaults to current year.
- `num_years::Int=5`: Number of operational years to simulate.
- `kwargs...`: Additional arguments passed to `make_schedule()` (e.g.,
  `report_interval`).

# Returns
- `dt`: Vector of time step sizes in seconds
- `forces`: Vector of forcing conditions for each time step
- `timestamps`: Vector of DateTime objects for temporal tracking

# Example
```julia
# Standard ATES schedule: charge Jun-Sep, discharge Dec-Mar
dt, forces, times = make_utes_schedule(
    charge_forces, discharge_forces, rest_forces;
    charge_period = ["June", "September"],
    discharge_period = ["December", "March"], 
    num_years = 5,
    report_interval = 7si_unit(:day)
)
```

# Notes
- Rest periods are automatically inserted between charge and discharge phases
- The function handles year transitions and ensures chronological ordering
- Periods with zero duration are automatically filtered out
"""
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
    ch_end < ch_start ? ch_end += Dates.Year(1) : nothing
    dch_start = process_time(start_year, discharge_period[1])
    dch_end = process_time(start_year, discharge_period[2], true)
    dch_end < dch_start ? dch_end += Dates.Year(1) : nothing

    periods = [ch_start, ch_end, dch_start, dch_end, ch_start + Dates.Year(1)]
    forces = [forces_charge; forces_rest; forces_discharge; forces_rest]
    remove, p0 = [], nothing
    for (pno, p) in enumerate(periods)
        if 1 < pno && p - p0 <= Dates.Second(0)
            push!(remove, pno)
        end
        p0 = p
    end
    deleteat!(periods, remove)
    deleteat!(forces, remove.-1)
    periods = process_periods(start_year, periods)

    periods = [(Dates.month(p), Dates.day(p), Dates.hour(p)) for p in periods]

    dt, forces, timestamps = make_schedule(forces, periods;
        start_year = start_year, num_years = num_years, kwargs...)

    return dt, forces, timestamps

end

function process_periods(year, periods)

    periods_proc = Vector{DateTime}(undef, length(periods))
    for (pno, period) in enumerate(periods)
        time = process_time(year, period)
        if pno > 1
            time <= periods_proc[pno-1] ? time += Dates.Year(1) : nothing
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

function topo_sort_well(cells, msh, N = missing, z = missing)

    N = ismissing(N) ? get_neighborship(UnstructuredMesh(msh)) : N
    z = ismissing(z) ? tpfv_geometry(msh).cell_centroids[3, :] : z

    sorted_cells = Int[]
    top_ix = last(findmin(z[cells]))
    push!(sorted_cells, popat!(cells, top_ix))

    current_cell = sorted_cells[end]
    while !isempty(cells)
        for r = 1:2
            face_ix = N[r,:] .== current_cell
            neighbor_cells = N[mod(r,2)+1, face_ix]
            is_neighbor = [n ∈ cells for n in neighbor_cells]
            wn = neighbor_cells[is_neighbor]
            isempty(wn) ? continue : nothing
            @assert length(wn) <= 1 "Branching wells not supported"

            push!(sorted_cells, wn[1])
            popat!(cells, findfirst(isequal(wn[1]), cells))
            current_cell = wn[1]
        end
    end

    return sorted_cells

end