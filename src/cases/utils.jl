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
    charge_months::Vector{String} = ["June", "July", "August", "September"],
    discharge_months::Vector{String} = ["December", "January", "February", "March"],
    start_month::Union{Missing, String} = missing,
    num_years = 5,
    years = missing,
    report_interval = 14day
    )

    # ## Process input
    # Validate months
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
        start_month = charge_months[1]
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