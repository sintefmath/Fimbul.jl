const DEFAULT_PHASE_PRESSURE_LIMITS = (1e5, 50e6)
const DEFAULT_PHASE_ENTHALPY_LIMITS = (500e3, 3500e3)

_plot_value(container, key::Symbol) = container[key]
_plot_value(container::NamedTuple, key::Symbol) = getproperty(container, key)

function plotting_timestep(case, res; time = nothing)
    plot_time = isnothing(time) ? get(case.input_data, :plot_time, last(res.time)) : time
    timestep = findfirst(t -> t >= plot_time, res.time)
    return isnothing(timestep) ? length(res.states) : timestep
end

phase_diagram_variable_key(variable::Symbol) = variable == :T ? :temperature : variable
phase_diagram_variable_key(variable) = variable

function phase_diagram_table(tables, variable)
    if variable isa Symbol
        key = phase_diagram_variable_key(variable)
        if haskey(tables, key)
            return tables[key]
        end
        error("Steam tables do not provide a lookup table for $(repr(variable)).")
    end
    return variable
end

default_phase_diagram_transform(variable::Symbol) = phase_diagram_variable_key(variable) == :temperature ? x -> x - 273.15 : identity
default_phase_diagram_transform(variable) = identity

function finite_extrema(values)
    vals = Float64[v for v in values if isfinite(v)]
    isempty(vals) && return nothing
    return extrema(vals)
end

function table_axes(table)
    if hasproperty(table, :X) && hasproperty(table, :Y)
        return getproperty(table, :X), getproperty(table, :Y)
    end
    return nothing
end

function default_phase_diagram_limits(table)
    axes = table_axes(table)
    if isnothing(axes)
        return DEFAULT_PHASE_PRESSURE_LIMITS, DEFAULT_PHASE_ENTHALPY_LIMITS
    end
    p_range, h_range = axes
    p_limits = finite_extrema(p_range)
    h_limits = finite_extrema(h_range)
    if isnothing(p_limits) || isnothing(h_limits)
        return DEFAULT_PHASE_PRESSURE_LIMITS, DEFAULT_PHASE_ENTHALPY_LIMITS
    end
    return p_limits, h_limits
end

function padded_limits(values; min_pad = 1e5)
    lo, hi = extrema(values)
    pad = max(0.05*(hi - lo), min_pad)
    return max(lo - pad, 0.0), hi + pad
end

function resolve_phase_diagram_limits(table; pressure_limits = nothing, enthalpy_limits = nothing)
    default_pressure_limits, default_enthalpy_limits = default_phase_diagram_limits(table)
    pressure_limits = isnothing(pressure_limits) ? default_pressure_limits : pressure_limits
    enthalpy_limits = isnothing(enthalpy_limits) ? default_enthalpy_limits : enthalpy_limits
    return pressure_limits, enthalpy_limits
end

function reservoir_state_ph(out; timestep = nothing, time = nothing)
    case = _plot_value(out, :case)
    res = _plot_value(out, :results)
    timestep = isnothing(timestep) ? plotting_timestep(case, res; time = time) : timestep
    state = res.states[timestep]
    return vec(state[:Pressure]), vec(state[:Enthalpy])
end

function evaluate_phase_diagram_grid(table, p_grid, h_grid, transform)
    values = Matrix{Float64}(undef, length(h_grid), length(p_grid))
    for (j, h) in enumerate(h_grid), (i, p) in enumerate(p_grid)
        value = table(p, h)
        value isa Number || error("Phase-diagram contours require a scalar table value, got $(typeof(value)).")
        values[j, i] = transform(value)
    end
    return values
end

function phase_diagram_envelope(tables, pressure_limits; n_pressure = 150)
    if !haskey(tables, :enthalpy_liquid_sat) || !haskey(tables, :enthalpy_vapor_sat)
        return nothing
    end
    p_min = max(first(pressure_limits), 0.0)
    p_max = min(last(pressure_limits), Fimbul.WATER_CRITICAL_PRESSURE)
    p_max < p_min && return nothing

    p_sat = collect(range(p_min, p_max; length = n_pressure))
    h_l_sat = tables[:enthalpy_liquid_sat].(p_sat)
    h_v_sat = tables[:enthalpy_vapor_sat].(p_sat)

    push!(p_sat, Fimbul.WATER_CRITICAL_PRESSURE)
    push!(h_l_sat, Fimbul.WATER_CRITICAL_ENTHALPY)
    push!(h_v_sat, Fimbul.WATER_CRITICAL_ENTHALPY)

    return h_l_sat, h_v_sat, p_sat
end

function phase_diagram_axis!(ax; axis_kwargs = (;))
    defaults = (
        xlabel = "Enthalpy [kJ/kg]",
        ylabel = "Pressure [MPa]",
    )
    for (key, value) in pairs((; defaults..., axis_kwargs...))
        setproperty!(ax, key, value)
    end
    return ax
end

function phase_diagram_figure(; figure_kwargs = (;), axis_kwargs = (;))
    fig = Figure(; figure_kwargs...)
    ax = Axis(fig[1, 1]; xlabel = "Enthalpy [kJ/kg]", ylabel = "Pressure [MPa]", axis_kwargs...)
    return fig, ax
end

function Fimbul.plot_phase_diagram_contours!(ax, tables;
    variable = :temperature,
    pressure_limits = nothing,
    enthalpy_limits = nothing,
    n_pressure = 150,
    n_enthalpy = 150,
    levels = 20,
    filled = true,
    lines = false,
    include_envelope = true,
    value_transform = nothing,
    contourf_kwargs = (;),
    contour_kwargs = (;),
    envelope_kwargs = (;))

    table = phase_diagram_table(tables, variable)
    pressure_limits, enthalpy_limits = resolve_phase_diagram_limits(table;
        pressure_limits = pressure_limits,
        enthalpy_limits = enthalpy_limits,
    )
    transform = isnothing(value_transform) ? default_phase_diagram_transform(variable) : value_transform

    p_grid = range(first(pressure_limits), last(pressure_limits); length = n_pressure)
    h_grid = range(first(enthalpy_limits), last(enthalpy_limits); length = n_enthalpy)
    values = evaluate_phase_diagram_grid(table, p_grid, h_grid, transform)

    if filled
        contourf_defaults = (; levels = levels, colormap = cgrad(:seaborn_icefire_gradient, levels, categorical = true, alpha = 1.0))
        hcf = contourf!(ax, h_grid ./ 1e3, p_grid ./ 1e6, values; contourf_defaults..., contourf_kwargs...)
    end

    if lines
        contour_defaults = (; levels = max(levels - 1, 1), color = (:white, 0.6), linewidth = 0.75)
        hcl = contour!(ax, h_grid ./ 1e3, p_grid ./ 1e6, values; contour_defaults..., contour_kwargs...)
    end

    if include_envelope
        envelope = phase_diagram_envelope(tables, pressure_limits; n_pressure = n_pressure)
        if !isnothing(envelope)
            h_l_sat, h_v_sat, p_sat = envelope
            keep = (p_sat .>= pressure_limits[1]) .& (p_sat .<= pressure_limits[2]) .& (h_l_sat .>= enthalpy_limits[1]) .& (h_v_sat .<= enthalpy_limits[2])
            h_l_sat = h_l_sat[keep]
            h_v_sat = h_v_sat[keep]
            p_sat = p_sat[keep]
            envelope_defaults = (; color = :white, linewidth = 1, linestyle = :solid)
            lines!(ax, h_l_sat ./ 1e3, p_sat ./ 1e6; envelope_defaults..., envelope_kwargs...)
            lines!(ax, h_v_sat ./ 1e3, p_sat ./ 1e6; envelope_defaults..., envelope_kwargs...)
        end
    end

    phase_diagram_axis!(ax)
    return ax, hcf, hcl
end

function Fimbul.plot_phase_diagram_contours(tables;
    figure_kwargs = (;),
    axis_kwargs = (;),
    kwargs...)
    fig, ax = phase_diagram_figure(; figure_kwargs = figure_kwargs, axis_kwargs = axis_kwargs)
    Fimbul.plot_phase_diagram_contours!(ax, tables; kwargs...)
    return fig, ax
end

function Fimbul.plot_reservoir_state_ph!(ax, pressure::AbstractVector, enthalpy::AbstractVector;
    line_kwargs = (;),
    kwargs...)
    defaults = (; color = :blue, linewidth = 2)
    lines!(ax, enthalpy ./ 1e3, pressure ./ 1e6; defaults..., line_kwargs..., kwargs...)
    phase_diagram_axis!(ax)
    return ax
end

function Fimbul.plot_reservoir_state_ph!(ax, out;
    timestep = nothing,
    time = nothing,
    kwargs...)
    pressure, enthalpy = reservoir_state_ph(out; timestep = timestep, time = time)
    return Fimbul.plot_reservoir_state_ph!(ax, pressure, enthalpy; kwargs...)
end

function Fimbul.plot_reservoir_state_ph(pressure::AbstractVector, enthalpy::AbstractVector;
    figure_kwargs = (;),
    axis_kwargs = (;),
    kwargs...)
    fig, ax = phase_diagram_figure(; figure_kwargs = figure_kwargs, axis_kwargs = axis_kwargs)
    Fimbul.plot_reservoir_state_ph!(ax, pressure, enthalpy; kwargs...)
    return fig, ax
end

function Fimbul.plot_reservoir_state_ph(out;
    figure_kwargs = (;),
    axis_kwargs = (;),
    kwargs...)
    fig, ax = phase_diagram_figure(; figure_kwargs = figure_kwargs, axis_kwargs = axis_kwargs)
    Fimbul.plot_reservoir_state_ph!(ax, out; kwargs...)
    return fig, ax
end

function Fimbul.plot_reservoir_state_phase_diagram!(ax, out, tables;
    timestep = nothing,
    time = nothing,
    pressure_limits = nothing,
    enthalpy_limits = nothing,
    contour_kwargs = (;),
    state_kwargs = (;),
    kwargs...)

    pressure, enthalpy = reservoir_state_ph(out; timestep = timestep, time = time)
    if isnothing(pressure_limits)
        pressure_limits = padded_limits(pressure)
    end
    if isnothing(enthalpy_limits)
        enthalpy_limits = padded_limits(enthalpy)
    end

    Fimbul.plot_phase_diagram_contours!(ax, tables;
        pressure_limits = pressure_limits,
        enthalpy_limits = enthalpy_limits,
        contour_kwargs...,
        kwargs...,
    )
    Fimbul.plot_reservoir_state_ph!(ax, pressure, enthalpy; state_kwargs...)
    return ax
end

function Fimbul.plot_reservoir_state_phase_diagram(out, tables;
    figure_kwargs = (;),
    axis_kwargs = (;),
    kwargs...)
    fig, ax = phase_diagram_figure(; figure_kwargs = figure_kwargs, axis_kwargs = axis_kwargs)
    Fimbul.plot_reservoir_state_phase_diagram!(ax, out, tables; kwargs...)
    return fig, ax
end