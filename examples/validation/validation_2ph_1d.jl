using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie
using NetworkLayout, LayeredLayouts, GraphMakie
# using HYPRE

##
n = 500
tabs = Fimbul.build_steam_tables_2ph(n_pressure = n, n_enthalpy = n)

## Case a, no gravity (default)
n_step = 999
cases = [:d]
results = Dict{Tuple{Symbol, Bool}, Any}()
for c in cases
    for gravity in (false, true)
        if c == :e && gravity
            @info "Skipping case :e with gravity (not defined)"
            continue
        end
        @info "Setting up case $c with gravity = $gravity"
        case = benchmark_2ph_1d(
            benchmark_case = c,
            nx = 100,
            cell_size = 10.0,
            enthalpy_tables = tabs,
            gravity = gravity)
        sim, cfg = setup_reservoir_simulator(case;
            # tol_mb = 1e-4,
            tol_cnv=1e-3,
            tol_mb=1e-7,
            info_level = 0,
            max_timestep = maximum(case.dt),
            # max_nonlinear_iterations = 1000,
        );
        result = simulate_reservoir(case[1:min(n_step, length(case.dt))], simulator=sim, config=cfg)
        out = Dict(:case => case, :results => result)
        results[(c, gravity)] = out
    end
end

##

function plotting_timestep(case, res)
    time = get(case.input_data, :plot_time, last(res.time))
    timestep = findfirst(t -> t >= time, res.time)
    return isnothing(timestep) ? length(res.states) : timestep
end

function temperature_lookup_table(tabs)
    if haskey(tabs, :temperature)
        return tabs[:temperature]
    elseif haskey(tabs, :T)
        return tabs[:T]
    else
        error("Steam tables do not provide a temperature lookup table.")
    end
end

function plot_reservoir_state_phase_diagram!(ax, out, tabs;
    n_pressure = 150, n_enthalpy = 150, p_min = 1e5, p_max = 50e6, h_min = 500e3, h_max = 3500e3)
    case, res = out[:case], out[:results]
    timestep = plotting_timestep(case, res)
    state = res.states[timestep]

    pressure = vec(state[:Pressure])
    enthalpy = vec(state[:Enthalpy])
    temperature = temperature_lookup_table(tabs)

    p_min, p_max = extrema(pressure)
    h_min, h_max = extrema(enthalpy)
    # p_max = ifelse(ismissing(p_max), p_max_data, p_max)
    # h_max = ifelse(ismissing(h_max), h_max_data, h_max)
    p_pad = max(0.05 * (p_max - p_min), 1e5)
    h_pad = max(0.05 * (h_max - h_min), 1e5)

    p_grid = range(max(p_min - p_pad, 0.0), p_max + p_pad; length = n_pressure)
    h_grid = range(max(h_min - h_pad, 0.0), h_max + h_pad; length = n_enthalpy)
    T_grid = [temperature(p, h) - 273.15 for h in h_grid, p in p_grid]

    p_sat = collect(range(max(p_min - p_pad, 0.0), min(p_max + p_pad, 22.064e6); length = n_pressure))
    p_crit = 22.064e6
    h_crit = 2085e3
    h_l_sat = tabs[:enthalpy_liquid_sat].(p_sat)
    push!(h_l_sat, h_crit)
    h_v_sat = tabs[:enthalpy_vapor_sat].(p_sat)
    push!(h_v_sat, h_crit)
    push!(p_sat, p_crit)

    # ax = Axis(fig[pos],
    #     title = "Reservoir state in P-H space",
    #     xlabel = "Enthalpy [kJ/kg]",
    #     ylabel = "Pressure [MPa]")

    levels = 20
    cmap = cgrad(:seaborn_icefire_gradient, levels, categorical = true, alpha=1.0)
    # α = 0.8
    # cmap = [c*α + 1.0*(1-α) for c in cmap]
    # cmap = [c.*0.1 + ]
    contourf!(ax, h_grid ./ 1e3, p_grid ./ 1e6, T_grid; levels = levels, colormap = cmap)
    # contour!(ax, h_grid ./ 1e3, p_grid ./ 1e6, T_grid; levels = levels-1, color = (:white, 0.6), linewidth = 0.75)
    lines!(ax, h_l_sat ./ 1e3, p_sat ./ 1e6, color = :white, linewidth = 1, linestyle = :solid)
    lines!(ax, h_v_sat ./ 1e3, p_sat ./ 1e6, color = :white, linewidth = 1, linestyle = :solid)
    lines!(ax, enthalpy ./ 1e3, pressure ./ 1e6, color = :black, linewidth = 2)
    # scatter!(ax, enthalpy ./ 1e3, pressure ./ 1e6, color = :black)
    # Colorbar(fig[pos[1], pos[2] + 1], ax.plots[1], label = "Temperature [C]")
    return ax
end

timestep = plotting_timestep(results[(:d, false)][:case], results[(:d, false)][:results])
state = results[(:d, false)][:results].states[timestep]
pressure = vec(state[:Pressure])
enthalpy = vec(state[:Enthalpy])

GLMakie.closeall()
fig = Figure()

p_min, p_max = extrema(tabs[:temperature].X[2:end])
h_min, h_max = extrema(tabs[:temperature].Y)

p_min = 1e5
p_max = 50e6
h_min = 500e3
h_max = 3000e3

ax = Axis(fig[1, 1],
# limits = ((h_min, h_max)./1e3, (p_min, p_max)./1e6)
)
plot_reservoir_state_phase_diagram!(ax, results[(:d, false)], tabs)#; p_min=p_min, p_max=p_max, h_min=h_min, h_max=h_max)
# T_grid = [temperature(p, h) - 273.15 for h in h_grid, p in p_grid]
fig


##

function plot_case!(fig, row, c, out, gravity)

    (case, res) = out[:case], out[:results]

    timestep = plotting_timestep(case, res)
    # timestep = length(results.states)

    x = tpfv_geometry(physical_representation(case.model.models[:Reservoir].data_domain)).cell_centroids
    if gravity
        x = x[3,:]  # z-axis for vertical flow
        do_reverse = true
    else
        x = x[1,:]  # x-axis for horizontal flow
        do_reverse = false
    end
    
    props = [:Enthalpy, :Pressure, :Temperature, :Saturations, :PhaseMassDensities,
        :PhaseViscosities, :FluidEnthalpy, :LiquidMassFractions, :VaporMassFractions]

    props = [:Enthalpy, :Pressure, :Temperature, :Saturations, :PhaseMassDensities]
    for (k, prop) in enumerate(props)
        ax = Axis(fig[row, k], title = string(prop))
        y = res.states[timestep][prop]
        if size(y, 1) == 2
            for ph in 1:2
                if do_reverse
                    y[ph,:] = reverse(y[ph,:])
                end
                lines!(ax, x, vec(y[ph,:]), label = "Phase $ph")
            end
        else
            if prop == :Temperature
                y = convert_from_si.(y, :Celsius)
            elseif prop == :Pressure
                y = convert_from_si.(y, "megapascal")
            end
            if do_reverse
                y = reverse(y)
            end
            lines!(ax, x, vec(y))
        end
    end
    return fig

end

GLMakie.closeall()
for ((c, gravity), out) in results
    fig = Figure(size=(2400, 500))
    plot_case!(fig, 1, c, out, gravity)
    # plot_reservoir_state_phase_diagram!(fig, (1, 6), out, tabs)
    display(GLMakie.Screen(), fig)
end

##
case = results[(:c, false)][:case]
st = deepcopy(case.state0[:Reservoir])
# st[:Pressure] = 1si"megapascal"
# st[:Enthalpy] .= range(100e3, 3000e3; length = length(st[:Enthalpy]))
st = Jutul.evaluate_all_secondary_variables(case.model.models[:Reservoir], st)

##
Jutul.plot_variable_graph(case.model.models[:Reservoir])

##
state = results[(:c, false)][:results].result.states[end][:Reservoir]
h = state[:Enthalpy]
h_ph = state[:FluidEnthalpy]

lines(isapprox.(h .- h_ph[2,:], 0.0), label = "Phase 1")