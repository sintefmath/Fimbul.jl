using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie
using NetworkLayout, LayeredLayouts, GraphMakie
# using HYPRE

##
n = 100
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

(fig, ax) = plot_reservoir_state_phase_diagram(results[(:d, false)], tabs;
    pressure_limits = (p_min, p_max),
    enthalpy_limits = (h_min, h_max),
)
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