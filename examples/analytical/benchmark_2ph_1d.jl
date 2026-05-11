using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie
using NetworkLayout, LayeredLayouts, GraphMakie
# using HYPRE

##
tabs = Fimbul.build_steam_tables_2ph()

## Case a, no gravity (default)
case_a = benchmark_2ph_1d(benchmark_case = :b, nx = 100, enthalpy_tables = tabs, gravity = false)

##
st = deepcopy(case_a.state0[:Reservoir])
# st[:Pressure] = 1si"megapascal"
# st[:Enthalpy] .= range(100e3, 3000e3; length = length(st[:Enthalpy]))
st = Jutul.evaluate_all_secondary_variables(case_a.model.models[:Reservoir], st)

# ## Case b, with gravity (flow in z direction)
# case_b = benchmark_2ph_1d(benchmark_case = :b, gravity = true, enthalpy_tables = tabs)
sim, cfg = setup_reservoir_simulator(case_a;
    tol_mb = 1e-5,
    info_level = 2,
    max_timestep = maximum(case_a.dt)*0.1,
    max_nonlinear_iterations = 1000,
);

n_step = 48
result = simulate_reservoir(case_a[1:min(n_step, length(case_a.dt))], simulator=sim, config=cfg)

##
plot_reservoir(case_a.model, result.states)

##
using GLMakie
using GLMakie
GLMakie.closeall()
timestep = length(result.states)÷2
time = 120si"year"
timestep = findfirst(t -> t >= time, result.time)
timestep = isnothing(timestep) ? length(result.states) : timestep
timestep = length(result.states)

x = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[1,:]
z = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[3,:]
# x = z

fig = Figure()
for (k, prop) in enumerate([:Temperature, :Pressure, :Saturations, :PhaseMassDensities, :WaterPhase, :PhaseViscosities])
    ax = Axis(fig[1, k], title = string(prop))
    y = result.states[timestep][prop]
    if size(y, 1) == 2
        for ph in 1:2
            lines!(ax, x, vec(y[ph,:]), label = "Phase $ph")
        end
    else
        if prop == :Temperature
            y = convert_from_si.(y, :Celsius)
        elseif prop == :Pressure
            y = convert_from_si.(y, "megapascal")
        end
        lines!(ax, x, y)
    end
end
display(fig)
