using Jutul, JutulDarcy, Fimbul, CoolProp, GLMakie
using NetworkLayout, LayeredLayouts, GraphMakie
# using HYPRE

##
tabs = Fimbul.build_steam_tables_2ph()

## Case a, no gravity (default)
case_a = benchmark_2ph_1d(benchmark_case = :a, nx = 100, enthalpy_tables = tabs, gravity = false)

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
);

result = simulate_reservoir(case_a, simulator=sim, config=cfg)

##
plot_reservoir(case_a.model, result.states)

##
using GLMakie
GLMakie.closeall()
timestep = length(result.states)÷2
time = 120si"year"
timestep = findfirst(t -> t >= time, result.time)
timestep = isnothing(timestep) ? length(result.states) : timestep

x = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[1,:]
z = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[3,:]

# x = z

fig = Figure()
for (k, prop) in enumerate([:Temperature, :Pressure, :Saturations, :TotalThermalEnergy, :PhaseMassDensities])
    ax = Axis(fig[1, k], title = string(prop))
    y = result.states[timestep][prop]
    if prop == :Saturations
        y = vec(y[1,:])
    elseif prop == :Temperature
        y = convert_from_si.(y, :Celsius)
    elseif prop == :Pressure
        y = convert_from_si.(y, "megapascal")
    end
    lines!(ax, x, y)
end
display(fig)

