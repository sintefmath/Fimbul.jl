using Jutul, JutulDarcy, Fimbul, SteamTables

tabs = build_steam_tables_2ph()

## Case a, no gravity (default)
case_a = benchmark_2ph_1d(benchmark_case = :a, enthalpy_tables = tabs, gravity = true)

# ## Case b, with gravity (flow in z direction)
# case_b = benchmark_2ph_1d(benchmark_case = :b, gravity = true, enthalpy_tables = tabs)

result = simulate_reservoir(case_a[1:48]; info_level = 1)

##
using GLMakie
GLMakie.closeall()
timestep = length(result.states)÷2
time = 120si"year"
timestep = findfirst(t -> t >= time, result.time)

x = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[1,:]
z = tpfv_geometry(physical_representation(case_a.model.models[:Reservoir].data_domain)).cell_centroids[3,:]

x = z

fig = Figure()
for (k, prop) in enumerate([:Temperature, :Pressure, :Saturations])
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