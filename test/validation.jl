using Jutul, JutulDarcy, Fimbul, Test, LinearAlgebra

@testset "Validation 1D" begin

    function compute_error(res, sol, x, t)
        dx = x[2] - x[1]
        dt = t[2] - t[1]
        ϵ = 0.0
        for (k, tk) in enumerate(t)
            ϵ += dt*norm(dx*(res.states[k][:Temperature] .- sol(x, tk))/L, 2)
        end
        ϵ /= sum(t)
        return ϵ
    end

    L = 100.0
    case, sol, x, t = analytical_1d(L = L, num_cells = 100, num_steps = 100)
    res = simulate_reservoir(case, info_level = -1)
    ϵ = compute_error(res, sol, x, t)
    @test ϵ < 1e-3

    T_b = convert_to_si(10.0, :Celsius)
    T_0 = x -> convert_to_si(100.0, :Celsius).*(40 <= x < 60) .+ T_b
    case, sol, x, t = analytical_1d(L = L, num_cells = 100, num_steps = 100,
        temperature_boundary = T_b, initial_condition = T_0)
    res = simulate_reservoir(case, info_level = -1)
    ϵ = compute_error(res, sol, x, t)
    @test ϵ < 1e-2

end