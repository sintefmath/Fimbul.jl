using Fimbul, JutulDarcy, Test

@testset "Egg" begin
    case = egg_geothermal()
    res = simulate_reservoir(case[1:1])
    @test true
    case = egg_geothermal_doublet()
    res = simulate_reservoir(case[1:1])
    @test true
    case = egg_ates()
    res = simulate_reservoir(case[1:1])
    @test true
end