using Fimbul
using Test

@testset "Fimbul.jl" begin
    include("closed_loop.jl")
    include("validation.jl")
    include("cases.jl")
    include("layered_domain.jl")
end
