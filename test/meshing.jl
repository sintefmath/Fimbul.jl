using Jutul, JutulDarcy, Fimbul, Test, LinearAlgebra

@testset "Extruded mesh" begin

    # Test functionality for interpolation in z-direction
    depths = [0.0, 10.0, 100.0, 101.0]
    hz = [2.5, 10.0, 0.1]
    z, layers = Fimbul.interpolate_z(depths, hz)
    ok = true
    for l = 1:length(depths)-1
        layer_ix = findall(layers .== l)
        zl = z[layer_ix[1]:layer_ix[end]+1]
        @test zl[1] == depths[l]
        @test zl[end] == depths[l+1]
        dz = diff(zl)
        if 1 < l < length(depths)-1
            @test isapprox(maximum(dz), hz[l])
        else
            @test all(isapprox.(dz, hz[l]))
        end
    end

    # Check extruded mesh with cell constraints
    x = fibonacci_pattern_2d(10, spacing = 5.0)
    cc = map(x -> [x], x)
    depths = [0.0, 10.0, 50.0, 100.0]
    hz = [2.5, 10.0, 25.0]
    mesh, layers, metrics = extruded_mesh(cc, depths; hz = hz)
    @test unique(layers) == collect(1:length(depths)-1)
    for xi in x
        xi = [xi[1], xi[2]]
        d = map(xn -> norm(xn[1:2] .- xi,2), mesh.node_points)
        ix = isapprox.(d, 0.0)
        @test sum(ix) == length(layers)+1
    end

end