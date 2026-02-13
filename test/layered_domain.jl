using Test
using Fimbul
using Jutul

@testset "Layered Domain" begin

    # Define simple inputs
    depths = [0.0, 10.0, 20.0, 30.0] # 3 layers
    # Simple constraints: a few points defining a region
    constraints = [[(0.0, 0.0)], [(10.0, 10.0)]]
    
    # Layer properties
    # prop_scalar: same for all layers
    # prop_vec: specific for each layer
    layer_props = (
        porosity = 0.1,
        permeability = [1.0, 2.0, 3.0].*si"darcy"
    )
    
    # Create domain
    domain, layers, metrics = Fimbul.layered_reservoir_domain(constraints, depths;
        layer_properties = layer_props,
        mesh_args = (hxy_min = 5.0, hxy_max = 20.0, hz = 10.0),
        general_ad = false # pass generic arg to reservoir_domain
    )
    
    @testset "Return types" begin
        @test domain isa DataDomain
        @test layers isa AbstractVector{Int}
        @test metrics isa NamedTuple
    end
    
    @testset "Layer Mapping" begin
        # Check number of cells matches
        n_cells = number_of_cells(domain)
        @test length(layers) == n_cells
        
        # Check layer indices are valid
        @test all(1 .<= layers .<= 3)
        @test minimum(layers) == 1
        @test maximum(layers) == 3
    end
    
    @testset "Property values" begin
        # Check scalar property expansion
        # porosity should be 0.1 everywhere
        @test haskey(domain, :porosity)
        @test all(domain[:porosity] .≈ 0.1)
        
        # Check vector property mapping
        # permeability should match layer index
        @test haskey(domain, :permeability)
        perm = domain[:permeability]
        
        for cell in 1:number_of_cells(domain)
            layer = layers[cell]
            expected_perm = layer_props.permeability[layer]
            @test perm[cell] ≈ expected_perm
        end
    end
    
    @testset "Error handling" begin
        # Property with wrong length (2 values for 3 layers)
        bad_props = (
            porosity = [0.13, 0.169],
        )
        
        @test_throws ErrorException Fimbul.layered_reservoir_domain(constraints, depths;
            layer_properties = bad_props,
            mesh_args = (hxy_min = 5.0, hxy_max = 20.0, hz = 10.0)
        )
    end
end
