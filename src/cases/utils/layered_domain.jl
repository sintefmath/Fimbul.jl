function layered_reservoir_domain(constraints_2d, depths, layer_properties::NamedTuple=NamedTuple();
    mesh_args = NamedTuple(),
    reservoir_domain_args...
    )
    
    # Process layer properties
    num_layers = length(depths) - 1
    processed_layer_properties = process_properties(layer_properties, num_layers)
    # Create extruded mesh
    mesh, layers, metrics = extruded_mesh(constraints_2d, depths; mesh_args...)
    # Expand layer properties to cells
    layer_properties = Dict{Symbol, Any}()
    for (name, value) in processed_layer_properties
        layer_properties[name] = value[layers]
    end
    layer_properties = (; layer_properties...) # Convert to NamedTuple
    # Create reservoir domain
    domain = reservoir_domain(mesh; layer_properties..., reservoir_domain_args...)
    # Return domain, layers and metrics
    return domain, layers, metrics

end

function layered_reservoir_domain(mesh::JutulMesh, layers, layer_properties::NamedTuple=NamedTuple();
    reservoir_domain_args...)

    processed_layer_properties = process_properties(layer_properties, maximum(layers))
    # Expand layer properties to cells
    layer_properties = Dict{Symbol, Any}()
    for (name, value) in processed_layer_properties
        layer_properties[name] = value[layers]
    end
    layer_properties = (; layer_properties...) # Convert to NamedTuple

    # Create reservoir domain
    domain = reservoir_domain(mesh; layer_properties..., reservoir_domain_args...)

    # Return domain, layers and metrics
    return domain
end

function process_properties(layer_properties::NamedTuple, num_layers::Int)
    processed_layer_properties = Dict{Symbol, Any}()
    for (name, value) in pairs(layer_properties)
        if length(value) == 1
            processed_layer_properties[name] = fill(value[1], num_layers)
        elseif length(value) == num_layers
            processed_layer_properties[name] = value
        else
            error("Length of property $name ($(length(value))) \
            does not match number of layers $num_layers")
        end
    end
    return processed_layer_properties

end