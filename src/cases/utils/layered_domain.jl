function layered_reservoir_domain(constraints_2d, depths;
    layer_properties = NamedTuple(),
    mesh_args = NamedTuple(),
    reservoir_domain_args...
    )
    
    # Process layer properties
    num_layers = length(depths) - 1
    function process_property(name, prop)
        if length(prop) == 1
            return fill(prop[1], num_layers)
        elseif length(prop) == num_layers
            return prop
        else
            error("Length of property $name ($(length(prop))) \
            does not match number of layers $num_layers")
        end
    end
    processed_layer_properties = Dict{Symbol, Any}()
    for (name, value) in pairs(layer_properties)
        value = process_property(name, value)
        processed_layer_properties[name] = value
    end

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