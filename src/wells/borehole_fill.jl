## Borehole fill materials for closed-loop wells

"""
Predefined borehole fill materials for closed-loop wells. Each entry holds the
thermal conductivity [W/(m K)], density [kg/m³] and heat capacity [J/(kg K)].
"""
const BOREHOLE_FILL_MATERIALS = Dict{Symbol, NamedTuple}(
    :grout => (thermal_conductivity = 2.3, density = 2000.0, heat_capacity = 1500.0),
)

"""
    borehole_fill_properties(fill = :grout;
        thermal_conductivity = missing, density = missing, heat_capacity = missing)

Thermal properties of the material filling the borehole around the pipes of a
closed-loop well, returned as a named tuple with fields `thermal_conductivity`
[W/(m K)], `density` [kg/m³] and `heat_capacity` [J/(kg K)].

`fill` is either the name of a predefined material (`:grout`) or a named tuple
with the three properties. Each property can be overridden by the
corresponding keyword argument.
"""
function borehole_fill_properties(fill::Symbol = :grout; kwarg...)
    haskey(BOREHOLE_FILL_MATERIALS, fill) || error("Unknown borehole fill \
        material :$fill. Available materials: \
        $(join(":" .* string.(collect(keys(BOREHOLE_FILL_MATERIALS))), ", ")).")
    return borehole_fill_properties(BOREHOLE_FILL_MATERIALS[fill]; kwarg...)
end

function borehole_fill_properties(fill::NamedTuple;
    thermal_conductivity = missing, density = missing, heat_capacity = missing)
    override(value, key) = ismissing(value) ? fill[key] : value
    for key in (:thermal_conductivity, :density, :heat_capacity)
        haskey(fill, key) || error("Borehole fill must define $key.")
    end
    return (
        thermal_conductivity = override(thermal_conductivity, :thermal_conductivity),
        density = override(density, :density),
        heat_capacity = override(heat_capacity, :heat_capacity),
    )
end

# Fill properties from the fill keyword arguments of the closed-loop setup
# functions, also accepting the deprecated grouting_* keyword arguments
function closed_loop_fill_properties(fill, fill_thermal_conductivity,
    fill_density, fill_heat_capacity;
    grouting_thermal_conductivity = missing,
    grouting_density = missing,
    grouting_heat_capacity = missing)

    fill_thermal_conductivity = deprecated_kwarg(fill_thermal_conductivity,
        grouting_thermal_conductivity,
        :grouting_thermal_conductivity, :fill_thermal_conductivity)
    fill_density = deprecated_kwarg(fill_density, grouting_density,
        :grouting_density, :fill_density)
    fill_heat_capacity = deprecated_kwarg(fill_heat_capacity,
        grouting_heat_capacity, :grouting_heat_capacity, :fill_heat_capacity)
    return borehole_fill_properties(fill;
        thermal_conductivity = fill_thermal_conductivity,
        density = fill_density,
        heat_capacity = fill_heat_capacity)

end

# Use the value of a deprecated keyword argument if it is given
function deprecated_kwarg(value, old_value, old_name, new_name)
    if !ismissing(old_value)
        Base.depwarn("Keyword argument `$old_name` is deprecated, use \
            `$new_name` instead.", :setup_closed_loop_well)
        ismissing(value) || error("Cannot give both `$old_name` and `$new_name`.")
        return old_value
    end
    return value
end
