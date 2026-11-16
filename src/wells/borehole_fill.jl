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

# Error for keyword arguments that were removed from the closed-loop setup
function check_removed_kwargs(; kwarg...)
    replacements = Dict(
        :radius_grout => "radius_borehole (1U) or grouting_thickness (coaxial)",
        :grouting_thermal_conductivity => "fill_thermal_conductivity",
        :grouting_density => "fill_density",
        :grouting_heat_capacity => "fill_heat_capacity",
    )
    for (k, v) in kwarg
        ismissing(v) && continue
        error("Keyword argument `$k` has been removed, use \
            $(replacements[k]) instead.")
    end
end
