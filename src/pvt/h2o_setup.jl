"""
    JutulDarcy.add_thermal_to_model!(model::GeothermalModel)

Activate the pressure-enthalpy thermal formulation for an `H2OSystem` model.

This replaces `Temperature` as a primary variable with `Enthalpy`, registers
temperature and phase enthalpy as secondary variables derived from the
pressure-enthalpy tables, and stores the `(P, T) -> H` lookup table in
`model.extra[:enthalpy]` for later use.
"""
function JutulDarcy.add_thermal_to_model!(model::GeothermalModel)

    invoke(JutulDarcy.add_thermal_to_model!, Tuple{SimulationModel}, model)
    Jutul.delete_variable!(model, :Temperature)

    pvt = model.system.pvt_tables

    set_secondary_variables!(model,
        Temperature = TemperatureFromEnthalpy(pvt[:temperature]),
        FluidEnthalpy = PhaseEnthalpyH2O(pvt[:enthalpy_liquid_ph], pvt[:enthalpy_vapor_ph]),
        FluidInternalEnergy = FluidInternalEnergyFromEnthalpy(),
    )
    set_primary_variables!(model, Enthalpy = Enthalpy())
    model.extra[:enthalpy] = pvt[:enthalpy]

    unique!(push!(model.output_variables, :Enthalpy, :Saturations))

    return model

end

JutulDarcy.system_uses_cnv_mb(system::H2OSystem) = true

using LazyArtifacts

"""
    steam_tables_h2o(; kwargs...) -> Dict{Symbol, Any}

Load the packaged H2O steam tables artifact. If the artifact cannot be loaded,
fall back to building the tables with CoolProp when that extension is available.
"""
function steam_tables_h2o(; kwargs...)
    try
        return steam_tables_h2o_from_artifact()
    catch
        if !isempty(methods(build_steam_tables_h2o))
            return build_steam_tables_h2o(; kwargs...)
        end
        throw(ErrorException(
            "Could not load the SteamTablesH2O artifact and CoolProp is not available to rebuild the steam tables."
        ))
    end
end

function steam_tables_h2o_from_artifact()
    tables_dir = artifact"SteamTablesH2O"
    table_files = filter(f -> endswith(lowercase(f), ".jld2"), readdir(tables_dir; join = true))
    isempty(table_files) && throw(ErrorException("SteamTablesH2O artifact does not contain a JLD2 table file."))
    tables = Jutul.JLD2.load(first(table_files))
    if tables isa Dict
        if haskey(tables, :data)
            return tables[:data]
        elseif haskey(tables, "data")
            return tables["data"]
        end
    end
    return tables
end
