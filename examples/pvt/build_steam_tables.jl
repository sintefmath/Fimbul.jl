# Build the pure-water steam tables used by `H2OSystem` and save them as a JLD2
# file in the format expected by the `SteamTablesH2O` artifact.
#
# The tables are generated with CoolProp through the `FimbulCoolPropExt`
# extension, using the default settings of `build_steam_tables_h2o`.
#
# Usage (from the Fimbul.jl root folder):
#
#     julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
#     julia --project=scripts scripts/build_steam_tables.jl [output_folder]
#
# The tables are saved as `steam_tables_h2o.jld2` in `output_folder` (default:
# `scripts/steam_tables_h2o`). To publish the tables as an artifact, compress
# the folder content and update `Artifacts.toml` accordingly. On macOS, set
# `COPYFILE_DISABLE=1` when running `tar` to avoid including `._*` metadata
# files in the archive.
using Fimbul, CoolProp
using Printf
const Jutul = Fimbul.Jutul

output_folder = length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "steam_tables_h2o")
mkpath(output_folder)
output_file = joinpath(output_folder, "steam_tables_h2o.jld2")

# ## Build tables
println("Building steam tables with CoolProp ...")
t = @elapsed tables = build_steam_tables_h2o()
T = tables[:temperature]
@printf("Built tables in %.1f s: %d pressure points (%.3g–%.3g Pa), %d enthalpy points (%.3g–%.3g J/kg)\n",
    t, length(T.X), extrema(T.X)..., length(T.Y), extrema(T.Y)...)

# ## Save tables
# The artifact loader expects a dictionary with the tables stored under "data".
Jutul.JLD2.save(output_file, "data", tables)
println("Saved tables to $output_file")

# ## Check saved tables
# Load the file back and check that the phase labelling is continuous across the
# critical pressure: vapor-like states (h > h_c) are labelled as vapor on both
# sides of p_c, and liquid-like states (h < h_c) as liquid.
loaded = Jutul.JLD2.load(output_file)["data"]
@assert keys(loaded) == keys(tables) "Saved tables do not match the generated tables"
Sv = loaded[:saturation_vapor_ph]
p_c, h_c = Fimbul.WATER_CRITICAL_PRESSURE, Fimbul.WATER_CRITICAL_ENTHALPY
for (h, Sv_expected) in ((h_c + 200e3, 1.0), (h_c - 200e3, 0.0))
    for p in (p_c - 50e3, p_c, p_c + 50e3, p_c + 1e6)
        Sv_ph = Sv(p, h)
        @assert isapprox(Sv_ph, Sv_expected; atol = 1e-8) "Unexpected vapor saturation $Sv_ph at p = $p Pa, h = $h J/kg"
    end
end
println("Checked phase labelling across the critical pressure.")
