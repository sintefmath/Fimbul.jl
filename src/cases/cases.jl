# Useful SI units
meter, litre, kilogram = si_units(:meter, :litre, :kilogram)
atm, = si_units(:atm)
second, hour, day, year = si_units(:second, :hour, :day, :year)
Kelvin, joule, watt, = si_units(:Kelvin, :joule, :watt)
darcy = si_unit(:darcy)

# Include utilities and setup functions
include("utils.jl")
include("analytical.jl")
include("doublet.jl")
include("egs.jl")
include("ags.jl")
include("ates.jl")
include("ftes.jl")
include("btes.jl")
include("egg_geothermal.jl")