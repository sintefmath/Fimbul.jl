module Fimbul

    using Jutul, JutulDarcy
    using LinearAlgebra

    using Gmsh
    include("meshing/extruded.jl")
    include("meshing/utils.jl")

    export fibonachi_pattern_2d, extruded_mesh

    using Dates
    include("cases/egg_geothermal.jl")

    using Integrals
    include("cases/analytical.jl")

    export egg_geothermal, egg_geothermal_doublet, egg_ates
    export analytical_1d

    using LBFGSB
    include("optimization/objectives.jl")
    include("optimization/utils.jl")

    export well_mismatch_thermal
    export calibrate_case

    using GLMakie
    include("visualization/well_data_plotting.jl")
    
    export plot_well_data

end
