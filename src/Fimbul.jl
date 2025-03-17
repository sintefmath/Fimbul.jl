module Fimbul

    using Jutul, JutulDarcy
    using LinearAlgebra

    using Gmsh
    include("meshing/extruded.jl")
    include("meshing/utils.jl")

    export fibonacci_pattern_2d, extruded_mesh

    # Cases
    using Dates
    using Integrals
    include("cases/utils.jl")
    include("cases/egg_geothermal.jl")
    include("cases/analytical.jl")
    include("cases/btes.jl")

    export make_utes_schedule
    export egg_geothermal, egg_geothermal_doublet, egg_ates
    export analytical_1d
    export btes

    using LBFGSB
    include("optimization/objectives.jl")
    include("optimization/utils.jl")

    export well_mismatch_thermal
    export calibrate_case

    using GLMakie
    include("visualization/well_data_plotting.jl")
    
    export plot_well_data

end
