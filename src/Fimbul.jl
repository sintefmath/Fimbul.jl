module Fimbul

    using Jutul, JutulDarcy
    using Dates
    include("cases/egg_geothermal.jl")

    export egg_geothermal, egg_geothermal_doublet, egg_ates

    using LBFGSB
    include("optimization/objectives.jl")
    include("optimization/utils.jl")

    export well_mismatch_thermal
    export calibrate_case

    using GLMakie
    include("visualization/well_data_plotting.jl")
    
    export plot_well_data

end
