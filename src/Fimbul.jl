module Fimbul

    # Meshing
    export fibonacci_pattern_2d
    export extruded_mesh, horizontal_fractured_mesh
    # Wells
    export setup_btes_well, setup_vertical_btes_well
    # Utils for setting up cases
    export make_schedule, make_utes_schedule
    export set_dirichlet_bcs
    # Visualization
    export plot_well_data!, plot_mswell_values!
    # Optimization
    export well_mismatch_thermal
    # Cases
    export analytical_1d
    export geothermal_doublet
    export egs
    export ags
    export ates
    export btes
    export egg_geothermal, egg_geothermal_doublet, egg_ates
    # Other utilities
    export thermal_radius_aquifer

    # Load dependencies into namespace
    using Jutul, JutulDarcy
    using LinearAlgebra, Statistics
    using Gmsh
    using Dates
    using Integrals

    # Meshing
    include("meshing/extruded.jl")
    include("meshing/fractured.jl")
    include("meshing/utils.jl")
    # Wells
    include("wells/btes.jl")
    # Cases
    include("cases/cases.jl")
    # Optimization
    include("optimization/objectives.jl")
    # Other utilities
    include("utils.jl")
    # Externals
    include("ext/ext.jl")

end
