module Fimbul

    # Meshing
    export fibonacci_pattern_2d
    export extruded_mesh, horizontal_fractured_mesh
    # Utils for setting up cases
    export make_utes_schedule
    export set_dirichlet_bcs
    # Visualization
    export plot_well_data, plot_mswell_values!
    # Optimization
    export well_mismatch_thermal
    export calibrate_case
    # Cases
    export analytical_1d
    export btes
    export egg_geothermal, egg_geothermal_doublet, egg_ates

    # Load dependencies into namespace
    using Jutul, JutulDarcy
    using LinearAlgebra
    using Gmsh
    using Meshes
    using Dates
    using Integrals
    using LBFGSB

    # Meshing
    include("meshing/extruded.jl")
    include("meshing/fractured.jl")
    include("meshing/utils.jl")
    # Cases
    include("cases/utils.jl")
    include("cases/analytical.jl")
    include("cases/doublet.jl")
    include("cases/btes.jl")
    include("cases/egg_geothermal.jl")
    # Optimization
    include("optimization/objectives.jl")
    include("optimization/utils.jl")
    # Externals
    include("ext/ext.jl")

end
