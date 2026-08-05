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
    export plot_phase_diagram_contours!, plot_phase_diagram_contours
    export plot_reservoir_state_ph!, plot_reservoir_state_ph
    export plot_reservoir_state_phase_diagram!, plot_reservoir_state_phase_diagram
    # Optimization
    export well_mismatch_thermal
    # Cases
    export analytical_1d
    export analytical_ates
    export benchmark_ht_1d
    export benchmark_ht_2d
    export geothermal_doublet
    export egs
    export ags
    export ates
    export btes
    export ftes
    export coaxial_bhe
    export egg_geothermal, egg_geothermal_doublet, egg_ates
    # Other utilities
    export thermal_radius_aquifer
    # Two-phase geothermal
    export H2OSystem
    export build_steam_tables_h2o
    export steam_tables_h2o

    # Load dependencies into namespace
    using Jutul, JutulDarcy
    using LinearAlgebra, Statistics
    using DelimitedFiles
    using Gmsh
    using Dates
    using Integrals

    # Meshing
    include("meshing/extruded.jl")
    include("meshing/fractured.jl")
    include("meshing/utils.jl")
    # H2O P-H formulation
    include("pvt/h2o_variables.jl")
    include("pvt/h2o_system.jl")
    include("pvt/h2o_setup.jl")
    # Wells
    include("wells/closed_loop.jl")
    include("wells/closed_loop_u1.jl")
    include("wells/closed_loop_coaxial.jl")
    include("wells/closed_loop_analytical.jl")
    # Cases
    include("cases/cases.jl")
    # Optimization
    include("optimization/objectives.jl")
    # Other utilities
    include("utils.jl")
    # Externals
    include("ext/ext.jl")

end
