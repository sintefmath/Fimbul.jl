# Fimbul -- geothermal energy simulation in Julia

Fimbul.jl is a Julia package for simulation of geothermal systems, written on
top of the porous media flow simulator
[JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl) developed by the
[Applied Computational Science
group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/)
at [SINTEF Digital](https://www.sintef.no/en/digital/).

To get started, clone this repository, start a julia REPL in the root folder of
the repository, and do
```julia
using Pkg; Pkg.instantiate()
```