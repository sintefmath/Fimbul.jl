[![Fimbul logo](https://github.com/sintefmath/Fimbul.jl/raw/main/docs/src/assets/logo_text_wide.png)](https://sintefmath.github.io/Fimbul.jl/dev/)

> [!TIP]
> Visit the docs at https://sintefmath.github.io/Fimbul.jl/dev/

# Geothermal simulation in Julia

Fimbul.jl is a [Julia](https://julialang.org/)-based toolbox for geothermal simulations based on [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl), a fully differentiable, high-performance porous media simulator toolbox. Fimbul and JutulDarcy are developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Getting started

To get started with Fimbul, install [Julia](https://julialang.org/) and follow these steps:
- Clone JutulDarcy and Fimbul
```bash
git clone https://github.com/sintefmath/JutulDarcy.jl.git
git clone https://github.com/sintefmath/Fimbul.jl.git
```
NOTE: Fimbul currently relies on the development version of JutulDarcy, and this repository therefore has to be cloned as well. This will likely change in a future release, so that the release version of JutulDarcy can be used instead.
- Make a project folder in a suitable location and navigate to it
```bash
mkdir fimbul-testing
cd fimbul-testing/
```
- Start a Julia REPL in the project folder, activate an environment, and add Fimbul and JutulDarcy. We will use the development versions that we just cloned.
```julia
using Pkg; Pkg.activate(".");
Pkg.develop(path="path/to/JutulDarcy/");
Pkg.develop(path="path/to/Fimbul/");
Pkg.instantiate()
```

You are now ready to run your first simulation! Fimbul comes with a number of example cases for geothermal energy applications. To check that everything works, you can run a small geothermal doublet case:
```julia
using Fimbul, JutulDarcy
using GLMakie
case = egg_geothermal_doublet()
result = simulate_reservoir(case)
plot_reservoir(case, result.states;
colormap = :seaborn_icefire_gradient, key = :Temperature)
```
Note that interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.
