[![Fimbul logo](https://github.com/sintefmath/Fimbul.jl/raw/main/docs/src/assets/logo_wide.png)](https://sintefmath.github.io/JutulDarcy.jl/dev/)

> [!TIP]
> Visit the docs at https://sintefmath.github.io/Fimbul.jl/dev/

# Geothermal simulation in Julia

Fimbul.jl is a [Julia](https://julialang.org/)-based toolbox for geothermal simulations based on [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl), a fully differentiable, high-performance porous media simulator toolbox. Fimbul and JutulDarcy are developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Getting started

To get started with Fimbul, install [Julia](https://julialang.org/) and follow these steps:
- Clone Fimbul
```bash
git clone https://github.com/sintefmath/Fimbul.jl.git
```
- Make a project folder in a suitable location and navigate to it
```bash
mkdir fimbul-testing
cd fimbul-testing/
```
- Start a Julia REPL in the project folder, activate an environment, and add Fimbul. As this is still in a pre-release state, you have to dev it,
```julia
using Pkg; Pkg.activate(".");
Pkg.dev("path/to/Fimbul/"); Pkg.instantiate()
```
You are now ready to run your first simulation! Fimbul comes with a number of example cases for geothermal energy applications. To check that everything works, you can run a small geothermal doublet case:
```julia
using Fimbul
using GLMakie
case = egg_geothermal_doublet()
results = simulate_reservoir(case)
plot_reservoir(case, result.states;
colormap = :seaborn_icefire_gradient, key = :Temperature, step = length(case.dt))
```
Note that interactive plotting requires GLMakie, which may not work if you are running Julia over SSH.
