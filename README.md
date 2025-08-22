[![Fimbul logo](https://github.com/sintefmath/Fimbul.jl/raw/main/docs/src/assets/logo_text_wide.png)](https://sintefmath.github.io/Fimbul.jl/dev/)

> [!TIP]
> Visit the docs at https://sintefmath.github.io/Fimbul.jl/dev/

# Geothermal simulation in Julia

Fimbul.jl is a [Julia](https://julialang.org/)-based toolbox for geothermal simulations based on [JutulDarcy.jl](https://github.com/sintefmath/JutulDarcy.jl), a fully differentiable, high-performance porous media simulator toolbox. Fimbul and JutulDarcy are developed by the [Applied Computational Science group](https://www.sintef.no/en/digital/departments-new/applied-mathematics/applied-computational-sciences/) at [SINTEF Digital](https://www.sintef.no/en/digital/).

## Getting started

To get started with Fimbul, install [Julia](https://julialang.org/) and follow these steps:
- Make a project folder in a suitable location and navigate to it
```bash
mkdir fimbul-testing
cd fimbul-testing/
```
- Start a Julia REPL in the project folder, activate an environment, and add Fimbul and JutulDarcy. We will also add GLMakie for plotting.
```julia
using Pkg; Pkg.activate(".");
Pkg.add("Fimbul");
Pkg.add("JutulDarcy");
Pkg.add("GLMakie");
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
The first time you run this code, Julia will compile the packages, which may take a few minutes. Subsequent runs will be much faster.

>[!NOTE]
>Interactive plotting requires `GLMakie`, which may not work if you are running Julia over SSH.

## License

Fimbul.jl is licensed under the MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2025 Øystein Klemetsdal, SINTEF Digital and Contributors
