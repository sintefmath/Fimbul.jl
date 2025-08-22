````@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: Fimbul
  text: Geothermal simulation in Julia
  image:
    src: logo_text_square.png
    alt: Fimbul.jl
  tagline: High-performance geothermal simulation toolbox based on automatic differentiation
  actions:
    - theme: brand
      text: Getting started
      link: /index
    - theme: alt
      text: View on Github
      link: https://github.com/sintefmath/Fimbul.jl
---
````

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

## Working with Fimbul
If you plan to use Fimbul extensively in your work, we strongly recommend that you read the documentation of [JutulDarcy](https://sintefmath.github.io/JutulDarcy.jl/dev/), in particular on [getting started](https://sintefmath.github.io/JutulDarcy.jl/dev/man/intro). This also covers the basics of installing Julia and creating a Julia environments, written for users who may not already be familiar with Julia package management.