# Closed-loop wells

Closed-loop wells (e.g. BTES wells and coaxial borehole heat exchangers) are
set up as a single merged multisegment well, where the supply and return parts
are connected by a regular well segment. The two parts are controlled as
separate wells, named by [`closed_loop_well_names`](@ref).

## Setup

```@docs
setup_closed_loop_well
setup_vertical_closed_loop_well
setup_btes_well
setup_vertical_btes_well
closed_loop_well_names
```

## Closed-loop well types

```@docs
AbstractClosedLoopWell
ClosedLoopGenericWell
ClosedLoop1UWell
ClosedLoopCoaxialWell
Fimbul.setup_closed_loop_supply
```

## Controls

```@docs
ClosedLoopRestControl
```

## Borehole fill

The material filling the borehole around the pipes (grout by default) is set
with the `fill_material` keyword argument of the closed-loop setup functions, either as
a predefined material or as a named tuple of properties. Each property can be
overridden by the `fill_thermal_conductivity`, `fill_density` and
`fill_heat_capacity` keyword arguments.

```@docs
Fimbul.borehole_fill_properties
```
