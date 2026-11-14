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
