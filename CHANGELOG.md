# Fimbul.jl Changelog

## [Unreleased]

### Breaking changes

- Closed-loop wells are now a single merged multisegment well instead of separate `<name>_supply` and `<name>_return` wells coupled by a cross term. The supply and return parts are merged with `JutulDarcy.WellMerging.merge_well_domains` and connected by a regular well segment. `setup_closed_loop_well` (and `setup_btes_well`, `setup_vertical_btes_well`) returns one well domain named `<name>`, and states and the model use this name. The two parts are controlled under names with a type-dependent suffix (see `closed_loop_well_names`): `_a`/`_b` for generic, `_supply`/`_return` for U1, and `_inner`/`_outer` for coaxial wells. For example, `coaxial_bhe` now uses `CoaxialWell_inner` and `CoaxialWell_outer`, and `ags` uses `AGS_a` and `AGS_b`. The return node comes after the supply nodes and has tag `:return`.
- `btes`: the sectors in `input_data[:sectors]` list the closed-loop wells instead of the supply and return names.

### New features

- New abstract type `AbstractClosedLoopWell` with the concrete types `GenericClosedLoopWell`, `U1ClosedLoopWell` and `CoaxialClosedLoopWell`. The type is stored as the `type` of the `MultiSegmentWell`, and `closed_loop_type` accepts either a type or `:generic` (alias `:simple`), `:u1` and `:coaxial`. New closed-loop types implement `setup_closed_loop_supply`.
- `ags` uses `GenericClosedLoopWell` with a general well topology.

### Compatibility

- Requires a JutulDarcy version where `setup_reservoir_model` adds the original names of merged wells to the facility.

## [v0.4.0]

### Breaking changes

- Fimbul now requires Jutul 0.5 and JutulDarcy 0.4, which are breaking releases. Cases return JutulDarcy objects, so code that uses them may need updates for the new JutulDarcy version.

### New features

**Magmatic intrusion case (`magmatic_intrusion`)**
- New high-enthalpy geothermal case with an ellipsoidal magmatic intrusion, based on the production and reinjection scenarios of Yapparova et al. (2023). Uses the pressure–enthalpy formulation (`H2OSystem`), so the reservoir may contain liquid, vapor, two-phase and supercritical water.
- The schedule has a natural-state period with no wells, followed by a production period with a producer at constant wellhead pressure and an optional cold-water injector. Wells are multisegment wells that are cased above the open interval. A natural-state result can be passed as `initial_state` to a production case.
- Tensor mesh with uniform cells across the intrusion and the plume region above it, and geometrically growing cells towards the boundaries.
- New example: `examples/production/magmatic_intrusion_demo.jl`.

**Pressure–enthalpy formulation**
- `set_reservoir_variable_defaults!` for `H2OSystem` models accepts `dh_max_abs` to limit the absolute enthalpy change per Newton iteration (not limited by default).
- `add_thermal_to_model!` for `H2OSystem` models now passes keyword arguments on to JutulDarcy. Previously, calls with keyword arguments silently bypassed the `H2OSystem` overload.

### Improvements

- `H2OSystem`: gravity uses a saturation-weighted face density, which varies continuously when a phase disappears in one of the cells. This avoids cycling Newton iterations close to the critical point.
- `H2OSystem`: adapted to the reworked compositional fluxes in JutulDarcy 0.4.
- Steam tables: supercritical states are labelled as liquid or vapor by the critical enthalpy, which removes a spurious two-phase band across the critical pressure. The `SteamTablesH2O` artifact has been updated accordingly.
- New script for building the steam tables artifact with CoolProp: `examples/pvt/build_steam_tables.jl`.
- Closed-loop wells (`ags`, `btes`, `coaxial_bhe`) can no longer be shut automatically, since shutting the supply or return gives a dead end.
- `coaxial_bhe`: uses the total energy formulation, and the grout cells have a small fluid volume to avoid a singular mass balance.
- `ags`: well sections are now connected also when the junction point lies on a cell face. With Jutul ≥ 0.5.3, the traced cells of the connected sections then do not share a cell, which made the case setup fail.

### Compatibility

- Requires Jutul ≥ 0.5.2 and JutulDarcy ≥ 0.4.0.

## [v0.3.6]

### New features

**Pressure–enthalpy formulation for high-temperature geothermal systems**
- New single-component, two-phase (liquid/vapor) H₂O system, `H2OSystem`, using pressure and specific enthalpy as primary variables. This makes it possible to simulate phase transitions such as boiling and condensation in high-enthalpy reservoirs. Temperature, phase enthalpies and internal energy are computed as secondary variables from tabulated steam properties.
- Steam tables are distributed as a lazy artifact and loaded with `steam_tables_h2o`. With the new `CoolProp` package extension, `build_steam_tables_h2o` can regenerate the tables, for example over a different pressure/enthalpy range.
- New governing-equations section in the documentation describing the formulation.

**High-temperature benchmarks**
- `benchmark_ht_1d`: 1D benchmark cases `:a`–`:e` from Weis et al. (2014), covering single-phase liquid, single-phase vapor and two-phase regimes, with optional vertical flow.
- `benchmark_ht_2d`: 2D upper-crust magmatic fluid-source plume benchmarks (`:single_phase_source`, `:two_phase_source`).
- HYDROTHERM reference solutions are distributed as a lazy artifact (`HYDROTHERMBenchmarks`).
- New validation examples: `examples/validation/validation_ht_1d.jl` and `examples/validation/validation_ht_2d.jl`. Analytical examples moved from `examples/analytical` to `examples/validation`.

**Plotting (GLMakie extension)**
- New phase-diagram plotting: `plot_phase_diagram_contours(!)`, `plot_reservoir_state_ph(!)` and `plot_reservoir_state_phase_diagram(!)`.

### Improvements

- `coaxial_bhe`: new `hz_min` keyword for vertical refinement near the well, and better default meshing. For vertical wells, perforation lengths are now taken from the cell `dz`, which avoids very short perforation segments in thin layers.
- `extruded_mesh`: now works with single-point cell constraints and no longer fails on degenerate refinement distances.
- Closed-loop BHE wells: avoid a singular grout system.
- `ftes`: added a full docstring.
- `ags`: fixed well coordinates.
- Updated plotting and resolution in the doublet, EGS, AGS, ATES, BTES, FTES, coaxial BHE and Egg HT-ATES examples.
- Documentation: new VitePress-based site with example tags, author badges and a version picker.

### Compatibility

- Requires Jutul ≥ 0.4.31 and JutulDarcy ≥ 0.3.12.
- New dependencies: `Artifacts`, `LazyArtifacts`, `DelimitedFiles`. New weak dependency: `CoolProp`.

## [v0.3.5]

### New features

**Fractured Thermal Energy Storage (`ftes`)**
- New `ftes` case for simulating Fractured Thermal Energy Storage systems. The setup places one central injector surrounded by multiple producer wells, with horizontal fracture planes connecting them to enable thermal transport through tight rock. Supports configurable fracture density, orientation, and radius.
- New example: `examples/storage/ftes_demo.jl`.

**EGS improvements**
- The `egs` case now uses cut_mesh with polygonal fracture cuts (`PolygonalSurface`), allowing more general fracture geometries.
- `examples/production/egs_demo.jl` now compares varying well spacing and effect of angled fractures

**Coaxial BHE**
- New deep `coaxial_bhe` case, with option to control whether fluid is injected into the inner pipe or the outer annulus.
- New example: `examples/production/coaxial_bhe_demo.jl`.

### New utilities

- `get_well_neighborship`: computes well cell connectivity for multi-branch wells. Used internally in `egs` and `ftes`.
- `scaled_rate`: computes a volumetric flow rate scaled to a given fraction of the pore volume in a region of interest over a time duration.
- Fracture utilities (`src/cases/utils/fractures.jl`): `add_fractures` and `strike_dip_to_normal`.
