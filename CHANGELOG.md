# Fimbul.jl Changelog

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
