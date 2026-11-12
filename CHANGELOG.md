# Fimbul.jl Changelog

## [Unreleased]

### New features

**ATES in the SPE10 model (`ates_spe10`)**
- Aquifer thermal energy storage with a hot and a cold well in the Upper Ness or Tarbert formation of the SPE10 model 2, in 3D or in a vertical slice along the wells (`slice`). The SPE10 block, optionally a lateral subset and coarsened, is surrounded by cap and base rock and by padding cells with hydrostatic and geothermal boundary conditions. The injector reinjects what the other well produces.
- `Fimbul.spe10_block` gives the (coarsened) SPE10 permeability and porosity of a set of layers and columns.
- Reads the SPE10 data with `JutulDarcy.SPE10.load_data`.

## [v0.5.0]

### Breaking changes

- Closed-loop wells are now a single merged multisegment well instead of separate `<name>_supply` and `<name>_return` wells coupled by a cross term. The supply and return parts are merged with `JutulDarcy.WellMerging.merge_well_domains` and connected by a regular well segment. `setup_closed_loop_well` (and `setup_btes_well`, `setup_vertical_btes_well`) returns one well domain named `<name>`, and states and the model use this name. The two parts are controlled under names with a type-dependent suffix (see `closed_loop_well_names`): `_supply`/`_return` for generic and 1U wells, and `_inner`/`_outer` for coaxial wells. For example, `coaxial_bhe` now uses `CoaxialWell_inner` and `CoaxialWell_outer`. The return node comes after the supply nodes and has tag `:return`.
- `ClosedLoopGenericWell` (formerly `closed_loop_type = :simple`) without a given topology is now a single pipe from the first to the last reservoir cell instead of a U-shaped pipe. The `:simple` alias is removed.
- By default, the return part is perforated in the reservoir cell of the supply end node (if perforated).
- Consistent keyword arguments for the closed-loop setup: `wall_thickness` (generic) is renamed to `wall_thickness_pipe`, and `inner_pipe_thermal_conductivity`/`outer_pipe_thermal_conductivity` (coaxial) to `pipe_thermal_conductivity_inner`/`pipe_thermal_conductivity_outer`. Topology arguments (cell groups, `well_cell_centers`, `end_nodes`) without `neighborship` now give an error for all types instead of being ignored.
- `setup_btes_well` and `setup_vertical_btes_well` default to 1U wells (`closed_loop_type = :u1`).
- `fibonacci_pattern_2d`: point `n` is placed at radius `sqrt((n - 1/2)/(N - 1/2))` times the outer radius instead of `sqrt(n/N)`, which closes the gap at the center of the pattern, and `spacing` is the median distance between neighboring points. The sunflower pattern of `btes` therefore places the wells slightly differently.
- `btes`: `num_wells` and `well_spacing` can no longer both be given to the pattern form, since one follows from the other and the field size. Without a field size, the missing one defaults to 48 wells or 5 m.
- `btes`: with `topology = :sectors_parallel`, each sector is by default a single closed-loop well `S<n>` with its boreholes in series (`well_grouping = :sector`), instead of one well `B<n>` per borehole chained with `ReinjectionTarget` controls. Use `well_grouping = :borehole` for the previous behavior, e.g. to control individual boreholes. With `reversed_discharge = true`, a sector well is discharged by injecting into its return part and producing from its supply part.
- `btes`: the sectors in `input_data[:sectors]` list the closed-loop wells instead of the supply and return names.
- `btes`: the wells are perforated in layer 2 only by default (`well_layers = [2]`), so the boreholes start below the insulation layer. During rest periods, the supply parts of the closed-loop wells are put under `ClosedLoopRestControl` and the return parts are disabled, instead of disabling all wells. By default, discharging runs through each sector in the same direction as charging (`reversed_discharge = false`), whereas it previously ran from the outermost to the innermost well.
- Borehole geometry and fill of closed-loop wells: 1U wells take the borehole radius `radius_borehole` (formerly `radius_grout`), and coaxial wells take the fill thickness outside the outer pipe `grouting_thickness` (formerly `radius_grout`), as generic wells. The fill properties are set with `fill_material`, `fill_thermal_conductivity`, `fill_density` and `fill_heat_capacity` (formerly `grouting_thermal_conductivity`, `grouting_density` and `grouting_heat_capacity`). The old keyword arguments are removed and give an error.
- Coaxial wells: the default fill is `:grout` (density 2000 kg/m³ and heat capacity 1500 J/(kg K), as for 1U wells, instead of 1500 kg/m³ and 1460 J/(kg K)), and the default fill thickness is 30 mm (borehole radius ≈ 100.7 mm instead of 100 mm).

### New features

- New abstract type `AbstractClosedLoopWell` with the concrete types `ClosedLoopGenericWell`, `ClosedLoop1UWell` and `ClosedLoopCoaxialWell`. The type is stored as the `type` of the `MultiSegmentWell`, and `closed_loop_type` accepts either a type or `:generic`, `:u1` and `:coaxial`. New closed-loop types implement `setup_closed_loop_supply`.
- `ags` uses `ClosedLoopGenericWell` with a general well topology.
- `btes`: configurable well fields. The wells are placed with a built-in pattern (`btes(pattern; ...)` with `:sunflower`, `:rectangular` or `:circular`), or given explicitly as sectors of well trajectories (`btes(field; ...)`). The sectors are coupled in parallel with the boreholes of each sector in series (`topology = :sectors_parallel`), or in series with the boreholes of each sector in parallel between a shared inlet and outlet (`topology = :sectors_series`). New keyword arguments `reversed_discharge`, `well_grouping` and `permeability`. A sector modelled as a single well needs far fewer timesteps than chained wells, since the coupling between its boreholes is resolved within the well.
- `btes`: the patterns take a field size (`radius` for `:sunflower` and `:circular`, `extent = (Lx, Ly)` for `:rectangular`) and either `num_wells` or `well_spacing`. The other one follows from the field size, and the number of wells is always honored exactly. The wells are grouped into sectors with `sector_grouping`: bundles of adjacent spiral spines (`:spiral`, default for `:sunflower`), circle sectors (`:angular`, default for `:circular`), consecutive wells from the center outwards (`:sequential`, for `:sunflower` and `:circular`), or blocks of the grid (`:cartesian`, for `:rectangular`: blocks of `sector_dims` wells, by default as square as possible, tile the full rows from the bottom up, and the remaining wells at the top form the last sectors; without `extent`, the field is a whole number of blocks wide). Spiral sectors are rebalanced by moving outermost wells between neighboring sectors, so that their sizes differ by at most one.
- `plot_btes_field` and `plot_btes_field!` (GLMakie extension): plot the wells of a BTES field colored by sector, with the flow path through each sector and the structure of the pattern (sunflower spines, circular rings or rectangular grid). `btes` stores the well field, pattern and spacing in `input_data`.
- `ClosedLoopRestControl`: control for closed-loop wells at rest, which holds the wellhead pressure so that the fluid in the loop can expand and contract with temperature.
- Borehole fill materials for closed-loop wells (`fill_material = :grout` by default, or a named tuple of properties), see `Fimbul.borehole_fill_properties`.
- Coaxial wells without grout (`grouting_thickness = 0`): the outer pipe exchanges heat directly with the reservoir through its wall.
- All closed-loop types accept `dir` (also per reservoir cell), generic wells tag their nodes `:pipe`, and `setup_vertical_closed_loop_well` is exported. The setup for each type is documented in `setup_closed_loop_supply`.

### Improvements

- `H2OSystem`: adapted to the upwinding on potential differences in the compositional fluxes of JutulDarcy 0.4.1, and accepts (and ignores) the flash-specific variable defaults that JutulDarcy passes on.
- `well_mismatch_thermal` works with merged closed-loop wells: a part (e.g. `CoaxialWell_inner`) is looked up in the merged well, and the temperature is taken at the top node of that part. This is needed to calibrate against the production temperature of a closed-loop well.

### Documentation

- New API page for closed-loop wells (setup functions, closed-loop well types, controls and borehole fill).
- Reworked coaxial BHE demo (`examples/production/coaxial_bhe_demo.jl`): a layered reservoir, a study of the inner pipe conductivity, and calibration of a homogeneous digital twin against the production temperature and against the temperature along the well.
- New example: design optimization of a system of tilted coaxial BHEs (`examples/production/coaxial_bhe_design_demo.jl`), which minimizes the drilling cost under a constraint on the thermal power, with the drilling cost correlation of Lukawski et al. (2014).
- Reworked BTES demo (`examples/storage/btes_demo.jl`): compares the sunflower, circular and rectangular patterns, with plots of the well fields, the thermal plume over the charge cycles, and the charged and discharged energy and recovery factor per cycle.
- Example notebooks are published with each documentation version and preview, and example scripts link to the documented revision instead of `main`.
- Versioned documentation: released versions get their own documentation, selectable alongside `dev`.

### Development

- The code is formatted with [Runic](https://github.com/fredrikekre/Runic.jl). CI checks the formatting and suggests fixes on pull requests, and a pre-commit hook is provided. The formatting commit is listed in `.git-blame-ignore-revs`.

### Compatibility

- Requires Jutul ≥ 0.5.5, whose `find_enclosing_cells` finds all cells along a trajectory. `btes` relies on it to give every borehole the same cells.
- Requires JutulDarcy ≥ 0.4.1 (and < 0.5), where `setup_reservoir_model` adds the original names of merged wells to the facility. Fimbul 0.4.0 does not work with JutulDarcy 0.4.1 (`H2OSystem` models fail to set up and simulate).

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
