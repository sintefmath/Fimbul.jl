function plot_well_data! end

function plot_mswell_values! end

"""
    plot_btes_field!(ax, case; kwargs...)
    plot_btes_field!(ax, field; pattern = missing, well_spacing = missing, kwargs...)

Plot the wells of a BTES field in `ax`, as seen from above, colored by sector.
The wells of each sector are connected in the order the fluid passes them
during charging, and the first well of each sector is marked with a larger
symbol. The field is either a case set up by [`btes`](@ref), or a well field
`[sector_1, sector_2, ...]` as accepted by `btes` (see also
`Fimbul.btes_well_field`).

With `structure = true` (default), the structure underlying the pattern is
drawn in `structure_color`: the spiral spines of a `:sunflower` pattern, the
rings of a `:circular` pattern, or the grid of a `:rectangular` pattern. For a
case set up with a pattern, the pattern is taken from the case.

# Keyword arguments
- `pattern = missing`: Pattern of the field (`:sunflower`, `:circular` or
  `:rectangular`), used to draw its structure.
- `well_spacing = missing`: Spacing between the rings of a `:circular` pattern.
  Defaults to the distance from the center to the innermost ring.
- `structure = true`: Draw the structure of the pattern.
- `structure_color = (:gray, 0.4)`: Color of the structure.
- `sector_colors = Makie.to_colormap(:tab10)`: Colors of the sectors.
- `linewidth = 2`, `markersize = 10`, `first_markersize = 20`: Line width of
  the flow paths, and marker sizes of the wells and of the first well of each
  sector.

Requires GLMakie to be loaded.
"""
function plot_btes_field! end

"""
    plot_btes_field(case_or_field; title = "", legend = true, kwargs...)

Plot the wells of a BTES field in a new figure, see [`plot_btes_field!`](@ref).
Returns the figure. Requires GLMakie to be loaded.
"""
function plot_btes_field end

function plot_phase_diagram_contours! end

function plot_phase_diagram_contours end

function plot_reservoir_state_ph! end

function plot_reservoir_state_ph end

function plot_reservoir_state_phase_diagram! end

function plot_reservoir_state_phase_diagram end

function build_steam_tables_h2o end
