# Governing equations in geothermal energy

Fimbul currently supports pure H<sub>2</sub>O as a single-component, two-phase
fluid system. Water may be present in a liquid phase ($L$) and/or a vapor phase
($V$).

## Conservation of mass
Mass conservation is written as

```math
    \partial_t\left([1-\phi]\rho_R + \phi[\rho_L s_L + \rho_V s_V]\right)
    + \nabla \cdot \left(\rho_L \vec{v}_L + \rho_V \vec{v}_V\right)
    = Q_M,
```

where $\phi$ is porosity, $\rho_\alpha$ is density for the rock ($\alpha = R$)
or fluid phase ($\alpha \in \{L, V\}$), and $s_\alpha$ is phase saturation
(the fraction of pore volume occupied by phase $\alpha$, with $s_L + s_V = 1$).
The term $Q_M$ denotes mass sources and sinks, including wells and boundary
fluxes. Phase velocities are given by Darcy's law,

```math
    \vec{v}_\alpha = - \lambda_\alpha \mathbf{K} \left(\nabla p - \rho_\alpha g \nabla z\right),
```

where $\lambda_\alpha = k_{r,\alpha}/\mu_\alpha$ is phase mobility (ratio of
relative permeability to viscosity), $\mathbf{K}$ is absolute permeability, $p$
is pressure, $g$ is gravitational acceleration, and $z$ is elevation.

> NOTE: In practice, the interface between liquid and vapor phases in an
> H<sub>2</sub>O system may exhibit non-negligible capillary pressure effects.
> This is currently not modeled in Fimbul; support is planned for a future
> release.

## Conservation of thermal energy
Thermal energy conservation is expressed as

```math
    \partial_t\left([1-\phi]\rho_R u_R + \phi[\rho_L s_L u_L + \rho_V s_V u_V]\right)
    + \nabla \cdot \left(\rho_L h_L \vec{v}_L + \rho_V h_V \vec{v}_V\right)
    - \nabla \cdot \left(\mathbf{\Lambda}\nabla T\right)
    = Q_H,
```

where $u_\alpha$ is specific internal energy, $h_\alpha = u_\alpha +
p/\rho_\alpha$ is specific enthalpy, $T$ is temperature, and $\mathbf{\Lambda}$
is the effective thermal conductivity tensor of the rock–fluid mixture. The
second term on the left represents advective heat transport by the flowing
phases; the third term is conductive heat transport via Fourier's law,

```math
    \vec{q}_c = - \mathbf{\Lambda}\nabla T.
```

The term $Q_H$ denotes thermal energy sources and sinks.

## Fluid properties
The fluid properties of H<sub>2</sub>O—such as density, viscosity, enthalpy,
and phase state—depend on both pressure and enthalpy and are provided through
precomputed lookup tables.

For the high-enthalpy (two-phase) regime, the steam tables are generated using
the [`CoolProp`](https://www.coolprop.org/) library [bell_pure_2014](@cite) via
the `FimbulCoolPropExt` extension. They are stored as `Artifacts` that are
loaded automatically when using the high-enthalpy formulation. The tables can
also be loaded directly using the `Artifacts` API. To generate your own tables,
see `build_steam_tables_h2o` in `FimbulCoolPropExt`.

For the low-enthalpy (single-phase liquid) regime, Fimbul uses PVT tables from
JutulDarcy generated using the NIST database.

The figure below shows the temperature, density, and vapor saturation of
H<sub>2</sub>O as functions of pressure and specific enthalpy, covering the
full range of both single-phase and two-phase states. The saturation line
separates the single-phase liquid region (left), two-phase region (center), and
single-phase vapor region (right). These phase diagram contours can be
visualized using `plot_phase_diagram_contours`.

## Single- and two-phase systems
Most Fimbul examples currently operate in the low-enthalpy regime, i.e.,
single-phase liquid conditions. The high-enthalpy examples use the two-phase
steam-table formulation described above.

## Multi-continuum and fractures
The formulation above assumes local thermal equilibrium between rock and fluid,
i.e., a single shared temperature field $T$. This is appropriate for many
porous-medium geothermal applications, but may not be adequate when flow is
dominated by fractures or other high-contrast heterogeneities.

To model fractured systems, Fimbul can be combined with the discrete fracture
modeling framework in JutulDarcy—see e.g., 
[ftes](@ref).
