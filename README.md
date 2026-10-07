# chtReactingNanoParcelFoam

`chtReactingNanoParcelFoam` is a custom multi-region OpenFOAM solver that combines conjugate heat transfer, compressible reacting fluid regions, solid heat conduction, surface-film coupling, and the custom reacting nano-parcel Lagrangian models from `reactingNanoParcelFoam`. It is useful for heated or cooled multi-region systems where reacting nano-droplets exchange mass, species, momentum, and enthalpy with the carrier gas while the fluid and solids are thermally coupled.

This solver has been checked only with OpenFOAM-v2406 through OpenFOAM-v2506.

## Features

- Extends the `chtMultiRegionFoam` style workflow with reacting fluid thermodynamics, solid regions, energy coupling, and optional reacting parcel clouds in each fluid region.
- Reuses the custom nano-parcel model stack from `reactingNanoParcelFoam`, including Fuchs-Knudsen phase change, electric-force support, droplet density/diameter update controls, and custom distribution/cloud-function models.
- Supports two workflows: steady carrier/species/CHT calculations without parcels, and transient carrier plus reacting nano parcels with two-way source coupling.

## Compilation

Load an OpenFOAM-v2406, v2412, or v2506 environment first:

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc
```

Build from this directory:

```bash
./Allwmake
```

The build script compiles the bundled custom Lagrangian libraries first, then:

```text
$FOAM_USER_APPBIN/chtReactingNanoParcelFoam
```

To clean generated build files:

```bash
./Allwclean
```

If your checkout loses executable bits, run:

```bash
chmod +x Allwmake Allwclean PREFLIGHT.sh
```

## Required Case Setup

Start from a valid multi-region CHT case. The case must include:

- `constant/regionProperties`
- fluid-region thermophysical and turbulence dictionaries
- solid-region thermophysical dictionaries
- region-specific `system/<region>/fvSchemes` and `system/<region>/fvSolution`
- region-specific initial fields for pressure, velocity, temperature, species, and solid temperature as required by the selected physics

For each fluid region, provide:

```text
constant/<fluidRegion>/surfaceFilmProperties
```

Use OpenFOAM's `none` surface-film model when no film is needed:

```text
surfaceFilmModel none;
```

For parcel calculations in a fluid region, provide:

```text
constant/<fluidRegion>/reactingCloud1Properties
```

Fields used by custom parcel forces, such as `E` for `CoulombForce`, must be placed in the corresponding fluid-region object registry/time directory.

## Usage

Run in serial:

```bash
chtReactingNanoParcelFoam
```

For parallel cases:

```bash
decomposePar -allRegions
mpirun -np <N> chtReactingNanoParcelFoam -parallel
reconstructPar -allRegions
```

Each fluid region reads split-solve controls under `PIMPLE` in that region's `fvSolution`:

```text
PIMPLE
{
    nOuterCorrectors          3;
    nCorrectors               1;
    nNonOrthogonalCorrectors  0;

    solveFlow       true;
    solveSpecies    true;
    solveEnergy     true;
    solveParcels    false;
}
```

`solveFlow false` freezes the carrier flow but still allows species and energy equations to advance when those switches are enabled. For energy-coupled CHT, keep `solveEnergy true`.

## Recommended Workflows

### Steady Carrier, Species, and CHT

Use `steadyState` as the default `ddtScheme` in every fluid and solid region and set:

```text
solveParcels false;
```

In this mode:

- parcel dictionaries and parcel clouds are not constructed;
- no physical-time parcel evolution is performed;
- the transient density equation is skipped;
- steady pressure correction is used for closed-volume CHT behavior.

This mode is recommended for obtaining a converged carrier and solid temperature field before injecting parcels.

### Transient Carrier with Reacting Nano Parcels

Use a transient `ddtScheme` in all regions and set:

```text
solveParcels true;
```

The solver advances reacting parcels once per physical time step before the PIMPLE/CHT outer iterations. Parcel and film source terms are then reused while the Eulerian flow, species, energy, and fluid-solid heat transfer converge at that time level.

### Two-Stage Procedure

1. Run the steady carrier/species/CHT case with `solveParcels false`.
2. Restart from the converged fields.
3. Change all regions to transient time schemes.
4. Set `solveParcels true` and run the coupled parcel calculation.

## Coupled Source Terms

When parcels are enabled, source terms are coupled into the carrier equations:

- continuity: parcel and film `Srho`
- pressure/continuity correction: parcel and film mass sources
- momentum: parcel `SU(U)`
- species: parcel `SYi` and film species sources
- energy: parcel `Sh(he)` and film heat sources

For closed-volume cases, the solver updates the target carrier mass using volume-integrated parcel and film mass transfer so that the pressure reference remains consistent with mass added to or removed from the gas.

## Notes

- All custom Lagrangian source directories needed by the nano-parcel model are bundled and built locally.
- Generated `lnInclude` directories and old object files are intentionally not stored in the repository; `wmake` recreates them.
- This solver uses a static multi-region mesh architecture. The dynamic-mesh operations from `reactingNanoParcelFoam` are not transplanted into the multi-region CHT workflow.
- The supplied `LiquidEvapFuchsKnudsen` model accepts either one volatile liquid entry, such as `solution (H2O);`, or two entries with the volatile liquid first and a non-volatile/solid component second.

