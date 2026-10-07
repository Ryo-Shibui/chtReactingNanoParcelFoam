# Full review notes — 2026-09-07

Reviewed against the uploaded full `chtReactingNanoParcelFoam` tree and the
OpenFOAM-v2406 `chtMultiRegionFoam` / `chtMultiRegionSimpleFoam` structure.

## Corrected issues

1. `parcelsFluid[i].evolve()` was unconditional in the main solver.
2. `SLGThermo` and the reacting cloud were always constructed, even when a
   steady carrier-only calculation was desired.
3. Parcel sources (`SU`, `SYi`, `Sh`, `Srho`) were unconditional in Eulerian
   equations.
4. The transient density equation `rhoEqn` was still solved with
   `default steadyState`; its `fvm::ddt(rho)` then has no transient diagonal.
5. The transient PIMPLE pressure equation was being reused blindly for a
   steady-state run. A separate steady pressure correction based on v2406
   `chtMultiRegionSimpleFoam/fluid/pEqn.H` was added.
6. `solveFlow`, `solveSpecies` and `solveEnergy` were only partly re-read from
   region `fvSolution`; the runtime controls are now consistent.
7. `solveParcels` is now read from each fluid region's `PIMPLE` dictionary and
   is intentionally startup-only.
8. Mixed default steady/transient ddt modes across CHT regions now stop with a
   clear error instead of producing an inconsistent multi-region solve.
9. Physical-time Courant/Diffusion timestep adjustment and Lagrangian evolve
   are skipped in steady carrier mode.
10. The uploaded archive had regressed to the original
    `LiquidEvapFuchsKnudsen` one-element `solution` bug. The pure-liquid fix was
    restored.

## Deliberate limitation

`solveParcels=true` with a `steadyState` carrier is rejected. OpenFOAM has a
steady cloud-tracking mode, but combining that with this custom reacting-nano
source coupling requires a separate, carefully defined steady Lagrangian
algorithm. The supported and recommended workflow is steady carrier/CHT first,
then a transient parcel restart.

## Files changed

- `chtReactingNanoParcelFoam.C`
- `fluid/createFluidFields.H`
- `fluid/setRegionFluidFields.H`
- `fluid/readFluidMultiRegionPIMPLEControls.H`
- `fluid/solveFluid.H`
- `fluid/UEqn.H`
- `fluid/YEqn.H`
- `fluid/EEqn.H`
- `fluid/rhoEqn.H`
- `fluid/pEqn.H`
- `fluid/pEqnSteady.H` (new)
- `lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvapFuchsKnudsen/LiquidEvapFuchsKnudsen.C`
- `verifyIntegration.py`
- `README.md`

## 2026-09-07 compile fix: surfaceFilm.Sh() type
OpenFOAM-v2406 `surfaceFilm.Sh()` returns `tmp<DimensionedField<scalar, volMesh>>`.
After parcel-source conditionalization, placing it directly after `rho*(U&g)` caused an invalid
`volScalarField + DimensionedField` expression. It is now removed from the inline RHS and added
as an RHS source after matrix construction using `EEqn -= surfaceFilm.Sh();`. The conditional
parcel heat source remains `EEqn -= parcelsPtr->Sh(he);`.
