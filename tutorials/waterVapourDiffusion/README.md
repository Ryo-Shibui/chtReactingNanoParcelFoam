# H2O vapour diffusion in nitrogen: legacy versus species-dependent D

Run `./Allrun` after building the solver and sourcing OpenFOAM. Both cases are
generated from identical physical inputs: N2 background, dilute H2O mass fraction
`1e-4 + 2e-5*cos(pi*y/0.020)`, 320 K, 101325 Pa, 20 mm height, zero velocity,
laminar flow, no reactions, no parcels and frozen energy. The native `blockMesh`
and `setExprFields` utilities construct the mesh and initial fields. `0.orig/`
holds the editable input fields. Use `./Allclean` before repeating a run.

The only case-input difference is `speciesDiffusionProperties`:

- `legacy`: explicit `model legacy`, equivalent to molecular Sc=1 here.
- `speciesDependent`: `model dilute`, using the measured H2O/N2 correlation of
  Marrero and Mason (1972), Table 13/Figure 76: `D=1.87e-10*T^2.072*(101325/p)`
  in m2/s for 282-373 K, estimated uncertainty about 4%.

At this state the reference D is 2.9007637639e-5 m2/s. The nitrogen-background
Sc=1 limit is about 1.75920e-5 m2/s, so the reference diffusion is 1.6489 times
faster. A steady one-dimensional linear profile would hide this difference;
this transient cosine mode exposes its decay directly.

`compare.py` writes profiles, amplitude history, analytical comparisons and a
figure to `postProcessing/comparison/`. It checks mass-fraction closure, actual
diffusion rates and a visible difference between the two solutions. The
analytical curves use the dilute-background limit; small mixture-property
changes and finite spatial/time discretization are included in tolerances.

Gas data use `sutherland`, `janaf`, and `perfectGas`; N2, O2, H2O and DEG have
distinct property entries. The active carrier species here are only N2/H2O;
additional entries in foam.dat do not silently activate extra species. Read
`../../docs/transportProperties.md` for sources, fit ranges and fit errors.
This is a laminar demonstration; turbulent mixing can reduce the visible
difference and should not be used to disguise an incorrect molecular D.

The fluid is a named region. Its dictionaries are in `constant/fluid/`,
including `speciesDiffusionProperties`; global region selection is in
`constant/regionProperties`. This diffusion example intentionally has no solid.
