# Prescribed-field drift of a charged 10 nm particle

Run `./Allrun`. Native `blockMesh` constructs a small fluid mesh. The gas is dry
N2 at 320 K/101325 Pa, using its own JANAF/Sutherland/perfectGas properties.
H2O is declared with zero mass fraction; species and energy equations are frozen.
The electric field is prescribed as (50000 0 0) V/m. This solver does not infer
electrode potentials or dielectric fields: use multiRegionChargedNanoParcelEFFoam
for that problem.

One initially stationary 10 nm particle has charge +1e and density 1000 kg/m3.
`nanoStokesDrag` uses mean free path 65 nm and `CoulombSep` uses the native E
field. Brownian motion and other forces are disabled for deterministic checking.
`check.py` computes the expected velocity and relaxation trajectory from the
actual case Sutherland coefficients, reads native Lagrangian output, and writes
a comparison to `postProcessing/drift/`. This is a short drift check; the particle
does not reach a wall during the 0.5 ms run.

Use `./Allclean` before a rerun. Native editable fields are in `0.orig/`.
