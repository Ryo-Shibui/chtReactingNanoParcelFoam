# Physical tutorials

Start with `waterVapourDiffusion/Allrun` for a controlled, side-by-side H2O/N2
comparison of the original Sc=1 transport and species-dependent diffusion.
All gas properties are species-specific JANAF/Sutherland data; the human-facing
examples no longer use the manufactured constant-property gases from unit
verification. Those numerical verification cases remain separately in `tests/`.

`chargedParticleDrift/Allrun` is a supplementary physical N2 prescribed-field
particle example. `tutorials/Allrun` runs both tutorials.
