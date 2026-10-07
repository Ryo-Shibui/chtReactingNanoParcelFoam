#!/usr/bin/env python3
from pathlib import Path
import sys

root = Path(__file__).resolve().parent

def text(rel):
    return (root/rel).read_text()

checks = {
    "solveParcels control declared": "solveParcelsFluid" in text("fluid/createFluidFields.H"),
    "cloud construction conditional": "if (solveParcelsFluid[i])" in text("fluid/createFluidFields.H"),
    "SLGThermo construction conditional": text("fluid/createFluidFields.H").find("if (solveParcelsFluid[i])") < text("fluid/createFluidFields.H").find("new SLGThermo"),
    "parcel pointer nullable": "basicReactingTypeCloud* parcelsPtr = nullptr" in text("fluid/setRegionFluidFields.H"),
    "parcel evolve guarded": "if (solveParcelsFluid[i])" in text("chtReactingNanoParcelFoam.C") and "parcelsFluid[i].evolve()" in text("chtReactingNanoParcelFoam.C"),
    "steady mode detected": "steadyRun" in text("chtReactingNanoParcelFoam.C"),
    "steady p equation present": (root/"fluid/pEqnSteady.H").is_file(),
    "transient rho equation guarded": "mesh.transient()" in text("fluid/solveFluid.H") and "#include \"rhoEqn.H\"" in text("fluid/solveFluid.H"),
    "momentum parcel source guarded": "if (solveParcels)" in text("fluid/UEqn.H") and "parcelsPtr->SU(U)" in text("fluid/UEqn.H"),
    "species parcel source guarded": "parcelsPtr->SYi" in text("fluid/YEqn.H"),
    "energy parcel source guarded": "parcelsPtr->Sh" in text("fluid/EEqn.H"),
    "rho parcel source guarded": "parcelsPtr->Srho(rho)" in text("fluid/rhoEqn.H"),
    "pressure parcel source guarded": "parcelsPtr->Srho()" in text("fluid/pEqn.H"),
    "CHT matrix assembly retained": "fvMatrixAssemblyPtr->addFvMatrix(EEqn)" in text("fluid/EEqn.H"),
    "pure liquid Fuchs-Knudsen supported": "solution_.size() < 1 || solution_.size() > 2" in text("lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvapFuchsKnudsen/LiquidEvapFuchsKnudsen.C") and "solToSolMap_ < 0" in text("lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvapFuchsKnudsen/LiquidEvapFuchsKnudsen.C"),
    "runtime solve controls reread": "solveSpeciesFluid" in text("fluid/readFluidMultiRegionPIMPLEControls.H") and "requestedSolveParcels" in text("fluid/readFluidMultiRegionPIMPLEControls.H"),
}

for name, ok in checks.items():
    print(f"{'PASS' if ok else 'FAIL'}  {name}")

# Catch accidental legacy unguarded reference names in the Eulerian equation files.
legacy = {
    "fluid/UEqn.H": "parcels.SU",
    "fluid/YEqn.H": "parcels.SYi",
    "fluid/EEqn.H": "parcels.Sh",
    "fluid/rhoEqn.H": "parcels.Srho",
    "fluid/pEqn.H": "parcels.Srho",
}
for rel, token in legacy.items():
    ok = token not in text(rel)
    checks[f"no legacy unguarded {token} in {rel}"] = ok
    print(f"{'PASS' if ok else 'FAIL'}  no legacy unguarded {token} in {rel}")

sys.exit(0 if all(checks.values()) else 1)
