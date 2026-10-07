#!/bin/sh
set -eu
cd "${0%/*}" || exit 1

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    echo "ERROR: OpenFOAM environment is not loaded."
    exit 1
fi

if ! command -v wmake >/dev/null 2>&1; then
    echo "ERROR: wmake is not on PATH."
    exit 1
fi

case "${WM_PROJECT_VERSION:-unknown}" in
    *2406*) echo "OpenFOAM version: ${WM_PROJECT_VERSION}" ;;
    *) echo "WARNING: designed from OpenFOAM-v2406 sources; current version is ${WM_PROJECT_VERSION:-unknown}." ;;
esac

for f in \
    lagrangian/intermediate/submodels/Reacting/PhaseChangeModel/LiquidEvapFuchsKnudsen/LiquidEvapFuchsKnudsen.C \
    lagrangian/intermediate/submodels/Kinematic/ParticleForces/Coulomb/CoulombForce.C \
    lagrangian/turbulence/submodels/Thermodynamic/ParticleForces/BrownianMotion/BrownianMotionForce.C
 do
    test -s "$f" || { echo "ERROR: missing/non-readable $f"; exit 1; }
 done

echo "Source-tree preflight passed."
