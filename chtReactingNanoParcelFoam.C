/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtReactingNanoParcelFoam

Group
    grpHeatTransferSolvers

Description
    Transient multi-region solver for compressible reacting fluid flow,
    solid heat conduction and conjugate heat transfer, augmented with the
    reactingNanoParcelFoam Lagrangian reacting-multiphase cloud models.

    Each fluid region owns an independent reacting cloud and optional surface
    film. Parcel-to-carrier mass, species, momentum and enthalpy source terms
    are coupled to the Eulerian equations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "surfaceFilmModel.H"
#include "SLGThermo.H"
#include "cloudMacros.H"

#ifndef CLOUD_BASE_TYPE
    #define CLOUD_BASE_TYPE ReactingMultiphase
    #define CLOUD_BASE_TYPE_NAME "reacting"
#endif

#include CLOUD_INCLUDE_FILE(CLOUD_BASE_TYPE)
#define basicReactingTypeCloud CLOUD_TYPE(CLOUD_BASE_TYPE)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Multi-region reacting CHT solver with optional reacting nano parcel"
        " clouds. Supports transient carrier/parcel operation and steady"
        " carrier/CHT operation with solveParcels=false."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"

    // All regions participating in one CHT solve must use the same default
    // time character (steadyState or transient).
    bool steadyRun = true;
    bool transientRun = true;

    forAll(fluidRegions, regionI)
    {
        steadyRun = steadyRun && fluidRegions[regionI].steady();
        transientRun = transientRun && fluidRegions[regionI].transient();
    }
    forAll(solidRegions, regionI)
    {
        steadyRun = steadyRun && solidRegions[regionI].steady();
        transientRun = transientRun && solidRegions[regionI].transient();
    }

    if (!steadyRun && !transientRun)
    {
        FatalErrorInFunction
            << "Mixed steady/transient default ddtSchemes across CHT regions "
            << "are not supported. Use steadyState for all regions or a "
            << "transient ddt scheme for all regions." << exit(FatalError);
    }

    if (transientRun)
    {
        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "surfaceFilmMultiRegionCourantNo.H"
        #include "setInitialMultiRegionDeltaT.H"
    }
    else
    {
        Info<< "Steady-state carrier/CHT mode detected (steadyState ddt)."
            << nl << "Courant/Diffusion based deltaT adjustment and physical-"
            << "time Lagrangian evolution are disabled." << nl << endl;
    }

    #include "createCoupledRegions.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        if (transientRun)
        {
            #include "compressibleMultiRegionCourantNo.H"
            #include "solidRegionDiffusionNo.H"
            #include "surfaceFilmMultiRegionCourantNo.H"
            #include "setMultiRegionDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Lagrangian and surface-film evolution is physical-time evolution.
        // It is therefore executed only in transient carrier mode.
        if (transientRun)
        {
            forAll(fluidRegions, i)
            {
                if (solveParcelsFluid[i])
                {
                    Info<< "\nEvolving reacting nano parcels in fluid region "
                        << fluidRegions[i].name() << endl;
                    parcelsFluid[i].evolve();
                }

                surfaceFilmFluid[i].evolve();

                scalar carrierMassRate =
                    fvc::domainIntegrate(surfaceFilmFluid[i].Srho()).value();

                if (solveParcelsFluid[i])
                {
                    carrierMassRate +=
                        fvc::domainIntegrate(parcelsFluid[i].Srho()).value();
                }

                // Closed-volume pressure reference must follow physical mass
                // transferred from parcels/film to the carrier.
                targetMassFluid[i] +=
                    runTime.deltaTValue()*carrierMassRate;
            }
        }

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = transientRun && (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }

            if (coupled)
            {
                Info<< "\nSolving energy coupled regions " << endl;
                fvMatrixAssemblyPtr->solve();
                #include "correctThermos.H"

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    if (!frozenFlow)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        // --- PISO loop
                        for (int corr=0; corr<nCorr; corr++)
                        {
                            if (mesh.steady())
                            {
                                #include "pEqnSteady.H"
                            }
                            else
                            {
                                #include "pEqn.H"
                            }
                        }
                        turbulence.correct();
                    }

                    if (transientRun)
                    {
                        rho = thermo.rho();
                    }
                    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                        << max(thermo.T()).value() << endl;
                }

                fvMatrixAssemblyPtr->clear();
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        fvMesh& mesh = fluidRegions[i];

                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        #include "readFluidMultiRegionPIMPLEControls.H"
                        #include "setRegionFluidFields.H"
                        // This loop is the stock CHT energy-only correction,
                        // not the user-requested frozen-flow scalar mode.
                        frozenFlow = true;
                        solveSpecies = false;
                        #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        fvMesh& mesh = solidRegions[i];

                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "setRegionSolidFields.H"
                        #include "solveSolid.H"
                    }

                    if (coupled)
                    {
                        Info<< "\nSolving energy coupled regions " << endl;
                        fvMatrixAssemblyPtr->solve();
                        #include "correctThermos.H"

                        forAll(fluidRegions, i)
                        {
                            #include "setRegionFluidFields.H"
                            rho = thermo.rho();
                        }

                        fvMatrixAssemblyPtr->clear();
                    }
                }
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
