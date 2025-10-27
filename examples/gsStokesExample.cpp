/** @file gsStokesExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/


#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    bool steady = false;
    bool unsteady = false;

    int deg = 1;
    int numRefine = 3;
    int maxIt = 10;
    real_t viscosity = 0.1;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;

    std::string outMode = "terminal"; // terminal/file/all/quiet
    bool plot = false;
    int plotPts = 10000;

    //command line
    gsCmdLine cmd("Solves the Stokes problem in a backward facing step (BFS) domain.");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    cmd.addInt("d", "deg", "B-spline degree for geometry representation", deg);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);

    cmd.addString("", "outMode", "Output mode (terminal/file/all/quiet)", outMode);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ( !(steady || unsteady) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --steady, --unsteady\n\n";

    gsFlowLogger logger(gsFlowLogger::parseOutputMode(outMode), "gsStokesExample.log");

    logger << "Solving Stokes problem in a backward facing step (BFS) domain.\n";
    logger << "viscosity = " << viscosity << "\n";

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 4;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(deg, a, b, a_in);

    logger << "domain: " << patches << "\n";


    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f("0", "0", 2); // external force

    defineBCs_step(bcInfo, 2);


    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    
    refineBasis_step(basis, numRefine, 0, 0, 0, 0, 2);
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)


    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases, &logger);

    if (steady)
    {
        gsINSSolverSteady<real_t> NSsolver(params);

        logger << "\nSolving the steady Stokes problem with direct linear solver.\n";
        logger << "numDofs: " << NSsolver.numDofs() << "\n";

        NSsolver.initialize();
        NSsolver.solveStokes();

        if (plot) 
        {
            gsField<> velocity = NSsolver.constructSolution(0);
            gsField<> pressure = NSsolver.constructSolution(1);
    
            logger << "Plotting in Paraview...";
            gsWriteParaview<>(velocity, "BFS_Stokes_st_velocity", plotPts, true);
            gsWriteParaview<>(pressure, "BFS_Stokes_st_pressure", plotPts);
            logger << " done.\n";
        }
    }

    if (unsteady)
    {
        logger << "\n";
        gsWarn << "Unsteady Stokes solver is currently not implemented.\n";

        // params.options().setReal("timeStep", timeStep);

        // gsINSSolverUnsteady<real_t> NSsolver(params);

        // logger << "\nSolving the unsteady Stokes problem with direct linear solver.\n";
        // logger << "numDofs: " << NSsolver.numDofs() << "\n";

        // NSsolver.initialize();
        // NSsolver.solveGeneralizedStokes(maxIt, tol);

        // if (plot) 
        // {
        //     gsField<> velocity = NSsolver.constructSolution(0);
        //     gsField<> pressure = NSsolver.constructSolution(1);
    
        //     logger << "Plotting in Paraview...";
        //     gsWriteParaview<>(velocity, "BFS_Stokes_unst_velocity", plotPts, true);
        //     gsWriteParaview<>(pressure, "BFS_Stokes_unst_pressure", plotPts);
        //     logger << " done.\n";
        // }
    }

    return 0; 
}
