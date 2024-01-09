/** @file gsStokesExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/


#include <gismo.h>

// #include <gsIncompressibleFlow/src/gsINSSolverSteady.h>
// #include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
// #include <gsIncompressibleFlow/src/gsINSUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsInfo << "Nothing happens here at the moment.n";

//     // ========================================= Settings ========================================= 

//     bool steady = true;
//     bool unsteady = true;

//     int deg = 1;
//     int numRefine = 3;
//     int maxIt = 10;
//     real_t viscosity = 0.1;
//     real_t timeStep = 0.1;
//     real_t tol = 1e-5;

//     bool plot = false;
//     int plotPts = 10000;
//     int numThreads = 1; 

//     //command line
//     gsCmdLine cmd("Solves the Stokes problem in a backward facing step (BFS) domain.");

//     cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
//     cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
//     cmd.addSwitch("plot", "Plot result in ParaView format", plot);

//     cmd.addInt("d", "deg", "B-spline degree for geometry representation", deg);
//     cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
//     cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
//     cmd.addInt("t", "nthreads", "Number of threads for parallel assembly", numThreads);
//     cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);

//     cmd.addReal("v", "visc", "Viscosity value", viscosity);
//     cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
//     cmd.addReal("", "tol", "Stopping tolerance", tol);

//     try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

//     gsInfo << "Solving Stokes problem in a backward facing step (BFS) domain.\n";
//     gsInfo << "viscosity = " << viscosity << "\n";

//     // ========================================= Define geometry ========================================= 
    
//     gsMultiPatch<> patches;

//     real_t a = 4;
//     real_t b = 2;
//     real_t a_in = 1;

//     patches = BSplineStep2D<real_t>(deg, a, b, a_in);

//     gsInfo << patches << "\n";


//     // ========================================= Define problem ========================================= 

//     gsBoundaryConditions<> bcInfo;
//     std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
//     gsFunctionExpr<> f("0", "0", 2); // external force

//     defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 2); // bcInfo, bndIn, bndOut, bndWall are defined here


//     // ========================================= Define basis ========================================= 

//     // Define discretization space by refining the basis of the geometry
//     gsMultiBasis<> basis(patches);
    
//     refineBasis_step(basis, numRefine, 0, 0, 0, 0, 2, a, b);
    
//     std::vector< gsMultiBasis<> >  discreteBases;
//     discreteBases.push_back(basis); // basis for velocity
//     discreteBases.push_back(basis); // basis for pressure
//     discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)


//     // ========================================= Solve ========================================= 

//     gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
//     gsINSSolverParams<real_t> params(NSpde, discreteBases);
//     params.options().setInt("numThreads",numThreads);

//     if (steady)
//     {
//         gsINSSolverSteady<real_t> NSsolver(params);

//         gsInfo << "\nSolving the steady Stokes problem with direct linear solver.\n";
//         gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

//         NSsolver.initialize();
//         NSsolver.solveStokes();

//         if (plot) 
//         {
//             gsField<> velocity = NSsolver.constructSolution(0);
//             gsField<> pressure = NSsolver.constructSolution(1);
    
//             gsInfo << "Plotting in Paraview...";
//             gsWriteParaview<>(velocity, "BFS_Stokes_st_velocity", plotPts, true);
//             gsWriteParaview<>(pressure, "BFS_Stokes_st_pressure", plotPts);
//             gsInfo << " done.\n";
//         }
//     }

//     if (unsteady)
//     {
//         params.options().setReal("timeStep", timeStep);

//         gsINSSolverUnsteady<real_t> NSsolver(params);

//         gsInfo << "\nSolving the unsteady Stokes problem with direct linear solver.\n";
//         gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

//         NSsolver.initialize();
//         NSsolver.solveGeneralizedStokes(maxIt, tol);

//         if (plot) 
//         {
//             gsField<> velocity = NSsolver.constructSolution(0);
//             gsField<> pressure = NSsolver.constructSolution(1);
    
//             gsInfo << "Plotting in Paraview...";
//             gsWriteParaview<>(velocity, "BFS_Stokes_unst_velocity", plotPts, true);
//             gsWriteParaview<>(pressure, "BFS_Stokes_unst_pressure", plotPts);
//             gsInfo << " done.\n";
//         }
//     }

    return 0; 
}
