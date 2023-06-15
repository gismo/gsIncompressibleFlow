/** @file gsINSSolversExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/


#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolverSteady.h>
#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsINSUtils.h>

using namespace gismo;

void solveProblem(gsINSSolverBase<real_t>& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, int id);

int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    bool steady = true;
    bool steadyIt = false;
    bool unsteady = false;
    bool unsteadyIt = false;

    int deg = 1;
    int numRefine = 3;
    int maxIt = 10;
    int picardIt = 5;
    int linIt = 100;
    real_t viscosity = 0.1;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    real_t linTol = 1e-6;
    std::string precond = "PCDmod_FdiagEqual";

    bool plot = false;
    int plotPts = 10000;
    int numThreads = 1; 

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem in a backward facing step (BFS) domain.");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    cmd.addSwitch("steadyIt", "Solve steady problem with preconditioned GMRES as linear solver", steadyIt);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    cmd.addSwitch("unsteadyIt", "Solve unsteady problem with preconditioned GMRES as linear solver", unsteadyIt);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

    cmd.addInt("d", "deg", "B-spline degree for geometry representation", deg);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addInt("t", "nthreads", "Number of threads for parallel assembly", numThreads);
    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);

    cmd.addString("p", "precond",
    "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", 
    precond);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "Solving Navier-Stokes problem in a backward facing step (BFS) domain.\n";
    gsInfo << "viscosity = " << viscosity << "\n";

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 8;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(deg, a, b, a_in);

    gsInfo << patches << "\n";


    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f("0", "0", 2); // external force

    defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 2); // bcInfo, bndIn, bndOut, bndWall are defined here


    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    
    refineBasis_step(basis, numRefine, 0, 0, 0, 0, 2, a, b);
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)


    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsINSSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setInt("numThreads",numThreads);

    int id = 1;

    if (steady)
    {
        gsINSSolverSteady<real_t> NSsolver(params);

        gsInfo << "\nSolving the steady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solveProblem(NSsolver, maxIt, tol, plot, plotPts, id);

        id++;
    }

    if (unsteady)
    {
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("maxIt_picard", picardIt);
        params.options().setReal("tol_picard", picardTol);

        gsINSSolverUnsteady<real_t> NSsolver(params);

        gsInfo << "\nSolving the unsteady problem with direct linear solver.\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        solveProblem(NSsolver, maxIt, tol, plot, plotPts, id);
        
        id++;
    }

    if (steadyIt)
    {
        params.options().setInt("maxIt_lin", linIt);
        params.options().setReal("tol_lin", linTol);
        params.options().setString("precType", precond);
        // params.precOptions().setReal("gamma", 1); // parameter for AL preconditioner

        gsINSSolverSteadyIter<real_t, LinSolver > NSsolver(params);

        gsInfo << "\nSolving the steady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("precType") << "\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        if (params.options().getString("precType").substr(0, 3) == "PCD")
            NSsolver.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall, params.precOptions().getInt("pcd_bcType"));

        solveProblem(NSsolver, maxIt, tol, plot, plotPts, id);

        std::vector<index_t> itVector = NSsolver.getLinIterVector();

        gsInfo << "Iterations of linear solver in each Picard iteration:\n";
        for (size_t i = 0; i < itVector.size(); i++)
            gsInfo << itVector[i] << ", ";

        gsInfo << "\nAverage number of linear solver iterations per Picard iteration: " << NSsolver.getAvgLinIterations() << "\n";
        
        id++;
    }

    if (unsteadyIt)
    {
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("maxIt_picard", picardIt);
        params.options().setReal("tol_picard", picardTol);
        params.options().setInt("maxIt_lin", linIt);
        params.options().setReal("tol_lin", linTol);
        params.options().setString("precType", precond);
        // params.precOptions().setReal("gamma", 10); // parameter for AL preconditioner

        gsINSSolverUnsteadyIter<real_t, LinSolver > NSsolver(params);

        gsInfo << "\nSolving the unsteady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("precType") << "\n";
        gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

        if (params.options().getString("precType").substr(0, 3) == "PCD")
            NSsolver.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall, params.precOptions().getInt("pcd_bcType"));

        solveProblem(NSsolver, maxIt, tol, plot, plotPts, id);
        
        std::vector<index_t> itVector = NSsolver.getLinIterVector();

        gsInfo << "Iterations of linear solver in each Picard iteration:\n";
        for (size_t i = 0; i < itVector.size(); i++)
            gsInfo << itVector[i] << ", ";

        gsInfo << "\nAverage number of linear solver iterations per Picard iteration: " << NSsolver.getAvgLinIterations() << "\n";

        id++;
    }

    return 0; 
}


void solveProblem(gsINSSolverBase<real_t>& NSsolver, int maxIt, real_t tol, bool plot, int plotPts, int id)
{
    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    NSsolver.solve(maxIt, tol);

    //NSsolver.solveStokes(); // solve the Stokes problem
    //NSsolver.solveGeneralizedStokes(maxIt, tol); // solve the time-dependent Stokes problem (only in unsteady solvers)

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";

    if (plot) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);
 
        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, "BFS_solver" + util::to_string(id) + "_velocity", plotPts, true);
        gsWriteParaview<>(pressure, "BFS_solver" + util::to_string(id) + "_pressure", plotPts);
        gsInfo << " done.\n";
    }
}