/** @file gsFlowExample.cpp
 
    In this example, the problem definition is read from a given xml file containing
    the computational domain geometry, boundary conditions and source function.
    The user can choose if steady and/or unsteady incompressible Navier-Stokes
    problem is computed. Several parameters for the computation can also be set (see help).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, gsFlowLogger& logger);

int main(int argc, char *argv[])
{

#if defined(gsPetsc_ENABLED) && defined(GISMO_WITH_MPI)

    // initialize MPI
    const gsMpi & mpi = gsMpi::init(argc, argv);
    gsMpiComm comm = mpi.worldComm();

    int rank  = comm.rank();

    PetscCall( PetscInitializeNoArguments() );

    // ========================================= Settings ========================================= 

    // solvers
    bool steady = false;
    bool unsteady = false;

    // file with problem definition
    std::string inputFile = "BFS2D_problem.xml";
    
    // discretization settings
    int numRefine = 3;
    int numElevate = 0; // number of degree elevations (before refinement)

    // problem parameters
    real_t viscosity = 0.1;
    
    // solver settings
    int maxIt = 10;
    int picardIt = 5;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;

    // linear solver settings
    std::string linSolver = "petsc"; // direct / iter / petsc
    int linIt = 50;
    real_t linTol = 1e-6;

    // output settings
    std::string outMode = "terminal"; // terminal/file/all/quiet
    bool plot = false;
    bool plotMesh = false;
    int plotPts = 10000;
    bool animation = false;
    int animStep = 5;

    // ---------------------------------------------------------------------------------

    //command line
    gsCmdLine cmd("Solves the incompressible Navier-Stokes problem in a given domain.");

    cmd.addSwitch("steady", "Solve steady INS problem", steady);
    cmd.addSwitch("unsteady", "Solve unsteady INS problem", unsteady);

    cmd.addString("", "input", "Path to the input xml file containing geometry, right-hand side functin and boundary conditions", inputFile);

    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);

    cmd.addInt("", "maxIt", "Max. number of Picard iterations (in steady case) or time steps (in unsteady case)", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    
    cmd.addString("", "linSolver", "Linear system solver (direct / iter / petsc)", linSolver);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);

    cmd.addString("o", "outMode", "Output mode (terminal/file/all/quiet)", outMode);
    cmd.addSwitch("plot", "Plot the final result in ParaView format", plot);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addSwitch("animation", "Plot animation of the unsteady problem", animation);
    cmd.addInt("", "animStep", "Number of iterations between screenshots for animation (used when animation = true)", animStep);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if ( (rank == 0) && !(steady || unsteady) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --steady, --unsteady\n\n";

    // ========================================= Define problem (geometry, BCs, rhs) ========================================= 
    
    gsMultiPatch<> patches;
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f; // external force
    std::string fn;

    gsFlowLogger logger(gsFlowLogger::parseOutputMode(outMode), "gsFlowExample.log", rank);

    logger << "Reading problem definition from file:\n" << inputFile << "\n\n";

    std::string path = gsFileManager::find(inputFile);
    if ( path.empty() )
    {
        if (rank == 0)
            gsWarn << "Input file not found, quitting.\n";

        return 1;
    }

    gsFileData<> fd(inputFile);
    fd.getId(0, patches);   // id=0: multipatch domain
    fd.getId(1, f);         // id=1: source function
    fd.getId(2, bcInfo);    // id=2: boundary conditions

    logger << "domain: " << patches << "\n";
    logger << "viscosity = " << viscosity << "\n";
    logger << "source function = " << f << "\n";

    // prepare ID string for the given problem
    size_t lastSlash = inputFile.find_last_of("/\\");
    size_t lastDot = inputFile.find_last_of('.');
    size_t start = (lastSlash == std::string::npos) ? 0 : lastSlash + 1;
    size_t count = (lastDot == std::string::npos || lastDot < start)
                     ? std::string::npos
                     : lastDot - start;
    std::string geoStr = inputFile.substr(start, count);

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    basis.degreeElevate(numElevate);

    for (int r = 0; r < numRefine; ++r)
        basis.uniformRefine();

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases, &logger, comm);
    params.options().setReal("timeStep", timeStep);
    params.options().setInt("nonlin.maxIt", picardIt);
    params.options().setReal("nonlin.tol", picardTol);
    params.options().setString("lin.solver", linSolver);
    params.options().setInt("lin.maxIt", linIt);
    params.options().setReal("lin.tol", linTol);

    gsOptionList solveOpt;
    solveOpt.addInt("maxIt", "", maxIt);
    solveOpt.addInt("plotPts", "", plotPts);
    solveOpt.addInt("animStep", "", animStep);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addSwitch("plot", "", plot);
    solveOpt.addSwitch("animation", "", animation);
    solveOpt.addSwitch("plotMesh", "", plotMesh);
    solveOpt.addString("id", "", "");
    
    if (steady)
    {
        solveOpt.setString("id", geoStr + "_steady");
        gsINSSolverSteady<real_t> NSsolver(params);

        logger << "\n----------\n";
        logger << "Solving the steady INS problem.\n";

        solveProblem(NSsolver, solveOpt, logger);

        if(linSolver != "direct")
            NSsolver.getLinSolver()->reportLinIterations();
    }

    if (unsteady)
    {
        solveOpt.setString("id", geoStr + "_unsteady");
        gsINSSolverUnsteady<real_t> NSsolver(params);

        logger << "\n----------\n";
        logger << "Solving the unsteady INS problem.\n";

        solveProblem(NSsolver, solveOpt, logger);

        if(linSolver != "direct")
            NSsolver.getLinSolver()->reportLinIterations();
    }

    PetscCall( PetscFinalize() );

#else
    gsWarn << "This version of INS solver requires MPI and PETSc, but some of them is not enabled!"
#endif

    return 0; 
}


template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, gsFlowLogger& logger)
{
    gsStopwatch clock;
    std::string id = opt.getString("id");

    int rank = NSsolver.getParams()->options().getInt("mpi.rank");

    // ------------------------------------
    // solve problem

    logger << "\ninitialization...\n";

    NSsolver.initialize();

    logger << "numDofs: " << NSsolver.numDofs() << "\n";

    // std::stringstream ss;
    // ss << "\nrank " << rank << "\nglobalStartEnd =\n" << NSsolver.getAssembler()->getGlobalStartEnd() << "\n" <<
    //     "\nblockUU:\n" << matStructureStr(NSsolver.getAssembler()->getBlockUU()) << "\n" <<
    //     "\nblockUP:\n" << matStructureStr(NSsolver.getAssembler()->getBlockUP()) << "\n" <<
    //     "\nblockPU:\n" << matStructureStr(NSsolver.getAssembler()->getBlockPU()) << "\n";
    // printOrderedOutput(ss.str(), NSsolver.getParams()->getMpiComm());

    gsINSSolverUnsteady<T, MatOrder>* pSolver = dynamic_cast<gsINSSolverUnsteady<T, MatOrder>* >(&NSsolver);

    if (pSolver) // start unsteady computation from Stokes solution
        pSolver->solveStokes();

    if (pSolver && opt.getSwitch("animation"))
        pSolver->solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), id, opt.getReal("tol"), opt.getInt("plotPts"));
    else
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0); // the last argument = min. number of iterations

    real_t totalT = clock.stop();

    logger << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    logger << "Solve time:" << NSsolver.getSolveTime() << "\n";
    logger << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    logger << "Total solveProblem time:" << totalT << "\n\n";

    // ------------------------------------
    // plot

    if ((rank == 0) && opt.getSwitch("plot")) 
    {
        logger << "Plotting in Paraview...\n";

        int plotPts = opt.getInt("plotPts");
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);
        gsWriteParaview<>(velocity, id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, id + "_pressure", plotPts);
        logger << "Done.\n";
    }
}