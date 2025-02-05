/** @file gsINSSolversExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/


#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo);
template<class T, int MatOrder, class LinSolver> void reportLinIterations(gsFlowLinSystSolver_iter<T, MatOrder, LinSolver>* linSolverPtr);
template<class T, int MatOrder> void markElimDof(gsINSSolver<T, MatOrder>& NSsolver);

int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    // solvers
    bool steady = false;
    bool steadyIt = false;
    bool unsteady = false;
    bool unsteadyIt = false;

    // domain definition
    int geo = 1; // 0 - custom input file, 1 - step, 2 - cavity, 3 - blade profile 2D
    int dim = 2; // relevant for step and cavity
    std::string inputFile = "";
    
    // discretization settings
    int numRefine = 3;
    int wallRefine = 0; // relevant for step, cavity, profile2D
    int leadRefine = 0; // relevant for profile2D
    int numElevate = 0; // number of degree elevations (before refinement)

    // problem parameters
    real_t viscosity = 0.1;
    real_t inVelX = 1; // inlet x-velocity for profile2D
    real_t inVelY = 0; // inlet y-velocity for profile2D
    
    // solver settings
    int maxIt = 10;
    int picardIt = 5;
    int linIt = 50;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    real_t linTol = 1e-6;
    std::string matFormation = "EbE";
    std::string precond = "MSIMPLER_FdiagEqual";
    bool stokesInit = false; // start unsteady problem from Stokes solution

    // output settings
    bool quiet = false;
    bool plot = false;
    bool plotMesh = false;
    int plotPts = 10000;
    bool animation = false;
    int animStep = 5;

    // ---------------------------------------------------------------------------------

    //command line
    gsCmdLine cmd("Solves the Navier-Stokes problem in a given domain (step, cavity, blade profile).");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    cmd.addSwitch("steadyIt", "Solve steady problem with preconditioned GMRES as linear solver", steadyIt);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    cmd.addSwitch("unsteadyIt", "Solve unsteady problem with preconditioned GMRES as linear solver", unsteadyIt);

    cmd.addInt("g", "geo", "Computational domain (0 - custom file, 1 - step, 2 - cavity, 3 - profile (only 2D))", geo);
    cmd.addInt("d", "dim", "Space dimension", dim);
    cmd.addString("", "input", "Full path to the input xml file containing geometry, right-hand side functin and boundary conditions", inputFile);

    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near step corner, cavity walls of blade profile", wallRefine);
    cmd.addInt("", "leadRefine", "Number of h-refinement steps near the beginning of the blade (for profile geometry)", leadRefine);
    cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "inVelX", "x-coordinate of inflow velocity (for profile geometry)", inVelX);
    cmd.addReal("", "inVelY", "y-coordinate of inflow velocity (for profile geometry)", inVelY);

    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);
    cmd.addString("", "loop", "Matrix formation method (EbE = element by element, RbR = row by row)", matFormation);
    cmd.addString("p", "precond", "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", precond);
    cmd.addSwitch("stokesInit", "Set Stokes initial condition", stokesInit);

    cmd.addSwitch("quiet", "Supress (some) terminal output", quiet);
    cmd.addSwitch("plot", "Plot the final result in ParaView format", plot);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addSwitch("animation", "Plot animation of the unsteady problem", animation);
    cmd.addInt("", "animStep", "Number of iterations between screenshots for animation (used when animation = true)", animStep);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ========================================= Define problem (geometry, BCs, rhs) ========================================= 
    
    gsMultiPatch<> patches;
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f; // external force

    std::string fn, geoStr;

    switch(geo)
    {
        case 0:
            break; // inputFile is given from cmd
        
        default:
            gsWarn << "Unknown geometry ID, using backward-facing step.\n";
            geo = 1;

        case 1:
            geoStr = "BFS" + util::to_string(dim) + "D";
            fn = geoStr + "_problem.xml";
            inputFile = FLOW_DATA_DIR + fn;
            break;

        case 2:
            geoStr = "LDC" + util::to_string(dim) + "D";
            fn = geoStr + "_problem.xml";
            inputFile = FLOW_DATA_DIR + fn;
            break;

        case 3:
            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";

            geoStr = "profile2D";
            inputFile = FLOW_DATA_DIR + geoStr + "_problem.xml";
            break;
    }

    gsInfo << "Reading problem definition from file:\n" << inputFile << "\n\n";

    gsFileData<> fd(inputFile);
    fd.getId(0, patches);   // id=0: multipatch domain
    fd.getId(1, f);         // id=1: source function
    fd.getId(2, bcInfo);    // id=2: boundary conditions

    gsInfo << "Solving Navier-Stokes problem in " << geoStr << " domain.\n";
    gsInfo << patches;
    gsInfo << "viscosity = " << viscosity << "\n";
    gsInfo << "source function = " << f << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);
    basis.degreeElevate(numElevate);

    switch(geo)
    {
        case 0:
            for (int r = 0; r < numRefine; ++r)
                basis.uniformRefine();
            break;

        case 1:
        default:
            refineBasis_step(basis, numRefine, 0, wallRefine, 0, 0, dim, 8.0, 2.0, 2.0); // 8, 2, 2 are dimensions of the domain in the input xml file
            break;

        case 2:
            refineBasis_cavity(basis, numRefine, wallRefine, dim);
            break;

        case 3:
            refineBasis_profile2D(basis, numRefine, wallRefine, leadRefine);
            break;
    }    

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("assemb.loop", matFormation);

    gsOptionList solveOpt;
    solveOpt.addInt("geo", "", geo);
    solveOpt.addInt("maxIt", "", maxIt);
    solveOpt.addInt("plotPts", "", plotPts);
    solveOpt.addInt("animStep", "", animStep);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addSwitch("plot", "", plot);
    solveOpt.addSwitch("animation", "", animation);
    solveOpt.addSwitch("plotMesh", "", plotMesh);
    solveOpt.addSwitch("stokesInit", "", stokesInit);
    solveOpt.addString("id", "", "");

    if (steady)
    {
        solveOpt.setString("id", "steady");
        params.options().setString("lin.solver", "direct");

        gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);

        // example of flow rate computation
        if (geo == 1)
        {
            std::vector<std::pair<int, boxSide> > bndOut;
            bndOut.push_back(std::make_pair(0, boundary::east));
            bndOut.push_back(std::make_pair(1, boundary::east));
            gsField<> velocity = NSsolver.constructSolution(0);
            gsFlowBndEvaluator_flowRate<real_t> flowRateEval(params, bndOut);
            flowRateEval.setVelocityField(velocity);
            flowRateEval.evaluate();
            gsInfo << "bndOut flow rate = " << flowRateEval.getValue() << "\n";
        }
    }

    if (steadyIt)
    {
        solveOpt.setString("id", "steadyIt");
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);

        gsINSSolverSteady<real_t, ColMajor > NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";

        solveProblem(NSsolver, solveOpt, geo);

        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver());
        reportLinIterations(linSolverPtr);
    }

    if (unsteady)
    {
        solveOpt.setString("id", "unsteady");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "direct");

        gsINSSolverUnsteady<real_t, RowMajor> NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the unsteady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);
    }

    if (unsteadyIt)
    {
        solveOpt.setString("id", "unsteadyIt");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);
        // params.precOptions().setReal("gamma", 10); // parameter for AL preconditioner

        gsINSSolverUnsteady<real_t, ColMajor > NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the unsteady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";

        solveProblem(NSsolver, solveOpt, geo);
        
        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver());
        reportLinIterations(linSolverPtr);
    }

    return 0; 
}


template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo)
{
    gsStopwatch clock;

    // ------------------------------------
    // prepare strings for output filenames

    bool plot = opt.getSwitch("plot");
    std::string geoStr = "";
    std::string id = opt.getString("id");
    if (plot)
    {
        index_t dim = NSsolver.getParams()->getPde().domain().geoDim();
        std::string dimStr = util::to_string(dim) + "D";

        switch(opt.getInt("geo"))
        {
            case 0:
                geoStr = "customGeo";
            case 1:
            default:
                geoStr = "BFS" + dimStr;
                break;

            case 2:
                geoStr = "LDC" + dimStr;
                break;

            case 3:
                geoStr = "profile2D";
                break;
        }
    }

    // ------------------------------------
    // solve problem

    if(geo == 2)
        markElimDof(NSsolver);

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

    gsINSSolverUnsteady<real_t>* pSolver = dynamic_cast<gsINSSolverUnsteady<real_t>* >(&NSsolver);

    if (pSolver)
        if (opt.getSwitch("stokesInit"))
            pSolver->solveStokes();

    if (pSolver && opt.getSwitch("animation"))
        pSolver->solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"));
    else
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0);    

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time:" << totalT << "\n\n";

    // ------------------------------------
    // plot

    if (opt.getSwitch("plot")) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);

        int plotPts = opt.getInt("plotPts");
 
        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, geoStr + "_" + id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_" + id + "_pressure", plotPts);
        // plotQuantityFromSolution("divergence", velocity, geoStr + "_" + id + "_velocityDivergence", plotPts);
        gsInfo << " done.\n";
    }
}

template<class T, int MatOrder, class LinSolver>
void reportLinIterations(gsFlowLinSystSolver_iter<T, MatOrder, LinSolver>* linSolverPtr)
{
    std::vector<index_t> itVector = linSolverPtr->getLinIterVector();

    gsInfo << "Iterations of linear solver in each Picard iteration:\n";
    for (size_t i = 0; i < itVector.size(); i++)
        gsInfo << itVector[i] << ", ";

    gsInfo << "\nAverage number of linear solver iterations per Picard iteration: " << linSolverPtr->getAvgLinIterations() << "\n";
}

// mark one pressure DoF in the domain corner as fixed zero (for the lid driven cavity problem)
template<class T, int MatOrder>
void markElimDof(gsINSSolver<T, MatOrder>& NSsolver)
{
    std::vector<gsMatrix<index_t> > elimDof(1);
    elimDof[0].setZero(1,1);

    NSsolver.markDofsAsEliminatedZeros(elimDof, 1);
}