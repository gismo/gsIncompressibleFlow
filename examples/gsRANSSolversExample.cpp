/** @file gsRANSSolversExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/


#include <gismo.h>

#include <gsIncompressibleFlow/src/gsRANSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsRANSSolverUnsteady<T, MatOrder>& NSsolver, gsOptionList opt, int geo);
template<class T, int MatOrder, class LinSolver> void reportLinIterations(gsFlowLinSystSolver_iter<T, MatOrder, LinSolver>* linSolverPtr);
template<class T, int MatOrder> void markElimDof(gsRANSSolverUnsteady<T, MatOrder>& NSsolver);

int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    // solvers
    bool steady = false;
    bool steadyIt = false;
    bool unsteady = true;
    bool unsteadyIt = false;

    // domain definition
    int geo = 1; // 1 - step, 2 - cavity, 3 - blade profile 2D
    int dim = 2;
    
    // discretization settings
    int deg = 1;
    int numRefine = 3;
    int wallRefine = 0;
    int leadRefine = 0; // for profile2D

    // problem parameters
    real_t viscosity = 0.01;
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

    cmd.addInt("g", "geo", "Computational domain (1 - step, 2 - cavity, 3 - profile (only 2D))", geo);
    cmd.addInt("d", "dim", "Space dimension", dim);

    cmd.addInt("", "deg", "B-spline degree for geometry representation", deg);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near step corner, cavity walls of blade profile", wallRefine);
    cmd.addInt("", "leadRefine", "Number of h-refinement steps near the beginning of the blade (for profile geometry)", leadRefine);

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


    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;
    real_t a, b, c;
    std::string geoStr;

    switch(geo)
    {
        case 1:
        {
            geoStr = util::to_string(dim) + "D backward-facing step";

            a = 8;
            b = 2;
            real_t a_in = 1;

            switch(dim)
            {
                case 2:
                default:
                    patches = BSplineStep2D<real_t>(deg, a, b, a_in);
                    break;

                case 3:
                    c = 2;
                    patches = BSplineStep3D<real_t>(deg, a, b, c, a_in);
                    break;
            }

            break;
        }
        case 2:
        {
            geoStr = util::to_string(dim) + "D lid-driven cavity";

            a = 1;
            b = 1; 

            switch(dim)
            {
                case 2:
                default:
                    patches = BSplineCavity2D<real_t>(deg, a, b);
                    break;

                case 3:
                    c = 1;
                    patches = BSplineCavity3D<real_t>(deg, a, b, c);
                    break;
            }

            break;
        }
        case 3:
        {
            geoStr = util::to_string(dim) + "D blade profile";

            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";

            //gsReadFile<>(FLOW_DATA_DIR "geo_profile2D.xml", patches);
            break;
        }
        default:
            GISMO_ERROR("Unknown domain.");
    }

    gsInfo << "Solving RANS problem with k-omega SST model in " << geoStr << " domain.\n";
    gsInfo << "viscosity = " << viscosity << "\n";
    gsInfo << patches;


    // ========================================= Define problem and basis ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f; // external force

    switch(dim)
    {
        case 2:
        default:
            f = gsFunctionExpr<>("0", "0", 2);
            break;

        case 3:
            f = gsFunctionExpr<>("0", "0", "0", 3);
            break;
    }

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);

    switch(geo)
    {
        case 1:
        {
            defineBCs_step(bcInfo, bndIn, bndOut, bndWall, dim); // bcInfo, bndIn, bndOut, bndWall are defined here
            refineBasis_step(basis, numRefine, 0, wallRefine, 0, 0, dim, a, b);
            break;
        }
        case 2:
        {
            defineBCs_cavity(bcInfo, bndWall, dim, 1); // bcInfo and bndWall are defined here, bndIn and bndOut remain empty
            refineBasis_cavity(basis, numRefine, wallRefine, dim);
            break;
        }
        case 3:
        {
            defineBCs_profile2D(bcInfo, bndIn, bndOut, bndWall, inVelX, inVelY);
            refineBasis_profile2D(basis, numRefine, wallRefine, leadRefine);
            break;
        }
        default:
            GISMO_ERROR("Unknown domain.");
    }
        
    std::vector< gsMultiBasis<> >  discreteBases, discreteBasesTM;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)

    discreteBasesTM.push_back(basis); // basis for k
    discreteBasesTM.push_back(basis); // basis for omega

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases, discreteBasesTM);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("assemb.loop", matFormation);
    params.setBndParts(bndIn, bndOut, bndWall);

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

    /*
    if (steady)
    {
        solveOpt.setString("id", "steady");
        params.options().setString("lin.solver", "direct");

        gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);

        // example of flow rate computation
        gsField<> velocity = NSsolver.constructSolution(0);
        gsFlowBndEvaluator_flowRate<real_t> flowRateEval(params, bndOut);
        flowRateEval.setVelocityField(velocity);
        flowRateEval.evaluate();
        gsInfo << "bndOut flow rate = " << flowRateEval.getValue() << "\n";
    }

    if (steadyIt)
    {
        solveOpt.setString("id", "steadyIt");
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);
        // params.precOptions().setReal("gamma", 1); // parameter for AL preconditioner

        gsINSSolverSteady<real_t, ColMajor > NSsolver(params);

        gsInfo << "\n----------\n";
        gsInfo << "Solving the steady problem with preconditioned GMRES as linear solver.\n";
        gsInfo << "Used preconditioner: " << params.options().getString("lin.precType") << "\n";

        solveProblem(NSsolver, solveOpt, geo);

        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver());
        reportLinIterations(linSolverPtr);
    }
    */

    if (unsteady)
    {
        solveOpt.setString("id", "unsteady");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "direct");

        gsRANSSolverUnsteady<real_t, RowMajor> NSsolver(params);

        gsInfo << "\n-----------------------------------------------------------\n";
        gsInfo << "Solving the unsteady RANS problem with direct linear solver\n";
        gsInfo << "-----------------------------------------------------------\n";

        solveProblem(NSsolver, solveOpt, geo);
    }

    /*
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
    */

    return 0; 
}


template<class T, int MatOrder>
void solveProblem(gsRANSSolverUnsteady<T, MatOrder>& NSsolver, gsOptionList opt, int geo)
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
            case 1:
            {
                geoStr = "BFS" + dimStr;
                break;
            }
            case 2:
            {
                geoStr = "LDC" + dimStr;
                break;
            }
            case 3:
            {
                geoStr = "profile" + dimStr;
                break;
            }
            default:
                GISMO_ERROR("Unknown domain.");
        }
    }

    // ------------------------------------
    // solve problem

    if(geo == 2)
        markElimDof(NSsolver);

    gsInfo << "\nInitialization...\n";
    NSsolver.initialize();

    gsInfo << "RANS numDofs: " << NSsolver.numDofs() << "\n";
    gsInfo << "TM numDofs: " << NSsolver.numDofsTM() << "\n\n";

    gsRANSSolverUnsteady<real_t>* pSolver = dynamic_cast<gsRANSSolverUnsteady<real_t>* >(&NSsolver);

    if (pSolver)
        if (opt.getSwitch("stokesInit"))
            pSolver->solveStokes();

    //if (pSolver && opt.getSwitch("animation"))
    //    pSolver->solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"), false, 1);
    //else
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
void markElimDof(gsRANSSolverUnsteady<T, MatOrder>& NSsolver)
{
    std::vector<gsMatrix<index_t> > elimDof(1);
    elimDof[0].setZero(1,1);

    NSsolver.markDofsAsEliminatedZeros(elimDof, 1);
}