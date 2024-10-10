/** @file gsINSSolversExample.cpp
 
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/


#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo);
template<class T, int MatOrder> void markElimDof(gsINSSolver<T, MatOrder>& NSsolver);

int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    bool steady = false;
    //bool steadyIt = false;
    bool unsteady = false;
    //bool unsteadyIt = false;
    //bool stokesInit = false;

    int deg = 1;
    int geo = 1; // 1 - step, 2 - cavity
    int numRefine = 3;
    int wallRefine = 0;
    int maxIt = 10;
    int picardIt = 5;
    //int linIt = 100;
    real_t viscosity = 0.1;
    real_t timeStep = 0.1;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    //real_t linTol = 1e-6;
    //std::string precond = "PCDmod_FdiagEqual";
    std::string matFormation = "EbE";

    bool plot = false;
    bool plotMesh = false;
    int plotPts = 10000;
    bool quiet = false;

    // for profile2D
    real_t inVelX = 1;
    real_t inVelY = 0;
    int leadRefine = 0;

    //command line
    gsCmdLine cmd("Solves the Navier-Stokes problem in a 2D domain (step, cavity).");

    cmd.addSwitch("steady", "Solve steady problem with direct linear solver", steady);
    //cmd.addSwitch("steadyIt", "Solve steady problem with preconditioned GMRES as linear solver", steadyIt);
    cmd.addSwitch("unsteady", "Solve unsteady problem with direct linear solver", unsteady);
    //cmd.addSwitch("unsteadyIt", "Solve unsteady problem with preconditioned GMRES as linear solver", unsteadyIt);
    //cmd.addSwitch("stokesInit", "Set Stokes initial condition", stokesInit);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addSwitch("quiet", "Display output in terminal", quiet);

    cmd.addInt("d", "deg", "B-spline degree for geometry representation", deg);
    cmd.addInt("g", "geo", "Computational domain (1 - step, 2 - cavity, 3 - profile)", geo);
    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near step corner, cavity walls of blade profile", wallRefine);
    cmd.addInt("", "leadRefine", "Number of h-refinement steps near the beginning of the blade (for profile geometry)", leadRefine);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    //cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "inVelX", "x-coordinate of inflow velocity (for profile geometry)", inVelX);
    cmd.addReal("", "inVelY", "y-coordinate of inflow velocity (for profile geometry)", inVelY);
    //cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);

    cmd.addString("", "matForm", "Matrix formation method (EbE = element by element, RbR = row by row)", matFormation);
    //cmd.addString("p", "precond", "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", precond);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "Solving Navier-Stokes problem in domain " << geo << ".\n";
    gsInfo << "viscosity = " << viscosity << "\n";

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a, b;

    switch(geo)
    {
        case 1:
        {
            a = 8;
            b = 2;
            real_t a_in = 1;

            patches = BSplineStep2D<real_t>(deg, a, b, a_in);

            break;
        }
        case 2:
        {
            a = 1;
            b = 1; 

            patches = BSplineCavity2D<real_t>(deg, a, b);

            break;
        }
        case 3:
        {
            gsReadFile<>(FLOW_DATA_DIR "geo_profile2D.xml", patches);
            break;
        }
        default:
            GISMO_ERROR("Unknown domain.");
    }

    gsInfo << patches << "\n";


    // ========================================= Define problem and basis ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall; // containers of patch sides corresponding to inflow, outflow and wall boundaries
    gsFunctionExpr<> f("0", "0", 2); // external force

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> basis(patches);

    switch(geo)
    {
        case 1:
        {
            defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 2); // bcInfo, bndIn, bndOut, bndWall are defined here
            refineBasis_step(basis, numRefine, 0, wallRefine, 0, 0, 2, a, b);
            break;
        }
        case 2:
        {
            defineBCs_cavity2D(bcInfo, 1, bndWall); // bcInfo and bndWall are defined here, bndIn and bndOut remain empty
            refineBasis_cavity2D(basis, numRefine, wallRefine);
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
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(basis); // basis for velocity
    discreteBases.push_back(basis); // basis for pressure
    discreteBases[0].degreeElevate(1); // elevate the velocity space (Taylor-Hood element type)


    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("matFormation", matFormation);

    // bool stokesInit, bool plot
    gsOptionList solveOpt;
    solveOpt.addInt("geo", "", geo);
    solveOpt.addInt("id", "", 0);
    solveOpt.addInt("maxIt", "", maxIt);
    solveOpt.addInt("plotPts", "", plotPts);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addSwitch("plot", "", plot);
    solveOpt.addSwitch("plotMesh", "", plotMesh);
    //solveOpt.addSwitch("stokesInit", "", stokesInit);

    int id = 1;

    if (steady)
    {
        solveOpt.setInt("id", id);

        gsINSSolverSteady<real_t, ColMajor> NSsolver(params);

        gsInfo << "\nSolving the steady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);

        // example of flow rate computation
        gsField<> velocity = NSsolver.constructSolution(0);
        gsFlowBndEvaluator_flowRate<real_t> flowRateEval(params, bndOut);
        flowRateEval.setVelocityField(velocity);
        flowRateEval.evaluate();
        gsInfo << "bndOut flow rate = " << flowRateEval.getValue() << "\n\n";

        id++;
    }

    if (unsteady)
    {
        solveOpt.setInt("id", id);

        params.options().setReal("timeStep", timeStep);
        params.options().setInt("maxIt_picard", picardIt);
        params.options().setReal("tol_picard", picardTol);

        gsINSSolverUnsteady<real_t, RowMajor> NSsolver(params);

        gsInfo << "\nSolving the unsteady problem with direct linear solver.\n";

        solveProblem(NSsolver, solveOpt, geo);
        
        id++;
    }

    // if (steadyIt)
    // {
    //     solveOpt.setInt("id", id);
    //     params.options().setInt("maxIt_lin", linIt);
    //     params.options().setReal("tol_lin", linTol);
    //     params.options().setString("precType", precond);
    //     // params.precOptions().setReal("gamma", 1); // parameter for AL preconditioner

    //     gsINSSolverSteadyIter<real_t, LinSolver > NSsolver(params);

    //     gsInfo << "\nSolving the steady problem with preconditioned GMRES as linear solver.\n";
    //     gsInfo << "Used preconditioner: " << params.options().getString("precType") << "\n";

    //     if (params.options().getString("precType").substr(0, 3) == "PCD")
    //         NSsolver.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall, params.precOptions().getInt("pcd_bcType"));

    //     solveProblem(NSsolver, solveOpt);

    //     std::vector<index_t> itVector = NSsolver.getLinIterVector();

    //     gsInfo << "Iterations of linear solver in each Picard iteration:\n";
    //     for (size_t i = 0; i < itVector.size(); i++)
    //         gsInfo << itVector[i] << ", ";

    //     gsInfo << "\nAverage number of linear solver iterations per Picard iteration: " << NSsolver.getAvgLinIterations() << "\n";
        
    //     id++;
    // }

    // if (unsteadyIt)
    // {
    //     solveOpt.setInt("id", id);
    //     params.options().setReal("timeStep", timeStep);
    //     params.options().setInt("maxIt_picard", picardIt);
    //     params.options().setReal("tol_picard", picardTol);
    //     params.options().setInt("maxIt_lin", linIt);
    //     params.options().setReal("tol_lin", linTol);
    //     params.options().setString("precType", precond);
    //     // params.precOptions().setReal("gamma", 10); // parameter for AL preconditioner

    //     gsINSSolverUnsteadyIter<real_t, LinSolver > NSsolver(params);

    //     gsInfo << "\nSolving the unsteady problem with preconditioned GMRES as linear solver.\n";
    //     gsInfo << "Used preconditioner: " << params.options().getString("precType") << "\n";

    //     if (params.options().getString("precType").substr(0, 3) == "PCD")
    //         NSsolver.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall, params.precOptions().getInt("pcd_bcType"));

    //     solveProblem(NSsolver, solveOpt);
        
    //     std::vector<index_t> itVector = NSsolver.getLinIterVector();

    //     gsInfo << "Iterations of linear solver in each Picard iteration:\n";
    //     for (size_t i = 0; i < itVector.size(); i++)
    //         gsInfo << itVector[i] << ", ";

    //     gsInfo << "\nAverage number of linear solver iterations per Picard iteration: " << NSsolver.getAvgLinIterations() << "\n";

    //     id++;
    // }

    return 0; 
}


template<class T, int MatOrder>
void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo)
{
    gsStopwatch clock;

    if(geo == 2)
        markElimDof(NSsolver);

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

    //gsINSSolverUnsteady<real_t>* pSolver = dynamic_cast<gsINSSolverUnsteady<real_t>* >(&NSsolver);

    // if (pSolver)
    // {
    //     if (opt.getSwitch("stokesInit"))
    //         pSolver->setStokesInitialCondition();

    //     pSolver->solveWithAnimation(opt.getInt("maxIt"), 10, opt.getReal("tol"), opt.getInt("plotPts"));
    // }
    // else
    // {
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"));
    //}      

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time:" << NSsolver.getInitAssemblyTime() + NSsolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time:" << NSsolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time:" << totalT << "\n";

    if (opt.getSwitch("plot")) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);

        std::string geoStr;
        switch(opt.getInt("geo"))
        {
            case 1:
            {
                geoStr = "BFS";
                break;
            }
            case 2:
            {
                geoStr = "LDC";
                break;
            }
            case 3:
            {
                geoStr = "profile";
                break;
            }
            default:
                GISMO_ERROR("Unknown domain.");
        }

        int id = opt.getInt("id");
        int plotPts = opt.getInt("plotPts");
 
        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, geoStr + "_solver" + util::to_string(id) + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_solver" + util::to_string(id) + "_pressure", plotPts);
        // plotQuantityFromSolution("divergence", velocity, geoStr + "_solver" + util::to_string(id) + "_velocityDivergence", plotPts);
        gsInfo << " done.\n";
    }
}

// mark one pressure DoF in the domain corner as fixed zero (for the lid driven cavity problem)
template<class T, int MatOrder>
void markElimDof(gsINSSolver<T, MatOrder>& NSsolver)
{
    std::vector<gsMatrix<index_t> > elimDof(1);
    elimDof[0].setZero(1,1);

    NSsolver.markDofsAsEliminatedZeros(elimDof, 1);
}