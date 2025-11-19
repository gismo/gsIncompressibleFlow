/** @file gsRANSSolversExample.cpp
 
    This example demonstrates capabilities of the RANS solver
    for several predefined domains: backward-facing step (2D and 3D), lid-driven
    cavity (2D and 3D) and a 2D blade profile. The user can choose parameters
    of the domain, problem definition and computation (see help).

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

template<class T, int MatOrder> void solveProblem(gsRANSSolverUnsteady<T, MatOrder>& NSsolver, gsOptionList opt, int geo, gsFlowLogger& logger);

int main(int argc, char *argv[])
{
    typedef gsGMRes<real_t> LinSolver;

    // ========================================= Settings ========================================= 

    // lin. solvers
    bool direct = false;
    bool iter = false;

    // domain definition
    index_t geo = 1; // 1 - step, 2 - cavity, 3 - blade profile 2D
    index_t dim = 2;
    real_t a = 16;   // length of backward facing step domain behind step
    real_t b = 2;   // height of backward facing step domain
    real_t c = 2;   // depth of backward facing step domain (3D)
    real_t h = 1;   // height of the backward facing step
    real_t a_in = 3; // length of the inflow part of backward facing step domain
    real_t aa = 1;  // width of cavity
    real_t bb = 1;  // height of cavity
    real_t cc = 1;  // depth of cavity (3D)
    
    // discretization settings
    index_t numRefine = 3;
    index_t wallRefine = 0;
    index_t leadRefine = 0; // for profile2D
    int numElevate = 0; // number of degree elevations (before refinement)

    // problem parameters
    real_t viscosity = 0.0001;
    std::string UinStr = "default"; // inlet x-velocity for step (default = -4*(y-1.5)^2 + 1)
    std::string lidVelX = "1"; // x-velocity of the cavity lid
    std::string lidVelZ = "0"; // z-velocity of the cavity lid
    real_t inVelX = 1; // inlet x-velocity for profile2D
    real_t inVelY = 0.53; // inlet y-velocity for profile2D
    
    // solver settings
    index_t maxIt = 20;
    index_t picardIt = 5;
    index_t maxItTM = 10;
    index_t linIt = 50;
    real_t timeStep = 0.05;
    real_t tol = 1e-5;
    real_t picardTol = 1e-4;
    real_t linTol = 1e-6;
    std::string matFormation = "EbE";
    std::string precond = "MSIMPLER_FdiagEqual";
    bool stokesInit = false; // start unsteady problem from Stokes solution
    bool TCSD_RANS_stab = false; // use T-CSD stabilization
    bool TCSD_TM_stab = false; // use T-CSD stabilization

    // output settings
    std::string outMode = "terminal"; // terminal/file/all/quiet
    bool plot = false;
    bool plotMesh = false;
    index_t plotPts = 20000;
    bool animation = false;
    index_t animStep = 10;

    // OpenMP parallelization
    int numThreads = 1;

    // ---------------------------------------------------------------------------------

    //command line
    gsCmdLine cmd("Solves the RANS equations with turbulence model in a given domain (step, cavity, blade profile).");

    cmd.addSwitch("direct", "Solve with direct linear solver", direct);
    cmd.addSwitch("iter", "Solve with preconditioned GMRES as linear solver", iter);

    cmd.addInt("g", "geo", "Computational domain (1 - step, 2 - cavity, 3 - profile (only 2D))", geo);
    cmd.addInt("d", "dim", "Space dimension", dim);
    cmd.addReal("a", "stepA", "Backward-facing step domain: lenght behind step", a);
    cmd.addReal("", "stepAin", "Backward-facing step domain: lenght before step", a_in);
    cmd.addReal("b", "stepB", "Backward-facing step domain: total height", b);
    cmd.addReal("c", "stepC", "Backward-facing step domain (3D): domain depth", c);
    cmd.addReal("", "stepH", "Backward-facing step domain: step height", h);
    cmd.addReal("", "cavityA", "Lid-driven cavity domain: width", aa);
    cmd.addReal("", "cavityB", "Lid-driven cavity domain: height", bb);
    cmd.addReal("", "cavityC", "Lid-driven cavity domain (3D): depth", cc);

    cmd.addInt("r", "uniformRefine", "Number of uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("", "wallRefine", "Number of h-refinement steps near step corner, cavity walls of blade profile", wallRefine);
    cmd.addInt("", "leadRefine", "Number of h-refinement steps near the beginning of the blade (for profile geometry)", leadRefine);
    cmd.addInt("e", "degElevate", "Number of degree elevations (performed before h-refinement)", numElevate);

    cmd.addReal("v", "visc", "Viscosity value", viscosity);
    cmd.addString("", "stepUinX", "Backward-facing step: x-coordinate of inflow velocity (string expression)", UinStr);
    cmd.addString("", "lidVelX", "Cavity: x-coordinate of lid velocity (string expression)", lidVelX);
    cmd.addString("", "lidVelZ", "Cavity: z-coordinate of lid velocity (string expression)", lidVelZ);
    cmd.addReal("", "profUinX", "Blade profile: x-coordinate of inflow velocity", inVelX);
    cmd.addReal("", "profUinY", "Blade profile: y-coordinate of inflow velocity", inVelY);

    cmd.addInt("", "maxIt", "Max. number of Picard iterations or time steps", maxIt);
    cmd.addInt("", "picardIt", "Max. number of inner Picard iterations for unsteady problem", picardIt);
    cmd.addInt("", "maxItTM", "Max. number of turbulence model iterations", maxItTM);
    cmd.addInt("", "linIt", "Max. number of GMRES iterations (if the lin. systems are solved iteratively)", linIt);
    cmd.addReal("", "timeStep", "Time discretization step for unsteady problem", timeStep);
    cmd.addReal("", "tol", "Stopping tolerance", tol);
    cmd.addReal("", "picardTol", "Tolerance for inner Picard iteration for unsteady problem", picardTol);
    cmd.addReal("", "linTol", "Tolerance for iterative linear solver", linTol);
    cmd.addString("", "loop", "Matrix formation method (EbE = element by element, RbR = row by row)", matFormation);
    cmd.addString("p", "precond", "Preconditioner type (format: PREC_Fstrategy, PREC = {PCD, PCDmod, LSC, AL, SIMPLE, SIMPLER, MSIMPLER}, Fstrategy = {FdiagEqual, Fdiag, Fmod, Fwhole})", precond);
    cmd.addSwitch("stokesInit", "Set Stokes initial condition", stokesInit);
    cmd.addSwitch("TCSD_RANS", "Use T-CSD stabilization of numerical solution of RANS", TCSD_RANS_stab);
    cmd.addSwitch("TCSD_TM", "Use T-CSD stabilization of numerical solution of TM", TCSD_TM_stab);

    cmd.addString("o", "outMode", "Output mode (terminal/file/all/quiet)", outMode);
    cmd.addSwitch("plot", "Plot the final result in ParaView format", plot);
    cmd.addSwitch("plotMesh", "Plot the computational mesh", plotMesh);
    cmd.addInt("", "plotPts", "Number of sample points for plotting", plotPts);
    cmd.addSwitch("animation", "Plot animation of the unsteady problem", animation);
    cmd.addInt("", "animStep", "Number of iterations between screenshots for animation (used when animation = true)", animStep);

    cmd.addInt("t", "nthreads", "Number of threads", numThreads);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    
    if ( !(direct || iter) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --direct, --iter\n\n";

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;
    std::string geoStr;

    switch(geo)
    {
        case 1:
        {
            geoStr = util::to_string(dim) + "D backward-facing step";

            switch(dim)
            {
                case 2:
                default:
                    patches = BSplineStep2D<real_t>(1, a, b, a_in, h);
                    break;

                case 3:
                    patches = BSplineStep3D<real_t>(1, a, b, c, a_in, h);
                    break;
            }

            break;
        }
        case 2:
        {
            geoStr = util::to_string(dim) + "D lid-driven cavity";

            switch(dim)
            {
                case 2:
                default:
                    patches = BSplineCavity2D<real_t>(1, aa, bb);
                    break;

                case 3:
                    c = 1;
                    patches = BSplineCavity3D<real_t>(1, aa, bb, cc);
                    break;
            }

            break;
        }
        case 3:
        {
            geoStr = util::to_string(dim) + "D blade profile";

            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";

            gsReadFile<>(FLOW_DATA_DIR "geo_profile2D.xml", patches);
            break;
        }
        default:
            GISMO_ERROR("Unknown domain.");
    }

    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;
    defineBndParts(geo, dim, bndIn, bndOut, bndWall);

    // ========================================= Define problem and basis ========================================= 

    gsBoundaryConditions<> bcInfo;
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
    basis.degreeElevate(numElevate);

    switch(geo)
    {
        case 1:
        {
            defineBCs_step(bcInfo, dim, false, UinStr);
            refineBasis_step(basis, numRefine, 0, wallRefine, 0, 0, dim);
            break;
        }
        case 2:
        {
            defineBCs_cavity(bcInfo, dim, 1, lidVelX, lidVelZ);
            refineBasis_cavity(basis, numRefine, wallRefine, dim);
            break;
        }
        case 3:
        {
            defineBCs_profile2D(bcInfo, inVelX, inVelY);
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

    discreteBases.push_back(basis); // basis for k
    discreteBases.push_back(basis); // basis for omega

    // ========================================= Solve ========================================= 

    gsFlowLogger logger(gsFlowLogger::parseOutputMode(outMode), "gsRANSSolversExample.log");
    logger << "Solving RANS problem with k-omega SST model in " << geoStr << " domain.\n";
    logger << "domain: " << patches;
    logger << "viscosity = " << viscosity << "\n";
    logger << "source function = " << f << "\n";

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    
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
    solveOpt.addSwitch("TCSD_RANS", "", TCSD_RANS_stab);
    solveOpt.addSwitch("TCSD_TM", "", TCSD_TM_stab);
    solveOpt.addString("id", "", "");

    if (direct)
    {
        logger << "\n-----------------------------------------------------------\n";
        logger << "Solving the unsteady RANS problem with direct linear solver";
        logger << "\n-----------------------------------------------------------\n\n";

        solveOpt.setString("id", "rans_direct");

        gsFlowSolverParams<real_t> paramsDir(NSpde, discreteBases, &logger);
        paramsDir.options().setString("assemb.loop", matFormation);
        paramsDir.options().setInt("TM.maxIt", maxItTM);
        paramsDir.options().setReal("timeStep", timeStep);
        paramsDir.options().setInt("nonlin.maxIt", picardIt);
        paramsDir.options().setReal("nonlin.tol", picardTol);
        paramsDir.setBndParts(bndIn, bndOut, bndWall);
        paramsDir.options().setString("lin.solver", "direct");
        paramsDir.options().setInt("numThreads", numThreads);
        paramsDir.options().setSwitch("TCSD_RANS", TCSD_RANS_stab);
        paramsDir.options().setSwitch("TCSD_TM", TCSD_TM_stab);

        gsRANSSolverUnsteady<real_t, RowMajor> NSsolver(paramsDir);
        solveProblem(NSsolver, solveOpt, geo, logger);
    }

    if (iter)
    {
        logger << "\n--------------------------------------------------------------\n";
        logger << "Solving the unsteady RANS problem with iterative linear solver";
        logger << "\n--------------------------------------------------------------\n\n";

        solveOpt.setString("id", "rans_iter");

        gsFlowSolverParams<real_t> paramsIter(NSpde, discreteBases, &logger);
        paramsIter.options().setString("assemb.loop", matFormation);
        paramsIter.options().setInt("TM.maxIt", maxItTM);
        paramsIter.options().setReal("timeStep", timeStep);
        paramsIter.options().setInt("nonlin.maxIt", picardIt);
        paramsIter.options().setReal("nonlin.tol", picardTol);
        paramsIter.setBndParts(bndIn, bndOut, bndWall);
        paramsIter.options().setString("lin.solver", "iter");
        paramsIter.options().setInt("lin.maxIt", linIt);
        paramsIter.options().setReal("lin.tol", linTol);
        paramsIter.options().setString("lin.precType", precond);
        paramsIter.options().setInt("numThreads", numThreads);
        paramsIter.options().setSwitch("TCSD_RANS", TCSD_RANS_stab);
        paramsIter.options().setSwitch("TCSD_TM", TCSD_TM_stab);

        gsRANSSolverUnsteady<real_t, RowMajor> NSsolver(paramsIter);
        solveProblem(NSsolver, solveOpt, geo, logger);
    }

    return 0; 
}


template<class T, int MatOrder>
void solveProblem(gsRANSSolverUnsteady<T, MatOrder>& NSsolver, gsOptionList opt, int geo, gsFlowLogger& logger)
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
                geoStr = "profile2D";
                break;
            }
            default:
                GISMO_ERROR("Unknown domain.");
        }
    }

    // ------------------------------------
    // solve problem

    logger << "\nInitialization...\n";
    NSsolver.initialize();

    logger << "RANS numDofs: " << NSsolver.numDofs() << "\n";
    logger << "TM numDofs: " << NSsolver.numDofsTM() << "\n\n";

    if (opt.getSwitch("stokesInit"))
        NSsolver.solveStokes();

    if (opt.getSwitch("animation"))
        NSsolver.solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"), true, 1);
    else
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0);    

    real_t totalT = clock.stop();

    logger << "\nAssembly time:" << NSsolver.getAssemblyTime() << "\n";
    logger << "Solve time:" << NSsolver.getSolveTime() << "\n";
    logger << "Solver setup time:" << NSsolver.getSolverSetupTime() << "\n";
    logger << "Total TM time: " << NSsolver.getTMtime() << "\n";
    logger << "Total solveProblem time:" << totalT << "\n\n";

    logger << "RANS lin. it.:\n";
    NSsolver.getLinSolver()->reportLinIterations();

    logger << "\nTM lin. it.:\n";
    NSsolver.getTMsolver().getLinSolver()->reportLinIterations();

    // ------------------------------------
    // plot

    if (opt.getSwitch("plot")) 
    {
        gsField<> velocity = NSsolver.constructSolution(0);
        gsField<> pressure = NSsolver.constructSolution(1);
        gsField<> ksol = NSsolver.constructSolution(2);
        gsField<> omegasol = NSsolver.constructSolution(3);

        int plotPts = opt.getInt("plotPts");
 
        logger << "Plotting in Paraview... ";
        gsWriteParaview<>(velocity, geoStr + "_" + id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_" + id + "_pressure", plotPts);
        gsWriteParaview<>(ksol, geoStr + "_" + id + "_k", plotPts);
        gsWriteParaview<>(omegasol, geoStr + "_" + id + "_omega", plotPts);
        logger << "Done.\n";
    }
}