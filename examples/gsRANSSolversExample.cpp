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
    std::string inputFile = "";
    
    // discretization settings
    index_t numRefine = 4;
    index_t wallRefine = 0;
    index_t leadRefine = 0; // for profile2D
    int numElevate = 0; // number of degree elevations (before refinement)

    // problem parameters
    real_t viscosity = 0.00001;
    real_t inVelX = 1; // inlet x-velocity for profile2D
    real_t inVelY = 0; // inlet y-velocity for profile2D
    
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

    // output settings
    bool quiet = false;
    bool plot = false;
    bool plotMesh = false;
    index_t plotPts = 10000;
    bool animation = false;
    index_t animStep = 5;

    // ---------------------------------------------------------------------------------

    //command line
    gsCmdLine cmd("Solves the RANS equations with turbulence model in a given domain (step, cavity, blade profile).");

    cmd.addSwitch("direct", "Solve with direct linear solver", direct);
    cmd.addSwitch("iter", "Solve with preconditioned GMRES as linear solver", iter);

    cmd.addInt("g", "geo", "Computational domain (1 - step, 2 - cavity, 3 - profile (only 2D))", geo);
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
    cmd.addInt("", "maxItTM", "Max. number of turbulence model iterations", maxItTM);
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
    
    if (!inputFile.empty())
        geo = 0;

    if ( !(direct || iter) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --direct, --iter\n\n";

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
            inputFile = fn;
            break;

        case 2:
            geoStr = "LDC" + util::to_string(dim) + "D";
            fn = geoStr + "_problem.xml";
            inputFile = fn;
            break;

        case 3:
            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";

            geoStr = "profile2D";
            inputFile = geoStr + "_problem.xml";
            break;
    }

    gsInfo << "Reading problem definition from file:\n" << inputFile << "\n\n";

    std::string path = gsFileManager::find(inputFile);
    if ( path.empty() )
    {
        gsWarn<<"Input file not found, quitting.\n";
        return 1;
    }

    gsFileData<> fd(inputFile);
    fd.getId(0, patches);   // id=0: multipatch domain
    fd.getId(1, f);         // id=1: source function
    fd.getId(2, bcInfo);    // id=2: boundary conditions

    gsInfo << "Solving RANS problem with k-omega SST model in " << geoStr << " domain.\n";
    gsInfo << patches;
    gsInfo << "viscosity = " << viscosity << "\n";
    gsInfo << "source function = " << f << "\n";

    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;
    defineBndParts(geo, dim, bndIn, bndOut, bndWall);


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

    discreteBases.push_back(basis); // basis for k
    discreteBases.push_back(basis); // basis for omega

    // ========================================= Solve ========================================= 

    gsNavStokesPde<real_t> NSpde(patches, bcInfo, &f, viscosity);
    gsFlowSolverParams<real_t> params(NSpde, discreteBases);
    params.options().setSwitch("quiet", quiet);
    params.options().setString("assemb.loop", matFormation);
    params.options().setInt("TM.maxIt", maxItTM);
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

    solveOpt.setString("id", "unsteady");
    params.options().setReal("timeStep", timeStep);
    params.options().setInt("nonlin.maxIt", picardIt);
    params.options().setReal("nonlin.tol", picardTol);

    if (direct)
    {
        solveOpt.setString("id", "direct");
        params.options().setString("lin.solver", "direct");

        gsRANSSolverUnsteady<real_t, RowMajor> NSsolver(params);

        gsInfo << "\n-----------------------------------------------------------\n";
        gsInfo << "Solving the unsteady RANS problem with direct linear solver\n";
        gsInfo << "-----------------------------------------------------------\n";

        solveProblem(NSsolver, solveOpt, geo);
    }

    if (iter)
    {
        solveOpt.setString("id", "iter");
        params.options().setString("lin.solver", "iter");
        params.options().setInt("lin.maxIt", linIt);
        params.options().setReal("lin.tol", linTol);
        params.options().setString("lin.precType", precond);

        gsRANSSolverUnsteady<real_t, RowMajor> NSsolver(params);

        gsInfo << "\n-----------------------------------------------------------\n";
        gsInfo << "Solving the unsteady RANS problem with iterative linear solver\n";
        gsInfo << "-----------------------------------------------------------\n";

        solveProblem(NSsolver, solveOpt, geo);
    }


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
            case 0:
                geoStr = "customGeo";
                break;
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

    gsInfo << "\nInitialization...\n";
    NSsolver.initialize();

    gsInfo << "RANS numDofs: " << NSsolver.numDofs() << "\n";
    gsInfo << "TM numDofs: " << NSsolver.numDofsTM() << "\n\n";

    if (opt.getSwitch("stokesInit"))
        NSsolver.solveStokes();

    if (opt.getSwitch("animation"))
        NSsolver.solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"), true, 1);
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
        gsField<> ksol = NSsolver.constructSolution(2);
        gsField<> omegasol = NSsolver.constructSolution(3);

        int plotPts = opt.getInt("plotPts");
 
        gsInfo << "Plotting in Paraview...";
        gsWriteParaview<>(velocity, geoStr + "_" + id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_" + id + "_pressure", plotPts);
        gsWriteParaview<>(ksol, geoStr + "_" + id + "_k", plotPts);
        gsWriteParaview<>(omegasol, geoStr + "_" + id + "_omega", plotPts);
        gsInfo << "Done.\n";
    }
}