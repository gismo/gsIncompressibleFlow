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
#include <gsHLBFGS/gsHLBFGS.h>

using namespace gismo;

template<class T, int MatOrder> void solveProblem(gsINSSolver<T, MatOrder>& NSsolver, gsOptionList opt, int geo);

// reset the solver with a new mesh
template<class T, int MatOrder>
void resetSolverWithNewMesh(gsINSSolverUnsteady<T, MatOrder>*& solver,
                            const gsMultiPatch<T>& newPatches);

// transfer the solution to a new mesh
template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField);

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
    int geo = 1; // 0 - custom input file, 1 - step, 2 - cavity, 3 - blade profile 2D, 4 - flapping beam
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
    int maxIt = 3;
    int picardIt = 5;
    int linIt = 50;
    real_t timeStep = 0.01;
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

    cmd.addInt("g", "geo", "Computational domain (0 - custom file, 1 - step, 2 - cavity, 3 - profile (only 2D), 4 - flapping beam)", geo);
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

    if (!inputFile.empty())
        geo = 0;

    if ( !(steady || steadyIt || unsteady || unsteadyIt) )
        gsWarn << "All computation flags set to false - nothing will be computed.\nPlease select at least one of the flags: --steady, --steadyIt, --unsteady, --unsteadyIt\n\n";

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
            
        case 4:
            if (dim == 3)
                gsWarn << "Flapping beam geometry is only 2D!\n";
                
            geoStr = "flappingBeam";
            inputFile = geoStr + "_flow.xml";
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
            
        case 4:
            // refine the flapping beam appropriately
            for (int r = 0; r < numRefine; ++r)
                basis.uniformRefine();
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
    solveOpt.addInt("nonlin.maxIt", "", picardIt);
    solveOpt.addReal("tol", "", tol);
    solveOpt.addReal("picardTol", "", picardTol);
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

        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver() );
        reportLinIterations(linSolverPtr);
    }

    if (unsteady)
    {
        solveOpt.setString("id", "unsteady");
        params.options().setReal("timeStep", timeStep);
        params.options().setInt("nonlin.maxIt", picardIt);
        params.options().setReal("nonlin.tol", picardTol);
        params.options().setString("lin.solver", "direct");

        gsINSSolverUnsteady<real_t, ColMajor> NSsolver(params);

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
        
        gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* linSolverPtr = dynamic_cast<gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<> >* >( NSsolver.getLinSolver() );
        reportLinIterations(linSolverPtr);
    }

    return 0; 
}

// ... 然后是 solveProblem 的定义 ...
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
                break;
            case 1:
                geoStr = "BFS" + dimStr;
                break;
            case 2:
                geoStr = "LDC" + dimStr;
                break;
            case 3:
                geoStr = "profile2D";
                break;
            case 4:
                geoStr = "flapping_beam" + dimStr;
                break;
            default:
                // 如果没有匹配的几何体ID，使用一个通用名称
                geoStr = "geo" + util::to_string(opt.getInt("geo")) + "_" + dimStr;
                break;
        }
    }

    // ------------------------------------
    // solve problem

    gsInfo << "\ninitialization...\n";
    NSsolver.initialize();

    gsInfo << "numDofs: " << NSsolver.numDofs() << "\n";

    gsINSSolverUnsteady<T, MatOrder>* pSolver = dynamic_cast<gsINSSolverUnsteady<T, MatOrder>*>(&NSsolver);

    
    if (pSolver)
    {
        if (opt.getSwitch("stokesInit"))
            pSolver->solveStokes();
        if (opt.getSwitch("animation"))
        {
            pSolver->solveWithAnimation(opt.getInt("maxIt"), opt.getInt("animStep"), 
                                        geoStr + "_" + id, opt.getReal("tol"), opt.getInt("plotPts"));
        }
        else
        {
            // create ParaView files for time step collection
            std::string fileNameU = geoStr + "_" + id + "_velocity_animation.pvd";
            std::ofstream fileU(fileNameU.c_str());
            GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

            std::string fileNameP = geoStr + "_" + id + "_pressure_animation.pvd";
            std::ofstream fileP(fileNameP.c_str());
            GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

            // use the existing functions in gsFlowUtils.h
            gismo::startAnimationFile(fileU);
            gismo::startAnimationFile(fileP);
            
            // determine the animation step
            int animStep = 5; 
            try {
                animStep = opt.getInt("animStep");
            } catch (...) {
                gsInfo << "Using default animation step: " << animStep << "\n";
            }
            
            // save the initial state
            if (opt.getSwitch("plot"))
            {
                gsField<T> velocity = pSolver->constructSolution(0);
                gsField<T> pressure = pSolver->constructSolution(1);
                
                std::string velocityBaseName = geoStr + "_" + id + "_velocity_step0";
                std::string pressureBaseName = geoStr + "_" + id + "_pressure_step0";
                
                gsWriteParaview<>(velocity, velocityBaseName, opt.getInt("plotPts"), opt.getSwitch("plotMesh"));
                gsWriteParaview<>(pressure, pressureBaseName, opt.getInt("plotPts"));
                
                // add references to the collection file
                int numPatches = NSsolver.getParams()->getPde().patches().nPatches();
                for (int p = 0; p < numPatches; p++)
                {
                    fileU << "<DataSet timestep=\"0\" part=\"" << p 
                          << "\" file=\"" << velocityBaseName << p << ".vts\"/>\n";
                    fileP << "<DataSet timestep=\"0\" part=\"" << p 
                          << "\" file=\"" << pressureBaseName << p << ".vts\"/>\n";
                }
            }
            
            // create temporary variables to store the previous solution
            gsVector<T> prevVelocityCoefs;
            gsVector<T> prevPressureCoefs;
            
            // get the initial coefficients
            if (opt.getSwitch("stokesInit")) {
                prevVelocityCoefs = pSolver->solutionCoefs(0);
                prevPressureCoefs = pSolver->solutionCoefs(1);
            }
            
            for (int i = 0; i < opt.getInt("maxIt"); ++i)
            {
                // 如果不是第一步，使用上一步的解作为初始条件
                if (i > 0)
                {
                    pSolver->setSolutionCoefs(prevVelocityCoefs, 0);
                    pSolver->setSolutionCoefs(prevPressureCoefs, 1);

                }
                
                // perform a time step
                pSolver->nextIteration();
                
                // get the coefficients of the current solution
                gsVector<T> currentVelocityCoefs = pSolver->solutionCoefs(0);
                gsVector<T> currentPressureCoefs = pSolver->solutionCoefs(1);

                gsField<T> velocity = pSolver->constructSolution(0);

                gsField<T> pressure = pSolver->constructSolution(1);
                // Change to new mesh
                const gsMultiPatch<T>& currentPatches = NSsolver.getParams()->getPde().patches();
                gsMultiPatch<T> newPatches = currentPatches;
                    
                newPatches.uniformRefine();
                gsInfo << "Degree of freedom last step " << currentPatches.coefs().dim() << "\n";
                gsInfo << "Degree of freedom new step " << newPatches.coefs().dim() << "\n";

                resetSolverWithNewMesh(pSolver, newPatches);
                
                // 计算解的变化 (仅在第二步及以后)
                T relChange = 0;
                if (i > 0)
                {
                    // 计算速度和压力变化的平方和
                    gsVector<T> velocityDiff = currentVelocityCoefs - prevVelocityCoefs;
                    gsVector<T> pressureDiff = currentPressureCoefs - prevPressureCoefs;
                    
                    T velocityNorm = velocityDiff.norm();
                    T pressureNorm = pressureDiff.norm();
                    T totalNorm = std::sqrt(velocityNorm*velocityNorm + pressureNorm*pressureNorm);
                    
                    T prevNorm = std::sqrt(prevVelocityCoefs.squaredNorm() + prevPressureCoefs.squaredNorm());
                    if (prevNorm > 1e-10)
                        relChange = totalNorm / prevNorm;
                    else
                        relChange = totalNorm;
                    
                    gsInfo << "Time step " << i+1 << " solution change: " << relChange << "\n";
                }
                
                prevVelocityCoefs = currentVelocityCoefs;
                prevPressureCoefs = currentPressureCoefs;
                
                if (opt.getSwitch("plot"))
                {
                    // 使用最新的解决方案字段
                    try {
                        gsField<T> velocity = pSolver->constructSolution(0);
                        gsField<T> pressure = pSolver->constructSolution(1);
                        
                        std::string velocityBaseName = geoStr + "_" + id + "_velocity_step" + util::to_string(i+1);
                        std::string pressureBaseName = geoStr + "_" + id + "_pressure_step" + util::to_string(i+1);
                        
                        gsWriteParaview<>(velocity, velocityBaseName, opt.getInt("plotPts"), opt.getSwitch("plotMesh"));
                        gsWriteParaview<>(pressure, pressureBaseName, opt.getInt("plotPts"));
                        
                        // 向集合文件添加引用
                        int numPatches = pSolver->getParams()->getPde().patches().nPatches();
                        for (int p = 0; p < numPatches; p++)
                        {
                            fileU << "<DataSet timestep=\"" << i+1 << "\" part=\"" << p 
                                  << "\" file=\"" << velocityBaseName << p << ".vts\"/>\n";
                            fileP << "<DataSet timestep=\"" << i+1 << "\" part=\"" << p 
                                  << "\" file=\"" << pressureBaseName << p << ".vts\"/>\n";
                        }
                    } catch (const std::exception& e) {
                        gsInfo << "Error " << e.what() << "\n";
                    }
                }
            }
            
            // 关闭ParaView集合文件
            gismo::endAnimationFile(fileU);
            gismo::endAnimationFile(fileP);
            
            gsInfo << "\nSimulation completed " << opt.getInt("maxIt") << " time steps.\n";
            gsInfo << "ParaView animation files created: \n";
            gsInfo << "  " << fileNameU << "\n";
            gsInfo << "  " << fileNameP << "\n";
        }
    }
    else
    {
        NSsolver.solve(opt.getInt("maxIt"), opt.getReal("tol"), 0);
    }

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
 
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, geoStr + "_" + id + "_velocity", plotPts, opt.getSwitch("plotMesh"));
        gsWriteParaview<>(pressure, geoStr + "_" + id + "_pressure", plotPts);
        // plotQuantityFromSolution("divergence", velocity, geoStr + "_" + id + "_velocityDivergence", plotPts);
    }
}
template<class T, int MatOrder>
void resetSolverWithNewMesh(gsINSSolverUnsteady<T, MatOrder>*& solver,
                            const gsMultiPatch<T>& newPatches)
{
    // 1. save the current solver parameters
    gsFlowSolverParams<T>* currentParams = solver->getParams().get();
    gsOptionList solverOptions = currentParams->options();
    T viscosity = currentParams->getPde().viscosity();
    const gsFunction<T>* f = currentParams->getPde().rhs();
    gsBoundaryConditions<> bcInfo = currentParams->getPde().bc();

    // 2. save the current solution (as fields)
    gsField<T> velocityFieldOld = solver->constructSolution(0);
    gsField<T> pressureFieldOld = solver->constructSolution(1);
    
    gsInfo << "Current dofs: " << solver->numDofs() << "\n";
    
    // 3. create new basis functions for the new geometry
    std::vector<gsMultiBasis<T>> currentBases = solver->getAssembler()->getBases();
    std::vector<gsMultiBasis<T>> newBases;
    
    for (size_t i = 0; i < currentBases.size(); ++i)
    {
        gsMultiBasis<T> newBasis(newPatches);
        
        // match the degree of the original basis functions
        for (size_t p = 0; p < newBasis.nBases(); ++p)
        {
            if (p < currentBases[i].nBases())
                newBasis.basis(p).setDegree(currentBases[i].basis(p).degree(0));
        }
        // apply refinement
        for (int r = 0; r < 3; ++r)
            newBasis.uniformRefine();
        
        newBases.push_back(newBasis);
    }

    // ensure Taylor-Hood elements (velocity one order higher than pressure)
    if (newBases.size() >= 2)
        newBases[0].degreeElevate(1);
    
    // 4. create a new PDE and solver parameters
    gsNavStokesPde<T> newPde(newPatches, bcInfo, f, viscosity);
    
    typename gsFlowSolverParams<T>::Ptr newParamsPtr =
        std::make_shared<gsFlowSolverParams<T>>(newPde, newBases);
    newParamsPtr->options() = solverOptions;
    
    // 5. create a new solver instance
    auto* newSolver = new gsINSSolverUnsteady<T, MatOrder>(newParamsPtr);
    
    // initialize the new solver
    newSolver->initialize();
    
    // print the number of dofs
    gsInfo << "New dofs: " << newSolver->numDofs() << "\n";
    
    // 6. project the old solution to the new mesh
    transferSolutionToNewMesh(newSolver, velocityFieldOld, pressureFieldOld);
    
    // safely replace the solver
    solver = newSolver;
}

template<class T, int MatOrder>
void transferSolutionToNewMesh(gsINSSolverUnsteady<T, MatOrder>* solver,
                              const gsField<T>& velocityField,
                              const gsField<T>& pressureField)
{
    // Get the new basis functions
    const std::vector<gsMultiBasis<T>>& newBases = solver->getAssembler()->getBases();
    const gsMultiPatch<T>& newPatches = solver->getParams()->getPde().patches();

    gsDebugVar(solver->numDofs());
    
    // get the target dimension and the number of dofs
    const index_t vDim = velocityField.function(0).targetDim();
    const index_t pDim = pressureField.function(0).targetDim();
    const index_t fullUdofs = solver->getAssembler()->getUdofs();
    const index_t tarDim = solver->getParams()->getPde().domain().geoDim();
    const index_t pShift = tarDim * fullUdofs;
    
    // 创建与setSolutionCoefs期望的尺寸匹配的矩阵
    gsMatrix<T> newVelocityCoefs(pShift, 1);
    gsMatrix<T> newPressureCoefs(solver->numDofs() - pShift, 1);
    
    newVelocityCoefs.setZero();
    newPressureCoefs.setZero();
    
    // 对速度场进行准插值
    gsQuasiInterpolate<T> quasiInterp;
    
    // get the mappers
    const std::vector<gsDofMapper>& mappers = solver->getAssembler()->getMappers();
    
    // interpolate for each patch
    for (size_t p = 0; p < newPatches.nPatches(); ++p)
    {
        // 获取当前patch的基函数
        const gsBasis<T>& vBasis = newBases[0].basis(p);
        
        for (size_t i = 0; i < vBasis.size(); ++i)
        {
            try {
                // Interpolate for each basis function
                gsMatrix<T> coef = quasiInterp.localIntpl(vBasis, velocityField.function(p), i);
                
                // Store the coefficients to the global coefficient matrix - use the mapper to get the global index
                if (mappers[0].is_free(i, p)) {
                    index_t globalIndex = mappers[0].index(i, p);
                    

                    for (index_t d = 0; d < vDim; ++d) {
                        newVelocityCoefs(globalIndex + d * fullUdofs, 0) = coef(0, d);
                    }
                }
            } catch (const std::exception& e) {
                gsInfo << "Velocity field interpolation error: " << e.what() << " on patch " << p << " basis function " << i << "\n";
            }
        }
    }
    
    // Interpolate the pressure field
    for (size_t p = 0; p < newPatches.nPatches(); ++p)
    {
        // Get the basis function of the current patch
        const gsBasis<T>& pBasis = newBases[1].basis(p);
        
        for (size_t i = 0; i < pBasis.size(); ++i)
        {
            try {
                // Interpolate for each basis function
                gsMatrix<T> coef = quasiInterp.localIntpl(pBasis, pressureField.function(p), i);
                
                // Store the coefficients to the global coefficient matrix - use the mapper to get the global index
                if (mappers[1].is_free(i, p)) {
                    index_t globalIndex = mappers[1].index(i, p);
                    newPressureCoefs(globalIndex, 0) = coef(0, 0);
                }
            } catch (const std::exception& e) {
                gsInfo << "Pressure field interpolation error: " << e.what() << " on patch " << p << " basis function " << i << "\n";
            }
        }
    }
    
    // set the projected coefficients to the new solution
    solver->setSolutionCoefs(newVelocityCoefs, 0);
    solver->setSolutionCoefs(newPressureCoefs, 1);

    gsDebugVar(newVelocityCoefs.size());
    gsDebugVar(newPressureCoefs.size());
    gsField<T> velocityResult = solver->constructSolution(0);
    gsField<T> pressureResult = solver->constructSolution(1);

    gsWriteParaview<>(velocityResult, "velocityResult", 1000, false);
    gsWriteParaview<>(pressureResult, "pressureResult", 1000, false);

    
    gsInfo << "Projected solution to new mesh\n";
}