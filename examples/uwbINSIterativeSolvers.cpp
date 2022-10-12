/** @file backwardStep.cpp

@brief Stokes example for uwb experimental purposes.

Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>

#include "uwbGeometryUtils.h"

#include "uwbINSSolverSteady.h"
#include "uwbINSSolverUnsteady.h"
#include "uwbRANSSolver.h"
#include "uwbTMSolverKOmega.h"

#include "uwbINSSolverSteadyIterative.h"
#include "uwbINSSolverUnsteadyIterative.h"
#include "uwbRANSSolverIterative.h"

#include "uwbLinSolvers.h"

using namespace gismo;

template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, const T kIn, const T kWall, const T oIn, const T oWall);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 

    bool plot = false;
    bool outFile = false;
    bool dg = false; // use DGFEM to connect patches
    int numRefine = 2;
    int plot_pts = 10000;

    int deg = 3;
    bool periodic = true;

    real_t viscosity = 0.01;
    real_t timeStep = 0.1; // time step for unsteady computation
    //int numThreads = 10; // number of threads for assembly
    int numIter = 5;

    real_t uMax = 2;
    real_t turbIntensity = 0.05;
    real_t kConst = 1e-6; // k on the wall
    real_t oConst = 35; // omega on the wall
    bool TMsupg = false;

    // preconditioner type options:
    // LSC_AdiagEqual, LSC_Adiag, LSC_Awhole, LSC_Amod
    // PCD_AdiagEqual, PCD_Adiag, PCD_Awhole, PCD_Amod
    // AL_Awhole, AL_Amod
    // SIMPLE_AdiagEqual, SIMPLE_Adiag, SIMPLE_Awhole, SIMPLE_Amod
    // SIMPLER_AdiagEqual, SIMPLER_Adiag, SIMPLER_Awhole, SIMPLER_Amod
    // MSIMPLER_AdiagEqual, MSIMPLER_Adiag, MSIMPLER_Awhole, MSIMPLER_Amod
    std::string precType = "LSC_AdiagEqual";
    bool precIter = false; // set direct/iterative method to solve subsystems in preconditioner

    int linMaxIt = 100;
    real_t gammaSt = 0.1; // parameter for AL precond, steady case
    real_t gammaUnst = 2.2; // parameter for AL precond, unsteady case
    real_t gammaRans = 2.2; // parameter for AL precond, RANS case
    real_t linTol = 1e-4;

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("file", "Generate output file", outFile);
    cmd.addSwitch("dg", "Use DGFEM to connect patches", dg);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("s", "plotSamples",
        "Number of sample points to use for plotting", plot_pts);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine < 0)
    {
        gsInfo << "Number of refinements must be non-negative, quitting.\n";
        return -1;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    if (dg)
        opt.intStrategy = iFace::dg;
    else
        opt.intStrategy = iFace::glue;

    gsInfo << "Solving the backward facing step flow example.\n";

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;

    defineBCs_step(bcInfo, bndIn, bndOut, bndWall, 3, periodic);

    gsFunctionExpr<> f("0", "0", "0", 3);
       
    // ========================================= Define geometry ========================================= 

    gsMultiPatch<> patches;
    real_t a = 6;
    real_t b = 2;
    real_t c = 2;
    real_t a_in = 1;
 

    patches = BSplineStep3D<real_t>(deg, a, b, c, a_in, periodic);
    gsInfo << patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry    

    gsMultiBasis<> tbasis(patches);

    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure      
    discreteBases[0].degreeElevate(1);//p-refine the velocity space (Taylor-Hood type)

    // ========================================= Turbulence solver ========================================= 
    gsBoundaryConditions<> bcInfoTurb;

    real_t kInConst = 1.5 * math::pow(uMax * turbIntensity, 2); // (3/2)*(UI)^2
    real_t oInConst = math::pow(0.09, -0.25) * math::sqrt(kInConst) / 0.07; // (C_mu)^(-1/4) * sqrt(k)/l, l = 0.07 * d
    defineBCs_TM(bcInfoTurb, kInConst, kConst, oInConst, oConst);

    gsDofMapper koMapper;
    discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, discreteBases, opt);

    koParams.settings().set(constantsINS::TMsupg, TMsupg);

    uwbTMSolverKOmega<real_t> koSolver(koParams);
    uwbTMSolverKOmega<real_t> koSolverIt(koParams);

    gsMatrix<> koInitial(koSolver.getAssembler()->numVarDofs(), 2);
    koInitial.col(0).setConstant(kConst);
    koInitial.col(1).setConstant(oConst);

    koSolver.setInitialCondition(koInitial);
    koSolverIt.setInitialCondition(koInitial);

    const gsMatrix<> tmp = koInitial;
    gsMatrix<> tmp1 = -1 * tmp;

    // ========================================= Define solver ========================================= 
  
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> paramsSt(NSpde, discreteBases, opt);
    
    paramsSt.settings().setPrecondType(precType);
    paramsSt.getPrecOptions().setReal("gamma",gammaSt);
    paramsSt.settings().set(constantsINS::iter_maxIt, linMaxIt);
    paramsSt.settings().set(constantsINS::iter_tol, linTol);

    uwbINSSolverParams<real_t> paramsUnst = paramsSt;
    paramsUnst.settings().set(constantsINS::timeStep, timeStep);
    paramsUnst.getPrecOptions().setReal("gamma", gammaUnst);

    uwbINSSolverParams<real_t> paramsRans = paramsUnst;
    paramsRans.getPrecOptions().setReal("gamma", gammaRans);
    paramsRans.settings().setTurbulenceSolver(&koSolver);

    uwbINSSolverParams<real_t> paramsRansIt = paramsUnst;
    paramsRansIt.getPrecOptions().setReal("gamma", gammaRans);
    paramsRansIt.settings().setTurbulenceSolver(&koSolverIt);
    paramsRansIt.getPrecOptions().setSwitch("iter", precIter);

    uwbINSSolverSteady<real_t> solverSteady(paramsSt);
    uwbINSSolverSteadyIterative<real_t, uwbGMResRight<real_t> > solverSteadyIt(paramsSt);

    uwbINSSolverUnsteady<real_t> solverUnsteady(paramsUnst);
    uwbINSSolverUnsteadyIterative<real_t, uwbGMResRight<real_t> > solverUnsteadyIt(paramsUnst);

    uwbRANSSolver<real_t> solverRans(paramsRans);
    uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> >  solverRansIt(paramsRansIt);

    // ========================================= Solving ========================================= 
    gsInfo << "numDofs: " << solverSteady.numDofs() << "\n";

    gsStopwatch time;

    gsInfo << "initialization...\n";
    time.restart();
    solverSteady.initialize();

    gsInfo << "Solving steady with direct method:\n";
    solverSteady.solve(numIter, 1e-8);
    real_t timeSt = time.stop();

    time.restart();

    if (precType.substr(0, 3) == "PCD")
        solverSteadyIt.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall);

    solverSteadyIt.initialize();

    gsInfo << "\nSolving steady with iterative method:\n";
    solverSteadyIt.solve(numIter, 1e-8);
    real_t timeStIt = time.stop();

    gsInfo << "Steady (direct) total time: " << timeSt << "\n";
    gsInfo << "Steady (iter) total time: " << timeStIt << "\n";
    gsInfo << "Assembly time (direct): " << solverSteady.getAssemblyTime() << "\n";
    gsInfo << "Assembly time (iter): " << solverSteadyIt.getAssemblyTime() << "\n";
    gsInfo << "Solver setup time (direct): " << solverSteady.getSolverSetupTime() << "\n";
    gsInfo << "Solver setup time (iter): " << solverSteadyIt.getSolverSetupTime() << "\n";
    gsInfo << "Solve time (direct): " << solverSteady.getSolveTime() << "\n";
    gsInfo << "Solve time (iter): " << solverSteadyIt.getSolveTime() << "\n\n";

    real_t errSt = (solverSteady.getSolution() - solverSteadyIt.getSolution()).norm() / solverSteady.getSolution().norm();
    gsInfo << "Steady solution error (direct vs. iterative): " << errSt << "\n\n";

    gsInfo << "Avg lin. iter. per Picard iter.: " << solverSteadyIt.getAvgLinIterations() << "\n\n";

    gsInfo << "initialization...\n";
    time.restart();
    solverUnsteady.initialize(); // unsteady solver

    gsInfo << "Solving unsteady with direct method:\n";
    solverUnsteady.solve(numIter, 1e-8);
    real_t timeUnst = time.stop();

    gsInfo << "initialization...\n";
    time.restart();

    if (precType.substr(0, 3) == "PCD")
        solverUnsteadyIt.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall);

    solverUnsteadyIt.initialize();

    gsInfo << "\nSolving unsteady with iterative method:\n";
    solverUnsteadyIt.solve(numIter, 1e-8);
    real_t timeUnstIt = time.stop();

    gsInfo << "Unsteady (direct) total time: " << timeUnst << "\n";
    gsInfo << "Unsteady (iter) total time: " << timeUnstIt << "\n";
    gsInfo << "Assembly time (direct): " << solverUnsteady.getAssemblyTime() << "\n";
    gsInfo << "Assembly time (iter): " << solverUnsteadyIt.getAssemblyTime() << "\n";
    gsInfo << "Solver setup time (direct): " << solverUnsteady.getSolverSetupTime() << "\n";
    gsInfo << "Solver setup time (iter): " << solverUnsteadyIt.getSolverSetupTime() << "\n";
    gsInfo << "Solve time (direct): " << solverUnsteady.getSolveTime() << "\n";
    gsInfo << "Solve time (iter): " << solverUnsteadyIt.getSolveTime() << "\n\n";

    real_t errUnst = (solverUnsteady.getSolution() - solverUnsteadyIt.getSolution()).norm() / solverUnsteady.getSolution().norm();
    gsInfo << "Unsteady solution error (direct vs. iterative): " << errUnst << "\n\n";

    gsInfo << "Avg lin. iter. per Picard iter.: " << solverUnsteadyIt.getAvgLinIterations() << "\n\n";

    gsInfo << "initialization...\n";
    time.restart();
    solverRans.initialize(); // rans solver

    gsInfo << "Solving rans with direct method:\n";
    solverRans.solve(numIter, 1e-8);
    real_t timeRans = time.stop();

    gsInfo << "initialization...\n";
    time.restart();

    if (precType.substr(0, 3) == "PCD")
        solverRansIt.getAssembler()->preparePCDboundary(bndIn, bndOut, bndWall);

    solverRansIt.initialize();

    gsInfo << "\nSolving rans with iterative method:\n";
    solverRansIt.solve(numIter, 1e-8);
    real_t timeRansIt = time.stop();

    gsInfo << "RANS total time (direct): " << timeRans << "\n";
    gsInfo << "RANS total time (iter): " << timeRansIt << "\n";
    gsInfo << "Assembly time (direct): " << solverRans.getAssemblyTime() << "\n";
    gsInfo << "Assembly time (iter): " << solverRansIt.getAssemblyTime() << "\n";
    gsInfo << "Solver setup time (direct): " << solverRans.getSolverSetupTime() << "\n";
    gsInfo << "Solver setup time (iter): " << solverRansIt.getSolverSetupTime() << "\n";
    gsInfo << "Solve time (direct): " << solverRans.getSolveTime() << "\n";
    gsInfo << "Solve time (iter): " << solverRansIt.getSolveTime() << "\n\n";

    gsInfo << "Turb. model time (direct, should be cca equal for iter): " << solverRans.getTurbModelTime() << "\n";
    gsInfo << "Assembly time (turb direct): " << solverRans.getTurbAssemblyTime() << "\n";
    gsInfo << "Solver setup time (turb direct): " << solverRans.getTurbSolverSetupTime() << "\n";
    gsInfo << "Solve time (turb direct): " << solverRans.getTurbSolveTime() << "\n\n";
    /*gsInfo << "Turb. model time (iter): " << solverRansIt.getTurbModelTime() << "\n";
    gsInfo << "Assembly time (turb iter): " << solverRansIt.getTurbAssemblyTime() << "\n";
    gsInfo << "Solver setup time (turb iter): " << solverRansIt.getTurbSolverSetupTime() << "\n";
    gsInfo << "Solve time (turb iter): " << solverRansIt.getTurbSolveTime() << "\n\n";*/

    real_t errRans = (solverRans.getSolution() - solverRansIt.getSolution()).norm() / solverRans.getSolution().norm();
    gsInfo << "Rans solution error (direct vs. iterative): " << errRans << "\n\n";

    gsInfo << "Avg lin. iter. per Picard iter.: " << solverRansIt.getAvgLinIterations() << "\n\n";


    if (outFile)
    {
        std::ofstream file;
        file.open("INSiterative_test.txt");

        file << "Solving 3D backward facing step with iterative solvers.\n";
        file << "numRefine = " << numRefine << "\n";
        file << "viscosity = " << viscosity << "\n";
        file << "timeStep = " << timeStep << "\n";
        file << "numIter = " << numIter << "\n";
        file << "precType = " << precType << "\n";
        file << "linTol = " << linTol << "\n";
        file << "TM SUPG = " << TMsupg << "\n\n";
        file << "numDofs = " << solverSteady.numDofs() << "\n";

        file << "\nSteady solver:\n";
        file << "Steady (direct) total time: " << timeSt << "\n";
        file << "Steady (iter) total time: " << timeStIt << "\n";
        file << "Assembly time (direct): " << solverSteady.getAssemblyTime() << "\n";
        file << "Assembly time (iter): " << solverSteadyIt.getAssemblyTime() << "\n";
        file << "Solver setup time (direct): " << solverSteady.getSolverSetupTime() << "\n";
        file << "Solver setup time (iter): " << solverSteadyIt.getSolverSetupTime() << "\n";
        file << "Solve time (direct): " << solverSteady.getSolveTime() << "\n";
        file << "Solve time (iter): " << solverSteadyIt.getSolveTime() << "\n\n";
        file << "Last solution change norm (direct): " << solverSteady.solutionChangeRelNorm() << "\n";
        file << "Last solution change norm (iterative): " << solverSteadyIt.solutionChangeRelNorm() << "\n";
        file << "Steady solution error (direct vs. iterative): " << errSt << "\n\n";

        file << "\nUnsteady solver:\n";
        file << "Unsteady (direct) total time: " << timeUnst << "\n";
        file << "Unsteady (iter) total time: " << timeUnstIt << "\n";
        file << "Assembly time (direct): " << solverUnsteady.getAssemblyTime() << "\n";
        file << "Assembly time (iter): " << solverUnsteadyIt.getAssemblyTime() << "\n";
        file << "Solver setup time (direct): " << solverUnsteady.getSolverSetupTime() << "\n";
        file << "Solver setup time (iter): " << solverUnsteadyIt.getSolverSetupTime() << "\n";
        file << "Solve time (direct): " << solverUnsteady.getSolveTime() << "\n";
        file << "Solve time (iter): " << solverUnsteadyIt.getSolveTime() << "\n\n";
        file << "Last solution change norm (direct): " << solverUnsteady.solutionChangeRelNorm() << "\n";
        file << "Last solution change norm (iterative): " << solverUnsteadyIt.solutionChangeRelNorm() << "\n";
        file << "Unsteady solution error (direct vs. iterative): " << errUnst << "\n\n";

        file << "\nRANS solver:\n";
        file << "RANS total time (direct): " << timeRans << "\n";
        file << "RANS total time (iter): " << timeRansIt << "\n";
        file << "Assembly time (direct): " << solverRans.getAssemblyTime() << "\n";
        file << "Assembly time (iter): " << solverRansIt.getAssemblyTime() << "\n";
        file << "Solver setup time (direct): " << solverRans.getSolverSetupTime() << "\n";
        file << "Solver setup time (iter): " << solverRansIt.getSolverSetupTime() << "\n";
        file << "Solve time (direct): " << solverRans.getSolveTime() << "\n";
        file << "Solve time (iter): " << solverRansIt.getSolveTime() << "\n";
        file << "Turb. model time (direct, should be cca equal for iter): " << solverRans.getTurbModelTime() << "\n";
        file << "Assembly time (turb): " << solverRans.getTurbAssemblyTime() << "\n";
        file << "Solver setup time (turb): " << solverRans.getTurbSolverSetupTime() << "\n";
        file << "Solve time (turb): " << solverRans.getTurbSolveTime() << "\n\n";
        file << "Last solution change norm (direct): " << solverRans.solutionChangeRelNorm() << "\n";
        file << "Last solution change norm (iterative): " << solverRansIt.solutionChangeRelNorm() << "\n";
        file << "Rans solution error (direct vs. iterative): " << errRans << "\n\n";

        file.close();
    }

    if (plot)
    {
        gsInfo << "Plotting in Paraview...\n";
        
        // steady 
        gsField<> velocity = solverSteady.constructSolution(0);
        gsField<> pressure = solverSteady.constructSolution(1);

        gsField<> velocityIt = solverSteadyIt.constructSolution(0);
        gsField<> pressureIt = solverSteadyIt.constructSolution(1);

        gsWriteParaview<>(velocity, "step_velocitySteady", plot_pts, true);
        gsWriteParaview<>(pressure, "step_pressureSteady", plot_pts);
        gsWriteParaview<>(velocityIt, "step_velocitySteadyIt", plot_pts, true);
        gsWriteParaview<>(pressureIt, "step_pressureSteadyIt", plot_pts);

        gsField<> velocityUn = solverUnsteady.constructSolution(0);
        gsField<> pressureUn = solverUnsteady.constructSolution(1);

        gsField<> velocityUnIt = solverUnsteadyIt.constructSolution(0);
        gsField<> pressureUnIt = solverUnsteadyIt.constructSolution(1);

        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocityUn, "step_velocityUnsteady", plot_pts, true);
        gsWriteParaview<>(pressureUn, "step_pressureUnsteady", plot_pts);
        gsWriteParaview<>(velocityUnIt, "step_velocityUnsteadyIt", plot_pts, true);
        gsWriteParaview<>(pressureUnIt, "step_pressureUnsteadyIt", plot_pts);

        gsField<> velocityRans = solverRans.constructSolution(0);
        gsField<> pressureRans = solverRans.constructSolution(1);

        gsField<> velocityRansIt = solverRansIt.constructSolution(0);
        gsField<> pressureRansIt = solverRansIt.constructSolution(1);

        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocityRans, "step_velocityRans", plot_pts, true);
        gsWriteParaview<>(pressureRans, "step_pressureRans", plot_pts);
        gsWriteParaview<>(velocityRansIt, "step_velocityRansIt", plot_pts, true);
        gsWriteParaview<>(pressureRansIt, "step_pressureRansIt", plot_pts);
    }

    return 0;
}

template<class T>
void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, const T kIn, const T kWall, const T oIn, const T oWall)
{
    // Boundary conditions
    gsFunctionExpr<T> Kin(util::to_string(kIn), 3);
    gsFunctionExpr<T> Oin(util::to_string(oIn), 3);
    gsFunctionExpr<T> K(util::to_string(kWall), 3);
    gsFunctionExpr<T> O(util::to_string(oWall), 3);

    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Kin, 0);

    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Oin, 1);
}