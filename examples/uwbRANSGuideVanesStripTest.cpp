/** @file uwbbl ades.cpp

    Author(s): B. Bastl, K. Michalkova
*/

#include <iostream>
#include <gismo.h>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <fstream>

//#include "../jku/gsMotorUtils.h"
#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"
#include "uwbDraftTube.h"
#include "uwbProfileOptimization.h"
#include "uwbHydraulicProfile.h"
#include "uwbKaplanTurbineRunnerBlade.h"
#include "uwbKaplanTurbineRunnerWheelDomain.h"
#include "uwbKaplanTurbineGuideVane.h"

// solvers
#include "uwbINSSolverSteady.h"
#include "uwbRANSSolver.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"

#include "uwbINSSolverSteadyIterative.h"
#include "uwbINSSolverUnsteadyIterative.h"
#include "uwbRANSSolverIterative.h"

#include "uwbLinSolvers.h"

using namespace gismo;

//template<class TT>
//gsMultiPatch<TT>  BSplineTurbine3D();
std::string omegaXr_component(real_t omega, std::string variable);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, gsFunctionExpr<> & Uin, real_t omega = 0);
template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, real_t kIn, real_t kWall, real_t kBlade, real_t oIn, real_t oWall, real_t oBlade);
template<class T> void refineBasis(gsMultiBasis<T>& tbasis, int numRefine, int numRefineLocal);

int main(int argc, char *argv[])
{
    // ========================================= Settings =========================================

    bool plot = false;
    bool outFile = false;
    int plot_pts = 30000;
    const int numThreads = 2;

    real_t timeStep = 0.01;

    int numIterSteadyNS = 20;
    int numIterKOmegaSteady = 50;
    int numIterRANS = 10;

    real_t omega = 0; // angular velocity for rotation

    int numRefine = 0;
    int numRefineLocal = 3;

    real_t viscosity = 0.001;
    real_t viscositySteady = 0.5;
    real_t turbIntensity = 0.01;

    // preconditioner type options:
    // LSC_AdiagEqual, LSC_Adiag, LSC_Awhole, LSC_Amod
    // AL_Awhole, AL_Amod
    // SIMPLE_AdiagEqual, SIMPLE_Adiag, SIMPLE_Awhole, SIMPLE_Amod
    // SIMPLER_AdiagEqual, SIMPLER_Adiag, SIMPLER_Awhole, SIMPLER_Amod
    // MSIMPLER_AdiagEqual, MSIMPLER_Adiag, MSIMPLER_Awhole, MSIMPLER_Amod
    std::string precType = "precond::LSC_Amod";
    bool precIter = true; // set direct/iterative method to solve subsystems in preconditioner
    int linMaxIt = 150;
    real_t gammaSt = 0.1; // parameter for AL precond, steady case
    real_t gammaUnst = 2.2; // parameter for AL precond, unsteady case
    real_t linTol = 1e-6;

    bool TMsupg = false;
    int tauSUPGType = 1;

    if (numRefine<0)
    {
        gsInfo << "Number of refinements must be non-negative, quitting.\n";
        return -1;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;

    gsInfo << "Solving turbulent flow in the Kaplan Turbine with guide vanes, but without runner blades.\n";

    // ========================================= Define problem =========================================

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> Uin(" 3*(-(100 / 9) * (sqrt(y^2+z^2)-0.7)^2 + 1)", "0", "0", 3);
    //gsFunctionExpr<> Uwall("0", "0", "0", 3);

    defineBCs_NS(bcInfo, Uin, omega);

    gsFunctionExpr<> f("0", "0", "0", 3);

    // ========================================= Define geometry =========================================

    std::string input(MOTOR_DATA_DIR "/uwb-pilsen/pruh_withoutTrailingEdge_Update.xml");

    gsFileData<> fileData(input);

    gsMultiPatch<> patches;
    if (fileData.has< gsMultiPatch<> >())
    {
        patches = *(fileData.getFirst< gsMultiPatch<> >());
    }
    else
    {
        gsWarn << "Input file doesn't have a gsMultiPatch inside.\n";
        return -1;
    }

    gsMultiPatch<> TMpatches = patches;

    TMpatches.addInterface(0, boundary::west, (size_t) 0, boundary::east);
    TMpatches.addInterface(2, boundary::west, 2, boundary::east);
    TMpatches.addAutoBoundaries();

    gsInfo << patches << "\n";
    gsInfo << TMpatches << "\n";

    // ========================================= Define basis =========================================
    gsMultiBasis<> tbasis(patches); // basis for RANS equations

    //-------- in case of refinement, 'refineBasis' function used only once - not finished!!!----------
    //refineBasis(tbasis, numRefine, numRefineLocal); // basis for TM model
    //gsMultiBasis<> TMtbasis = tbasis;
    //TMtbasis.addInterface(&patches.basis(0), boundary::front, &patches.basis((0)), boundary::back);
    //TMtbasis.addInterface(&patches.basis(2), boundary::front, &patches.basis((2)), boundary::back);
    //-------------------------------------------------------------------------------------------------

    //-------- in case of refinement, the same 'refineBasis' function must be used twice --------------
    gsMultiBasis<> TMtbasis(TMpatches); // basis for TM model
    refineBasis(tbasis, numRefine, numRefineLocal);
    refineBasis(TMtbasis, numRefine, numRefineLocal);
    //-------------------------------------------------------------------------------------------------

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    std::vector< gsMultiBasis<> >  TMdiscreteBases;
    TMdiscreteBases.push_back(TMtbasis);
    TMdiscreteBases.push_back(TMtbasis);
    TMdiscreteBases[0].degreeElevate(1);

    // ========================================= Compute Steady NS =========================================
    gsInfo << "Solving Steady case: \n";
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt, numThreads);

    if (omega)
        params.settings().set(constantsINS::omega, omega); // set rotation (assumed around x-axis)

    params.settings().setPrecondType(precType);
    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::iter_maxIt, linMaxIt);
    params.settings().set(constantsINS::iter_tol, linTol);
    params.getPrecOptions().setSwitch("iter", precIter);
    params.getPrecOptions().setReal("gamma", gammaSt);

    // solvers
    //uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    uwbINSSolverSteadyIterative<real_t, uwbGMResRight<real_t> > navStokes(params); // steady coupled solver

    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    navStokes.initialize(); // steady solver
    navStokes.solve(numIterSteadyNS, 1e-5);
    //gsMatrix<> steadySolution = navStokes.getSolution();

    gsField<> absVelocitySteady = navStokes.constructSolution(0);
    gsWriteParaview<>(absVelocitySteady, "guidevanesstrip_absSteadyVelocity", plot_pts, true);

    gsFileData<> fd;
    fd << navStokes.getSolution();
    fd.save("guidevanesstrip_steadySolution_nu" + util::to_string(viscositySteady) + ".xml");
    //---------------------

    /*gsInfo << "Computing initial velocity by solving Stokes problem...\n";

    uwbINSPde<real_t> NSpde(*patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    uwbINSSolverSteady<real_t> stokesSolver(params);

    gsInfo << "numDofs: " << stokesSolver.numDofs() << "\n";

    stokesSolver.initialize();

    stokesSolver.setStokesSolution();*/

    //================================= wall distance estimation ==================================================
    real_t inletWidth = 0.6170365; //
    //gsVector<real_t> Uin(3);
    real_t Umax_inlet = 3;
    real_t uFreeStream = Umax_inlet;
    real_t Re = uFreeStream * inletWidth / viscosity;

    //vector of the sides of the patches from which the wall distance is computed
    std::vector<boxSide> distanceSidesWalls;
    distanceSidesWalls.push_back(boundary::north);
    distanceSidesWalls.push_back(boundary::south);
    distanceSidesWalls.push_back(boundary::north);
    distanceSidesWalls.push_back(boundary::south);
    distanceSidesWalls.push_back(boundary::north);
    distanceSidesWalls.push_back(boundary::south);
    std::vector<boxSide> distanceSidesBlade;
    distanceSidesBlade.push_back(boundary::east);
    distanceSidesBlade.push_back(boundary::west);

    //vector of indexes of the patches corresponding to distanceSides
    //length of the vector distancePatches must be equal to the length of vector distanceSides
    gsVector<int> distancePatchesWalls(6);
    distancePatchesWalls << 0, 0, 1, 1, 2, 2;
    gsVector<int> distancePatchesBlade(2);
    distancePatchesBlade << 1, 1;

    int numSamplePts = 500; //number of sample points for which the distance to the boundary is computed
    real_t maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

    //table with wall distance information for every chosen side is printed, if the last input parameter is set as true
    //real_t wallDistanceWalls = navStokes.getAssembler()->getBlockAssembler().estimateDimensionlessWallDistance(distancePatchesWalls, distanceSidesWalls, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true);
    real_t wallDistanceWalls = navStokes.computeDimensionlessWallDistance(distancePatchesWalls, distanceSidesWalls, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
    gsInfo << "\nminimum wallDistance at walls = " << wallDistanceWalls << "\n";
    //real_t wallDistanceBlade = navStokes.getAssembler()->getBlockAssembler().estimateDimensionlessWallDistance(distancePatchesBlade, distanceSidesBlade, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true);
    real_t wallDistanceBlade = navStokes.computeDimensionlessWallDistance(distancePatchesBlade, distanceSidesBlade, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
    gsInfo << "\nminimum wallDistance at blade = " << wallDistanceBlade << "\n";

    // ========================================= Define turbulence solver =========================================
    gsInfo << "Define turbulence solver: \n";

    gsBoundaryConditions<> bcInfoTurb;

    real_t kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
    gsInfo << "\nkInConst = " << kInConst << "\n";
    real_t viscosityRatio = 10.;
    real_t oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
    gsInfo << "oInConst = " << oInConst << "\n";

    real_t kWall = 0.0;//1e-10;
    real_t beta = 0.0708;
    real_t oWall = 6 * viscosity / (beta * math::pow(wallDistanceWalls, 2));
    gsInfo << "kWall = " << kWall << "\n";
    gsInfo << "oWall = " << oWall << "\n";
    real_t kBlade = 1.5 * math::pow(omega * turbIntensity, 2);
    real_t oBlade;
    if (omega)
        oBlade = kBlade / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at the wall
    else
        oBlade = 6 * viscosity / (beta * math::pow(wallDistanceBlade, 2));

    gsInfo << "kBlade = " << kBlade << "\n";
    gsInfo << "oBlade = " << oBlade << "\n\n";
    defineBCs_TM(bcInfoTurb, kInConst, kWall, kBlade, oInConst, oWall, oBlade);

    gsDofMapper koMapper;
    TMdiscreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPdeSteady(TMpatches, bcInfoTurb, f, viscositySteady);
    uwbINSSolverParams<real_t> koParamsSteady(koPdeSteady, TMdiscreteBases, opt, numThreads);

    gsMatrix<> koInitial(koMapper.freeSize(), 2);
    koInitial.col(0).setConstant(kWall);
    koInitial.col(1).setConstant(oWall);

    uwbTMSolverKOmegaLinSteady<real_t> turbSolver(koParamsSteady);
    turbSolver.setInitialCondition(koInitial);

    // ========================================= Solving k-omega =========================================
    gsInfo << "\nComputing steady linearized k-omega...\n";

    gsInfo << "initialization...\n";

    gsField<> velocity = navStokes.constructSolution(0);

    //gsStopwatch time;
    turbSolver.initialize(velocity);
    //real_t Tassembly = time.stop();
    //gsInfo << "Assembly time:" << Tassembly << "\n";

    //time.restart();
    turbSolver.solve(numIterKOmegaSteady, 1e-5); // solution change norm tol = 10^(-5)
    //real_t Tsolve = time.stop();

    //gsInfo << "Solve time:" << Tsolve << "\n";

    gsField<> kOmegaSteady = turbSolver.constructSolution();
    gsWriteParaview<>(kOmegaSteady, "guidevanesstrip_kOmegaSteady", plot_pts, true);

    gsFileData<> fdTM;
    fdTM << turbSolver.getSolution();
    fdTM.save("guidevanesstrip_TMsteadySolution_nu" + util::to_string(viscositySteady) + ".xml");

    // ========================================= Solving RANS =========================================
    gsInfo << "\nSolving RANS...\n";

    uwbINSPde<real_t> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> RANSparams(RANSpde, discreteBases, opt, numThreads);

    uwbINSPde<real_t> koPde(TMpatches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, TMdiscreteBases, opt, numThreads);

    if (TMsupg) {
        koParams.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType); //set formula for stabilization parameter tau
    }

    uwbTMSolverKOmega<real_t> turbSolver_unsteady(koParams);
    turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());

    if (omega)
        RANSparams.settings().set(constantsINS::omega, omega); // set rotation (assumed around x-axis)

    RANSparams.settings().setPrecondType(precType);
    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().set(constantsINS::iter_maxIt, linMaxIt);
    RANSparams.settings().set(constantsINS::iter_tol, linTol);
    RANSparams.getPrecOptions().setSwitch("iter", precIter);
    RANSparams.getPrecOptions().setReal("gamma", gammaUnst);

    RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);
    //uwbRANSSolver<real_t> ransSolver(RANSparams);
    uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> > ransSolver(RANSparams);

    gsInfo << "initialization...\n";

    ransSolver.setInitialCondition(navStokes.getSolution());
    ransSolver.initialize();



    ransSolver.solve(numIterRANS, 1e-5); // solution change norm tol = 10^(-5)

    real_t Tassembly = ransSolver.getAssemblyTime();
    real_t Tsolve = ransSolver.getSolverSetupTime();
    real_t Tsetupsolve = ransSolver.getSolveTime();
    real_t Tturbmodel = ransSolver.getTurbModelTime();


    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
    gsInfo << "Turbulent model time:" << Tturbmodel << "\n";


    if (outFile)
    {
        time_t t = std::time(0);   // get time now
        struct tm * now = std::localtime(&t);
        std::stringstream filename;
        filename << "guidevanesstrip_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
            << now->tm_hour << "-" << now->tm_min << ".txt";

        std::ofstream ofile;
        ofile.open(filename.str());
        ofile << "Solving RANS problem in Kaplan turbine with guide vanes, but without runner blades:\n";
        ofile << "viscosity = " << viscosity << "\n";
        ofile << "angular velocity = " << omega << "\n";
        ofile << "time step = " << timeStep << "\n";
        ofile << "number of iterations = " << ransSolver.getIterationNumber() << "\n";
        ofile << "last relative solution change = " << ransSolver.solutionChangeRelNorm() << "\n";
        ofile << "Assembly time:" << Tassembly << "\n";
        ofile << "Solve time:" << Tsolve << "\n";
        ofile << "Solver setup time:" << Tsetupsolve << "\n";
        ofile << "Turbulent model time:" << Tturbmodel << "\n";
        ofile.close();
    }

    // Optionally plot solution in paraview
    if (plot)
    {
        gsInfo << "Plotting in Paraview...\n";

        //gsVector<size_t> relPatches(1);
        //relPatches << 1;

        gsField<> absVelocity = ransSolver.constructSolution(0);
        gsField<> relVelocity = ransSolver.constructSolution(0, true);
        //gsField<> combVelocity = ransSolver.constructSolutionCombined(relPatches);
        gsField<> pressure = ransSolver.constructSolution(1);
        gsField<> kOmega = turbSolver_unsteady.constructSolution();


        // Write solution to paraview files
        gsWriteParaview<>(absVelocity, "guidevanesstrip_RANS_absVelocity", plot_pts, true);
        gsWriteParaview<>(relVelocity, "guidevanesstrip_RANS_relVelocity", plot_pts);
        //gsWriteParaview<>(combVelocity, "turbine_RANS_combVelocity", plot_pts);
        gsWriteParaview<>(pressure, "guidevanesstrip_RANS_pressure", plot_pts);;
        gsWriteParaview<>(kOmega, "guidevanesstrip_TM_komega", plot_pts);
        ransSolver.plotTurbulentViscosity("guidevanesstrip_turb_viscosity", plot_pts);

        // Run paraview
        //system("paraview guidevanesstrip_RANS_absVelocity.pvd");
    }

    //getchar();

    return 0;
}

std::string omegaXr_component(real_t omega, std::string variable)
{
    std::ostringstream s;

    s << omega << " * " << variable;

    return s.str();
}

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, gsFunctionExpr<> & Uin, real_t omega)
{
    gsMatrix<real_t> transformMatrix(3, 3);

    const real_t phi = -(1. / 7)*EIGEN_PI;
    const real_t cos = math::cos(phi);
    const real_t sin = math::sin(phi);
    transformMatrix(0, 0) = 1;
    transformMatrix(0, 1) = 0;
    transformMatrix(0, 2) = 0;
    transformMatrix(1, 0) = 0;
    transformMatrix(1, 1) = cos;
    transformMatrix(1, 2) = sin;
    transformMatrix(2, 0) = 0;
    transformMatrix(2, 1) = -sin;
    transformMatrix(2, 2) = cos;

    gsFunctionExpr<> Uwall("0", "0", "0", 3); // wall velocity
    //gsFunctionExpr<> Ublade("0", omegaXr_component(-omega, "z"), omegaXr_component(omega, "y"), 3);

    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);

    bcInfo.setIdentityMatrix(3);
    bcInfo.addPeriodic(0, boundary::west, 0, boundary::east, 3);
    bcInfo.addPeriodic(2, boundary::west, 2, boundary::east, 3);

    bcInfo.setTransformMatrix(transformMatrix);
}

template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfoTurb, real_t kIn, real_t kWall, real_t kBlade, real_t oIn, real_t oWall, real_t oBlade)
{
    // Boundary conditions
    gsFunctionExpr<> Kin(util::to_string(kIn), 3);
    gsFunctionExpr<> Oin(util::to_string(oIn), 3);
    gsFunctionExpr<> Kwall(util::to_string(kWall), 3);
    gsFunctionExpr<> Kblade(util::to_string(kBlade), 3);
    gsFunctionExpr<> Owall(util::to_string(oWall), 3);
    gsFunctionExpr<> Oblade(util::to_string(oBlade), 3);

    bcInfoTurb.addCondition(0, boundary::front, condition_type::dirichlet, Kin, 0);
    bcInfoTurb.addCondition(0, boundary::north, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, Kwall, 0);

    bcInfoTurb.addCondition(2, boundary::south, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(2, boundary::north, condition_type::dirichlet, Kwall, 0);

    bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Kblade, 0);
    bcInfoTurb.addCondition(1, boundary::east, condition_type::dirichlet, Kblade, 0);
    //---

    bcInfoTurb.addCondition(0, boundary::front, condition_type::dirichlet, Oin, 1);
    bcInfoTurb.addCondition(0, boundary::north, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, Owall, 1);

    bcInfoTurb.addCondition(2, boundary::south, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(2, boundary::north, condition_type::dirichlet, Owall, 1);

    bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Oblade, 1);
    bcInfoTurb.addCondition(1, boundary::east, condition_type::dirichlet, Oblade, 1);
}

template<class T> void refineBasis(gsMultiBasis<T>& tbasis, int numRefine, int numRefineLocal)
{
    //
    // Implemented for "pruh_withoutTrailingEdge.xml
    //
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    real_t uRefineKnot;
    real_t vRefineKnot;
    //real_t wallRefineKnot = 0.1;
    gsMatrix<> box_u0(3, 2);
    gsMatrix<> box_v0(3, 2);
    gsMatrix<> box_w0(3, 2);

    // Initial refinemenet in u
    for (int i = 0; i < 1; i++)
    {
        int k = static_cast<int>(pow(2, i+1) + 1);
        for (int j = 0; j < k; j++)
        {
            box_u0 << j/(k-1), (j+1)/(k-1), 0, 0, 0, 0;
            tbasis.refine(0, box_u0);
            tbasis.refine(1, box_u0);
            tbasis.refine(2, box_u0);
        }
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";
    // Initial refinemenet in v
    for (int i = 0; i < 2; i++)
    {
        int k = static_cast<int>(pow(2, i) + 1);
        for (int j = 0; j < k; j++)
        {
            box_u0 << 0, 0, j/(k-1), (j+1)/(k-1), 0, 0;
            tbasis.refine(0, box_u0);
            tbasis.refine(1, box_u0);
            tbasis.refine(2, box_u0);
        }
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";
    // Initial refinemenet in w
    const gsTensorBSplineBasis<3, real_t>*  basis_p0 = dynamic_cast<const gsTensorBSplineBasis<3, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
    box_w0 << 0, 0, 0, 0, basis_p0->knot(2, (basis_p0->degree(2))+2), basis_p0->knot(2, (basis_p0->degree(2))+2+(basis_p0->degree(2)));
    tbasis.refine(0, box_w0);
    const gsTensorBSplineBasis<3, real_t>*  basis_p1 = dynamic_cast<const gsTensorBSplineBasis<3, real_t>*>(&(tbasis.basis(1))); //basis of the second patch
    box_w0 << 0, 0, 0, 0, 0, basis_p1->knot(2, (basis_p1->degree(2))+1);
    tbasis.refine(1, box_w0);
    box_w0 << 0, 0, 0, 0, 0.8, 1;
    tbasis.refine(1, box_w0);
    //const gsTensorBSplineBasis<3, real_t>*  basis_p2 = dynamic_cast<const gsTensorBSplineBasis<3, real_t>*>(&(tbasis.basis(2))); //basis of the third patch
    /*for (int i = 0; i < 2; i++)
    {
        int k = static_cast<int>(pow(2, i) + 1);
        for (int j = 0; j < k; j++)
        {
            box_u0 << 0, 0, 0, 0, , ;
            tbasis.refine(2, box_w0);
        }
    }*/
    for (int i = 0; i < 2; i++)
    {
        box_w0 << 0, 0, 0, 0, 0.25, 1;
        tbasis.refine(2, box_w0);
    }

    // Uniform refinement the whole multi-patch
    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();

    //refine in u, near bottom wall
    const gsTensorBSplineBasis<3, real_t>*  basis_u0 = dynamic_cast<const gsTensorBSplineBasis<3, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
    uRefineKnot = basis_u0->knot(0, (basis_u0->degree(0))+1);
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_u0 << 0, uRefineKnot / math::pow(2, i), 0, 0, 0, 0;
        gsInfo << box_u0 << "\n\n";

        tbasis.refine(0, box_u0);
        tbasis.refine(1, box_u0);
        tbasis.refine(2, box_u0);

        box_u0 << 1 - (uRefineKnot / math::pow(2, i)), 1, 0, 0, 0, 0;
        gsInfo << box_u0 << "\n\n";

        tbasis.refine(0, box_u0);
        tbasis.refine(1, box_u0);
        tbasis.refine(2, box_u0);
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    //refine in v, near blade
    vRefineKnot = basis_u0->knot(1, (basis_u0->degree(1))+1);
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_v0 << 0, 0, 0, vRefineKnot / math::pow(2, i), 0, 0;
        gsInfo << box_v0 << "\n\n";

        tbasis.refine(0, box_v0);
        tbasis.refine(1, box_v0);
        tbasis.refine(2, box_v0);

        box_v0 << 0, 0, 1 - (vRefineKnot / math::pow(2, i)), 1, 0, 0;
        gsInfo << box_v0 << "\n\n";

        tbasis.refine(0, box_v0);
        tbasis.refine(1, box_v0);
        tbasis.refine(2, box_v0);
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    gsInfo << "Knot vectors (0): \n" << tbasis.piece(0).component(0).detail() << tbasis.piece(0).component(1).detail() << tbasis.piece(0).component(2).detail() << "\n";
    gsInfo << "Knot vectors (1): \n" << tbasis.piece(1).component(0).detail() << tbasis.piece(1).component(1).detail() << tbasis.piece(1).component(2).detail() << "\n";
    gsInfo << "Knot vectors (2): \n" << tbasis.piece(2).component(0).detail() << tbasis.piece(2).component(1).detail() << tbasis.piece(2).component(2).detail() << "\n";


    /*
    //refine in u after blade
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_u1 << 0, uRefineKnot / math::pow(2, i + 1), 0, 0;

        tbasis.refine(4, box_u1);
        tbasis.refine(5, box_u1);
    }

    ////refine in v
    //for (unsigned k = 0; k < patches.nPatches(); k++) {
    //    for (int i = 0; i < numRefineLocal; i++)
    //    {
    //        box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
    //        box_v1 << 0, 0, 1 - (wallRefineKnot / math::pow(2, i)), 1;
    //        tbasis.refine(k, box_v0);
    //        tbasis.refine(k, box_v1);
    //    }
    //}

    //wall part
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
        box_v1 << 0, 0, 1 - (wallRefineKnot / math::pow(2, i)), 1;
        tbasis.refine(0, box_v1);
        tbasis.refine(2, box_v1);
        tbasis.refine(4, box_v1);
        tbasis.refine(1, box_v0);
        tbasis.refine(3, box_v0);
        tbasis.refine(5, box_v0);
    }

    //blade part
    wallRefineKnot = wallRefineKnot + 0.1;
    for (int i = 0; i < numRefineLocal + 1; i++)
    {
        box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
        box_v1 << 0, 0, 1 - (wallRefineKnot / math::pow(2, i)), 1;
        tbasis.refine(0, box_v0);
        tbasis.refine(2, box_v0);
        tbasis.refine(4, box_v0);
        tbasis.refine(1, box_v1);
        tbasis.refine(3, box_v1);
        tbasis.refine(5, box_v1);
    }
    */
}

