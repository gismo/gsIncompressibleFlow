/** @file uwbRANSTests.cpp

Author(s): E. Turnerova, K. Michalkova
*/
#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>



// solvers
#include "uwbINSSolverSteady.h"
#include "uwbRANSSolver.h"
#include "uwbRANSSolverIterative.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"

//problem settings
#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"

#include "uwbADRNonconstCoefsProfile.h"

#include <math.h>
#include <string.h>
#include <sstream>

using namespace gismo;

int main(int argc, char *argv[])
{
    const real_t PI = 3.141592653589793238463;


    // ========================================= Settings =========================================
    /*bool plot = true;
    bool plotMeshes = true;
    int plot_pts = 30000;

    real_t timeStep = 0.001;

    int numIterSteadyNS = 50;
    int numIterKOmegaSteady = 50;
    int maxNumIterRANS = 2;
    int numIterRANS = maxNumIterRANS;
    int minNumIterRANS = 1;
    real_t tolRelNorm = 1e-2;
    int maxRANSPicardIt = 2;
    int maxTMPicardFirstIt = 50;

    real_t viscosity = 0.000001;
    real_t viscositySteady = 0.01;
    real_t turbIntensity = 0.05;
    real_t viscosityRatio = 500.;

    unsigned index_of_profile = 6; // i_blade for all profiles; 0-6
    //-----------------------------------------------
    std::string tmEvaluator = "koWilcoxLRN";
    //std::string tmEvaluator = "koSST";
    //-----------------------------------------------

    int refineType = 1;

    bool TMsupg = false;
    int tauSUPGType = 0; // 0 'h/(2*deg*norm(u)) * (cotgh(Pe)) - 1/Pe'; 1 'h/(2*deg*norm(u)'; 2 '((2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
                         // 3 '((2/timeStep)^2 + (2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
    */

    //========================================== load data from file =======================================================
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(9);
    gsVector<real_t> inRealT(6);
    std::string inString;
    gsVector<bool> inBool(3);

    readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettings.txt", inInt, inRealT, inString, inBool);

    // ========================================= Settings =========================================
    bool plot = inBool(0);
    bool plotMeshes = inBool(1);
    int plot_pts = inInt(0);

    real_t timeStep = inRealT(0);

    int numIterSteadyNS = inInt(1);
    int numIterKOmegaSteady = inInt(2);
    int maxNumIterRANS = inInt(3);
    int numIterRANS = inInt(4);
    int minNumIterRANS = inInt(5);
    real_t tolRelNorm = inRealT(1);
    int maxRANSPicardIt = inInt(6);
    int maxTMPicardFirstIt = inInt(7);

    real_t viscosity = inRealT(2);
    real_t viscositySteady = inRealT(3);
    real_t turbIntensity = inRealT(4);
    real_t viscosityRatio = inRealT(5);

    std::string tmEvaluator = inString;

    bool TMsupg = inBool(2);
    int tauSUPGType = inInt(8);

    // ========================================= Geometry =========================================
    gsInfo << "Reading file geometry. \n";

    gsVector<int> inGeoInt(6);
    gsMatrix<int> refineUniformSettings(0,0);
    gsMatrix<int> refineLocalSettings(0,0);
    gsVector<real_t> geomParams(0);
    gsVector<bool> inGeoBool(2);
    std::vector<real_t> kvfit_knots;

    readInitialGeometry(MOTOR_DATA_DIR "uwb-pilsen/initialGeometry.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);

    int geomChoice = inGeoInt(0);
    int index_of_profile = inGeoInt(1);
    int uniformRefine = inGeoInt(2);

    bool uniform_knot = inGeoBool(0);
    bool coarse = inGeoBool(1);

    //gsInfo << refineUniformSettings << "\n";
    //gsInfo << refineLocalSettings << "\n";

    if (uniform_knot)
        kvfit_knots = {0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1};
    else
        kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};

    //=========================================================================================================================
    bool AFC = false;

    real_t omega = 0.;//2*PI*538.0/60.0;   // value 538 cycles per minute is given

    int setting_of_blade_rotation = 0; //0-optimal, 1-maximal, 2-minimal
    int setting_of_guide_vanes = 0; //0-compute with the formula, 1- 56° GV for optimal and 60° GV for maximal, 2-68° GV for optimal a 72° GV for maximal

    bool loadIC = false;

    int plot_pts_atProfile = 173;

    //------------direct solver----------------------
    typedef uwbRANSSolver<real_t> SolverType;

    //------------linear solvers---------------------
    bool iterative = false;
    //typedef uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> > SolverType;

    // preconditioner type options:
    // LSC_AdiagEqual, LSC_Adiag, LSC_Awhole, LSC_Amod
    // AL_Awhole, AL_Amod
    // SIMPLE_AdiagEqual, SIMPLE_Adiag, SIMPLE_Awhole, SIMPLE_Amod
    // SIMPLER_AdiagEqual, SIMPLER_Adiag, SIMPLER_Awhole, SIMPLER_Amod
    // MSIMPLER_AdiagEqual, MSIMPLER_Adiag, MSIMPLER_Awhole, MSIMPLER_Amod
    //std::string precType = "LSC_AdiagEqual";
    std::string precType = "MSIMPLER_AdiagEqual";

    real_t gamma = 2.5; // parameter for AL precond.

    int linMaxIt = 100;
    real_t linTol = 1e-7;

    bool printInfo = true;

    //int numThreads = 10; // number of threads for assembly

    gsInfo << "Solving turbulent flow for 2d blade profile.\n";

    //==================================================================================================
    //==================================================================================================
    //==================================================================================================
    //==================================================================================================
    // ========================================= Define problem =========================================
    gsInfo << "===================================================\n";
    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;
    opt.dirValues = dirichlet::interpolation;
    //opt.dirValues = dirichlet::l2Projection; //default

    std::string profile = "_blade" + std::to_string(index_of_profile);//strs_profile.str();
    gsInfo << "profil" + profile;
    gsInfo << "\n";

    const unsigned num_blades = 4;
    const unsigned num_bladeprofiles = 7;

    uwbRANSProfileExample<real_t> problemSettings(viscosity, profile, index_of_profile, num_bladeprofiles, num_blades, geomChoice, kvfit_knots, coarse, geomParams, plot_pts);
    problemSettings.setParameters(setting_of_blade_rotation, setting_of_guide_vanes);

    //=======================Input velocities==================================================================
    problemSettings.setProfileInputVelocities(omega);
    gsVector<> velocity_absolute_x = problemSettings.getVelocityAbsoluteX();
    gsVector<> velocity_absolute_y = problemSettings.getVelocityAbsoluteY();
    gsVector<> velocity_blade = problemSettings.getVelocityBlade();
    //========================================= Define geometry =========================================
    gsBoundaryConditions<> bcInfo;
    problemSettings.defineBCs_NS(bcInfo, velocity_absolute_x(index_of_profile), velocity_absolute_y(index_of_profile), velocity_blade(index_of_profile));
    //problemSettings.defineBCs_NS(bcInfo, 0.5, 0.0, 0.0);
    gsFunctionExpr<> f("0", "0", 2);


    gsMultiPatch<real_t> patches_old = problemSettings.makeMultiPatch();
    gsFileData<> fdpatches;
    fdpatches << patches_old;
    fdpatches.save("filePatches_created.xml");

     gsFileData<> fdp("filePatches_created.xml");
     gsMultiPatch<> * patches_file = fdp.getAnyFirst< gsMultiPatch<> >().release();

//      gsMultiPatch<> patches = *patches_file;
//              gsInfo << patches << "\n";
//      gsWriteParaview( patches, "patches2", 50000, true);

      gsMultiPatch<> patches = * patches_file;


    // ========================================= Define basis and refine =========================================
    gsMultiBasis<> tbasis(patches); // basis for RANS equations

    for(int i = 0; i < uniformRefine; i++)
        tbasis.uniformRefine();

    for(int i = 0; i < refineUniformSettings.rows(); i++)
        problemSettings.refineBasisUniformZones(tbasis, refineUniformSettings(i,0), refineUniformSettings(i,1), refineUniformSettings(i,2), refineUniformSettings(i,3), refineUniformSettings(i,4));

    for(int i = 0; i < refineLocalSettings.rows(); i++)
        problemSettings.refineBasisLocalZones(tbasis, refineLocalSettings(i,0), refineLocalSettings(i,1), refineLocalSettings(i,2), refineLocalSettings(i,3), refineLocalSettings(i,4));

    //-------------------------------------------------------------------------------------------------
    //==================================================================================================

    std::ofstream file;
    file.open("ProfileOutputInfo.txt");
    if (printInfo)
    {
        file << "============parameters===================" << "\n";
        file << "index_of_profile: " << index_of_profile << "\n";
        file << "timeStep: " << timeStep << "\n";
        file << "viscosity: " << viscosity << "\n";
        file << "viscositySteady: " << viscositySteady << "\n";
        file << "turbIntensity: " << turbIntensity << "\n";
        file << "numIterSteadyNS : " << numIterSteadyNS << "\n";
        file << "numIterKOmegaSteady : " << numIterKOmegaSteady << "\n";
        file << "maxNumIterRANS : " << maxNumIterRANS << "\n";
        file << "minNumIterRANS : " << minNumIterRANS << "\n";
        file << "maxRANSPicardIt : " << maxRANSPicardIt << "\n";
        file << "tolRelNorm : " << tolRelNorm << "\n";
        file << "tmEvaluator = " << tmEvaluator << "\n";
        file << "maxTMPicardFirstIt : " << maxTMPicardFirstIt << "\n";
        file << "iterative : " << iterative << "\n";
        file << "lin. max iter. : " << linMaxIt << "\n";
        file << "lin. tol : " << linTol << "\n";
        file << "geomChoice : " <<geomChoice<< "\n";
        file << "uniformRefine : " <<uniformRefine<< "\n";
        file << "uniform knot : " << uniform_knot << "\n";
        file << "coarse : " <<coarse<< "\n";
        file << "refineBasisUniformZones : " <<refineUniformSettings<< "\n";
        file << "refineBasisLocalZones : " <<refineLocalSettings<< "\n";
        file << "=======================================\n";
    }

    int num_meshes;
    num_meshes = patches.nPatches();

    if (plotMeshes)
    {
        for (int index_of_mesh = 0; index_of_mesh < num_meshes; index_of_mesh++)
        {
            gsMesh<> mesh;
            std::ostringstream strs_mesh;
            strs_mesh << index_of_mesh;
            makeMesh(tbasis.at(index_of_mesh), mesh, 10);
            patches.patch(index_of_mesh).evaluateMesh(mesh);
            gsWriteParaview(mesh, "meshPatch" + strs_mesh.str() + profile);
        }
    }

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n\n";
    file << "Velocity basis degree: " << discreteBases[0].degree() << "\n";

    // ========================================= Compute Steady NS =========================================
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver

    if (!loadIC)
        problemSettings.computeSteadyNS(navStokes, numIterSteadyNS);

    //gsField<> pressureSteadyAtProfile = navStokes.constructSolution(1, problemSettings.getWallBoundaryPatches(), problemSettings.getWallBoundarySides());
    //gsWriteParaview<>(pressureSteadyAtProfile, problemSettings.getGeometryName() + "pressureSteadyAtProfile" + profile, plot_pts_atProfile);

    //================================= wall distance estimation ==================================================
    real_t wallDistance = problemSettings.computeWallDistance(navStokes);

    gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";
    file << "\nminimum wallDistance = " << wallDistance << "\n";

    //==================================== compute aspect ratio  ==================================================
    real_t maxAspectRatio = navStokes.computeAspectRatio();
    gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";
    file << "maxAspectRatio = " << maxAspectRatio << "\n";
    // ========================================= Define turbulence solver =========================================
    gsBoundaryConditions<> bcInfoTurb;

    real_t uFreeStream = problemSettings.getFreeStream();
    real_t kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
    gsInfo << "\nkInConst = " << kInConst << "\n";
    file << "\nkInConst = " << kInConst << "\n";
    real_t oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
    gsInfo << "oInConst = " << oInConst << "\n";
    file << "oInConst = " << oInConst << "\n";

    real_t kBlade = 1.5 * math::pow(velocity_blade(index_of_profile) * turbIntensity, 2);
    real_t oBlade;
    if (omega)
        oBlade = kBlade / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at the wall
    else
    {
        real_t beta = 0.0708;
        oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2));
    }
    gsInfo << "kBlade = " << kBlade << "\n";
    gsInfo << "oBlade = " << oBlade << "\n\n";
    file << "kBlade = " << kBlade << "\n";
    file << "oBlade = " << oBlade << "\n\n";

    problemSettings.defineBCs_TM(bcInfoTurb, kInConst, kBlade, oInConst, oBlade);

    gsDofMapper koMapper;
    discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPdeSteady(patches, bcInfoTurb, f, viscositySteady);
    uwbINSSolverParams<real_t> koParamsSteady(koPdeSteady, discreteBases, opt);

    gsMatrix<> koInitial(koMapper.freeSize(), 2);
    koInitial.col(0).setConstant(kBlade);
    koInitial.col(1).setConstant(oBlade);

    uwbTMSolverKOmegaLinSteady<real_t> turbSolver(koParamsSteady);
    turbSolver.setInitialCondition(koInitial);
    // ========================================= Solving k-omega =========================================

    if (!loadIC)
        problemSettings.computeSteadyKOmegaTM(navStokes, turbSolver, numIterKOmegaSteady, koMapper.freeSize());

    // ========================================= Solving RANS =========================================
    gsInfo << "\nSolving RANS...\n";

    if (iterative)
    {
        gsInfo << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";
        file << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";

        if (precType == "AL_Awhole" || precType == "AL_Amod")
            file << "gamma: " << gamma << "\n\n";
    }
    else
    {
        gsInfo << "SOLVING with direct solver.\n\n";
        file << "SOLVING with direct solver.\n\n";
    }

    uwbINSPde<real_t> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> RANSparams(RANSpde, discreteBases, opt);
    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);
    //RANSparams.settings().set(constantsINS::unst_innerTol, picardTol);

    if (AFC)
        discreteBases[1].getMapper(dirichlet::none, opt.intStrategy, bcInfoTurb, koMapper, 0);
    else
        discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, discreteBases, opt);
    koParams.settings().set(constantsINS::timeStep, timeStep);
    koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
    koParams.settings().setTurbulenceEvaluator(tmEvaluator);

    if (TMsupg) {
        gsWarn << "SUPG stabilization only for turbulent model not for RANS equations.\n";
        koParams.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
        koParams.settings().set(constantsINS::tauStabType, tauSUPGType); //set formula for stabilization parameter tau
    }
    if (AFC)
        koParams.settings().set(constantsINS::TMafc, AFC);

    gsMatrix<real_t> RANSSolution, TMSolution;
    //------------- read from file ---------------------
    if (loadIC)
    {
        gsFileData<> fdRead("RANS_IC_solution.xml");
        RANSSolution = *(fdRead.getFirst< gsMatrix<real_t> >());

        gsFileData<> fdReadTM("TM_IC_solution.xml");
        TMSolution = *(fdReadTM.getFirst< gsMatrix<real_t> >());
    }
    //------------------------------------------------

    uwbTMSolverKOmega<real_t> turbSolver_unsteady(koParams);
    if (loadIC)
        turbSolver_unsteady.setInitialCondition(TMSolution);
    else
        turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());//(koInitial);

    if (tmEvaluator == "koSST")
        problemSettings.solvePoissonEquation(turbSolver_unsteady, geomChoice);

    if (iterative)
    {
        RANSparams.settings().setPrecondType(precType);
        RANSparams.getPrecOptions().setReal("gamma", gamma);
        RANSparams.settings().set(constantsINS::iter_maxIt, linMaxIt);
        RANSparams.settings().set(constantsINS::iter_tol, linTol);
    }
    RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);

    //--------
    //uwbRANSSolver<real_t> ransSolver(RANSparams);
    //uwbRANSSolverIterative<real_t, uwbGMResRight<real_t> > ransSolver(RANSparams);
    SolverType ransSolver(RANSparams);
    //--------


    gsInfo << "initialization...\n";

    if (loadIC)
        ransSolver.setInitialCondition(RANSSolution);
    else
        ransSolver.setInitialCondition(navStokes.getSolution());
    ransSolver.initialize();


    gsInfo << "numDofs: " << ransSolver.numDofs() << "\n";
    file << "numDofs: " << ransSolver.numDofs() << "\n";



    real_t pitch = ((2 * PI*problemSettings.getRr(index_of_profile)) / num_blades);
    gsInfo << "pitch = " << pitch << "\n";
    file << "pitch = " << pitch << "\n";

    //=== solving with solution change relative norm as stopping criterion ============
    ransSolver.solve(numIterRANS, tolRelNorm); // solution change norm tol = 10^(-5)
    //ransSolver.solveWithAnimation(numIterRANS, 2, tolRelNorm, plot_pts, true);

    real_t Tassembly = ransSolver.getAssemblyTime();
    real_t Tsolve = ransSolver.getSolverSetupTime();
    real_t Tsetupsolve = ransSolver.getSolveTime();
    real_t Tturbmodel = ransSolver.getTurbModelTime();


    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
    gsInfo << "Turbulent model time:" << Tturbmodel << "\n";
    file << "Assembly time:" << Tassembly << "\n";
    file << "Solve time:" << Tsolve << "\n";
    file << "Solver setup time:" << Tsetupsolve << "\n";
    file << "Turbulent model time:" << Tturbmodel << "\n";
    file.close();

    //--- save solution into file ---
    gsFileData<> fd, fd_TM;
    fd << ransSolver.getSolution();
    fd_TM << turbSolver_unsteady.getSolution();
    fd.save(problemSettings.getGeometryName() + "RANS_solution"  + profile + ".xml");
    fd_TM.save(problemSettings.getGeometryName() + "TM_solution"  + profile + ".xml");
    //=====================================================================================================

    gsInfo << "Plotting in Paraview...\n";

    gsField<> velocity = ransSolver.constructSolution(0);
    gsField<> pressure = ransSolver.constructSolution(1);
    gsField<> kOmega = turbSolver_unsteady.constructSolution();

    problemSettings.plotSolutionField(velocity, "velocity");
    problemSettings.plotSolutionField(pressure, "pressure");
    problemSettings.plotSolutionField(kOmega, "komega");
    ransSolver.plotTurbulentViscosity(problemSettings.getGeometryName() + "turb_viscosity" + profile, plot_pts);

    gsField<> pressureAtProfile = ransSolver.constructSolution(1, problemSettings.getWallBoundaryPatches(), problemSettings.getWallBoundarySides());
    gsWriteParaview<>(pressureAtProfile, problemSettings.getGeometryName() + "pressureAtProfile" + profile, plot_pts_atProfile);

    gsField<> kDiffusionCoeff = ransSolver.constructTMCoefficientSol("kDiffusionCoeff");
    problemSettings.plotSolutionField(kDiffusionCoeff, "kDiffusionCoeff");
    gsField<> oDiffusionCoeff = ransSolver.constructTMCoefficientSol("oDiffusionCoeff");
    problemSettings.plotSolutionField(oDiffusionCoeff, "oDiffusionCoeff");

    std::vector<gsField<> > advectionDiffusionCoeffs;
    advectionDiffusionCoeffs.push_back(velocity);
    advectionDiffusionCoeffs.push_back(kDiffusionCoeff);

    solveADRprofile(patches, advectionDiffusionCoeffs, problemSettings);

    return 0;
}
