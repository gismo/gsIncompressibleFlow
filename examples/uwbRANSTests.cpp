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
//#include "uwbRANSSolverIterative.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"

//problem settings
#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"


#include <math.h>
#include <string.h>
#include <sstream>

using namespace gismo;

int main(int argc, char *argv[])
{
    const real_t PI = 3.141592653589793238463;


    // ========================================= Settings =========================================
    /*int plot_pts = 10000;

    int numIterSteadyNS = 10;
    int numIterKOmegaSteady = 20;
    int numIterRANS = 10;
    int minNumIterRANS = 1;
    int maxRANSPicardIt = 3;
    int maxTMPicardFirstIt = 10;
    int maxTMPicardIt = 5;

    int tauStabTypeSUPG = 2; // 0 'h/(2*deg*norm(u)) * (cotgh(Pe)) - 1/Pe'; 1 'h/(2*deg*norm(u)'; 2 '((2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
                             // 3 '((2/timeStep)^2 + (2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
    int tauStabTypeCW = 2;
    int tauStabTypeAD = 2;

    int crosswindType = 0;

    int tauTMStabTypeSUPG = 3;
    int tauTMStabTypeCW = 3;
    int tauTMStabTypeAD = 3;

    int TMcrosswindType = 0;

    int diamPowerSUPG = 2;
    int diamPowerCW = 2;
    int diamPowerAD = 2;
    int TMdiamPowerSUPG = 2;
    int TMdiamPowerCW = 2;
    int TMdiamPowerAD = 2;

    int CWresidualType = 0;
    int TMCWresidualType = 0;

    int animateStep = 5;

    real_t timeStep = 0.001;
    real_t tolRelNorm = 1e-10;

    real_t viscosity = 0.0001;
    real_t viscositySteady = 0.05;

    real_t turbIntensity = 0.03;
    real_t viscosityRatio = 10.;

    real_t productionXPoint = -0.16;

    std::string tmEvaluator = "koWilcoxLRN";
    //std::string tmEvaluator = "koSST";

    bool plot = false;
    bool plotMeshes = false;

    bool SUPG = false;
    bool INScrosswind = false;
    bool RANScrosswind = false;
    bool RANSad = false;

    bool TMsupg = false;
    bool TMcrosswind = false;
    bool TMad = false;
    bool TMisoAD = false;
    bool TMfctLowOrder = false;
    bool TMafc = false;
    bool TMafcHO = false;

    bool animate = false;

    bool limitTMProduction = false;

    bool loadIC = false;
    bool NSsteadyInitial = false;
    */

    //========================================== load data from file =======================================================
    gsInfo << "Reading file settings. \n";
    gsVector<int> inInt(27);
    gsVector<real_t> inRealT(7);
    std::string inString;
    gsVector<bool> inBool(22);

    readInitialSetting(MOTOR_DATA_DIR "uwb-pilsen/initialSettings.txt", inInt, inRealT, inString, inBool);
    //readInitialSetting("initialSettings.txt", inInt, inRealT, inString, inBool);

    // ========================================= Settings =========================================
    int plot_pts = inInt(0);

    int numIterSteadyNS = inInt(1);
    int numIterKOmegaSteady = inInt(2);
    int numIterRANS = inInt(3);
    int minNumIterRANS = inInt(4);
    int maxRANSPicardIt = inInt(5);
    int maxTMPicardFirstIt = inInt(6);
    int maxTMPicardIt = inInt(7);

    int tauStabTypeSUPG = inInt(8);
    int tauStabTypeCW = inInt(9);
    int tauStabTypeAD = inInt(10);

    int crosswindType = inInt(11);

    int tauTMStabTypeSUPG = inInt(12);
    int tauTMStabTypeCW = inInt(13);
    int tauTMStabTypeAD = inInt(14);

    int TMcrosswindType = inInt(15);

    int diamPowerSUPG = inInt(16);
    int diamPowerCW = inInt(17);
    int diamPowerAD = inInt(18);
    int TMdiamPowerSUPG = inInt(19);
    int TMdiamPowerCW = inInt(20);
    int TMdiamPowerAD = inInt(21);

    int CWresidualType = inInt(22);
    int TMCWresidualType = inInt(23);

    int isoADtype = inInt(24);

    int animateStep = inInt(25);

    int plot_pts_atProfile = inInt(26);

    real_t timeStep = inRealT(0);
    real_t tolRelNorm = inRealT(1);

    real_t viscosity = inRealT(2);
    real_t viscositySteady = inRealT(3);

    real_t turbIntensity = inRealT(4);
    real_t viscosityRatio = inRealT(5);

    real_t productionXPoint = inRealT(6);

    std::string tmEvaluator = inString;

    bool plot = inBool(0);
    bool plotMeshes = inBool(1);

    bool SUPG = inBool(2);
    bool TCSD = inBool(3);
    //bool INScrosswind = inBool(4);
    bool RANScrosswind = inBool(5);
    bool RANSad = inBool(6);

    bool TMsupg = inBool(7);
    bool TMcrosswind = inBool(8);
    bool TMad = inBool(9);
    bool TMisoAD = inBool(10);
    bool TMfctLowOrder = inBool(11);
    bool TMafc = inBool(12);
    bool TMafcHO = inBool(13);

    bool animate = inBool(14);

    bool limitTMProduction = inBool(15);

    bool loadIC = inBool(16);
    bool NSsteadyInitial = inBool(17);

    bool UyInZero = inBool(18);
    bool computeTMfirst = inBool(19);
    bool finer_grid = inBool(20);
    bool reSaveSolution = inBool(21);

    //-----------
    bool printInfo = true;
    if (reSaveSolution)
    {
        loadIC = true;
        numIterRANS = 0;
        plotMeshes = false;
        animate = false;
        printInfo = false;
    }
    //-----------

    // ========================================= Geometry =========================================
    gsInfo << "Reading file geometry. \n";

    gsVector<int> inGeoInt(7);
    gsMatrix<int> refineUniformSettings(0,0);
    gsMatrix<int> refineLocalSettings(0,0);
    gsVector<real_t> geomParams(0);
    gsVector<bool> inGeoBool(2);
    std::vector<real_t> kvfit_knots;

    readInitialGeometry(MOTOR_DATA_DIR "uwb-pilsen/initialGeometry.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);
    //readInitialGeometry("initialGeometry.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool, kvfit_knots);

    if (!reSaveSolution)
    {
        gsInfo << "inGeoInt = \n" << inGeoInt << "\n";
        gsInfo << "inGeoBool = \n" << inGeoBool << "\n";
    }

    int geomChoice = inGeoInt(0);
    int index_of_profile = inGeoInt(1);
    int uniformRefine = inGeoInt(3);


    bool coarse = inGeoBool(0);
    bool makeLinearBool = inGeoBool(1);

    int makeLinearNumOfRefine = inGeoInt(2);

    gsVector<int> kvfit_print(kvfit_knots.size());
    for (size_t i = 0; i < kvfit_knots.size(); i++)
    {
        kvfit_print(i) = kvfit_knots[i];
    }

    //gsInfo << refineUniformSettings << "\n";
    //gsInfo << refineLocalSettings << "\n";

//    if (uniform_knot)
//        kvfit_knots = {0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1};
//    else
//        kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};

    //========================================================================================================================
    real_t omega = 0.;//2*PI*538.0/60.0;   // value 538 cycles per minute is given

    int setting_of_blade_rotation = 0; //0-optimal, 1-maximal, 2-minimal
    int setting_of_guide_vanes = 0; //0-compute with the formula, 1- 56° GV for optimal and 60° GV for maximal, 2-68° GV for optimal a 72° GV for maximal

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
    real_t linTol = 1e-9;

    int numThreads = 1; // number of threads for assembly

    if (!reSaveSolution)
        gsInfo << "Solving turbulent flow for 2d blade profile.\n";
    else
        gsInfo << "\nResave solution from loaded data.\n\n";

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
    problemSettings.setUyInZero(UyInZero);

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


    gsMultiPatch<real_t> patches_start = problemSettings.makeMultiPatch();
    gsFileData<> fdpatches_start;
    fdpatches_start << patches_start;
    fdpatches_start.save("filePatches_created.xml");
     gsWriteParaview( patches_start, "patches_start", 50000, true);

    //gsFileData<> fdp("filePatches_created.xml");
    //gsMultiPatch<> * patches_file = fdp.getAnyFirst< gsMultiPatch<> >().release();

    //      gsMultiPatch<> patches = *patches_file;
    //              gsInfo << patches << "\n";
    //      gsWriteParaview( patches, "patches2", 50000, true);

    //  gsMultiPatch<> patches = * patches_file;



//     if (makeLinearBool)
//         patches = makeLinear(patches_start, makeLinearNumOfRefine);
//     else
//         patches = patches_start;

     gsWriteParaview(patches_start, "patches_beforeBasisRefine", 50000, true);
     gsFileData<> fdpatches;
     fdpatches << patches_start;
     fdpatches.save("filePatches_beforeBasisRefine.xml");

    // ========================================= Define basis and refine =========================================
    gsMultiBasis<> tbasis_old(patches_start); // basis for RANS equations

    for(int i = 0; i < uniformRefine; i++)
        tbasis_old.uniformRefine();

    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;

    if (finer_grid)
    {
        tbasis_old.refine(0, box);
        tbasis_old.refine(1, box);
        tbasis_old.refine(2, box);

        tbasis_old.refine(1, box);
        tbasis_old.refine(2, box);

        box << 0, 0, 0, 0.25;
        tbasis_old.refine(0, box);
        tbasis_old.refine(1, box);
        tbasis_old.refine(2, box);
        box << 0, 0, 0.75, 1;
        tbasis_old.refine(0, box);
        tbasis_old.refine(1, box);
        tbasis_old.refine(2, box);
    }
    else
    {
        tbasis_old.refine(1, box);
    }

    for(int i = 0; i < refineUniformSettings.rows(); i++)
        problemSettings.refineBasisUniformZones(tbasis_old, refineUniformSettings(i,0), refineUniformSettings(i,1), refineUniformSettings(i,2), refineUniformSettings(i,3), refineUniformSettings(i,4));


    for(int i = 0; i < refineLocalSettings.rows(); i++)
        problemSettings.refineBasisLocalZones(tbasis_old, refineLocalSettings(i,0), refineLocalSettings(i,1), refineLocalSettings(i,2), refineLocalSettings(i,3), refineLocalSettings(i,4));



    //-------------------------------------------------------------------------------------------------
     gsMultiPatch<> patches;
      gsMultiBasis<> tbasis;

        if (makeLinearBool)
        {
           patches = makeLinearBasis(patches_start, tbasis_old, makeLinearNumOfRefine);
           gsMultiBasis<> tbasis_new(patches);
           tbasis = tbasis_new;

        }
        else
        {
           patches =  patches_start;
           tbasis = tbasis_old;
        }

    if (!reSaveSolution)
        gsInfo << "tbasis" << tbasis;


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
        file << "numIterRANS : " << numIterRANS << "\n";
        file << "minNumIterRANS : " << minNumIterRANS << "\n";
        file << "maxRANSPicardIt : " << maxRANSPicardIt << "\n";
        file << "tolRelNorm : " << tolRelNorm << "\n";
        file << "tmEvaluator = " << tmEvaluator << "\n";
        file << "maxTMPicardFirstIt : " << maxTMPicardFirstIt << "\n";
        file << "iterative : " << iterative << "\n";
        file << "lin. max iter. : " << linMaxIt << "\n";
        file << "lin. tol : " << linTol << "\n";
        file << "geomChoice : " <<geomChoice<< "\n";
        file << "knot vector : " << kvfit_print << "\n";
        file << "makeLinear : " << makeLinearBool << "\n";
        file << "makeLinearNumRef : " << makeLinearNumOfRefine << "\n";
        file << "uniformRefine : " <<uniformRefine<< "\n";
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

    if (!reSaveSolution)
        gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n\n";
    if (printInfo)
        file << "Velocity basis degree: " << discreteBases[0].degree() << "\n";

    // ========================================= Compute Steady NS =========================================
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.setNumThreads(numThreads); // set the number of threads for assembly

    /*if (SUPG)
    {
        params.settings().set(constantsINS::SUPG, SUPG);
        params.settings().set(constantsINS::tauStabTypeSUPG, tauStabTypeSUPG);
    }*/
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    if (!reSaveSolution)
        gsInfo << "numDofs = " << navStokes.numDofs() << "\n";

    if (!loadIC && NSsteadyInitial)
        problemSettings.computeSteadyNS(navStokes, numIterSteadyNS);

    //gsField<> pressureSteadyAtProfile = navStokes.constructSolution(1, problemSettings.getWallBoundaryPatches(), problemSettings.getWallBoundarySides());
    //gsWriteParaview<>(pressureSteadyAtProfile, problemSettings.getGeometryName() + "pressureSteadyAtProfile" + profile, plot_pts_atProfile);

    //================================= wall distance estimation ==================================================
    real_t wallDistance = problemSettings.computeWallDistance(navStokes);
     if (!reSaveSolution)
     {
        gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";
        if (printInfo)
            file << "\nminimum wallDistance = " << wallDistance << "\n";

        //==================================== compute aspect ratio  ==================================================
        real_t maxAspectRatio = navStokes.computeAspectRatio();
        gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";
        if (printInfo)
            file << "maxAspectRatio = " << maxAspectRatio << "\n";
     }
    // ========================================= Define turbulence solver =========================================
    gsBoundaryConditions<> bcInfoTurb;

    real_t uFreeStream = problemSettings.getFreeStream();
    real_t kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
    real_t oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet

    real_t kBlade = 1.5 * math::pow(velocity_blade(index_of_profile) * turbIntensity, 2);
    real_t oBlade;
    if (omega)
        oBlade = kBlade / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at the wall
    else
    {
        real_t beta = 0.0708;
        oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2));
    }

    if (!reSaveSolution)
    {
        gsInfo << "\nkInConst = " << kInConst << "\n";
        gsInfo << "oInConst = " << oInConst << "\n";
        gsInfo << "kBlade = " << kBlade << "\n";
        gsInfo << "oBlade = " << oBlade << "\n\n";
    }
    if (printInfo)
    {
        file << "\nkInConst = " << kInConst << "\n";
        file << "oInConst = " << oInConst << "\n";
        file << "kBlade = " << kBlade << "\n";
        file << "oBlade = " << oBlade << "\n\n";
    }

    problemSettings.defineBCs_TM(bcInfoTurb, kInConst, kBlade, oInConst, oBlade);

    gsDofMapper koMapper;
    discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPdeSteady(patches, bcInfoTurb, f, viscositySteady);
    uwbINSSolverParams<real_t> koParamsSteady(koPdeSteady, discreteBases, opt);

    koParamsSteady.setNumThreads(numThreads); // set the number of threads for assembly

    //if (TMsupg)
    //    koParamsSteady.settings().set(constantsINS::TMsupg, TMsupg);

    gsMatrix<> koInitial(koMapper.freeSize(), 2);
    koInitial.col(0).setConstant(kInConst);//(kBlade);
    koInitial.col(1).setConstant(oInConst);//(oBlade);

    uwbTMSolverKOmegaLinSteady<real_t> turbSolver(koParamsSteady);
    turbSolver.setInitialCondition(koInitial);
    // ========================================= Solving k-omega =========================================

    if (!loadIC && NSsteadyInitial)
        problemSettings.computeSteadyKOmegaTM(navStokes, turbSolver, numIterKOmegaSteady, koMapper.freeSize());

    // ========================================= Solving RANS =========================================
    if (!reSaveSolution)
    {
        gsInfo << "\nSolving RANS...\n";

        if (iterative)
        {
            gsInfo << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";
            if (printInfo)
            {
                file << "SOLVING with iterative solver, preconditioner type: " << precType << "\n\n";
                if (precType == "AL_Awhole" || precType == "AL_Amod")
                    file << "gamma: " << gamma << "\n\n";
            }
        }
        else
        {
            gsInfo << "SOLVING with direct solver.\n\n";
            if (printInfo)
                file << "SOLVING with direct solver.\n\n";
        }
    }

    uwbINSPde<real_t> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> RANSparams(RANSpde, discreteBases, opt);

    RANSparams.setNumThreads(numThreads); // set the number of threads for assembly

    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);
    //RANSparams.settings().set(constantsINS::unst_innerTol, picardTol);
    if (SUPG)
    {
        RANSparams.settings().set(constantsINS::SUPG, SUPG);
        RANSparams.settings().set(constantsINS::tauStabTypeSUPG, tauStabTypeSUPG);
        RANSparams.settings().set(constantsINS::diamPowerSUPG, diamPowerSUPG);
    }
    if (TCSD)
    {
        RANSparams.settings().set(constantsINS::TCSD, TCSD);
        RANSparams.settings().set(constantsINS::tauStabTypeSUPG, tauStabTypeSUPG);
    }
    if (RANScrosswind)
    {
        RANSparams.settings().set(constantsINS::RANScrosswind, RANScrosswind);
        RANSparams.settings().set(constantsINS::crosswindType, crosswindType);
        RANSparams.settings().set(constantsINS::tauStabTypeCW, tauStabTypeCW);
        RANSparams.settings().set(constantsINS::diamPowerCW, diamPowerCW);
        RANSparams.settings().set(constantsINS::CWresidualType, CWresidualType);
    }
    /*else if (INScrosswind)
    {
        RANSparams.settings().set(constantsINS::CROSSWIND, INScrosswind);
        RANSparams.settings().set(constantsINS::crosswindType, crosswindType);
        RANSparams.settings().set(constantsINS::tauStabTypeCW, tauStabTypeCW);
        RANSparams.settings().set(constantsINS::diamPowerCW, diamPowerCW);
        RANSparams.settings().set(constantsINS::CWresidualType, CWresidualType);
    }*/
    if (RANSad)
    {
        RANSparams.settings().set(constantsINS::RANSad, RANSad);
        RANSparams.settings().set(constantsINS::tauStabTypeAD, tauStabTypeAD);
        RANSparams.settings().set(constantsINS::diamPowerAD, diamPowerAD);
    }

    if (TMafc || TMafcHO)
        discreteBases[1].getMapper(dirichlet::none, opt.intStrategy, bcInfoTurb, koMapper, 0);
    else
        discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

    uwbINSPde<real_t> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, discreteBases, opt);

    koParams.setNumThreads(numThreads); // set the number of threads for assembly

    koParams.settings().set(constantsINS::timeStep, timeStep);
    koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
    koParams.settings().set(constantsINS::turb_innerIt, maxTMPicardIt);
    koParams.settings().setTurbulenceEvaluator(tmEvaluator);

    if (TMcrosswind)
    {
        koParams.settings().set(constantsINS::TMcrosswind, TMcrosswind);
        koParams.settings().set(constantsINS::crosswindType, TMcrosswindType);
        koParams.settings().set(constantsINS::tauStabTypeCW, tauTMStabTypeCW);
        koParams.settings().set(constantsINS::diamPowerCW, TMdiamPowerCW);
        koParams.settings().set(constantsINS::CWresidualType, TMCWresidualType);
    }
    if (TMsupg)
    {
        koParams.settings().set(constantsINS::TMsupg, TMsupg);
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauTMStabTypeSUPG); //set formula for stabilization parameter tau
        koParams.settings().set(constantsINS::diamPowerSUPG, TMdiamPowerSUPG);
    }
    else if (TMafc)
        koParams.settings().set(constantsINS::TMafc, TMafc);
    else if (TMafcHO)
        koParams.settings().set(constantsINS::TMafcHO, TMafcHO);
    else if (TMad)
    {
        koParams.settings().set(constantsINS::TMad, TMad);
        koParams.settings().set(constantsINS::tauStabTypeAD, tauTMStabTypeAD);
        koParams.settings().set(constantsINS::diamPowerAD, TMdiamPowerAD);
    }
    else if (TMisoAD)
    {
        koParams.settings().set(constantsINS::TMisoAD, TMisoAD);
        koParams.settings().set(constantsINS::isoADtype, isoADtype);
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauTMStabTypeSUPG);
    }
    else if (TMfctLowOrder)
        koParams.settings().set(constantsINS::TMfct_lowOrder, TMfctLowOrder);


    if (limitTMProduction)
    {
        koParams.settings().set(constantsINS::limitTMProduction, limitTMProduction);
        koParams.settings().set(constantsINS::productionXPoint, productionXPoint);
    }

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

    gsMatrix<> koInitialUnsteady(koMapper.freeSize(), 2);
    koInitialUnsteady.col(0).setConstant(kInConst);//(kBlade);
    koInitialUnsteady.col(1).setConstant(oInConst);//(oBlade);
    if (loadIC)
        turbSolver_unsteady.setInitialCondition(TMSolution);
    else if (NSsteadyInitial)
        turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());
    else
        turbSolver_unsteady.setInitialCondition(koInitialUnsteady);

    if (tmEvaluator != "koWilcoxLRN")
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
    if (!reSaveSolution)
        gsInfo << "initialization...\n";
    ransSolver.setComputationSequence(computeTMfirst);

    if (loadIC)
        ransSolver.setInitialCondition(RANSSolution);
    else if (NSsteadyInitial)
        ransSolver.setInitialCondition(navStokes.getSolution());

    ransSolver.initialize();

    //ransSolver.setStokesInitialCondition();*/

    real_t pitch = ((2 * PI*problemSettings.getRr(index_of_profile)) / num_blades);
    if (!reSaveSolution)
    {
        gsInfo << "numDofs: " << ransSolver.numDofs() << "\n";
        gsInfo << "pitch = " << pitch << "\n";
    }
    if (printInfo)
    {
        file << "numDofs: " << ransSolver.numDofs() << "\n";
        file << "pitch = " << pitch << "\n";
    }

    //=== solving with solution change relative norm as stopping criterion ============
    if (!reSaveSolution)
    {
        if (animate)
            ransSolver.solveWithAnimation(numIterRANS, animateStep, tolRelNorm, plot_pts, true);
        else
            ransSolver.solve(numIterRANS, tolRelNorm); // solution change norm tol = 10^(-5)

        real_t Tassembly = ransSolver.getAssemblyTime();
        real_t Tsetupsolve = ransSolver.getSolverSetupTime();
        real_t Tsolve = ransSolver.getSolveTime();
        real_t Tturbassembly = ransSolver.getTurbAssemblyTime();
        real_t Tturbsetupsolve = ransSolver.getTurbSolverSetupTime();
        real_t Tturbsolve = ransSolver.getTurbSolveTime();

        gsInfo << "Assembly time:" << Tassembly << "\n";
        gsInfo << "Solve time:" << Tsolve << "\n";
        gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
        gsInfo << "Turbulent assembly time:" << Tturbassembly << "\n";
        gsInfo << "Turbulent solver setup time:" << Tturbsetupsolve << "\n";
        gsInfo << "Turbulent solve time:" << Tturbsolve << "\n";

        if (printInfo)
        {
            file << "Assembly time:" << Tassembly << "\n";
            file << "Solve time:" << Tsolve << "\n";
            file << "Solver setup time:" << Tsetupsolve << "\n";
            file << "Turbulent assembly time:" << Tturbassembly << "\n";
            file << "Turbulent solver setup time:" << Tturbsetupsolve << "\n";
            file << "Turbulent solve time:" << Tturbsolve << "\n";
            file.close();
        }

        //--- save solution into file ---
        gsFileData<> fd, fd_TM;
        fd << ransSolver.getSolution();
        fd_TM << turbSolver_unsteady.getSolution();
        fd.save(problemSettings.getGeometryName() + "RANS_solution"  + profile + ".xml");
        fd_TM.save(problemSettings.getGeometryName() + "TM_solution"  + profile + ".xml");
    }

    //=====================================================================================================

    // Optionally plot solution in paraview
    if (plot)
    {
        gsInfo << "Plotting in Paraview...\n";

        gsField<> velocity = ransSolver.constructSolution(0);
        gsField<> pressure = ransSolver.constructSolution(1);

        const gsGeometry<> * geoResult = & velocity.igaFunction();
        gsWriteParaview(* geoResult,"geoResult",10000,true,true);

//        gsMatrix<> IGAvelsol =  velocity.igaFunction().coefs();
//        gsMatrix<> IGApresssol =  pressure.igaFunction().coefs();

//        gsInfo <<  "coefsSize" << velocity.igaFunction().coefsSize();

//        const gsBasis<real_t> & IGAvelbasis = velocity.igaFunction().basis();
//        const gsBasis<real_t> & IGApressbasis = pressure.igaFunction().basis();

//        gsMesh<> velCN;
//        gsMesh<> pressCN;
//        velocity.igaFunction().controlNet(velCN);
//        pressure.igaFunction().controlNet(pressCN);

//        gsInfo << "IGAvelsol" << IGAvelsol;
//        gsInfo << "IGApresssol" << IGApresssol;

//        gsWrite(IGAvelsol,"IGAvelsol.xml");
//        gsWrite(IGApresssol,"IGApresssol.xml");
//        //gsWrite(velCN,"IGAvelCN.xml");
//        //gsWrite(pressCN,"IGApressCN.xml");
//        gsWrite(IGAvelbasis,"IGAvelbasis.xml");
//        gsWrite(IGApressbasis,"IGApressbasis.xml");


        gsField<> kOmega = turbSolver_unsteady.constructSolution();

        problemSettings.plotSolutionField(velocity, problemSettings.getGeometryName() + "velocity");
        problemSettings.plotSolutionField(velocity, problemSettings.getGeometryName() + "velocity", true);
        problemSettings.plotSolutionField(pressure, problemSettings.getGeometryName() + "pressure");
        problemSettings.plotSolutionField(kOmega, problemSettings.getGeometryName() + "komega");
        ransSolver.plotTurbulentViscosity(problemSettings.getGeometryName() + "turb_viscosity" + profile, plot_pts);

        gsField<> pressureAtProfile = ransSolver.constructSolution(1, problemSettings.getWallBoundaryPatches(), problemSettings.getWallBoundarySides());
        gsWriteParaview<>(pressureAtProfile, problemSettings.getGeometryName() + "pressureAtProfile" + profile, plot_pts_atProfile);

        gsInfo << "Plotting done.\n";
    }

    return 0;
}
