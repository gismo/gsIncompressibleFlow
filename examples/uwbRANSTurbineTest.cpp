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

using namespace gismo;

//template<class TT>
//gsMultiPatch<TT>  BSplineTurbine3D();
std::string omegaXr_component(real_t omega, std::string variable);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, real_t uMax = 1, real_t omega = 0);
template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, real_t kIn, real_t kWall, real_t kBlade, real_t oIn, real_t oWall, real_t oBlade);
template<class T> void refineBasis(gsMultiBasis<T>& tbasis, int numRefine, int numRefineLocal);

int main(int argc, char *argv[])
{
    const real_t PI = 3.141592653589793238463;

    // ========================================= Settings =========================================

    bool plot2 = true;
    bool outFile = true;
    int plot_pts = 30000;
    const int numThreads = 4;

    real_t timeStep = 0.002;

    int numIterSteadyNS = 20;
    int numIterKOmegaSteady = 50;
    int numIterRANS = 50;

    real_t omega = 0; // angular velocity for rotation

    int numRefine = 0;
    int numRefineLocal = 0;

    real_t viscosity = 0.1;
    real_t viscositySteady = 0.5;
    real_t uMax = 1;
    real_t turbIntensity = 0.01;

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

    gsInfo << "Solving turbulent flow in the Kaplan Turbine.\n";

    // ========================================= Define problem =========================================

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, uMax, omega);

    gsFunctionExpr<> f("0", "0", "0", 3);

    // ========================================= Define geometry =========================================

    bool plot = false; // If set to true, paraview file is generated and launched on exit
    bool simplified_HP = true;

    // READING THE PARAMETERS FOR HYDRAULIC PROFILE DEFINITION FROM THE FILE
    std::string inputFile = MOTOR_DATA_DIR "uwb-pilsen/Parameters_HydraulicProfile.txt";

    std::map <std::string, real_t> param_list;

    std::ifstream cFile;
    cFile.open(inputFile.c_str());
        if (cFile.is_open())
        {
            std::string line;
            while(getline(cFile, line)){
                line.erase(remove_if(line.begin(), line.end(), isspace),
                                     line.end());
                if(line[0] == '#' || line.empty())
                    continue;
                char delimiterPos = line.find("=");
                std::string name = line.substr(0, delimiterPos);
                std::string value_read = line.substr(delimiterPos + 1);
                real_t value;
                if (value_read[0]=='-'){
                    std::string value_read_num = value_read.substr(1);
                value = - atof(value_read_num.c_str());}
                else {value = atof(value_read.c_str());}
                //value = atof(value_read.c_str());
                param_list[name] = value;
            }
            cFile.close();
        }
        else
           { gsWarn << "Unable to open param file." << '\n'; }

    // SETTING PARAMETERS FOR INNER AND OUTER SURFACE OF KAPLAN TURBINE
    // matrix with parameters has the following structure
    // first line: start_x, start_y, 0.0, 0.0
    // line : length, height_start, height_end, 1.0
    // circular arc: centreX, centreY, alpha, 2.0   ...//radius is computed from the previous section
    gsMatrix<real_t> param_inner(12,4);
    param_inner << -1.2, param_list["v1"], 0.0, 0.0,
            param_list["l1"], param_list["v1"], param_list["v1"], 1.0,
            param_list["d1"] - param_list["l1"], param_list["v1"], param_list["v1"], 1.0,
            param_list["d2"], param_list["v1"], param_list["v2"], 1.0,
            param_list["Sx"],0.0,  param_list["alfa02"],2.0,
            param_list["d3"],0.0,param_list["v3"],1.0,
            param_list["d4"],param_list["v3"],param_list["v3"],1.0,
            0.0, 0.0, param_list["alfa031"], 2.0,
            0.0, 0.0, param_list["alfa032"], 2.0,
            param_list["Sx6"],  param_list["Sy6"], param_list["alfa6"], 2.0,
            param_list["d6"], 0.0, param_list["v6"], 1.0,
            param_list["Sx7"], param_list["Sy7"], param_list["alfa7"], 2.0;

    gsMatrix<real_t> param_outer(10,4);
    param_outer << - 1.2, param_list["h1"], 0.0, 0.0,
            param_list["l1"], param_list["h1"], param_list["h1"], 1.0,
            0.75*param_list["l2"], param_list["h1"], 0.25*param_list["h1"]+0.75*param_list["h2"], 1.0,
            0.25*param_list["l2"], 0.75*param_list["h2"], param_list["h2"], 1.0,
            param_list["Sx"],0.0, param_list["alfa"],2.0,
            param_list["Sx2"],param_list["Sy2"], param_list["alfa2"],2.0,
            param_list["l3"], param_list["h3"], param_list["h3"],1.0,
            0.0,0.0, param_list["alfa3"],2.0,
            param_list["Sx4"],param_list["Sy4"],  param_list["alfa4"],2.0,
            param_list["Sx5"],param_list["Sy5"],  param_list["alfa5"],2.0;

    // DEFINITION AND COMPUTATION OF HYDRAULIC PROFILE OF KAPLAN TURBINE
    HydraulicProfile<real_t> mHydraulicProfile(param_inner, param_outer, simplified_HP);
    mHydraulicProfile.computeHydraulicProfile();

    gsVector<real_t> rotation_axis(3);
    rotation_axis << 0, 1, 0;
    gsTensorBSpline<2, real_t> surfaceInner;
    gsTensorBSpline<2, real_t> surfaceInner_rotated;
    surfaceInner = mHydraulicProfile.getHydraulicProfileInner();
    surfaceInner_rotated = surfaceInner;
    surfaceInner_rotated.rotate(PI/2, rotation_axis);

    gsTensorBSpline<2, real_t> surfaceOuter;
    gsTensorBSpline<2, real_t> surfaceOuter_rotated;
    surfaceOuter = mHydraulicProfile.getHydraulicProfileOuter();
    surfaceOuter_rotated = surfaceOuter;
    surfaceOuter_rotated.rotate(PI/2, rotation_axis);
    gsInfo << "HydraulicProfileComputed \n";

    int num_bladeprofiles = 7;

    gsVector<> camber_x(num_bladeprofiles);
    gsVector<> camber_y(num_bladeprofiles);
    gsVector<> leading_angle(num_bladeprofiles);
    gsVector<> trailing_angle(num_bladeprofiles);
    gsVector<> thickness_x(num_bladeprofiles);
    gsVector<> thickness_y(num_bladeprofiles);
    gsVector<> ending_offset(num_bladeprofiles);
    gsVector<> output_angle(num_bladeprofiles);
    gsVector<> radius(num_bladeprofiles);
    gsVector<> chord_length(num_bladeprofiles);
    gsVector<> angle(num_bladeprofiles);
    gsVector<> rotation_center_x(num_bladeprofiles);
    gsVector<> rotation_center_y(num_bladeprofiles);
    gsVector<> rr(num_bladeprofiles);

    // SETTING PARAMETERS FOR RUNNER WHEEL BLADE OF KAPLAN TURBINE
    rr << 0.175, 0.229, 0.283, 0.338, 0.392, 0.446, 0.5;
    camber_x << 0.425908805, 0.419, 0.416039882, 0.416576484, 0.416716376, 0.419203131, 0.43280567; // 0.411886743 na originalne druhe pozici
    //gsInfo << camber_x << "\n";
    camber_y << 0.095896144, 0.064997436, 0.037083707, 0.02724709, 0.024356984, 0.023262639, 0.019802704;
    leading_angle << 38.47692, 21.399859, 12.641776, 9.28275, 8.38282, 8.338553, 8.446091;
    leading_angle = leading_angle*PI/180;
    trailing_angle << 20.717182, 10.721469, 6.371868, 4.702573, 4.217487, 4.164196, 4.233241;
    trailing_angle = trailing_angle*pi/180;
    thickness_x << 0.281408739, 0.275195139, 0.273161229, 0.272459059, 0.272310247, 0.272053329, 0.271637628;
    thickness_y << 0.064950942, 0.044591162, 0.030435219, 0.020799775, 0.01464564, 0.011033004, 0.009989223;
    //ending_offset << 0.001184, 0.000812, 0.000555, 0.000379, 0.000267, 0.000201, 0.000182;
    ending_offset << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    output_angle << 8.381088, 5.891258, 4.042351, 2.765094, 1.950139, 1.466739, 1.329705;
    output_angle = output_angle*PI/180;
    radius << 0.020637, 0.014575, 0.007271, 0.003687, 0.001955, 0.001012, 0.000828;
    chord_length << 0.36686418, 0.42813715,	0.486298848, 0.545457321, 0.595371262, 0.615928719, 0.584588959;
    angle << 53.696487,	41.265848, 33.007703, 27.603276, 24.437586, 22.893162, 21.162381;
    angle = angle*PI/180;
    rotation_center_x << 0.494758396, 0.469497406, 0.444542963, 0.417724545, 0.390108787, 0.361175154, 0.330805204;
    rotation_center_y << 0.060495569, 0.028225794, 0.00125711, -0.006884641, -0.010228889, -0.010435203, -0.00079539;


    // DEFINITION AND COMMPUTATION OF RUNNER WHEEL BLADE OF KAPLAN TURBINE
    KaplanTurbineRunnerBlade<real_t> runnerBlade(num_bladeprofiles, camber_x, camber_y, leading_angle, trailing_angle, thickness_x, thickness_y, ending_offset,
                                         output_angle, radius, chord_length, angle, rotation_center_x, rotation_center_y, rr);
    runnerBlade.compute(plot, true);
    runnerBlade.trimming(surfaceInner, plot, true);


    // DEFINITION AND COMPUTATION OF RUNNER WHEEL DOMAIN (DOMAIN BETWEEN TWO RUNNER WHEEL BLADES) OF KAPLAN TURBINE
    KaplanTurbineRunnerWheelDomain<real_t> runnerWheelDomain(mHydraulicProfile, runnerBlade, 4);
    runnerWheelDomain.compute(plot, true);

    gsMultiPatch<> patches;
    patches = runnerWheelDomain.getRunnerDomain();

    gsMultiPatch<> TMpatches = patches;

    TMpatches.addInterface(0, boundary::front, (size_t) 0, boundary::back);
    TMpatches.addInterface(2, boundary::front, 2, boundary::back);
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

    // solvers
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver

    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    navStokes.initialize(); // steady solver
    navStokes.solve(numIterSteadyNS, 1e-5);
    //gsMatrix<> steadySolution = navStokes.getSolution();

    gsField<> absVelocitySteady = navStokes.constructSolution(0);
    gsWriteParaview<>(absVelocitySteady, "turbine_absSteadyVelocity", plot_pts, true);

    gsFileData<> fd;
    fd << navStokes.getSolution();
    fd.save("turbine_steadySolution_nu" + util::to_string(viscositySteady) + ".xml");
    //---------------------

    /*gsInfo << "Computing initial velocity by solving Stokes problem...\n";

    uwbINSPde<real_t> NSpde(*patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    uwbINSSolverSteady<real_t> stokesSolver(params);

    gsInfo << "numDofs: " << stokesSolver.numDofs() << "\n";

    stokesSolver.initialize();

    stokesSolver.setStokesSolution();*/

    //================================= wall distance estimation ==================================================
    real_t inletWidth = 0.355046; //
    gsVector<real_t> Uin(3);
    Uin << uMax, 0, 0;
    real_t uFreeStream = Uin.norm();
    real_t Re = uFreeStream * inletWidth / viscosity;

    //vector of the sides of the patches from which the wall distance is computed
    std::vector<boxSide> distanceSidesWalls;
    distanceSidesWalls.push_back(boundary::east);
    distanceSidesWalls.push_back(boundary::west);
    distanceSidesWalls.push_back(boundary::east);
    distanceSidesWalls.push_back(boundary::back);
    distanceSidesWalls.push_back(boundary::east);
    distanceSidesWalls.push_back(boundary::west);
    std::vector<boxSide> distanceSidesBlade;
    distanceSidesBlade.push_back(boundary::front);
    distanceSidesBlade.push_back(boundary::back);

    //vector of indexes of the patches corresponding to distanceSides
    //length of the vector distancePatches must be equal to the length of vector distanceSides
    gsVector<int> distancePatchesWalls(6);
    distancePatchesWalls << 0, 0, 1, 1, 2, 2;
    gsVector<int> distancePatchesBlade(2);
    distancePatchesBlade << 1, 1;

    int numSamplePts = 500; //number of sample points for which the distance to the boundary is computed
    real_t maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

    //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
    //estimation of the wall distance will be computed, if the last input parameter is set as true
    real_t wallDistanceWalls = navStokes.computeDimensionlessWallDistance(distancePatchesWalls, distanceSidesWalls, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
    gsInfo << "\nminimum wallDistance at walls = " << wallDistanceWalls << "\n";
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

    //gsStopwatch time; -> turbModelTime
    turbSolver.initialize(velocity);
    //real_t Tassembly = time.stop();
    //gsInfo << "Assembly time:" << Tassembly << "\n";

    //time.restart();
    turbSolver.solve(numIterKOmegaSteady, 1e-5); // solution change norm tol = 10^(-5)
    //real_t Tsolve = time.stop();

    //gsInfo << "Solve time:" << Tsolve << "\n";

    gsField<> kOmegaSteady = turbSolver.constructSolution();
    gsWriteParaview<>(kOmegaSteady, "turbine_kOmegaSteady", plot_pts, true);

    gsFileData<> fdTM;
    fdTM << turbSolver.getSolution();
    fdTM.save("turbine_TMsteadySolution_nu" + util::to_string(viscositySteady) + ".xml");

    // ========================================= Solving RANS =========================================
    gsInfo << "\nSolving RANS...\n";

    uwbINSPde<real_t> RANSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> RANSparams(RANSpde, discreteBases, opt, numThreads);

    uwbINSPde<real_t> koPde(TMpatches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, TMdiscreteBases, opt, numThreads);

    if (TMsupg) {
        gsWarn << "SUPG stabilization only for turbulent model not for RANS equations.\n";
        koParams.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType); //set formula for stabilization parameter tau
    }

    uwbTMSolverKOmega<real_t> turbSolver_unsteady(koParams);
    turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());

    if (omega)
        RANSparams.settings().set(constantsINS::omega, omega); // set rotation (assumed around x-axis)

    RANSparams.settings().set(constantsINS::timeStep, timeStep);
    RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);
    uwbRANSSolver<real_t> ransSolver(RANSparams);

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
        filename << "turbine_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
            << now->tm_hour << "-" << now->tm_min << ".txt";

        std::ofstream ofile;
        ofile.open(filename.str());
        ofile << "Solving RANS problem in the runner wheel:\n";
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
    if (plot2)
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
        gsWriteParaview<>(absVelocity, "turbine_RANS_absVelocity", plot_pts, true);
        gsWriteParaview<>(relVelocity, "turbine_RANS_relVelocity", plot_pts);
        //gsWriteParaview<>(combVelocity, "turbine_RANS_combVelocity", plot_pts);
        gsWriteParaview<>(pressure, "turbine_RANS_pressure", plot_pts);;
        gsWriteParaview<>(kOmega, "turbine_TM_komega", plot_pts);
        ransSolver.plotTurbulentViscosity("turbine_turb_viscosity", plot_pts);
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

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, real_t uMax, real_t omega)
{
    gsMatrix<real_t> transformMatrix(3, 3);

    const real_t phi = (1. / 2)*EIGEN_PI;
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

    gsFunctionExpr<> Uin(util::to_string(uMax), "0", "0", 3);
    gsFunctionExpr<> Uwall("0", "0", "0", 3); // wall velocity
    gsFunctionExpr<>Ublade("0", omegaXr_component(-omega, "z"), omegaXr_component(omega, "y"), 3);

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);

    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::front, condition_type::dirichlet, Ublade, 0);
    bcInfo.addCondition(1, boundary::back, condition_type::dirichlet, Ublade, 0);

    bcInfo.setIdentityMatrix(3);
    bcInfo.addPeriodic(0, boundary::front, 0, boundary::back, 3);
    bcInfo.addPeriodic(2, boundary::front, 2, boundary::back, 3);

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

    bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, Kin, 0);
    bcInfoTurb.addCondition(0, boundary::east, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Kwall, 0);

    bcInfoTurb.addCondition(2, boundary::east, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(2, boundary::west, condition_type::dirichlet, Kwall, 0);

    bcInfoTurb.addCondition(1, boundary::east, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Kwall, 0);
    bcInfoTurb.addCondition(1, boundary::front, condition_type::dirichlet, Kblade, 0);
    bcInfoTurb.addCondition(1, boundary::back, condition_type::dirichlet, Kblade, 0);
    //---

    bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, Oin, 1);
    bcInfoTurb.addCondition(0, boundary::east, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Owall, 1);

    bcInfoTurb.addCondition(2, boundary::east, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(2, boundary::west, condition_type::dirichlet, Owall, 1);

    bcInfoTurb.addCondition(1, boundary::east, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Owall, 1);
    bcInfoTurb.addCondition(1, boundary::front, condition_type::dirichlet, Oblade, 1);
    bcInfoTurb.addCondition(1, boundary::back, condition_type::dirichlet, Oblade, 1);
}

template<class T> void refineBasis(gsMultiBasis<T>& tbasis, int numRefine, int numRefineLocal)
{
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    real_t uRefineKnot = 0.2;
    //real_t wallRefineKnot = 0.1;

    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();

    gsMatrix<> box_u0(3, 2);
    //gsMatrix<> box_u1(2, 2);
    //gsMatrix<> box_v0(2, 2);
    //gsMatrix<> box_v1(2, 2);
    gsMatrix<> box_w0(3, 2);
    //gsMatrix<> box_w1(2, 2);

    //box_u0 << 0.8, 1, 0, 0, 0, 0;
    //tbasis.refine(0, box_u0);
    //tbasis.refine(1, box_u0);
    //tbasis.refine(2, box_u0);

    //refine in u, near bottom wall
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_u0 << 1 - (uRefineKnot / math::pow(2, i)), 1, 0, 0, 0, 0;
        gsInfo << box_u0 << "\n\n";

        tbasis.refine(0, box_u0);
        tbasis.refine(1, box_u0);
        tbasis.refine(2, box_u0);
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    //refine in u, near upper wall
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_u0 << 0, uRefineKnot / math::pow(2, i), 0, 0, 0, 0;
        gsInfo << box_u0 << "\n\n";

        tbasis.refine(0, box_u0);
        tbasis.refine(1, box_u0);
        tbasis.refine(2, box_u0);
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    //refine in w, near blade
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_w0 << 0, 0, 0, 0, 1 - (uRefineKnot / math::pow(2, i)), 1;
        gsInfo << box_w0 << "\n\n";

        tbasis.refine(0, box_w0);
        tbasis.refine(1, box_w0);
        tbasis.refine(2, box_w0);
    }
    gsInfo << "Number of basis functions: " << tbasis.piece(0).size() << ", " << tbasis.piece(1).size() << ", " << tbasis.piece(2).size() << "\n";

    //refine in w, near blade
    for (int i = 0; i < numRefineLocal; i++)
    {
        box_w0 << 0, 0, 0, 0, 0, uRefineKnot / math::pow(2, i);
        gsInfo << box_w0 << "\n\n";

        tbasis.refine(0, box_w0);
        tbasis.refine(1, box_w0);
        tbasis.refine(2, box_w0);
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
