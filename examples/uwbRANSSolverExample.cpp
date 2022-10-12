/** @file uwbINSSolversExample.cpp

Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers
#include <gsIncompressibleFlow/src/uwbINSSolverSteady.h>
#include <gsIncompressibleFlow/src/uwbRANSSolver.h>
#include <gsIncompressibleFlow/src/uwbTMSolverKOmega.h>
#include <gsIncompressibleFlow/src/uwbTMSolverKOmegaLinSteady.h>

using namespace gismo;

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a = 1, T const & b = 2, T const & a_in = 1);
template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax = 1);
template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, T kIn, T kWall, T oIn, T oWall);
template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal);
template<class T> void computeSteadyNS(uwbINSSolverSteady<T>& navStokesSteady, int numIterNSSteady, int plot_pts = 10000);
template<class T> T computeWallDistance(uwbINSSolverSteady<T>& navStokesSteady, T Re, T viscosity, T uMax);
template<class T> void solvePoissonEquation(gsMultiPatch<T> patches, uwbTMSolverKOmega<T>& turbSolver, int plot_pts = 10000);

int main(int argc, char *argv[])
{
    // ========================================= Settings ========================================= 
    //---------------------------------------------------------------------------------------------
    /*gsVector<int> inInt(11);
    gsVector<real_t> inRealT(6);
    std::string inString;
    gsVector<bool> inBool(4);

    std::ifstream inFile;
    inFile.open("initialSettings.txt");
    if (!inFile) {
        gsInfo << "Unable to open file";
        exit(1); // terminate with error
    }

    for(int i = 0; i < inInt.rows(); i++)
        inFile >> inInt(i);
    for(int i = 0; i < inRealT.rows(); i++)
        inFile >> inRealT(i);
    inFile >> inString;
    for(int i = 0; i < inBool.rows(); i++)
        inFile >> inBool(i);

    inFile.close();

    // ========================================= Settings =========================================
    bool plot = inBool(0);
    bool plotMeshes = inBool(1);
    int plot_pts = inInt(0);

    real_t timeStep = inRealT(0);

    int numIterNSSteady = inInt(1);
    int numIterKOmegaSteady = inInt(2);
    //int maxNumIterRANS = inInt(3);
    int numIter = inInt(4);
    //int minNumIterRANS = inInt(5);
    real_t tol = inRealT(1);
    int maxRANSPicardIt = inInt(6);
    int maxTMPicardFirstIt = inInt(7);

    real_t viscosity = inRealT(2);
    real_t viscositySteady = inRealT(3);
    real_t turbIntensity = inRealT(4);
    real_t viscosityRatio = inRealT(5);

    std::string tmEvaluator = inString;

    //int refineType = inInt(9);

    bool TMsupg = inBool(2);
    int tauSUPGType = inInt(10);
    //bool supg = inBool(3);*/

    //----------------------------------------------------

    bool plot = true;
    int plot_pts = 10000;

    real_t viscositySteady = 0.1;
    int numIterNSSteady = 10; // max number of time steps

    real_t viscosity = 0.0001;
    real_t timeStep = 0.001; // time step for unsteady computation
    real_t tol = 1e-5; // stopping tolerance
    int numIter = 1; // max number of time steps
    int maxRANSPicardIt = 5;

    real_t turbIntensity = 0.01;
    //real_t tmTfheta = 1; // theta for time discretization of TM (0 - forward Euler, 0.5 - Crank-Nicolson, 1 - backward Euler)
    //int turbInnerIt = 5; // max iterations of TM (turb. model)
    int maxTMPicardFirstIt = 50;
    //real_t turbInnerTol = 1e-4; // stopping tolerance for TM
    real_t viscosityRatio = 50.; // nu_T / nu

    bool TMsupg = false;
    int tauSUPGType = 0; // 0 'h/(2*deg*norm(u)) * (cotgh(Pe)) - 1/Pe'; 1 'h/(2*deg*norm(u)'; 2 '((2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
                         // 3 '((2/timeStep)^2 + (2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
    std::string tmEvaluator = "koSST"; // choose "koWilcoxLRN" or "koSST"

    int numRefine = 2;
    int numRefineLocal = 2;
    bool solveSteadyNS = true;
    real_t uMax = 2; // inlet velocity maximum
    int numThreads = 1; // number of threads for assembly
    bool afc = false; // AFC stabilization of TM

    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("l", "localRefine",
        "Number of local h-refinement steps to perform before solving", numRefineLocal);
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
    opt.intStrategy = iFace::glue;

    gsInfo << "Solving the backward step RANS example.\n";

    // ========================================= Define problem ========================================= 

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, uMax);

    gsFunctionExpr<> f("0", "0", 2);

    // ========================================= Define geometry ========================================= 
    
    gsMultiPatch<> patches;

    real_t a = 8;
    real_t b = 2;
    real_t a_in = 1;

    patches = BSplineStep2D<real_t>(a, b, a_in);

    gsInfo << patches << "\n";

    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);
    
    refineBasis_NS(tbasis, numRefine, numRefineLocal);
    
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    // ========================================= Compute Steady NS =========================================
    uwbINSPde<real_t> NSpdeSteady(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> paramsSteady(NSpdeSteady, discreteBases, opt);
    uwbINSSolverSteady<real_t> navStokesSteady(paramsSteady); // steady coupled solver

    if (solveSteadyNS)
        computeSteadyNS(navStokesSteady, numIterNSSteady, plot_pts);

    //================================= wall distance estimation ==================================================
    real_t inletWidth = b/2;
    real_t Re = uMax * inletWidth / viscosity;
    real_t wallDistance = computeWallDistance(navStokesSteady, Re, viscosity, uMax);

    //==================================== compute aspect ratio  ==================================================
    real_t maxAspectRatio = navStokesSteady.computeAspectRatio();
    gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";

    // ========================================= Define turbulence solver =========================================
    gsBoundaryConditions<> bcInfoTurb;

    // k and omega inlet values
    real_t kInConst = 1.5 * math::pow(uMax * turbIntensity, 2); // (3/2)*(UI)^2
    //real_t oInConst = math::pow(0.09, -0.25) * math::sqrt(kInConst) / 0.07; // (C_mu)^(-1/4) * sqrt(k)/l, l = 0.07 * d
    real_t oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
    real_t kWall = 1e-10; // k on the wall
    real_t oWall;
    if (wallDistance > 0)
    {
        real_t beta = 0.0708;
        oWall = 6 * viscosity / (beta * math::pow(wallDistance, 2)); // omega on the wall
    }
    else
        oWall = 15;

    gsInfo << "wallDistance = " << wallDistance << "\n";
    gsInfo << "kInConst = " << kInConst << "\n";
    gsInfo << "oInConst = " << oInConst << "\n";
    gsInfo << "kWall = " << kWall << "\n";
    gsInfo << "oWall = " << oWall << "\n";

    defineBCs_TM(bcInfoTurb, kInConst, kWall, oInConst, oWall);

    uwbINSPde<real_t> koPde(patches, bcInfoTurb, f, viscosity);
    uwbINSSolverParams<real_t> koParams(koPde, discreteBases, opt);

    if (TMsupg) {
        koParams.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType); //set formula for stabilization parameter tau
    }

    if (afc)
        koParams.settings().set(constantsINS::TMafc, afc);

    //koParams.settings().set(constantsINS::theta, tmTheta);
    //koParams.settings().set(constantsINS::turb_innerIt, turbInnerIt);
    //koParams.settings().set(constantsINS::turb_innerTol, turbInnerTol);
    koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
    koParams.settings().setTurbulenceEvaluator(tmEvaluator);

    uwbTMSolverKOmega<real_t> turbSolver(koParams);

    gsMatrix<> koInitial(turbSolver.getAssembler()->numVarDofs(), 2);
    koInitial.col(0).setConstant(kWall);
    koInitial.col(1).setConstant(oWall);

    turbSolver.setInitialCondition(koInitial);

    if (tmEvaluator == "koSST")
        solvePoissonEquation(patches, turbSolver, plot_pts);

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);
    params.settings().setTurbulenceSolver(&turbSolver);

    // --------- other options ---------
    params.setNumThreads(numThreads); // set the number of threads for assembly
 //   params.settings().set(constantsINS::unst_innerIt, 5); // max number of inner Picard iter. for unstready
 //   params.settings().set(constantsINS::unst_innerTol, 1e-4); // stopping tolerance for inner Picard iter. in unstready
 // ---------------------------------

    // solver
    uwbRANSSolver<real_t> navStokes(params);

    // ========================================= Solving ========================================= 
    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    gsInfo << "initialization...\n";
    navStokes.initialize();


    if (solveSteadyNS)
        navStokes.setInitialCondition(navStokesSteady.getSolution());
    else
        navStokes.setStokesInitialCondition();
 

    navStokes.solve(numIter, tol); 
    //navStokes.solveWithAnimation(numIter, 5, tol, plot_pts, plotTurb); // plot solution every 5 time steps (plotTurb = true -> plot also TM solution)

    real_t Tassembly = navStokes.getAssemblyTime();
    real_t Tsolve = navStokes.getSolverSetupTime();
    real_t Tsetupsolve = navStokes.getSolveTime();
    real_t Tturbmodel = navStokes.getTurbModelTime();
    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";
    gsInfo << "Turbulent model time:" << Tturbmodel << "\n";

    // Optionally plot solution in paraview
    if (plot)
    {
        gsField<> velocity = navStokes.constructSolution(0);
        gsField<> pressure = navStokes.constructSolution(1);
        gsField<> kOmega = turbSolver.constructSolution(); // k and omega plotted as two components of a vector fcn

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "step_velocity", plot_pts);
        gsWriteParaview<>(pressure, "step_pressure", plot_pts);
        gsWriteParaview<>(kOmega, "step_komega", plot_pts);
        navStokes.plotTurbulentViscosity("step_turb_viscosity", plot_pts);
    }


    return 0; 
}

template<class T> gsMultiPatch<T> BSplineStep2D(T const & a, T const & b, T const & a_in)
{
    gsMultiPatch<T> mp;

    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, a, b / 2));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, b / 2, a, b));
    mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-a_in, b / 2, 0.0, b));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);
    mp.addAutoBoundaries();

    return mp;
}

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, T uMax)
{
    std::string xVel = util::to_string(uMax) + " * (-4*(y-1.5)^2 + 1)";
    gsFunctionExpr<T> Uin(xVel, "0", 2); // inlet velocity
    gsFunctionExpr<T> Uwall("0", "0", 2); // wall velocity

    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);
}

template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, T kIn, T kWall, T oIn, T oWall)
{
    // Boundary conditions
    gsFunctionExpr<T> Kin(util::to_string(kIn), 2);
    gsFunctionExpr<T> Oin(util::to_string(oIn), 2);
    gsFunctionExpr<T> K(util::to_string(kWall), 2);
    gsFunctionExpr<T> O(util::to_string(oWall), 2);

    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Kin, 0);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, K, 0);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, K, 0);

    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Oin, 1);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, O, 1);
    bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, O, 1);
}

template<class T>
void computeSteadyNS(uwbINSSolverSteady<T>& navStokesSteady, int numIterNSSteady, int plot_pts)
{
    gsInfo << "Solving Steady case: \n";
    gsInfo << "numDofs: " << navStokesSteady.numDofs() << "\n";

    navStokesSteady.initialize(); // steady solver
    navStokesSteady.solve(numIterNSSteady, 1e-5);

    gsField<> velocitySteady = navStokesSteady.constructSolution(0);
    gsField<> pressureSteady = navStokesSteady.constructSolution(1);
    gsWriteParaview<>(velocitySteady, "step_SteadyVelocity", plot_pts, true);
    gsWriteParaview<>(pressureSteady, "step_SteadyPressure", plot_pts);

    gsFileData<real_t> fd;
    fd << navStokesSteady.getSolution();
    fd.save("step_NSsteadySolution.xml");
}

template<class T>
T computeWallDistance(uwbINSSolverSteady<T>& navStokesSteady, T Re, T viscosity, T uMax)
{
    //vector of the sides of the patches from which the wall distance is computed
    std::vector<boxSide> distanceSides;
    distanceSides.push_back(boundary::west);
    distanceSides.push_back(boundary::south);
    distanceSides.push_back(boundary::north);
    distanceSides.push_back(boundary::south);
    distanceSides.push_back(boundary::north);

    //vector of indexes of the patches corresponding to distanceSides
    //length of the vector distancePatches must be equal to the length of vector distanceSides
    gsVector<int> distancePatches(5);
    distancePatches << 0, 0, 1, 2, 2;

    int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
    real_t maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

    //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
    //estimation of the wall distance will be computed, if the last input parameter is set as true
    return navStokesSteady.computeDimensionlessWallDistance(distancePatches, distanceSides, viscosity, Re, uMax, maxYplus, numSamplePts, true, true);
}

template<class T>
void solvePoissonEquation(gsMultiPatch<T> patches, uwbTMSolverKOmega<T>& turbSolver, int plot_pts)
{
    int numRefinePoisson = 4;

    gsMultiBasis<> tbasisPoisson(patches); // basis for RANS equations
    for (int i = 0; i < numRefinePoisson; ++i)
        tbasisPoisson.uniformRefine();

    gsFunctionExpr<real_t> fw("1", 2);
    gsFunctionExpr<real_t> gw("0", 2);
    gsFunctionExpr<real_t> wallw("0.0", 2);
    gsBoundaryConditions<real_t> bcInfow;
    bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
    bcInfow.addCondition(0, boundary::east, condition_type::neumann, gw, 0);
    bcInfow.addCondition(1, boundary::east, condition_type::neumann, gw, 0);
    bcInfow.addCondition(0, boundary::west, condition_type::dirichlet, wallw, 0);
    bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
    bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw, 0);
    bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw, 0);
    bcInfow.addCondition(2, boundary::south, condition_type::dirichlet, wallw, 0);

    gsInfo << "\nSolving Poisson equation.\n";

    turbSolver.setPoissonSolution(patches, tbasisPoisson, bcInfow, fw, true, plot_pts);

    gsInfo << "Poisson equation resolved.\n\n";
    turbSolver.plotWallDistance("step_wall_distance", plot_pts);
}

template<class T> void refineBasis_NS(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal)
{
    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;

    for (int i = 0; i < 3; i++)
    {
        basis.refine(0, box);
        basis.refine(1, box);
    }

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    // refinement near wal
    int numRefineLocalWal = 1;

    real_t parArea = 0.2;

    for (int j = 0; j < numRefineLocal; j++)
    {
        box << 0, 0, 0, parArea;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
        }

        box << 0, 0, 1 - parArea, 1;
        for (int i = 0; i < numRefineLocalWal; i++)
        {
            basis.refine(0, box);
            basis.refine(1, box);
            basis.refine(2, box);
        }

        parArea = parArea / 2;
    }
}
