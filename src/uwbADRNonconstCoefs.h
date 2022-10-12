#pragma once

#include "uwbINSSolverSteady.h"
#include "uwbRANSSolver.h"
#include "uwbTMSolverKOmegaLinSteady.h"
#include "uwbTMSolverKOmega.h"

namespace gismo
{

template<class T>
std::vector<gsField<T> > computeADCoeffs(const gsMultiPatch<T>& patches, const gsMultiBasis<T>& tbasis, T viscosity)
{
    int plot_pts = 30000;

    T timeStep = 0.001;

    int numIterSteadyNS = 50;
    int numIterKOmegaSteady = 50;
    int numIterRANS = 200;
    T tolRelNorm = 1e-3;
    int maxRANSPicardIt = 3;
    int maxTMPicardFirstIt = 50;

    //T viscosity = 0.0001;
    T viscositySteady = 0.01;
    T turbIntensity = 0.04;
    T viscosityRatio = 100.;

    //-----------------------------------------------
    std::string tmEvaluator = "koWilcoxLRN";
    /*std::string tmEvaluator = "koWilcoxLRN";
    std::string tmEvaluator = "koSST";
    std::string tmEvaluator = "koSSTMenter2009";
    std::string tmEvaluator = "koSAS";
    std::string tmEvaluator = "koSAS_SS";
    std::string tmEvaluator = "koSAS_SO";
    std::string tmEvaluator = "koSAS_OO";*/
    //-----------------------------------------------

    bool TMsupg = false;
    int tauSUPGType = 1; // 0 'h/(2*deg*norm(u)) * (cotgh(Pe)) - 1/Pe'; 1 'h/(2*deg*norm(u)'; 2 '((2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
                             // 3 '((2/timeStep)^2 + (2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'

    //------------direct solver----------------------
    typedef uwbRANSSolver<T> SolverType;

    gsInfo << "Solving turbulent flow for 2d Lshape.\n";

    gsInfo << "===================================================\n";
    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;
    opt.dirValues = dirichlet::interpolation;
    //opt.dirValues = dirichlet::l2Projection; //default

    gsBoundaryConditions<T> bcInfo;

    std::string str_x = std::to_string(3/sqrt(13));
    std::string str_y = std::to_string(-2/sqrt(13));
    gsFunctionExpr<T> Uin(str_x, str_y, 2);
    gsFunctionExpr<T> Uwall("0", "0", 2);

    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uin, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
    //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);


    std::vector< gsMultiBasis<T> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n\n";

    gsFunctionExpr<T> f("0", "0", 2);

    uwbINSPde<T> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<T> params(NSpde, discreteBases, opt);
    uwbINSSolverSteady<T> navStokes(params); // steady coupled solver

    gsInfo << "Solving Steady NS case: \n";
    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    navStokes.initialize(); // steady solver
    navStokes.solve(numIterSteadyNS, 1e-5);

    gsField<T> velocitySteady = navStokes.constructSolution(0);
    gsField<T> pressureSteady = navStokes.constructSolution(1);

    gsWriteParaview<T>(velocitySteady, "Lshape_steadyVelocity", plot_pts);
    gsWriteParaview<T>(pressureSteady, "Lshape_steadyPressure", plot_pts);

    gsFileData<T> fd;
    fd << navStokes.getSolution();
    fd.save("Lshape_NSsteadySolution.xml");


    T inletWidth = 1.;
    gsVector<T> Uinlet(2);
    Uinlet << 3/sqrt(13), -2/sqrt(13);
    T uFreeStream = Uinlet.norm();
    T Re = uFreeStream * inletWidth / viscosity;
    std::vector<boxSide> distanceSides;
     gsVector<int> distancePatches;

     //vector of the sides of the patches from which the wall distance is computed
     distanceSides.push_back(boundary::north);
     //distanceSides.push_back(boundary::east);
     distanceSides.push_back(boundary::south);
     distanceSides.push_back(boundary::south);
     distanceSides.push_back(boundary::west);
     distanceSides.push_back(boundary::north);
     distanceSides.push_back(boundary::east);
     //vector of indexes of the patches corresponding to distanceSides
     //length of the vector distancePatches must be equal to the length of vector distanceSides
     distancePatches.resize(6);
     distancePatches << 0, 0, 1, 1, 2, 2;

     int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
     T maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

     //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
     //estimation of the wall distance will be computed, if the last input parameter is set as true
     T wallDistance = navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
     gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";

     T maxAspectRatio = navStokes.computeAspectRatio();
     gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";

     // ========================================= Define turbulence solver =========================================
     gsBoundaryConditions<T> bcInfoTurb;

     T kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
     gsInfo << "\nkInConst = " << kInConst << "\n";
     T oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
     gsInfo << "oInConst = " << oInConst << "\n";

     T kWall = 0.;
     T oWall;
     T beta = 0.0708;
     oWall = 6 * viscosity / (beta * math::pow(wallDistance, 2));
     gsInfo << "kWall = " << kWall << "\n";
     gsInfo << "oWall = " << oWall << "\n\n";

     gsFunctionExpr<T> Kin(util::to_string(kInConst), 2);
     gsFunctionExpr<T> Oin(util::to_string(oInConst), 2);
     gsFunctionExpr<T> KWall(util::to_string(kWall), 2);
     gsFunctionExpr<T> OWall(util::to_string(oWall), 2);

     bcInfoTurb.addCondition(2, boundary::west, condition_type::dirichlet, Kin, 0);
     bcInfoTurb.addCondition(2, boundary::north, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(2, boundary::east, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(0, boundary::north, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, KWall, 0);
     //bcInfoTurb.addCondition(0, boundary::east, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, KWall, 0);
     bcInfoTurb.addCondition(2, boundary::west, condition_type::dirichlet, Oin, 1);
     bcInfoTurb.addCondition(2, boundary::north, condition_type::dirichlet, OWall, 1);
     bcInfoTurb.addCondition(2, boundary::east, condition_type::dirichlet, OWall, 1);
     bcInfoTurb.addCondition(0, boundary::north, condition_type::dirichlet, OWall, 1);
     bcInfoTurb.addCondition(0, boundary::south, condition_type::dirichlet, OWall, 1);
     //bcInfoTurb.addCondition(0, boundary::east, condition_type::dirichlet, OWall, 1);
     bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, OWall, 1);
     bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, OWall, 1);

     gsDofMapper koMapper;
     discreteBases[1].getMapper(opt.dirStrategy, opt.intStrategy, bcInfoTurb, koMapper, 0);

     uwbINSPde<T> koPdeSteady(patches, bcInfoTurb, f, viscositySteady);
     uwbINSSolverParams<T> koParamsSteady(koPdeSteady, discreteBases, opt);

     gsMatrix<T> koInitial(koMapper.freeSize(), 2);
     koInitial.col(0).setConstant(kWall);
     koInitial.col(1).setConstant(oWall);

     uwbTMSolverKOmegaLinSteady<T> turbSolver(koParamsSteady);

     //koParamsSteady.settings().set(constantsINS::timeDerTerm, false);
     //koParamsSteady.settings().set(constantsINS::turb_innerFirstIt, numIterKOmegaSteady);
     //koParamsSteady.settings().setTurbulenceEvaluator(tmEvaluator);
     if (TMsupg)
     {
         koParamsSteady.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
         koParamsSteady.settings().set(constantsINS::tauStabType, tauSUPGType);
     }
     //uwbTMSolverKOmega<T> turbSolver(koParamsSteady);

     turbSolver.setInitialCondition(koInitial);

     if (tmEvaluator == "koSST" || tmEvaluator == "koSSTMenter2009" || tmEvaluator == "koSAS"
             || tmEvaluator == "koSAS_SS" || tmEvaluator == "koSAS_SO" || tmEvaluator == "koSAS_OO")
     {
         int numRefinePoisson = 4;

         gsMultiBasis<T> tbasisPoisson(patches); // basis for RANS equations
         for (int i = 0; i < numRefinePoisson; ++i)
             tbasisPoisson.uniformRefine();

         gsFunctionExpr<T> fw("1", 2);
         gsFunctionExpr<T> gw("0", 2);
         gsFunctionExpr<T> wallw("0.0", 2);
         gsBoundaryConditions<T> bcInfow;

         bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
         bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(2, boundary::east, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(0, boundary::north, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
         //bcInfow.addCondition(0, boundary::east, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(1, boundary::west, condition_type::dirichlet, wallw, 0);

         gsInfo << "\nSolving Poisson equation.\n";
         turbSolver.setPoissonSolution(patches, tbasisPoisson, bcInfow, fw, true, plot_pts);

         gsInfo << "Poisson equation resolved.\n\n";
         turbSolver.plotWallDistance("Lshape_wall_distance", plot_pts);
     }

     gsInfo << "\nComputing steady linearized k-omega...\n";
     gsInfo << "initialization...\n";
     gsField<T> velocitySol = navStokes.constructSolution(0);
     gsInfo << "TMnumDofs = " << koMapper.freeSize() << "\n";
     turbSolver.initialize(velocitySol);
     turbSolver.solve(numIterKOmegaSteady, 1e-5); // solution change norm tol = 10^(-5)
     //turbSolver.solve(1, 1e-5); // solution change norm tol = 10^(-5)
     gsField<T> kOmegaSteady = turbSolver.constructSolution();
     gsWriteParaview<T>(kOmegaSteady, "Lshape_steadykOmega", plot_pts, true);
     gsFileData<T> fdTM;
     fdTM << turbSolver.getSolution();
     fdTM.save("Lshape_TMsteadySolution.xml");

     // ========================================= Solving RANS =========================================
     uwbINSPde<T> RANSpde(patches, bcInfo, f, viscosity);
     uwbINSSolverParams<T> RANSparams(RANSpde, discreteBases, opt);
     RANSparams.settings().set(constantsINS::timeStep, timeStep);
     RANSparams.settings().set(constantsINS::unst_innerIt, maxRANSPicardIt);
     //RANSparams.settings().set(constantsINS::unst_innerTol, picardTol);

     uwbINSPde<T> koPde(patches, bcInfoTurb, f, viscosity);
     uwbINSSolverParams<T> koParams(koPde, discreteBases, opt);
     koParams.settings().set(constantsINS::timeStep, timeStep);
     //koParams.settings().set(constantsINS::turb_innerIt, maxTMPicardIt);
     koParams.settings().set(constantsINS::turb_innerFirstIt, maxTMPicardFirstIt);
     koParams.settings().setTurbulenceEvaluator(tmEvaluator);

     if (TMsupg) {
         koParams.settings().set(constantsINS::TMsupg, TMsupg); // set SUPG
         koParams.settings().set(constantsINS::tauStabType, tauSUPGType); //set formula for stabilization parameter tau
     }

     uwbTMSolverKOmega<T> turbSolver_unsteady(koParams);
     turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());//(koInitial);

     if (tmEvaluator == "koSST" || tmEvaluator == "koSSTMenter2009" || tmEvaluator == "koSAS"
             || tmEvaluator == "koSAS_SS" || tmEvaluator == "koSAS_SO" || tmEvaluator == "koSAS_OO")
     {
         int numRefinePoisson = 4;

         gsMultiBasis<T> tbasisPoisson(patches); // basis for RANS equations
         for (int i = 0; i < numRefinePoisson; ++i)
             tbasisPoisson.uniformRefine();

         gsFunctionExpr<T> fw("1", 2);
         gsFunctionExpr<T> gw("0", 2);
         gsFunctionExpr<T> wallw("0.0", 2);
         gsBoundaryConditions<T> bcInfow;

         bcInfow.addCondition(2, boundary::west, condition_type::neumann, gw, 0);
         bcInfow.addCondition(2, boundary::north, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(2, boundary::east, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(0, boundary::north, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(0, boundary::south, condition_type::dirichlet, wallw, 0);
         //bcInfow.addCondition(0, boundary::east, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw, 0);
         bcInfow.addCondition(1, boundary::west, condition_type::dirichlet, wallw, 0);

         gsInfo << "\nSolving Poisson equation.\n";
         turbSolver_unsteady.setPoissonSolution(patches, tbasisPoisson, bcInfow, fw, true, plot_pts);

         gsInfo << "Poisson equation resolved.\n\n";
         turbSolver_unsteady.plotWallDistance("Lshape_wall_distance", plot_pts);
     }

     RANSparams.settings().setTurbulenceSolver(&turbSolver_unsteady);

     //--------
     //uwbRANSSolver<T> ransSolver(RANSparams);
     //uwbRANSSolverIterative<T, uwbGMResRight<T> > ransSolver(RANSparams);
     SolverType ransSolver(RANSparams);
     //--------

     gsStopwatch time;
     T Tassembly;
     T Tsolve;

     gsInfo << "initialization...\n";
     time.restart();
     ransSolver.setInitialCondition(navStokes.getSolution());
     ransSolver.initialize();
     Tassembly = time.stop();

     gsInfo << "numDofs: " << ransSolver.numDofs() << "\n";

     time.restart();

     //=== solving with solution change relative norm as stopping criterion ============
     //ransSolver.solveCoupled(numIterRANS, tolRelNorm);
     ransSolver.solve(numIterRANS, tolRelNorm);
     //ransSolver.solveWithAnimation(numIterRANS, 2, tolRelNorm, plot_pts, true);

     Tsolve = time.stop();

     gsInfo << "Assembly time:" << Tassembly << "\n";
     gsInfo << "Solve time:" << Tsolve << "\n";


     //--- save solution into file ---
     gsFileData<T> fdu, fd_TMu;
     fdu << ransSolver.getSolution();
     fd_TMu << turbSolver_unsteady.getSolution();
     fdu.save("Lshape_RANS_solution.xml");
     fd_TMu.save("Lshape_TM_solution.xml");

     gsInfo << "Plotting in Paraview...\n";

     gsField<T> velocity = ransSolver.constructSolution(0);
     gsField<T> pressure = ransSolver.constructSolution(1);
     gsField<T> kOmega = turbSolver_unsteady.constructSolution();

     gsWriteParaview<T>(velocity, "Lshape_velocity", plot_pts);
     gsWriteParaview<T>(pressure, "Lshape_pressure", plot_pts);
     gsWriteParaview<T>(kOmega, "Lshape_komega", plot_pts);
     ransSolver.plotTurbulentViscosity("Lshape_turb_viscosity", plot_pts);

     gsField<T> kDiffusionCoeff = ransSolver.constructTMCoefficientSol("kDiffusionCoeff");
     gsWriteParaview<T>(kDiffusionCoeff, "Lshape_kDiffusionCoeff", plot_pts);
     gsField<T> oDiffusionCoeff = ransSolver.constructTMCoefficientSol("oDiffusionCoeff");
     gsWriteParaview<T>(oDiffusionCoeff, "Lshape_oDiffusionCoeff", plot_pts);

     std::vector<gsField<T> > advectionDiffusionCoeffs;
     advectionDiffusionCoeffs.push_back(velocity);
     advectionDiffusionCoeffs.push_back(kDiffusionCoeff);

     return advectionDiffusionCoeffs;
}

} //namespace gismo
