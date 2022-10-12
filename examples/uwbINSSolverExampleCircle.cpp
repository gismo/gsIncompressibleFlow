/** @file uwbINSSolversExampleCircle.cpp

Author(s): H. Hornikova K. Michalkova E.Turnerova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>
#include "gsModeling/gsCoonsPatch.h"
#include <../../motor/jku/gsMotorUtils.h>


// solvers
#include "uwbINSSolverSteady.h"
//#include "uwbINSSolverDecoupledInterface.h"
//#include "uwbINSSolverUnsteady.h"

using namespace gismo;

template<class TT> gsMultiPatch<TT>  BSplineCircle2D(TT const & length = 2.2, TT const & width = 0.4, TT const & widthExpand = 0.01, TT const & radius = 0.05, TT const & centreX = 0.0, TT const & centreY = 0.0, bool const & prepatch = true, TT const & prepatchWidth = 0.2);
template<class TT> gsMultiPatch<TT>  BSplineCircle3D(TT const & length = 2.2, TT const & width = 0.4, TT const & widthExpand = 0.01, TT const & radius = 0.05, TT const & centreX = 0.0, TT const & centreY = 0.0, TT const & depth = 0.4 , bool const & prepatch = true, TT const & prepatchWidth = 0.2);
template<class TT> gsMultiPatch<TT>*  BSplineCircle2D4Patches(TT const & length = 2.2, TT const & width = 0.4, TT const & widthExpand = 0.01, TT const & radius = 0.05, TT const & centreX = 0.0, TT const & centreY = 0.0, bool const & prepatch = true, TT const & prepatchWidth = 0.2);
template<class TT> gsMultiPatch<TT>*  BSplineCircle3D4Patches(TT const & length = 2.2, TT const & width = 0.4, TT const & widthExpand = 0.01, TT const & radius = 0.05, TT const & centreX = 0.0, TT const & centreY = 0.0, TT const & depth = 0.4 , bool const & prepatch = true, TT const & prepatchWidth = 0.2);


int main(int argc, char *argv[])
{
    // ========================================= Settings =========================================

    bool plot = false;
    bool outFile = false;
    bool dg = false; // use DGFEM to connect patches
    int numRefine = 2;
    int numRefineLocal = 2;   //local refine 2D - 8patches around circle
    real_t wallRefineKnot = 0.3;    //local refine 2D - 8patches around circle
    int plot_pts = 10000;
    bool prepatch = true;

    real_t viscosity = 0.1;
    char method = 'C'; // C - coupled, P - projection (decoupled solvers)
    decoupled::method decMethod = decoupled::iterative;
    decoupled::projection projVersion = decoupled::proj1;
    real_t alpha_u = 0.5; // relaxation parameter for velocity in projection method
    real_t alpha_p = 1; // relaxation parameter for pressure in projection method
    real_t timeStep = 0.1; // time step for unsteady computation
    //real_t omega = 0; // angular velocity for rotation
    //int numThreads = 10; // number of threads for assembly

    bool SUPG = false;
    int tauSUPGType = 1;


    //command line
    gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addSwitch("file", "Write info into .txt file", outFile);
    cmd.addSwitch("dg", "Use DGFEM to connect patches", dg);
    cmd.addInt("r", "uniformRefine",
        "Number of Uniform h-refinement steps to perform before solving", numRefine);
    cmd.addInt("s", "plotSamples",
        "Number of sample points to use for plotting", plot_pts);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine<0)
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

    gsInfo << "Solving the Circle example.\n";

   // ========================================= Define problem =========================================

    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f("0", "0", 2); // external force
    gsFunctionExpr<> P("0", "0", 2); // pressure boundary condition (decoupled case)
    gsFunctionExpr<> Uin("-7.139  (y^2) + 0.07139 y + 0.2998", "0", 2);
    gsFunctionExpr<> Uwall("0", "0", 2); // relative velocity BCs for inlet/wall
/*
    //----------------prepare for 3D tests------------------

    f = new gsFunctionExpr<>("0", "0", "0", 3);

    // Boundary conditions
    Uin = new gsFunctionExpr<>("-7.139  (y^2) + 0.07139 y + 0.2998", "0", "0", 3);
    Uwall = new gsFunctionExpr<>("0", "0", "0", 3);
    P = new gsFunctionExpr<>("0", "0", "0", 3);

    //------------------------------------------------------
*/
    //----boundary conditions for 4 patches around circle----
/*
    if(prepatch){
        bcInfo.addCondition(5, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(5, boundary::north, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, Uwall, 0);
    }
    else{

        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
    }


    bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(3, boundary::east, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(4, boundary::south, condition_type::dirichlet, Uwall, 0);
    if (method == 'P')
    {
        bcInfo.addCondition(4, boundary::east, condition_type::dirichlet, P, 1);

    }

    //---------3D periodic conditions 4 patches-------
    bcInfo.setIdentityMatrix(3);
    bcInfo.addPeriodic(0, boundary::front, 0, boundary::back, 3);
    bcInfo.addPeriodic(1, boundary::front, 1, boundary::back, 3);
    bcInfo.addPeriodic(2, boundary::front, 2, boundary::back, 3);
    bcInfo.addPeriodic(3, boundary::front, 3, boundary::back, 3);
    bcInfo.addPeriodic(4, boundary::front, 4, boundary::back, 3);
    if(prepatch){
        bcInfo.addPeriodic(5, boundary::front, 5, boundary::back, 3);
    }
*/

    //------------------------------------------------------

    //----boundary conditions for 8 patches around circle----

    if(prepatch){
        bcInfo.addCondition(10, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(10, boundary::north, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(11, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(11, boundary::south, condition_type::dirichlet, Uwall, 0);
    }
    else{

        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uin, 0);
    }


    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Uwall, 0);
     bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(4, boundary::north, condition_type::dirichlet, Uwall, 0);
     bcInfo.addCondition(4, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(6, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(8, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(3, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(5, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(5, boundary::north, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(7, boundary::south, condition_type::dirichlet, Uwall, 0);
    bcInfo.addCondition(9, boundary::south, condition_type::dirichlet, Uwall, 0);


    if (method == 'P')
    {
        bcInfo.addCondition(8, boundary::east, condition_type::dirichlet, P, 1);
        bcInfo.addCondition(9, boundary::east, condition_type::dirichlet, P, 1);

    }

    //------------------------------------------------------

/*
    //---------3D periodic conditions 8 patches-------
    bcInfo.setIdentityMatrix(3);
    bcInfo.addPeriodic(0, boundary::front, 0, boundary::back, 3);
    bcInfo.addPeriodic(1, boundary::front, 1, boundary::back, 3);
    bcInfo.addPeriodic(2, boundary::front, 2, boundary::back, 3);
    bcInfo.addPeriodic(3, boundary::front, 3, boundary::back, 3);
    bcInfo.addPeriodic(4, boundary::front, 4, boundary::back, 3);
    bcInfo.addPeriodic(5, boundary::front, 5, boundary::back, 3);
    bcInfo.addPeriodic(6, boundary::front, 6, boundary::back, 3);
    bcInfo.addPeriodic(7, boundary::front, 7, boundary::back, 3);
    bcInfo.addPeriodic(8, boundary::front, 8, boundary::back, 3);
    bcInfo.addPeriodic(9, boundary::front, 9, boundary::back, 3);
    if(prepatch){
        bcInfo.addPeriodic(10, boundary::front, 10, boundary::back, 3);
        bcInfo.addPeriodic(11, boundary::front, 11, boundary::back, 3);
    }
*/

// ========================================= Define geometry =========================================

    gsMultiPatch<> patches;

    real_t radius = 0.05;
    real_t centreX = 0.0;
    real_t centreY = 0.0;
    real_t width = 0.4;
    real_t widthExpand = 0.01;
    real_t length = 2.2;
    real_t prepatchWidth = width;
    //real_t depth = width;

    patches = BSplineCircle2D<real_t>(length, width, widthExpand, radius, centreX, centreY, prepatch, prepatchWidth);
    //patches = BSplineCircle3D<real_t>(length, width, widthExpand, radius, centreX, centreY,depth, prepatch, prepatchWidth);

    gsInfo << patches << "\n";

    // ========================================= Define basis =========================================

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<> tbasis(patches);
    for (int i = 0; i < numRefine; ++i)
     tbasis.uniformRefine();

      // local refine for 2D case
    gsMatrix<> box_v0(2, 2);
    gsMatrix<> box_v1(2, 2);

    for (size_t k = 0; k < patches.nPatches(); k++){
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
            box_v1 << 0, 0, 1 - wallRefineKnot / math::pow(2, i), 1;
            tbasis.refine(k, box_v0);
            tbasis.refine(k, box_v1);
        }
    }

    //3D refinement
    /*
    gsMatrix<> box_v0(3, 2);
     gsMatrix<> box_v1(3, 2);


       for (int k = 0; k < patches->nPatches(); k++){
           for (int i = 0; i < numRefineLocal; i++)
           {
               box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i), 0, 0;
               box_v1 << 0, 0, 1 - wallRefineKnot / math::pow(2, i), 1,0,0;
               tbasis.refine(k, box_v0);
               tbasis.refine(k, box_v1);
           }
       }
*/
    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space


    // ========================================= Define solver =========================================

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    // --------- other options ---------
    //params.setNumThreads(1); // set the number of threads for assembly

    //if (omega)
    //    params.settings().set(constantsINS::omega, omega); // set rotation (assumed around x-axis)

    if (SUPG) {
        params.settings().set(constantsINS::SUPG, SUPG); // set SUPG
        params.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType); //set formula for stabilization parameter tau
    }
    params.settings().set(constantsINS::timeStep, timeStep);
    params.settings().set(constantsINS::alpha_u, alpha_u);
    params.settings().set(constantsINS::alpha_p, alpha_p);
    params.settings().setDecoupledMethod(decMethod);
    params.settings().setProjVersion(projVersion);
    // ---------------------------------

    // solvers
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    //uwbINSSolverDecoupledInterface<real_t> navStokes(params); // interface for decoupled solvers
    //uwbINSSolverUnsteady<real_t> navStokes(params); // unsteady coupled solver
    // decoupled solvers can be used as unsteady with time step dt = alpha_u / (1 - alpha_u)

    // ========================================= Solving =========================================
    gsInfo << "numDofs: " << navStokes.numDofs() << "\n";

    gsInfo << "initialization...\n";


    navStokes.initialize(); // steady solver


    navStokes.solve(50, 1e-5); // max. 50 iterations, solution change norm tol = 10^(-5)


    real_t Tassembly = navStokes.getAssemblyTime();
    real_t Tsolve = navStokes.getSolverSetupTime();
    real_t Tsetupsolve = navStokes.getSolveTime();


    gsInfo << "Assembly time:" << Tassembly << "\n";
    gsInfo << "Solve time:" << Tsolve << "\n";
    gsInfo << "Solver setup time:" << Tsetupsolve << "\n";

    if (outFile)
    {
        time_t t = std::time(0);   // get time now
        struct tm * now = std::localtime(&t);
        std::stringstream filename;
        filename << "step_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
            << now->tm_hour << "-" << now->tm_min << ".txt";

        std::ofstream ofile;
        ofile.open(filename.str().c_str());
        ofile << "Solving Circle example:\n";
        ofile << "DG = " << dg << "\n";
        ofile << "viscosity = " << viscosity << "\n";
        ofile << "method = " << method << "\n";
        ofile << "inlet velocity = " << Uin << "\n";
        //ofile << "angular velocity = " << omega << "\n";
        /*if (method == 'P')
            ofile << "alpha_u = " << alpha_u << "\n";*/
        ofile << "number of iterations = " << navStokes.getIterationNumber() << "\n";
        ofile << "last relative solution change = " << navStokes.solutionChangeRelNorm() << "\n";
        ofile << "Assembly time:" << Tassembly << "\n";
        ofile << "Solve time:" << Tsolve << "\n";
        ofile << "Solver setup time:" << Tsetupsolve << "\n";
        ofile.close();
    }

    // Optionally plot solution in paraview
    if (plot)
    {
        gsField<> velocity = navStokes.constructSolution(0);
        //gsField<> *relVelocity = navStokes.constructSolution(0, true); // in case of rotation
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
           gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "circle_velocity", plot_pts, true);
        gsWriteParaview<>(pressure, "circle_pressure", plot_pts);


        // Run paraview
        //system("paraview circle_velocity.pvd &");

    }

    //system("pause");

    return 0;
}

//8 patches around circle
template<class TT>
gsMultiPatch<TT>  BSplineCircle2D(const TT & length, TT const & width, TT const & widthExpand, TT const & radius, TT const & centreX, TT const & centreY, bool const &prepatch, TT const & prepatchWidth)
{
    const real_t PI = 3.141592653589793238463;
    unsigned degree = 3;
    real_t alpha = PI/2.0;
    real_t const K = (4.0/3.0)*(tan(alpha/4.0)); //constant for circle
    bool plotMeshes = false;


    gsMultiPatch<TT> mp;

    //initial data for patches
    gsKnotVector<TT> kv_cub(0, 1, 0, 4);
    gsKnotVector<TT> kv_lin(0, 1, 0, 2);
    gsTensorBSplineBasis<2, TT> basis_lc(kv_lin, kv_cub);
    gsTensorBSplineBasis<2, TT> basis_lin(kv_lin, kv_lin);

    gsMatrix<TT> boundary_02(4,2);
    gsMatrix<TT> boundary_24(4,2);
    gsMatrix<TT> boundary_46(4,2);
    boundary_02 << centreX - radius, centreY,
                     centreX - (width/2 - radius)/2, centreY + radius, //2*(centreY + radius)/3 + (width/2 + widthExpand)/3,
            - width/4, (centreY + radius)/2 + (width/2 + widthExpand)/2, //(4*widthExpand + width)/6,
                  - width/4, width/2 + widthExpand;
    boundary_24 << centreX, centreY + radius,
                    2*centreX/3, 2*(centreY + radius)/3 + (width/2 + widthExpand)/3,
                    centreX/3, (centreY + radius)/3 + 2*(width/2 + widthExpand)/3,
                    0.0, width/2 + widthExpand;
    boundary_46 << centreX + radius, centreY,
                    centreX + (width/2 - radius)/2, centreY + radius,
                    width/4, (centreY + radius)/2 + (width/2 + widthExpand)/2,
                    width/4, width/2 + widthExpand;
    gsMatrix<TT> boundary_13(4,2);
    gsMatrix<TT> boundary_35(4,2);
    gsMatrix<TT> boundary_57(4,2);
    boundary_13 << - width/4, -width/2,
                   - width/4, (centreY - radius)/2 - width/4,
                   centreX - (width/2 - radius)/2, centreY - radius,
                   centreX - radius, centreY;

    boundary_35 << 0.0, -width/2,
                    centreX/3, (centreY - radius)/3 + 2*(-width/2)/3,
                     2*centreX/3, 2*(centreY - radius)/3 + (-width/2)/3,
                    centreX, centreY - radius;
    boundary_57 << width/4, -width/2,
                   width/4, (centreY - radius)/2 + (-width/2)/2,
                   centreX + (width/2 - radius)/2, centreY - radius,
                   centreX + radius, centreY;

    //--------------------------------patch 0-------------------------------------------


        gsMatrix<TT> coef_patch0(8, 2);
        coef_patch0 <<  - width/2, widthExpand/2,
                        boundary_02(0,0), boundary_02(0,1),
                        - width/2, (4*widthExpand + width)/6,
                        boundary_02(1,0), boundary_02(1,1),
                        - width/2, (5*widthExpand + 2*width)/6,
                        boundary_02(2,0), boundary_02(2,1),
                        - width/2, width/2 + widthExpand,
                        boundary_02(3,0), boundary_02(3,1);
        gsTensorBSpline<2, TT> patch0(basis_lc, coef_patch0);
        patch0.degreeElevate(degree - 1, 0);


    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(8, 2);
    coef_patch1 <<  - width/2, - width/2,
                     boundary_13(0,0), boundary_13(0,1),
                    - width/2, (widthExpand - 2*width)/6,
            boundary_13(1,0), boundary_13(1,1),
                    -width/2, (2*widthExpand - width)/6,
            boundary_13(2,0), boundary_13(2,1),
                    - width/2, widthExpand/2,
            boundary_13(3,0), boundary_13(3,1);

    gsTensorBSpline<2, TT> patch1(basis_lc, coef_patch1);
    patch1.degreeElevate(degree - 1, 0);

    //--------------------------------patch 2-------------------------------------------

    gsMultiPatch<TT> boundaries2;
    gsMatrix<TT> boundary_2s(4,2);
    gsMatrix<TT> boundary_2n(4,2);

    boundary_2s << centreX - radius, centreY,
            centreX - radius, centreY + K*radius,
            centreX - K*radius, centreY + radius,
            centreX, centreY + radius;
    boundary_2n << -width/4, width/2 + widthExpand,
            -width/6, width/2 + widthExpand,
            -width/12, width/2 + widthExpand,
            0.0, width/2 + widthExpand;


    boundaries2.addPatch(gsBSpline<>( kv_cub, boundary_02));
    boundaries2.addPatch(gsBSpline<>( kv_cub, boundary_2n));
    boundaries2.addPatch(gsBSpline<>( kv_cub, boundary_24));
    boundaries2.addPatch(gsBSpline<>( kv_cub, boundary_2s));

    gsCoonsPatch<real_t> patch2 = coonsPatch(boundaries2);
    patch2.compute();

    //--------------------------------patch 3-------------------------------------------

    gsMultiPatch<TT> boundaries3;
    gsMatrix<TT> boundary_3s(4,2);
    gsMatrix<TT> boundary_3n(4,2);

    boundary_3n << centreX - radius, centreY,
            centreX - radius, centreY - K*radius,
            centreX - K*radius, centreY - radius,
            centreX, centreY - radius;
    boundary_3s << -width/4, -width/2,
            -width/6, -width/2,
            -width/12, -width/2,
            0.0, -width/2;


    boundaries3.addPatch(gsBSpline<>( kv_cub, boundary_13));
    boundaries3.addPatch(gsBSpline<>( kv_cub, boundary_3n));
    boundaries3.addPatch(gsBSpline<>( kv_cub, boundary_35));
    boundaries3.addPatch(gsBSpline<>( kv_cub, boundary_3s));

    gsCoonsPatch<real_t> patch3 = coonsPatch(boundaries3);
    patch3.compute();


    //--------------------------------patch 4-------------------------------------------
    gsMultiPatch<TT> boundaries4;
    gsMatrix<TT> boundary_4s(4,2);
    gsMatrix<TT> boundary_4n(4,2);

    boundary_4s << centreX, centreY + radius,
            centreX + K*radius, centreY + radius,
            centreX + radius, centreY + K*radius,
            centreX + radius, centreY;
    boundary_4n << 0.0, width/2 + widthExpand,
            width/12, width/2 + widthExpand,
            width/6, width/2 + widthExpand,
            width/4, width/2 + widthExpand;


    boundaries4.addPatch(gsBSpline<>( kv_cub, boundary_24));
    boundaries4.addPatch(gsBSpline<>( kv_cub, boundary_4n));
    boundaries4.addPatch(gsBSpline<>( kv_cub, boundary_46));
    boundaries4.addPatch(gsBSpline<>( kv_cub, boundary_4s));

    gsCoonsPatch<real_t> patch4 = coonsPatch(boundaries4);
    patch4.compute();

    //--------------------------------patch 5-------------------------------------------
    gsMultiPatch<TT> boundaries5;
    gsMatrix<TT> boundary_5s(4,2);
    gsMatrix<TT> boundary_5n(4,2);

    boundary_5n << centreX, centreY - radius,
            centreX + K*radius, centreY - radius,
            centreX + radius, centreY - K* radius,
            centreX + radius, centreY;
    boundary_5s << 0.0, -width/2,
            width/12, -width/2,
            width/6, -width/2,
            width/4, -width/2;


    boundaries5.addPatch(gsBSpline<>( kv_cub, boundary_35));
    boundaries5.addPatch(gsBSpline<>( kv_cub, boundary_5n));
    boundaries5.addPatch(gsBSpline<>( kv_cub, boundary_57));
    boundaries5.addPatch(gsBSpline<>( kv_cub, boundary_5s));

    gsCoonsPatch<real_t> patch5 = coonsPatch(boundaries5);
    patch5.compute();


    //--------------------------------patch 6-------------------------------------------
    gsMatrix<TT> coef_patch6(8, 2);
    coef_patch6 << boundary_46(0,0), boundary_46(0,1),
                    width/2, centreY + widthExpand/2,
                     boundary_46(1,0), boundary_46(1,1),
                    width/2, (4*widthExpand + width)/6,
                   boundary_46(2,0), boundary_46(2,1),
                    width/2, (5*widthExpand + 2*width)/6,
                    boundary_46(3,0), boundary_46(3,1),
                    width/2, width/2 + widthExpand;

    gsTensorBSpline<2, TT> patch6(basis_lc, coef_patch6);
    patch6.degreeElevate(degree - 1, 0);


    //--------------------------------patch 7-------------------------------------------
    gsMatrix<TT> coef_patch7(8, 2);
    coef_patch7 <<  boundary_57(0,0), boundary_57(0,1),
                    + width/2, - width/2,
                    boundary_57(1,0), boundary_57(1,1),
                    + width/2, (widthExpand - 2*width)/6,
                    boundary_57(2,0), boundary_57(2,1),
                    + width/2, (2*widthExpand - width)/6,
                    boundary_57(3,0), boundary_57(3,1),
                    + width/2, widthExpand/2;

    gsTensorBSpline<2, TT> patch7(basis_lc, coef_patch7);
    patch7.degreeElevate(degree - 1, 0);


    //--------------------------------patch 8-------------------------------------------
    gsMatrix<TT> coef_patch8(4, 2);
    coef_patch8 <<  width/2, widthExpand/2,
                    length - width/2, widthExpand/2,
                    width/2, width/2 + widthExpand,
                    length - width/2, width/2 + widthExpand;
    gsTensorBSpline<2, TT> patch8(basis_lin, coef_patch8);
    patch8.degreeElevate(degree - 1, -1);



    //--------------------------------patch 9-------------------------------------------
    gsMatrix<TT> coef_patch9(4, 2);
    coef_patch9 <<  width/2, -width/2,
                    length - width/2, -width/2,
                    width/2, widthExpand/2,
                    length - width/2, widthExpand/2;
    gsTensorBSpline<2, TT> patch9(basis_lin, coef_patch9);
    patch9.degreeElevate(degree - 1, -1);

    //refining last patches according to length/width
    real_t insertKnotNum = trunc(length/width)-1;
    if(insertKnotNum > 0.00001){
        for(real_t i = 1/(3*insertKnotNum); i < 1.0; i+=1/(3*insertKnotNum)){
        patch8.insertKnot(i,0,1);
        patch9.insertKnot(i,0,1);
        }
    }

    mp.addPatch(patch0);
    mp.addPatch(patch1);
    mp.addPatch(patch2.result());
    mp.addPatch(patch3.result());
    mp.addPatch(patch4.result());
    mp.addPatch(patch5.result());
    mp.addPatch(patch6);
    mp.addPatch(patch7);
    mp.addPatch(patch8);
    mp.addPatch(patch9);

    mp.addInterface(0, boundary::south, 1, boundary::north);
    mp.addInterface(0, boundary::east, 2, boundary::west);
    mp.addInterface(1, boundary::east, 3, boundary::west);
    mp.addInterface(2, boundary::east, 4, boundary::west);
    mp.addInterface(3, boundary::east, 5, boundary::west);
    mp.addInterface(4, boundary::east, 6, boundary::west);
    mp.addInterface(5, boundary::east, 7, boundary::west);
    mp.addInterface(6, boundary::east, 8, boundary::west);
    mp.addInterface(6, boundary::south, 7, boundary::north);
    mp.addInterface(7, boundary::east, 9, boundary::west);
    mp.addInterface(8, boundary::south, 9, boundary::north);


    //add start patch
    if(prepatch){
        //--------------------------------patch 10-------------------------------------------
        gsMatrix<TT> coef_patch10(4, 2);
        coef_patch10 << - prepatchWidth - width/2, widthExpand/2,
                - width/2, widthExpand/2,
                - prepatchWidth - width/2, width/2 + widthExpand,
                 - width/2, width/2 + widthExpand;

        gsTensorBSpline<2, TT> patch10(basis_lin, coef_patch10);
        patch10.degreeElevate(degree - 1, -1);



        //--------------------------------patch 11-------------------------------------------
        gsMatrix<TT> coef_patch11(4, 2);
        coef_patch11 << - prepatchWidth - width/2, -width/2,
                - width/2,-width/2,
                - prepatchWidth - width/2, widthExpand/2,
                 - width/2, widthExpand/2;

        gsTensorBSpline<2, TT> patch11(basis_lin, coef_patch11);
        patch11.degreeElevate(degree - 1, -1);


        real_t insertKnotNumPrepatch = trunc(prepatchWidth/width);

        if(insertKnotNumPrepatch > 0.00001){
            for(real_t i = 1/(3*insertKnotNumPrepatch); i < 1.0; i+=1/(3*insertKnotNumPrepatch)){

            patch10.insertKnot(i,0,1);
            patch11.insertKnot(i,0,1);
            }
        }

        mp.addPatch(patch10);
        mp.addPatch(patch11);

        mp.addInterface(10, boundary::east, 0, boundary::west);
        mp.addInterface(11, boundary::east, 1, boundary::west);
        mp.addInterface(10, boundary::south, 11, boundary::north);

        if (plotMeshes){
            gsMesh<> mesh10;
            patch10.controlNet(mesh10);
            gsWriteParaview(mesh10,"cpCirclePatch10");

            gsMesh<> mesh11;
            patch11.controlNet(mesh11);
            gsWriteParaview(mesh11,"cpCirclePatch11");
        }

    }


    mp.addAutoBoundaries();

    if (plotMeshes) {
        gsMesh<> mesh0;
        patch0.controlNet(mesh0);

        gsWriteParaview(mesh0,"cpCirclePatch0");

        gsMesh<> mesh1;
        patch1.controlNet(mesh1);

        gsWriteParaview(mesh1,"cpCirclePatch1");

        gsMesh<> mesh2;
        (patch2.result()).controlNet(mesh2);

        gsWriteParaview(mesh2,"cpCirclePatch2");

        gsMesh<> mesh3;
        (patch3.result()).controlNet(mesh3);

        gsWriteParaview(mesh3,"cpCirclePatch3");

        gsMesh<> mesh4;
        (patch4.result()).controlNet(mesh4);
        gsWriteParaview(mesh4,"cpCirclePatch4");

        gsMesh<> mesh5;
        (patch5.result()).controlNet(mesh5);
        gsWriteParaview(mesh5,"cpCirclePatch5");

        gsMesh<> mesh6;
        patch6.controlNet(mesh6);
        gsWriteParaview(mesh6,"cpCirclePatch6");

        gsMesh<> mesh7;
        patch7.controlNet(mesh7);
        gsWriteParaview(mesh7,"cpCirclePatch7");

        gsMesh<> mesh8;
        patch8.controlNet(mesh8);
        gsWriteParaview(mesh8,"cpCirclePatch8");

        gsMesh<> mesh9;
        patch9.controlNet(mesh9);
        gsWriteParaview(mesh9,"cpCirclePatch9");


   }

    return mp;
}

//8 patches around circle
template<class TT>

gsMultiPatch<TT>  BSplineCircle3D(const TT & length, TT const & width, TT const & widthExpand, TT const & radius, TT const & centreX, TT const & centreY, const TT &depth, bool const &prepatch, TT const & prepatchWidth)
{

    unsigned degree = 3;
    bool plotMeshes = false;



    gsMultiPatch<TT> patches2D;
    patches2D = BSplineCircle2D<real_t>(length, width, widthExpand, radius, centreX, centreY, prepatch, prepatchWidth);

    gsMultiPatch<TT> mp;

    for(unsigned n=0; n<patches2D.nPatches();n++) {
        gsMatrix<TT> coef_patch(2*((patches2D.patch(n)).coefsSize()), 3);
        for(unsigned i=0;i<2*((patches2D.patch(n)).coefsSize());i++)
        {
            if (i < ((patches2D.patch(n)).coefsSize())) {
                coef_patch(i,0) = (patches2D.patch(n)).coef(i,0);
                coef_patch(i,1) = (patches2D.patch(n)).coef(i,1);
                coef_patch(i,2) = -depth/2;
            } else {
                coef_patch(i,0) = (patches2D.patch(n)).coef(i - (patches2D.patch(n)).coefsSize(),0);
                coef_patch(i,1) = (patches2D.patch(n)).coef(i - (patches2D.patch(n)).coefsSize(),1);
                coef_patch(i,2) = depth/2;
            }
        }


        const gsTensorBSplineBasis<2, real_t>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(patches2D->basis(n)));

           gsKnotVector<TT> kv(0,1,0,2);
           gsKnotVector<TT> kv1;
           kv1 =basis->knots(0);
           gsKnotVector<TT> kv2;
           kv2=basis->knots(1);


        gsTensorBSplineBasis<3, TT> patch_basis(kv1,kv2,kv);

        gsTensorBSpline<3, TT> patch(patch_basis, coef_patch);
        patch.degreeElevate(degree - 1, 2);

        mp.addPatch(patch);

   }

    mp.addInterface(0, boundary::south, 1, boundary::north);
    mp.addInterface(0, boundary::east, 2, boundary::west);
    mp.addInterface(1, boundary::east, 3, boundary::west);
    mp.addInterface(2, boundary::east, 4, boundary::west);
    mp.addInterface(3, boundary::east, 5, boundary::west);
    mp.addInterface(4, boundary::east, 6, boundary::west);
    mp.addInterface(5, boundary::east, 7, boundary::west);
    mp.addInterface(6, boundary::east, 8, boundary::west);
    mp.addInterface(6, boundary::south, 7, boundary::north);
    mp.addInterface(7, boundary::east, 9, boundary::west);
    mp.addInterface(8, boundary::south, 9, boundary::north);


    //add start patch
    if(prepatch){
        mp.addInterface(10, boundary::east, 0, boundary::west);
        mp.addInterface(11, boundary::east, 1, boundary::west);
        mp.addInterface(10, boundary::south, 11, boundary::north);

        if (plotMeshes){
            gsMesh<> mesh10;
            (mp.patch(10)).controlNet(mesh10);
            gsWriteParaview(mesh10,"cpCirclePatch10");

            gsMesh<> mesh11;
            (mp.patch(11)).controlNet(mesh11);
            gsWriteParaview(mesh11,"cpCirclePatch11");
        }

    }


    mp.addAutoBoundaries();

    if (plotMeshes) {
        gsMesh<> mesh0;
        (mp.patch(0)).controlNet(mesh0);
        gsWriteParaview(mesh0,"cpCirclePatch0");

        gsMesh<> mesh1;
        (mp.patch(1)).controlNet(mesh1);
        gsWriteParaview(mesh1,"cpCirclePatch1");

        gsMesh<> mesh2;
       (mp.patch(2)).controlNet(mesh2);
        gsWriteParaview(mesh2,"cpCirclePatch2");

        gsMesh<> mesh3;
       (mp.patch(3)).controlNet(mesh3);
        gsWriteParaview(mesh3,"cpCirclePatch3");

        gsMesh<> mesh4;
        (mp.patch(4)).controlNet(mesh4);
        gsWriteParaview(mesh4,"cpCirclePatch4");

        gsMesh<> mesh5;
        (mp.patch(5)).controlNet(mesh5);
        gsWriteParaview(mesh5,"cpCirclePatch5");

        gsMesh<> mesh6;
        (mp.patch(6)).controlNet(mesh6);
        gsWriteParaview(mesh6,"cpCirclePatch6");

        gsMesh<> mesh7;
        (mp.patch(7)).controlNet(mesh7);
        gsWriteParaview(mesh7,"cpCirclePatch7");

        gsMesh<> mesh8;
        (mp.patch(8)).controlNet(mesh8);
        gsWriteParaview(mesh8,"cpCirclePatch8");

        gsMesh<> mesh9;
        (mp.patch(9)).controlNet(mesh9);
        gsWriteParaview(mesh9,"cpCirclePatch9");


   }

    return mp;
}


//4 patches around circle
template<class TT>
gsMultiPatch<TT>  *BSplineCircle2D4Patches(const TT & length, TT const & width, TT const & widthExpand, TT const & radius, TT const & centreX, TT const & centreY, bool const &prepatch, TT const & prepatchWidth)
{
    const real_t PI = 3.141592653589793238463;
    unsigned degree = 3;
    real_t alpha = PI/2.0;
    real_t const K = (4.0/3.0)*(tan(alpha/4.0)); //constant for circle
    bool plotMeshes = false;


    gsMultiPatch<TT> * mp = new gsMultiPatch<TT>;




    //initial data for patches
    gsKnotVector<TT> kv_cub(0, 1, 0, 4);
    gsKnotVector<TT> kv_lin(0, 1, 0, 2);
    gsTensorBSplineBasis<2, TT> basis(kv_lin, kv_cub);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(8, 2);
    coef_patch0 <<  - width/2, - width/2,
                    - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                     - width/2, - (width + widthExpand)/4,
                    - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius * sin(alpha/2)),
                      - width/2, (width + widthExpand)/4,
                    - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius *  sin(alpha/2)),
                    - width/2,  width/2 + widthExpand,
                    - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    patch0.degreeElevate(degree - 1, 0);


    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(8, 2);
    coef_patch1 <<  - width/2, width/2 + widthExpand,
            - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY,
                    - width/4, width/2 + widthExpand,
            - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                     width/4, width/2 + widthExpand,
                        radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                     + width/2, width/2 + widthExpand,

                                radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch1(basis, coef_patch1);
    patch1.degreeElevate(degree - 1, 0);

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2(8, 2);
    coef_patch2 <<  - width/2,- width/2,
                   - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
             - width/4, - width/2,
            - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
             width/4, - width/2,
            radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
              + width/2,- width/2,
            radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY;

    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    patch2.degreeElevate(degree - 1, 0);

    //--------------------------------patch 3-------------------------------------------
    gsMatrix<TT> coef_patch3(8, 2);
    coef_patch3 << width/2, - width/2,
                   radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                      width/2, -(width + widthExpand)/4,
                    radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius*sin(alpha/2)),
                   width/2,  (width + widthExpand)/4,
                    radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius * sin(alpha/2)),
                    width/2,  (width/2 + widthExpand),
                     radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
    gsTensorBSpline<2, TT> patch3(basis, coef_patch3);
    patch3.degreeElevate(degree - 1, 0);

    //--------------------------------patch 4-------------------------------------------
    gsMatrix<TT> coef_patch4(8, 2);
    coef_patch4 << length - width/2, - width/2,
                   width/2, - width/2,
                      length - width/2, - (width + widthExpand)/4,
                    width/2, - (width + widthExpand)/4,
                   length - width/2,  (width + widthExpand)/4,
                    width/2, + (width + widthExpand)/4,
                    length - width/2,  (width/2 + widthExpand),
                   width/2,  (width/2 + widthExpand);
    gsTensorBSpline<2, TT> patch4(basis, coef_patch4);
    patch4.degreeElevate(degree - 1, 0);

    //refining last patch according to length/width
    real_t insertKnotNum = trunc(length/width)-1;
    if(insertKnotNum>1.0){
        for(real_t i = 1/(2*insertKnotNum); i < 1.0; i+=1/(2*insertKnotNum)){
        patch4.insertKnot(i,0,1);
        }
    }


    mp->addPatch(patch0);
    mp->addPatch(patch1);
    mp->addPatch(patch2);
    mp->addPatch(patch3);
    mp->addPatch(patch4);

    mp->addInterface(0, boundary::north, 1, boundary::south);
    mp->addInterface(0, boundary::south, 2, boundary::south);
    mp->addInterface(3, boundary::north, 1, boundary::north);
    mp->addInterface(3, boundary::south, 2, boundary::north);
    mp->addInterface(4, boundary::east, 3, boundary::west);

    //add start patch
    if(prepatch){
        //--------------------------------patch 5-------------------------------------------
        gsMatrix<TT> coef_patch5(8, 2);
        coef_patch5 << -prepatchWidth - width/2, - width/2,
                - width/2, - width/2,
                 - prepatchWidth - width/2, - (width + widthExpand)/4,
                - width/2 , - (width + widthExpand)/4,
                  - prepatchWidth - width/2, (width + widthExpand)/4,
                - width/2 , (width + widthExpand)/4,
                - prepatchWidth - width/2,  width/2 + widthExpand,
                - width/2,  width/2 + widthExpand;
        gsTensorBSpline<2, TT> patch5(basis, coef_patch5);
        patch5.degreeElevate(degree - 1, 0);

        real_t insertKnotNumPrepatch = trunc(prepatchWidth/width);

        if(insertKnotNumPrepatch > 0.00001){
            for(real_t i = 1/(2*insertKnotNumPrepatch); i < 1.0; i+=1/(2*insertKnotNumPrepatch)){

            patch5.insertKnot(i,0,1);
            }
        }


        mp->addPatch(patch5);
        mp->addInterface(0, boundary::west, 5, boundary::east);

        if (plotMeshes){
            gsMesh<> mesh5;
            patch5.controlNet(mesh5);
            gsWriteParaview(mesh5,"cpCirclePatch5");
        }
    }




    mp->addAutoBoundaries();

    if (plotMeshes) {
        gsMesh<> mesh0;
        patch0.controlNet(mesh0);
        gsWriteParaview(mesh0,"cpCirclePatch0");

        gsMesh<> mesh1;
        patch1.controlNet(mesh1);
        gsWriteParaview(mesh1,"cpCirclePatch1");

        gsMesh<> mesh2;
        patch2.controlNet(mesh2);
        gsWriteParaview(mesh2,"cpCirclePatch2");

        gsMesh<> mesh3;
        patch3.controlNet(mesh3);
        gsWriteParaview(mesh3,"cpCirclePatch3");

        gsMesh<> mesh4;
        patch4.controlNet(mesh4);
        gsWriteParaview(mesh4,"cpCirclePatch4");


   }

    return mp;
}



//4 patches around circle
template<class TT>
gsMultiPatch<TT>  *BSplineCircle3D4Patches(const TT & length, TT const & width, TT const & widthExpand, TT const & radius, TT const & centreX, TT const & centreY, const TT &depth, bool const &prepatch, const TT &prepatchWidth)
{


    unsigned degree = 3;
    bool plotMeshes = false;



    gsMultiPatch<TT> *patches2D = new gsMultiPatch<TT>;
    patches2D = BSplineCircle2D4Patches<real_t>(length, width, widthExpand, radius, centreX, centreY, prepatch, prepatchWidth);

    gsMultiPatch<TT> * mp = new gsMultiPatch<TT>;

    for(unsigned n=0; n<patches2D->nPatches();n++) {
        gsMatrix<TT> coef_patch(2*((patches2D->patch(n)).coefsSize()), 3);
        for(unsigned i=0;i<2*((patches2D->patch(n)).coefsSize());i++)
        {
            if (i < ((patches2D->patch(n)).coefsSize())) {
                coef_patch(i,0) = (patches2D->patch(n)).coef(i,0);
                coef_patch(i,1) = (patches2D->patch(n)).coef(i,1);
                coef_patch(i,2) = -depth/2;
            } else {
                coef_patch(i,0) = (patches2D->patch(n)).coef(i - (patches2D->patch(n)).coefsSize(),0);
                coef_patch(i,1) = (patches2D->patch(n)).coef(i - (patches2D->patch(n)).coefsSize(),1);
                coef_patch(i,2) = depth/2;
            }
        }


        const gsTensorBSplineBasis<2, real_t>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(patches2D->basis(n)));

           gsKnotVector<TT> kv(0,1,0,2);
           gsKnotVector<TT> kv1;
           kv1 =basis->knots(0);
           gsKnotVector<TT> kv2;
           kv2=basis->knots(1);

           gsTensorBSplineBasis<3, TT> patch_basis(kv1,kv2,kv);
        gsTensorBSpline<3, TT> patch(patch_basis, coef_patch);
        patch.degreeElevate(degree - 1, 2);

        mp->addPatch(patch);

   }

    delete patches2D;



    mp->addInterface(0, boundary::north, 1, boundary::south);
    mp->addInterface(0, boundary::south, 2, boundary::south);
    mp->addInterface(3, boundary::north, 1, boundary::north);
    mp->addInterface(3, boundary::south, 2, boundary::north);
    mp->addInterface(4, boundary::east, 3, boundary::west);

    //add start patch
    if(prepatch){

       mp->addInterface(0, boundary::west, 5, boundary::east);

        if (plotMeshes){
            gsMesh<> mesh5;
             (mp->patch(5)).controlNet(mesh5);
            gsWriteParaview(mesh5,"cpCirclePatch5");
        }
    }


    mp->addAutoBoundaries();

    if (plotMeshes) {
        gsMesh<> mesh0;
         (mp->patch(0)).controlNet(mesh0);
        gsWriteParaview(mesh0,"cpCirclePatch0");

        gsMesh<> mesh1;
        (mp->patch(1)).controlNet(mesh1);
        gsWriteParaview(mesh1,"cpCirclePatch1");

        gsMesh<> mesh2;
         (mp->patch(2)).controlNet(mesh2);
        gsWriteParaview(mesh2,"cpCirclePatch2");

        gsMesh<> mesh3;
         (mp->patch(3)).controlNet(mesh3);
        gsWriteParaview(mesh3,"cpCirclePatch3");

        gsMesh<> mesh4;
         (mp->patch(4)).controlNet(mesh4);
        gsWriteParaview(mesh4,"cpCirclePatch4");


   }

    return mp;
}

