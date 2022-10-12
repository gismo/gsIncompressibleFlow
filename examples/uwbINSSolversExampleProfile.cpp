/** @file uwbINSSolversExampleProfile.cpp

Author(s): K. Michalkova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include <gismo.h>

// solvers

#include <gismo.h>
#include "gsModeling/gsCoonsPatch.h"
#include <../../motor/jku/gsMotorUtils.h>


#include "uwbINSSolverSteady.h"
//#include "uwbINSSolverDecoupledInterface.h"
#include "uwbINSSolverUnsteady.h"

#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"

#include <sstream>

#include <gsOptimizer/gsQualityMeasure.h>

#include <math.h>
#include <string.h>
#include <sstream>
#include "../jku/gsMotorUtils.h"
//#include "../jku/gsQualityMeasure2.h"

using namespace gismo;

template<class TT>
gsMultiPatch<TT> BSplineProfile2DBetweenPatch(TT const & index,
                                       TT const & length_x1,
                                       TT const & length_x2,
                                       TT const & pitch,
                                       TT const & camberX,
                                       TT const & camberY,
                                       TT const & leadingAngle,
                                       TT const & trailingAngle,
                                       TT const & thicknessX,
                                       TT const & thicknessY,
                                       TT const & endingOffset,
                                       TT const & outputAngle,
                                       TT const & radius,
                                       TT const & chordLength,
                                       TT const & Angle,
                                       TT const & rotationCenterX,
                                       TT const & rotationCenterY);

template<class TT>  // should be rewritten to gsMultiPatch<TT> instead of pointer
gsMultiPatch<TT>* BSplineProfile2DPatch(TT const & index, TT const & length_x1, TT const & length_x2, TT const & pitch, TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
                                                     TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY);


//change of continuity between patch 0 and patch 1
//template<class TT>
//gsMultiPatch<TT>* BSplineProfile2DPatchC1(TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
//                                              TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY);

template<class TT>
gsMultiPatch<TT>* BSplineProfile3DPatch(TT const & index, TT const & length_x1, TT const & length_x2,TT const & pitch, TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
                                                     TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY);

int main(int argc, char *argv[])
{
    const real_t PI = 3.141592653589793238463;
    
    // ========================================= Settings ========================================= 
   
    bool plot = false;
    bool plotGeometry = false;
    bool outFile = false;
    bool dg = false; // use DGFEM to connect patches
    bool relative = true; //relative velocity is used
    int numRefine = 2;
    real_t wallRefineKnot = 0.1;
    real_t numRefineLocal = 3;
    real_t uRefineKnot = 0.1;
    int plot_pts = 10000;

    int numIter = 20;
    real_t viscosity = 0.001;
    char method = 'C'; // C - coupled, P - projection (decoupled solvers)
    decoupled::method decMethod = decoupled::iterative;
    decoupled::projection projVersion = decoupled::proj2;
    real_t alpha_u = 0.5; // relaxation parameter for velocity in projection method
    real_t alpha_p = 1; // relaxation parameter for velocity in projection method
    real_t timeStep = 0.005; // time step for unsteady computation
    //real_t omega = 0; // angular velocity for rotation
    //int numThreads = 10; // number of threads for assembly

	bool SUPG = false;
    int tauSUPGType = 1;


    // ========================================= Define geometry =========================================

    /*
    //------parametry definice------
    //parametry, ktere lze menit
    real_t camberX;    //x-ova souradnice max polohy strednice
    real_t camberY;    //y-ova souradnice max polohy strednice
    real_t leadingAngle;      //vstupni uhel strednice
    real_t trailingAngle;        //vystupni uhel strednice
    real_t thicknessX;         //x-ova souradnice max prohnuti tloustky
    real_t thicknessY;            //y-ova souranice max prohnuti tloustky
    real_t outputAngle;         //vystupni uhel tloustky
    real_t radius;            //polomer osk. kruz v nabehu

    //parametry, ktere "nelze" menit -> potreba upravit soubor uwbProfileOptimization.h
            real_t endingOffset;     //odsazeni tloustky na konci
            real_t chordLength;       //delka lopatky -> v pripade zmeny potreba upravit ohranicujici kanal
            real_t angle;   //otoceni do mrize
            real_t rotationCenterX;   //zde neni potreba, slouzi pro otoceni profilu do lop. mrize
            real_t rotationCenterY;
    */


    unsigned num_blades = 4;
        unsigned num_bladeprofiles = 7;

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
     gsVector<> angle_input(num_bladeprofiles);

    // SETTING PARAMETERS
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

//    gsVector<> length_x1(num_bladeprofiles);
//    length_x1 << 0.0, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
//    real_t length_x2 = 0.433544;

    real_t length_x1_fixed = -0.8;
    real_t length_x2_fixed = 2;
//    real_t rotate_angle = 0.0 * PI/180;


    //=======================Input velocities==================================================================
      real_t omega = 0.;//2*PI*538.0/60.0;   // value 538 cycles per minute is given

      gsVector<> velocity_blade(num_bladeprofiles);    //angular velocity
      for(unsigned i=0; i<num_bladeprofiles; i++){
          velocity_blade(i) = -rr(i)*omega;
      }

      real_t flow_rate = 5.73;
      real_t D_out = 0.5*2;
      real_t d_out = 0.122926*2;
      real_t vel_mer_out = flow_rate/(PI*(D_out*D_out-d_out*d_out)/4.0);
      real_t vel_mer_in = vel_mer_out;                //meridial velocity

      gsVector<> velocity_absolute_x(num_bladeprofiles);
      gsVector<> velocity_absolute_y(num_bladeprofiles);
      gsVector<> velocity_relative_x(num_bladeprofiles);
      gsVector<> velocity_relative_y(num_bladeprofiles);

      angle_input = angle-leading_angle;
    //   gsInfo<< "\n angle between input velocity and x-axis (degrees) \n";
    //  gsInfo << angle_input*180/PI;
      for(unsigned i=0; i<num_bladeprofiles; i++){
          velocity_relative_x(i) = vel_mer_in;
          velocity_relative_y(i) = vel_mer_in*math::tan(angle_input(i));
      }

      for(unsigned i=0; i<num_bladeprofiles; i++){
         velocity_absolute_x(i) = velocity_relative_x(i);
         velocity_absolute_y(i) = velocity_relative_y(i)+velocity_blade(i);
      }

  /*
          gsInfo<< "\n velocity_blade\n";
          gsInfo<< velocity_blade;
          gsInfo<< "\n";

          gsInfo<< "velocity_absolute_x\n";
          gsInfo<< velocity_absolute_x;
          gsInfo<< "\n";

          gsInfo<< "velocity_absolute_y\n";
          gsInfo<< velocity_absolute_y;
          gsInfo<< "\n";

          gsInfo<< "velocity_relative_x\n";
          gsInfo<< velocity_relative_x;
          gsInfo<< "\n";

          gsInfo<< "velocity_relative_y\n";
          gsInfo<< velocity_relative_y;
          gsInfo<< "\n";

  */

    //for(unsigned i_blade=1; i_blade<7; i_blade++){   //to compute profiles 1-6
        gsInfo << "===================================================\n";

        unsigned index_of_profile = 3; //  1-6, do not use 0

        std::ostringstream strs_profile;
        std::string profile;

            strs_profile << index_of_profile;
            if(relative == true){
            profile = "_relative_blade" + strs_profile.str();}
            else{
            profile = "_absolute_blade" + strs_profile.str();
            }

             gsInfo << "profil" + profile;
             gsInfo << "\n";




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

    gsInfo << "Solving the  bladeProfile example.\n";

    // ========================================= Compute geometry =========================================


    gsMultiPatch<real_t> patches;
//    patches = *BSplineProfile2DPatch<real_t>(index_of_profile,length_x1_fixed,length_x2_fixed,((2*PI*rr[index_of_profile])/num_blades), camber_x[index_of_profile], camber_y[index_of_profile],  leading_angle[index_of_profile],  trailing_angle[index_of_profile],  thickness_x[index_of_profile],  thickness_y[index_of_profile],  ending_offset[index_of_profile],  output_angle[index_of_profile],  radius[index_of_profile],  chord_length[index_of_profile],
//                                             angle[index_of_profile], rotation_center_x[index_of_profile],  rotation_center_y[index_of_profile]);

    patches = BSplineProfile2DBetweenPatch<real_t>(index_of_profile,length_x1_fixed,length_x2_fixed,((2*PI*rr[index_of_profile])/num_blades), camber_x[index_of_profile], camber_y[index_of_profile],  leading_angle[index_of_profile],  trailing_angle[index_of_profile],  thickness_x[index_of_profile],  thickness_y[index_of_profile],  ending_offset[index_of_profile],  output_angle[index_of_profile],  radius[index_of_profile],  chord_length[index_of_profile],
                                             angle[index_of_profile], rotation_center_x[index_of_profile],  rotation_center_y[index_of_profile]);

    patches.addInterface(0, boundary::north, (size_t) 0, boundary::south);
    patches.addInterface(2, boundary::north, 2, boundary::south);
   // patches.addAutoBoundaries();
    gsInfo << patches << "\n";


    if (plotGeometry){
        gsWriteParaview( patches, "ProfileGeometryPatches"+profile,10000);
        gsFileManager::open("ProfileGeometryPatches.pvd");
        return 0;
    }



    // ========================================= Define problem ========================================= 
    
    gsBoundaryConditions<> bcInfo;
    gsFunctionExpr<> f("0", "0", 2); // external force
    gsFunctionExpr<> *P; // pressure boundary condition (decoupled case)
    gsFunctionExpr<> *Uin, *Ublade; // relative velocity BCs for inlet/wall

    std::ostringstream strs_vel_x;
    std::ostringstream strs_vel_y;
    std::ostringstream strs_velblade;

    // Boundary conditions - old approach

    if(relative == true){
        strs_vel_x << velocity_relative_x(index_of_profile);
        std::string str_x = strs_vel_x.str();
        strs_vel_y << velocity_relative_y(index_of_profile);
        std::string str_y = strs_vel_y.str();

        Uin = new gsFunctionExpr<>(str_x, str_y, 2);
        Ublade = new gsFunctionExpr<>("0", "0", 2);
        P = new gsFunctionExpr<>("0", "0", 2);
    } else
    {
        strs_vel_x << velocity_absolute_x(index_of_profile);
        std::string str_x = strs_vel_x.str();
        strs_vel_y << velocity_absolute_y(index_of_profile);
        std::string str_y = strs_vel_y.str();
        strs_velblade << velocity_blade(index_of_profile);
        std::string str2 = strs_velblade.str();
        Uin = new gsFunctionExpr<>(str_x, str_y, 2);
        Ublade = new gsFunctionExpr<>("0",str2, 2);
        P = new gsFunctionExpr<>("0", "0", 2);
    }

    gsInfo << "\n";
    gsInfo << "Uin:" << *Uin;
    gsInfo << "\n";


    gsInfo << "\n";
    gsInfo << "Ublade:" << *Ublade;
    gsInfo << "\n";



/*
    //----------------prepare for 3D tests------------------

    f = new gsFunctionExpr<>("0", "0", "0", 3);

    // Boundary conditions
    Uin = new gsFunctionExpr<>("1", "0", "0", 3);
    Uwall = new gsFunctionExpr<>("0", "0", "0", 3);
    P = new gsFunctionExpr<>("0", "0", "0", 3);
*/



       bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Ublade, 0);
       bcInfo.addCondition(2, boundary::north, condition_type::dirichlet, Ublade, 0);

       bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
    if (method == 'P')
    {
        bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, P, 1);
        bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, P, 1);
    }

    //---------2D periodic conditions 8 patches-------
//       bcInfo.setIdentityMatrix(2);
//       bcInfo.addPeriodic(0, boundary::north, 0, boundary::south, 2);
//       bcInfo.addPeriodic(2, boundary::north, 2, boundary::south, 2);


/*
    //---------3D periodic conditions 8 patches-------
       bcInfo.setIdentityMatrix(3);
       bcInfo.addPeriodic(0, boundary::front, 0, boundary::back, 3);
       bcInfo.addPeriodic(1, boundary::front, 1, boundary::back, 3);
       bcInfo.addPeriodic(2, boundary::front, 2, boundary::back, 3);
       bcInfo.addPeriodic(3, boundary::front, 3, boundary::back, 3);
       bcInfo.addPeriodic(4, boundary::front, 4, boundary::back, 3);
       bcInfo.addPeriodic(5, boundary::front, 5, boundary::back, 3);
*/









    // ========================================= Define basis ========================================= 

    // Define discretization space by refining the basis of the geometry
    gsMultiBasis<real_t> tbasis(patches);
    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();


//     //-------------refine for 2D-------------------
//    gsMatrix<> box_u0(2, 2);
//    gsMatrix<> box_u1(2, 2);
//    gsMatrix<> box_v0(2, 2);
//    gsMatrix<> box_v1(2, 2);

//    //refine in u before blade
//    for (int i = 0; i < numRefineLocal+1; i++)
//    {
//        box_u0 << 1 - (uRefineKnot/ math::pow(2, i)), 1, 0, 0;

//        tbasis.refine(0, box_u0);
//        tbasis.refine(1, box_u0);
//    }

//    //refine in u at the end of blade
//    for (int i = 0; i < numRefineLocal - 1; i++)
//    {
//        box_u0 << 1 - (uRefineKnot / math::pow(2, i + 2)), 1, 0, 0;

//        tbasis.refine(2, box_u0);
//        tbasis.refine(3, box_u0);
//    }

//    //refine in u after blade
//    for (int i = 0; i < numRefineLocal; i++)
//    {
//        box_u1 << 0, uRefineKnot / math::pow(2, i + 1), 0, 0;

//        tbasis.refine(4, box_u1);
//        tbasis.refine(5, box_u1);
//    }

//    //refine in v
//    for (unsigned k = 0; k < patches.nPatches(); k++) {
//        for (int i = 0; i < numRefineLocal; i++)
//        {
//            box_v0 << 0, 0, 0, wallRefineKnot / math::pow(2, i);
//            box_v1 << 0, 0, 1 - (wallRefineKnot / math::pow(2, i)), 1;
//            tbasis.refine(k, box_v0);
//            tbasis.refine(k, box_v1);
//        }
//    }



    std::vector< gsMultiBasis<real_t> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

    // ========================================= Define solver ========================================= 

    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);

    // --------- other options ---------
    //params.setNumThreads(1); // set the number of threads for assembly

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

 

    navStokes.solve(numIter, 1e-5); // max. 50 iterations, solution change norm tol = 10^(-5)

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
        filename << "bladeProfile_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
            << now->tm_hour << "-" << now->tm_min << ".txt";

        std::ofstream ofile;
        ofile.open(filename.str().c_str());
        ofile << "Solving  bladeProfile example:\n";
        ofile << "DG = " << dg << "\n";
        ofile << "viscosity = " << viscosity << "\n";
        ofile << "method = " << method << "\n";
        ofile << "inlet velocity = " << *Uin << "\n";
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
        //gsField<> relVelocity = navStokes.constructSolution(0, true); // in case of rotation
        gsField<> pressure = navStokes.constructSolution(1);

        // Write solution to paraview files
        gsInfo << "Plotting in Paraview...\n";
        gsWriteParaview<>(velocity, "bladeProfile_velocity"+profile, plot_pts);
        gsWriteParaview<>(pressure, "bladeProfile_pressure"+profile, plot_pts);

        //system("paraview bladeProfile_velocity.pvd");
    }

    delete Uin;
    delete Ublade;
    delete P;

    //system("pause");

//} for compute all profiles

    return 0; 
}

template<class TT>  // should be rewritten to gsMultiPatch<TT> instead of raw pointer
gsMultiPatch<TT>* BSplineProfile2DPatch(TT const & index, TT const & length_x1, TT const & length_x2, TT const & pitch, TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
                                                     TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY)

{
    //----------------set parameters for blade profile----------------
    int num_samples = 30;
    int num_samples3d = 30;
    gsVector<TT> vec(2);
    /*vec(0) = rotationCenterX*(chordLength);
    vec(1) = rotationCenterY*(chordLength);*/
    vec(0) = 0.0;
    vec(1) = 0.0;
    gsMatrix<TT> mat(2,2);
    mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
           chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
    gsBSpline<TT> suction_side_curve;
    gsBSpline<TT> pressure_side_curve;
    gsBSpline<TT> suction_side_curve_transf;
    gsBSpline<TT> pressure_side_curve_transf;
    gsKnotVector<TT> kvfit(0, 1, 4, 4);
    unsigned num_cpblade = 8;
    BladeProfile<TT> * pBladeProfile = 0;


    unsigned degree = 3;

    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<TT>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
                                       thicknessY, endingOffset, outputAngle, radius, chordLength,
                                       Angle, rotationCenterX, rotationCenterY, 0.0);



    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

    vec (0) = rotationCenterX;
    vec (1) = rotationCenterY;
    suction_side_curve.translate(-vec);
    pressure_side_curve.translate(-vec);
    pBladeProfile->setSuctionSide(suction_side_curve);
    pBladeProfile->setPressureSide(pressure_side_curve);

    pressure_side_curve_transf = pBladeProfile->getPressureSide();
    suction_side_curve_transf = pBladeProfile->getSuctionSide();

    pressure_side_curve_transf.linearTransform(mat);
    suction_side_curve_transf.linearTransform(mat);

    pBladeProfile->setPressureSide(pressure_side_curve_transf);
    pBladeProfile->setSuctionSide(suction_side_curve_transf);



    gsBSpline < TT > bs = pBladeProfile -> getPressureSide ();
    gsBSpline < TT > bp = pBladeProfile -> getSuctionSide ();
    gsMatrix < TT > cp_bp (num_cpblade, 2);
    gsMatrix < TT > cp_bs (num_cpblade, 2);

/*
     gsBSpline<TT> pressure_side_curve_transf = pBladeProfile->getPressureSide();
    gsBSpline<TT> suction_side_curve_transf = pBladeProfile->getSuctionSide();

    vec(0) = suction_side_curve_transf.coef(0, 0);
    vec(1) = pressure_side_curve_transf.coef(0,1);


   suction_side_curve_transf.translate(-vec);
   pressure_side_curve_transf.translate(-vec);

   pBladeProfile->setPressureSide(pressure_side_curve_transf);
   pBladeProfile->setSuctionSide(suction_side_curve_transf);
*/


    //---------------set parameters for boundary of patches-----------------------
    //point [length_x1, width_y1] is the lower left point of boundary rectangle
    //point [length_x2 + length_end, width_y2] is the upper right point of boundary rectangle
    real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;
    real_t width_y1 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ - pitch/2.0;
//    real_t length_x1 = -0.8;
//    real_t length_x2 = 2.0;
    real_t length_end = 0.0;
    real_t length = math::abs(length_x1) + math::abs(length_x2);
    real_t point_x = bp.coef(bp.coefsSize()-2,0)+(((-bp.coef(bp.coefsSize()-2,0)+bp.coef(bp.coefsSize()-1,0))*(bp.coef(bp.coefsSize()-2,1)-width_y2))/(bp.coef(bp.coefsSize()-2,1)-bp.coef(bp.coefsSize()-1,1)));
     real_t point_y = (bp.coef(0,0)+(((bp.coef(0,0)-bp.coef(bp.coefsSize()-1,0))*(-bp.coef(0,1)+width_y1))/(bp.coef(0,1)-bp.coef(bp.coefsSize()-1,1))) + bp.coef(0,0))/2.0;

    //control points of pressure side
    for(unsigned i = 0; i < bp.coefsSize(); i++){
        cp_bp(i, 0) = bp.coef(i, 0);
        cp_bp(i, 1) = bp.coef(i, 1);
    }
    //control points of suction side
    for(unsigned i = 0; i < bs.coefsSize(); i++){
        cp_bs(i, 0) = bs.coef(i, 0);
        cp_bs(i, 1) = bs.coef(i, 1);
    }

/*
real_t alpha = math::atan(bs.coef( bs.coefsSize()-1, 1)/bs.coef( bs.coefsSize()-1, 0));

    //control points transformed of pressure side
    for(unsigned i = 0; i < bp.coefsSize(); i++){
        cp_bs(i, 0) = bp.coef( i, 0)*math::cos(2*alpha)+bp.coef( i, 1)*math::sin(2*alpha);
        cp_bs(i, 1) = -bp.coef( i, 1)*math::cos(2*alpha)+bp.coef( i, 0)*math::sin(2*alpha);
    }
    //control points transformed of suction side
    for(unsigned i = 0; i < bs.coefsSize(); i++){
        cp_bp(i, 0) = bs.coef( i, 0)*math::cos(2*alpha)+bs.coef( i, 1)*math::sin(2*alpha);
        cp_bp(i, 1) = -bs.coef( i, 1)*math::cos(2*alpha)+bs.coef( i, 0)*math::sin(2*alpha);
    }

    */
  gsMultiPatch<TT> * mp = new gsMultiPatch<TT>;



    //initial data for patches without blade profile - patches are linear -> elevate
    gsKnotVector<TT> kv_u(0, 1, 0, 2);
    gsKnotVector<TT> kv_v(0, 1, 0, 2);
    gsTensorBSplineBasis<2, real_t> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, real_t> basis_blade(kvfit, kv_v);



    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(4, 2);
    coef_patch0 <<  length_x1, (width_y1+width_y2)/2.0,
                    bp.coef(0, 0), bp.coef(0, 1),
                    length_x1, width_y2,
                    point_y, width_y2;


    gsTensorBSpline<2, real_t> patch0(basis, coef_patch0);
    patch0.degreeElevate(degree - 1, -1);




    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1 (4, 2);
    coef_patch1 <<  length_x1, width_y1,
                    point_y, width_y1,
                    length_x1, (width_y1+width_y2)/2.0,
                    bs.coef(0, 0), bs.coef(0, 1);


    gsTensorBSpline<2, real_t> patch1(basis, coef_patch1);
    patch1.degreeElevate(degree - 1, -1);
    //--------------------------------patch 4-------------------------------------------
    gsMatrix<TT> coef_patch4(4, 2);
    coef_patch4 <<  bp.coef(num_cpblade-1, 0), bp.coef(num_cpblade-1, 1),
                   // length_x2 + length_end, bp.coef(num_cpblade-1, 1),
                    length_x2 + length_end, (width_y1+width_y2)/2.0,
                    point_x, width_y2,
                    length_x2 + length_end, width_y2;

    gsTensorBSpline<2, real_t> patch4(basis, coef_patch4);
    patch4.degreeElevate(degree - 1, -1);


    //--------------------------------patch 5-------------------------------------------
    gsMatrix<TT> coef_patch5 (4, 2);
    coef_patch5 <<  point_x, width_y1,
                    length_x2 + length_end, width_y1,
                    bs.coef(num_cpblade - 1, 0), bs.coef(num_cpblade - 1, 1),
                  //  length_x2 + length_end, bs.coef(num_cpblade - 1, 1),
                    length_x2 + length_end, (width_y1+width_y2)/2.0;

    gsTensorBSpline<2, real_t> patch5(basis, coef_patch5);
    patch5.degreeElevate(degree - 1, -1);


    real_t insertKnotNum = trunc(length_end)+1;
    if(insertKnotNum > 1.0){
        for(real_t i = 1/(insertKnotNum); i < 1.0; i+=1/(insertKnotNum)){
        patch4.insertKnot(i,0,1);
        patch5.insertKnot(i,0,1);
        }
    }
    gsInfo << "patches 0,1,4,5 was created \n";
    //--------------------------------patch 2-------------------------------------------

    gsMatrix<TT> chordalbs(1,num_cpblade);
    chordalbs = chordalParameterization(cp_bs);

    gsMatrix<TT> coef_patch2(2 * num_cpblade, 2);
    for(unsigned i = 0; i < num_cpblade; i++){
          coef_patch2(i,0)= cp_bs(i, 0);
          coef_patch2(i,1)= cp_bs(i, 1);
                }


    for(unsigned i = num_cpblade; i < 2 * num_cpblade; i++){
            coef_patch2(i, 0)= point_y +   chordalbs(i-num_cpblade) * (point_x-point_y);
            coef_patch2(i, 1)= width_y2 ;
        }


    gsTensorBSpline<2, real_t> patch2(basis_blade, coef_patch2);
    patch2.degreeElevate(degree - 1, 1);

    if(index==1){
    patch0.coef(7,0)=bs.coef(0,0)+(((-bs.coef(0,0)+bs.coef(1,0))*(bs.coef(0,1)-patch0.coef(7,1)))/(bs.coef(0,1)-bs.coef(1,1)))-chordLength*0.01;
    patch0.coef(11,0)=patch0.coef(7,0);
    patch2.coef(8,0)=patch0.coef(7,0);
    patch2.coef(16,0)=patch0.coef(11,0);};


    //--------------------------------patch 3-------------------------------------------
    gsMatrix<TT> chordalbp = chordalParameterization(cp_bp);
    gsMatrix<TT> coef_patch3 (2 * num_cpblade, 2);
    for(unsigned i = num_cpblade; i < 2 * num_cpblade; i++){
            coef_patch3(i, 0)= bp.coef(i -num_cpblade, 0);
            coef_patch3(i, 1)= bp.coef(i -num_cpblade, 1) ;
        }
    for(unsigned i = 0; i < num_cpblade; i++){
            coef_patch3(i, 0)= point_y +   chordalbp(i) * (point_x-point_y);
            coef_patch3(i, 1)= width_y1 ;
        }


    gsTensorBSpline<2, real_t> patch3(basis_blade, coef_patch3);
    patch3.degreeElevate(degree - 1, 1);

    mp->addPatch(patch0);
    mp->addPatch(patch1);
    mp->addPatch(patch2);
    mp->addPatch(patch3);
    mp->addPatch(patch4);
    mp->addPatch(patch5);




    mp->addInterface(0, boundary::south, 1, boundary::north);
    mp->addInterface(0, boundary::east, 2, boundary::west);
    mp->addInterface(1, boundary::east, 3, boundary::west);
    mp->addInterface(2, boundary::east, 4, boundary::west);
    mp->addInterface(3, boundary::east, 5, boundary::west);
    mp->addInterface(4, boundary::south, 5, boundary::north);
    mp->addAutoBoundaries();


    bool plotMeshes = false;
    if (plotMeshes) {
        gsMesh<> mesh0;
        patch0.controlNet(mesh0);

        gsWriteParaview(mesh0,"cpPatch0");

        gsMesh<> mesh1;
        patch1.controlNet(mesh1);

        gsWriteParaview(mesh1,"cpPatch1");

        gsMesh<> mesh2;
        patch2.controlNet(mesh2);

        gsWriteParaview(mesh2,"cpPatch2");

        gsMesh<> mesh3;
        patch3.controlNet(mesh3);

        gsWriteParaview(mesh3,"cpPatch3");

        gsMesh<> mesh4;
       patch4.controlNet(mesh4);
        gsWriteParaview(mesh4,"cpPatch4");

        gsMesh<> mesh5;
       patch5.controlNet(mesh5);
        gsWriteParaview(mesh5,"cpPatch5");



   }
    return mp;
}

//template<class TT>
//gsMultiPatch<TT>*  BSplineProfile2DPatchC1(TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
//                                                     TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY)

//{
//    //----------------set parameters for blade profile----------------
//    int num_samples = 30;
//    gsBSpline<TT> suction_side_curve;
//    gsBSpline<TT> pressure_side_curve;
//    gsKnotVector<TT> kvfit(0, 1, 4, 4);
//    unsigned num_cpblade = 8;
//    BladeProfile<TT> * pBladeProfile = 0;

//    //---------------set parameters for boundary of patches-----------------------
//    //point [length_x1, width_y1] is the lower left point of boundary rectangle
//    //point [length_x2 + length_end, width_y2] is the upper right point of boundary rectangle
//    real_t width_y2 = 0.91412/2 + 0.1;
//    real_t width_y1 = -0.91412/2 - 0.25;
//    real_t length_x1 = -0.8;
//    real_t length_x2 = chordLength + 1.0;
//    real_t length_end = 2.0;
//    real_t length = math::abs(length_x1) + math::abs(length_x2);
//    unsigned degree = 3;

//    //---------------compute blade profile for given parameters----------------------
//    pBladeProfile = new BladeProfile<TT>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
//                                       thicknessY, endingOffset, outputAngle, radius, chordLength,
//                                       Angle, rotationCenterX, rotationCenterY, 0.0);



//    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 1);

//    suction_side_curve.rotate(-Angle);
//    pressure_side_curve.rotate(-Angle);
//    pBladeProfile->setSuctionSide(suction_side_curve);
//    pBladeProfile->setPressureSide(pressure_side_curve);

//    gsBSpline<TT> bp = pBladeProfile->getPressureSide();
//    gsBSpline<TT> bs = pBladeProfile->getSuctionSide();
//    gsMatrix<TT> cp_bp (num_cpblade,2);
//    gsMatrix<TT> cp_bs (num_cpblade,2);

//    //control points of pressure side
//    for(unsigned i = 0; i < bp.coefsSize(); i++){
//        cp_bp(i, 0) = bp.coef(i, 0);
//        cp_bp(i, 1) = bp.coef(i, 1);
//    }
//    //control points of suction side
//    for(unsigned i = 0; i < bs.coefsSize(); i++){
//        cp_bs(i, 0) = bs.coef(i, 0);
//        cp_bs(i, 1) = bs.coef(i, 1);
//    }

//   gsMultiPatch<TT> * mp = new gsMultiPatch<TT>;



//    //initial data for patches without blade profile - patches are linear -> elevate
//    gsKnotVector<TT> kv_u(0, 1, 0, 2);
//    gsKnotVector<TT> kv_v(0, 1, 0, 2);
//    gsTensorBSplineBasis<2, real_t> basis(kv_u, kv_v);
//    gsTensorBSplineBasis<2, real_t> basis_blade(kvfit, kv_v);



//    //--------------------------------patch 0-------------------------------------------
//    gsMatrix<TT> coef_patch0(4, 2);
//    coef_patch0 <<  length_x1, 0,
//                    bp.coef(0, 0), bp.coef(0, 1),
//                    length_x1, width_y2,
//                    length_x1 + length/4, width_y2;


//    gsTensorBSpline<2, real_t> patch0(basis, coef_patch0);
//    patch0.degreeElevate(degree - 1, -1);




//    //--------------------------------patch 1-------------------------------------------
//    gsMatrix<TT> coef_patch1pom(4, 2);
//    coef_patch1pom <<  length_x1, width_y1,
//                    length_x1 + length/4, width_y1,
//                    length_x1, 0,
//                    bs.coef(0, 0), bs.coef(0, 1);



//    gsTensorBSpline<2, real_t> patch1pom(basis, coef_patch1pom);
//    patch1pom.degreeElevate(degree - 1, -1);

//    gsMatrix<TT> coef_patch1(9, 2);
//    coef_patch1=patch1pom.coefs();

//    coef_patch1(4,0)=2*patch0.coef(1,0)-patch0.coef(4,0);
//    coef_patch1(4,1)=2*patch0.coef(1,1)-patch0.coef(4,1);



//    gsTensorBSpline<2, real_t> patch1(patch1pom.basis(), coef_patch1);

//    //--------------------------------patch 4-------------------------------------------
//    gsMatrix<TT> coef_patch4(4, 2);
//    coef_patch4 <<  bp.coef(num_cpblade-1, 0), bp.coef(num_cpblade-1, 1),
//                    length_x2 + length_end, bp.coef(num_cpblade-1, 1),
//                    length_x2 - length/4, width_y2,
//                    length_x2 + length_end, width_y2;

//    gsTensorBSpline<2, real_t> patch4(basis, coef_patch4);
//    patch4.degreeElevate(degree - 1, -1);


//    //--------------------------------patch 5-------------------------------------------
//    gsMatrix<TT> coef_patch5 (4, 2);
//    coef_patch5 <<  length_x2 - length/4, width_y1,
//                    length_x2 + length_end, width_y1,
//                    bs.coef(num_cpblade - 1, 0), bs.coef(num_cpblade - 1, 1),
//                    length_x2 + length_end, bs.coef(num_cpblade - 1, 1);

//    gsTensorBSpline<2, real_t> patch5(basis, coef_patch5);
//    patch5.degreeElevate(degree - 1, -1);


//    real_t insertKnotNum = trunc(length_end)+1;
//    if(insertKnotNum > 1.0){
//        for(real_t i = 1/(insertKnotNum); i < 1.0; i+=1/(insertKnotNum)){
//        patch4.insertKnot(i,0,1);
//        patch5.insertKnot(i,0,1);
//        }
//    }
//    gsInfo << "patches 0,1,4,5 was created \n";
//    //--------------------------------patch 2-------------------------------------------

//    gsMatrix<TT> chordalbs(1,num_cpblade);
//    chordalbs = chordalParameterization(cp_bs);

//    gsMatrix<TT> coef_patch2(2 * num_cpblade, 2);
//    for(unsigned i = 0; i < num_cpblade; i++){
//          coef_patch2(i,0)= cp_bs(i, 0);
//          coef_patch2(i,1)= cp_bs(i, 1);
//                }


//    for(unsigned i = num_cpblade; i < 2 * num_cpblade; i++){
//            coef_patch2(i, 0)= length_x1 + (length/4.0) +   chordalbs(i-num_cpblade) * (length/2);
//            coef_patch2(i, 1)= width_y2 ;
//        }


//    gsTensorBSpline<2, real_t> patch2(basis_blade, coef_patch2);
//    patch2.degreeElevate(degree - 1, 1);

//    //--------------------------------patch 3-------------------------------------------
//    gsMatrix<TT> chordalbp = chordalParameterization(cp_bp);
//    gsMatrix<TT> coef_patch3 (2 * num_cpblade, 2);
//    for(unsigned i = num_cpblade; i < 2 * num_cpblade; i++){
//            coef_patch3(i, 0)= bp.coef(i -num_cpblade, 0);
//            coef_patch3(i, 1)= bp.coef(i -num_cpblade, 1) ;
//        }
//    for(unsigned i = 0; i < num_cpblade; i++){
//            coef_patch3(i, 0)= length_x1 + length/4 +   chordalbp(i) * (length/2);
//            coef_patch3(i, 1)= width_y1 ;
//        }


//    gsTensorBSpline<2, real_t> patch3(basis_blade, coef_patch3);
//    patch3.degreeElevate(degree - 1, 1);

//    mp->addPatch(patch0);
//    mp->addPatch(patch1);
//    mp->addPatch(patch2);
//    mp->addPatch(patch3);
//    mp->addPatch(patch4);
//    mp->addPatch(patch5);




//    mp->addInterface(0, boundary::south, 1, boundary::north);
//    mp->addInterface(0, boundary::east, 2, boundary::west);
//    mp->addInterface(1, boundary::east, 3, boundary::west);
//    mp->addInterface(2, boundary::east, 4, boundary::west);
//    mp->addInterface(3, boundary::east, 5, boundary::west);
//    mp->addInterface(4, boundary::south, 5, boundary::north);
//    mp->addAutoBoundaries();


//    bool plotMeshes = true;
//    if (plotMeshes) {
//        gsMesh<> mesh0;
//        patch0.controlNet(mesh0);

//        gsWriteParaview(mesh0,"cpCirclePatch0");

//        gsMesh<> mesh1;
//        patch1.controlNet(mesh1);

//        gsWriteParaview(mesh1,"cpCirclePatch1");

//        gsMesh<> mesh2;
//        patch2.controlNet(mesh2);

//        gsWriteParaview(mesh2,"cpCirclePatch2");

//        gsMesh<> mesh3;
//        patch3.controlNet(mesh3);

//        gsWriteParaview(mesh3,"cpCirclePatch3");

//        gsMesh<> mesh4;
//       patch4.controlNet(mesh4);
//        gsWriteParaview(mesh4,"cpCirclePatch4");

//        gsMesh<> mesh5;
//       patch5.controlNet(mesh5);
//        gsWriteParaview(mesh5,"cpCirclePatch5");



//   }
//    return mp;
//}

template<class TT>
gsMultiPatch<TT>* BSplineProfile3DPatch(TT const & index, TT const & length_x1, TT const & length_x2, TT const & pitch, TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
                                        TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY)
{

    unsigned degree = 3;
    real_t depth = 1.8;

      gsMultiPatch<TT> *patches2D = new gsMultiPatch<TT>;
    patches2D = BSplineProfile2DPatch<real_t>(index, length_x1, length_x2, pitch, camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength,
                                       Angle, rotationCenterX, rotationCenterY);

 gsMultiPatch<TT> * mpatches = new gsMultiPatch<TT>;

 for(unsigned n=0; n<patches2D->nPatches();n++) {
     gsMatrix<TT> coef_patch(2*((patches2D->patch(n)).coefsSize()), 3);
     coef_patch.setZero(2*((patches2D->patch(n)).coefsSize()), 3);
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

        gsTensorBSpline<3,TT> patch3D(patch_basis,coef_patch);


      patch3D.degreeElevate(degree - 1, 2);
      mpatches->addPatch(patch3D);

   }

    mpatches->addInterface(0, boundary::south, 1, boundary::north);
    mpatches->addInterface(0, boundary::east, 2, boundary::west);
    mpatches->addInterface(1, boundary::east, 3, boundary::west);
    mpatches->addInterface(2, boundary::east, 4, boundary::west);
    mpatches->addInterface(3, boundary::east, 5, boundary::west);
    mpatches->addInterface(4, boundary::south, 5, boundary::north);


    mpatches->addAutoBoundaries();


    return mpatches;
}

template<class TT>
gsMultiPatch<TT> BSplineProfile2DBetweenPatch(TT const & index,
                                       TT const & length_x1,
                                       TT const & length_x2,
                                       TT const & pitch,
                                       TT const & camberX,
                                       TT const & camberY,
                                       TT const & leadingAngle,
                                       TT const & trailingAngle,
                                       TT const & thicknessX,
                                       TT const & thicknessY,
                                       TT const & endingOffset,
                                       TT const & outputAngle,
                                       TT const & radius,
                                       TT const & chordLength,
                                       TT const & Angle,
                                       TT const & rotationCenterX,
                                       TT const & rotationCenterY)
{


    //----------------set parameters for blade profile----------------
    bool plot = false;
    int num_samples = 30;
    gsVector<TT> vec(2);
    //gsInfo << pitch << "\n ";
    vec (0) = rotationCenterX;
    vec (1) = rotationCenterY;
    gsMatrix<TT> mat(2,2);
    mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
           chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
    //gsInfo << vec << "\n \n";
    //gsInfo << mat << "\n \n";
    gsBSpline<TT> suction_side_curve;
    gsBSpline<TT> pressure_side_curve;
    gsBSpline<TT> suction_side_curve_transf;
    gsBSpline<TT> pressure_side_curve_transf;
    gsKnotVector<TT> kvfit(0, 1, 4, 4);
    unsigned num_cpblade = 8;
    BladeProfile<TT> * pBladeProfile = 0;
    //unsigned degree = 3;
    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<TT>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
                                       thicknessY, endingOffset, outputAngle, radius, chordLength,
                                       Angle, rotationCenterX, rotationCenterY, 0.0);
    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);
     //---------------transform given profile----------------------
    //gsInfo << suction_side_curve.coefs();
    //gsInfo << pressure_side_curve.coefs();

    suction_side_curve.translate(-vec);
    pressure_side_curve.translate(-vec);
    pBladeProfile->setSuctionSide(suction_side_curve);
    pBladeProfile->setPressureSide(pressure_side_curve);
    pressure_side_curve_transf = pBladeProfile->getPressureSide();
    suction_side_curve_transf = pBladeProfile->getSuctionSide();
    pressure_side_curve_transf.linearTransform(mat);
    suction_side_curve_transf.linearTransform(mat);

    //gsInfo << suction_side_curve_transf.coefs();
    //gsInfo << pressure_side_curve_transf.coefs();

    pBladeProfile->setPressureSide(pressure_side_curve_transf);
    pBladeProfile->setSuctionSide(suction_side_curve_transf);
    gsBSpline < TT > bs = pBladeProfile -> getPressureSide ();
    gsBSpline < TT > bp = pBladeProfile -> getSuctionSide ();
    gsMatrix < TT > cp_bp (num_cpblade, 2);
    gsMatrix < TT > cp_bp_pom (num_cpblade, 2);
    gsMatrix < TT > cp_bs (num_cpblade, 2);

    //---------------set parameters for boundary of patches-----------------------
    //real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

    //control points of pressure side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bp(i, 0) = bp.coef(i, 0) ;
        cp_bp(i, 1) = bp.coef(i, 1) + pitch;
        cp_bp_pom(i, 0) = bp.coef(i, 0) ;
        cp_bp_pom(i, 1) = bp.coef(i, 1) ;
    }
    //control points of suction side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bs(i, 0) = bs.coef(i, 0);
        cp_bs(i, 1) = bs.coef(i, 1);
    }

    /*gsInfo << "\n pressure side \n";
    gsInfo << cp_bp;
    gsInfo << "\n pressure side pom \n";
    gsInfo << cp_bp_pom;
    gsInfo << "\n suction side \n";
    gsInfo << cp_bs;
    */

    gsMultiPatch<TT> mp;
    // compute discrete coons patch to optimize
    gsMatrix<TT> coef_patchAll (56, 2);
    coef_patchAll.setZero(56, 2);
    gsMatrix<TT> a_cp(14,2);
    gsMatrix<TT> b_cp(14,2);
    gsMatrix<TT> c_cp(4,2);
    gsMatrix<TT> d_cp(4,2);

    gsMatrix<TT> coef_patchStart (4, 2);
    coef_patchStart << length_x1, cp_bs(0,1),
                         3.0*length_x1/4.0 + cp_bs(0,0)/4.0, cp_bs(0,1),
                         length_x1/4.0 + 3.0*cp_bs(0,0)/4.0, cp_bs(0,1),
                         cp_bs(0,0),cp_bs(0,1);

     for (unsigned int i = 0; i<4;i++)
     {
        a_cp(i,0)=coef_patchStart(i,0);
        a_cp(i,1)=coef_patchStart(i,1);
        b_cp(i,0)=coef_patchStart(i,0);
        b_cp(i,1)=coef_patchStart(i,1)+pitch;
     }

    for (unsigned int i = 4; i<10;i++)
    {
        a_cp(i,0)=cp_bs(i-3,0);
        a_cp(i,1)=cp_bs(i-3,1);
        b_cp(i,0)=cp_bp(i-3,0);
        b_cp(i,1)=cp_bp(i-3,1);
    }

    real_t yend_coor = -((cp_bs(0,1)*cp_bs(7,0) - cp_bs(0,0)*cp_bs(7,1) - cp_bs(0,1)*length_x2 + cp_bs(7,1)*length_x2)/(
                          cp_bs(0,0) - cp_bs(7,0)));
    gsMatrix<TT> coef_patchEnd (4, 2);
    coef_patchEnd << cp_bs(7,0),cp_bs(7,1),
                     length_x2/4.0 + 3.0*cp_bs(7,0)/4.0, yend_coor/4.0 + 3.0*cp_bs(7,1)/4.0,
                     3.0*length_x2/4.0 + cp_bs(7,0)/4.0, 3.0*yend_coor/4.0 + cp_bs(7,1)/4.0,
                     length_x2, yend_coor;
    for (unsigned int i = 10; i<14;i++)
    {
        a_cp(i,0)=coef_patchEnd(i-10,0);
        a_cp(i,1)=coef_patchEnd(i-10,1);
        b_cp(i,0)=coef_patchEnd(i-10,0);
        b_cp(i,1)=coef_patchEnd(i-10,1)+pitch;
    }
        c_cp << length_x1,a_cp(0,1),
                length_x1,1.0*b_cp(0,1)/4.0 + 3.0*a_cp(0,1)/4.0,
                length_x1,3.0*b_cp(0,1)/4.0 + 1.0*a_cp(0,1)/4.0,
                length_x1, b_cp(0,1);
        d_cp << length_x2,a_cp(13,1),
                length_x2,1.0*b_cp(13,1)/4.0 + 3.0*a_cp(13,1)/4.0,
                length_x2,3.0*b_cp(13,1)/4.0 + 1.0*a_cp(13,1)/4.0,
                length_x2, b_cp(13,1);
    /*gsInfo << length_x1;
    gsInfo << "\n";
    gsInfo << length_x2;
    gsInfo << "\n";
    //gsInfo << width_y1;
    gsInfo << "\n";
    gsInfo << width_y2;
    gsInfo << "\n";
    gsInfo<< "\n-------------------------\n";
    gsInfo << a_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << b_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << c_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << d_cp;
    */
    gsKnotVector<TT> kv_uu(0, 1, 4, 4);
    gsKnotVector<TT> kv_vv(0, 1, 0, 4);
    kv_uu.insert(0.2/3.0,3);
    kv_uu.insert(1-(0.2/3.0),3);
    gsInfo<< kv_uu;
    gsMultiPatch<TT> * boundaries4 = new gsMultiPatch<TT>;
    boundaries4->addPatch(gsBSpline<TT>( kv_vv, c_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_uu, a_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_vv, d_cp));
    boundaries4->addPatch(gsBSpline<TT>( kv_uu, b_cp));
    gsCoonsPatch<TT> patchAll = coonsPatch(*boundaries4);
    patchAll.compute();
    mp.addPatch(patchAll.result());

   //=================================optimization===========================================

    gsVector<TT> area_vec(7);
    area_vec.setZero(7);
    area_vec << 1,1,0.25,0.25,0,0,0;

    real_t orthogonality = 0.0;
    real_t skewness = 0.0;
    real_t eccentricity = 0.0;
    real_t intersection = 0.0;
    real_t uniformity = 0.25;
    real_t area = area_vec(index);
    real_t length = 0;
    real_t epsilon = 1e-7;

    gsQualityMeasure<TT> optimization(mp.patch(0));
    real_t opt_val = optimization.functional(orthogonality, skewness,
                                             eccentricity, uniformity,
                                             length, area,
                                             intersection, epsilon);
    optimization.optimize(orthogonality, skewness,
                          eccentricity, uniformity,
                          length, area,
                          intersection, epsilon);
    gsInfo << "Value of functional: "
           << opt_val
           << "\n";

    if(plot)
    {
        gsFileData<> fileData;
        fileData << mp.patch(0);
        std::string out;
        out = "optimize_blade"+ util::to_string(index) +".xml";
        fileData.dump(out);
        gsMesh<TT> mesh;
        makeMesh(mp.patch(0).basis(), mesh, 10);
        mp.patch(0).evaluateMesh(mesh);
        out = "optimize_bladeMesh" +  util::to_string(index) ;
        gsWriteParaview(mesh, out);
        gsMesh<TT> mesh2;
        mp.patch(0).controlNet(mesh2);
        out = "optimize_bladeControlNet" +  util::to_string(index);
        gsWriteParaview(mesh2,out);
    }

    //=================================divide gsGeometry into three patches===========================================

    //initial data for patches
    gsMatrix<TT> coefs = mp.patch(0).coefs();
    gsMultiPatch<TT> mpFinal;
    gsKnotVector<TT> kv_u(0, 1, 0, 4);
    gsKnotVector<TT> kv_v(0, 1, 0, 4);
    gsTensorBSplineBasis<2, TT> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, TT> basis_blade(kvfit, kv_v);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(16, 2);
    coef_patch0.setZero(16,2);

    for(int i = 0; i < 4; i++)
    {
       for(int j = 0; j < 4; j++)
       {
          coef_patch0(i*4+j,0) = coefs(i*4+j+i*10,0);
          coef_patch0(i*4+j,1) = coefs(i*4+j+i*10,1);
       }
    }
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    for(real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch0.insertKnot(knot,0);
        patch0.insertKnot(knot,1);
    }

    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1 (num_cpblade*4, 2);
    coef_patch1.setZero(num_cpblade*4,2);

    for(int i = 0; i < 4; i++)
    {
       for(unsigned j = 0; j < num_cpblade; j++)
       {
          coef_patch1(i*num_cpblade+j,0) = coefs(i*4+j+3+i*10,0);
          coef_patch1(i*num_cpblade+j,1) = coefs(i*4+j+3+i*10,1);
       }
    }
    //gsInfo << "\n pressure side \n";
    //gsInfo << coef_patch1<< "\n";
    gsTensorBSpline<2, TT> patch1(basis_blade, coef_patch1);
    for(real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch1.insertKnot(knot,1);
    }

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2(16, 2);
    coef_patch2.setZero(16,2);

    for(int i = 0; i < 4; i++)
    {
       for(int j = 0; j < 4; j++)
       {
          coef_patch2(i*4+j,0) = coefs(i*4+j+10+i*10,0);
          coef_patch2(i*4+j,1) = coefs(i*4+j+10+i*10,1);
       }
    }
    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    for(real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch2.insertKnot(knot,0);
        patch2.insertKnot(knot,1);
    }



    mpFinal.addPatch(patch0);
    mpFinal.addPatch(patch1);
    mpFinal.addPatch(patch2);

    mpFinal.addInterface(0, boundary::east, 1, boundary::west);
    mpFinal.addInterface(1, boundary::east, 2, boundary::west);
    mpFinal.addAutoBoundaries();

    return mpFinal;
}



