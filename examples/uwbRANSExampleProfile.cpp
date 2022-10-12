/** @file uwbRANSExampleProfile.cpp

Author(s): H. Hornikova, E. Turnerova
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

#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"

//#include <sstream>


//#include <gsIO/gsIOUtils.h>
#include <gsOptimizer/gsQualityMeasure.h>

#include <math.h>
#include <string.h>
#include <sstream>
#include "../jku/gsMotorUtils.h"
//#include "../jku/gsQualityMeasure2.h"

//#include "uwbProfileOptimization.h"
#include "gsModeling/gsCoonsPatch.h"
//#include <../../motor/jku/gsMotorUtils.h>

#include "uwbObjectiveFunctionEvaluator.h"

using namespace gismo;




template<class TT>
gsMultiPatch<TT> BSplineProfile2DPatch(TT const & index,
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
    TT const & rotationCenterY,
    TT const & uniformity_param);

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, real_t velocity_x = 1.0, real_t velocity_y = 0.0, real_t velocity_blade = 0.0, bool geometry = false);
template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfo, real_t kIn, real_t kBlade, real_t oIn, real_t oBlade, bool geometry = false);
template<class T> void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, bool geometry = false);
template<class T> void refineBasisUniformZones(gsMultiBasis<T>& tbasis, int numRefineUniformLocal_v,
    int numUniformKnot_v); //in v direction
template<class T> void refineBasisZones(gsMultiBasis<T>& tbasis,
    int numRefineLocal_v, int numKnot_v,
    int numRefineLocal_u0, int numRefineLocal_u1_s, int numRefineLocal_u1_e, int numRefineLocal_u2, int numKnot_u0, int numKnot_u1_s, int numKnot_u1_e, int numKnot_u2);
    template<class T> void solvePoissonEquation(gsMultiPatch<T> patches, uwbTMSolverKOmega<T>& turbSolver_unsteady, std::string profile, bool geometry, int plot_pts = 10000);

int main(int argc, char *argv[])
{
    const real_t PI = 3.141592653589793238463;

    // ========================================= Settings =========================================
    bool plot = true;
    int plot_pts = 30000;
    bool geometry_between = true;

    real_t timeStep = 0.001;

    int numIterSteadyNS = 50;
    int numIterKOmegaSteady = 50;
    int maxNumIterRANS = 2;
    int minNumIterRANS = 1;
    real_t tolObjRelVal = 1e-3;
    real_t tolRelNorm = 1e-2;
    int maxRANSPicardIt = 2;
    int maxTMPicardFirstIt = 50;

    real_t viscosity = 0.000001;
    real_t viscositySteady = 0.01;
    real_t turbIntensity = 0.05;
    real_t viscosityRatio = 500.;

    unsigned index_of_profile = 6; // i_blade for all profiles; 0-6, do not use 0 for geometry_between = false
    //-----------------------------------------------
    std::string tmEvaluator = "koWilcoxLRN";
    //std::string tmEvaluator = "koSST";

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
    //-----------------------------------------------

    // geometry/mesh
    int setting_of_blade_rotation = 0; //0-optimal, 1-maximal, 2-minimal
    int setting_of_guide_vanes = 0; //0-compute with the formula, 1- 56° GV for optimal and 60° GV for maximal, 2-68° GV for optimal a 72° GV for maximal

    int numRefine;
    if (geometry_between)
    {
        numRefine = 1;
    }
    else
    {
        numRefine = 3;
    }

    //local refinement -- only for geometry_between = true (false - parameters in refineBasis())
    //int numRefineLocal = 6; //blade and periodic wall + before blade

    //uniform refine in v direction from the knot numUniformKnot_v
    int numRefineUniformLocal_v = 0;
    int numUniformKnot_v = 1;

    int numRefineLocal_v = 4;
    int numKnot_v = 1;
    int numRefineLocal_u0 = 4;
    int numKnot_u0 = 1;
    int numRefineLocal_u1_s = 4;
    int numKnot_u1_s = 1;
    int numRefineLocal_u1_e = 0;
    int numKnot_u1_e = 1;
    int numRefineLocal_u2 = 0;
    int numKnot_u2 = 1;

    int numRefineLocalFirstKnot_v = 0;

    /*int numRefineUniformLocal_v = 0;
    int numUniformKnot_v = 1;
    //refine zones in all patches numRefineLocal_v in v direction from the numKnot_v,
    //in the first patch in u direction  numRefineLocal_u0 from the numKnot_u0 at the end
    //in the second patch in u direction  numRefineLocal_u1_s from the numKnot_u1_s at the begin.
    //in the second patch in u direction  numRefineLocal_u1_e from the numKnot_u1_e at the end
    //in the second patch in u direction  numRefineLocal_u2 from the numKnot_u2 at the begin.

    int numRefineLocal_v = 4;
    int numKnot_v = 1;
    int numRefineLocal_u0 = 4;
    int numKnot_u0 = 1;
    int numRefineLocal_u1_s = 4;
    int numKnot_u1_s = 1;
    int numRefineLocal_u1_e = 0;
    int numKnot_u1_e = 1;
    int numRefineLocal_u2 = 0;
    int numKnot_u2 = 1;

    int numRefineLocalFirstKnot_v = 0; //after refinement one more time from the first knot (refined net) in v*/

    bool TMsupg = false;
    int tauSUPGType = 0; // 0 'h/(2*deg*norm(u)) * (cotgh(Pe)) - 1/Pe'; 1 'h/(2*deg*norm(u)'; 2 '((2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
                         // 3 '((2/timeStep)^2 + (2norm(u)/h)^2 + 9(4nu/h^2)^2)^-0.5'
    bool AFC = false;

    bool loadIC = false;
    gsMatrix<real_t> RANSSolution, TMSolution;
    //------------- read from file ---------------------
    if (loadIC)
    {
        gsFileData<> fdRead("RANS_IC_solution.xml");
        RANSSolution = *(fdRead.getFirst< gsMatrix<real_t> >());

        gsFileData<> fdReadTM("TM_IC_solution.xml");
        TMSolution = *(fdReadTM.getFirst< gsMatrix<real_t> >());
    }
    //-------------------------------------------------

    bool printInfo = true;
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
        file << "tolObjRelVal : " << tolObjRelVal << "\n";
        file << "tolRelNorm : " << tolRelNorm << "\n";
        file << "tmEvaluator = " << tmEvaluator << "\n";
        file << "maxTMPicardFirstIt : " << maxTMPicardFirstIt << "\n";
        file << "iterative : " << iterative << "\n";
        file << "lin. max iter. : " << linMaxIt << "\n";
        file << "lin. tol : " << linTol << "\n";
        file << "numRefine: " << numRefine << "\n";
        //file << "refineLocal: " << numRefineLocal << "\n";
        file << "numRefineUniformLocal_v: " << numRefineUniformLocal_v << "\n";
        file << "numUniformKnot_v: " << numUniformKnot_v << "\n";
        file << "numRefineLocal_v: " << numRefineLocal_v << "\n";
        file << "numKnot_v: " << numKnot_v << "\n";
        file << "numRefineLocal_u0: " << numRefineLocal_u0 << "\n";
        file << "numKnot_u0: " << numKnot_u0 << "\n";
        file << "numRefineLocal_u1: " << numRefineLocal_u1_s << "\n";
        file << "numKnot_u1: " << numKnot_u1_s << "\n";
        file << "numRefineLocal_u1: " << numRefineLocal_u1_e << "\n";
        file << "numKnot_u1: " << numKnot_u1_e << "\n";
        file << "numRefineLocal_u1: " << numRefineLocal_u2 << "\n";
        file << "numKnot_u1: " << numKnot_u2 << "\n";

        file << "=======================================\n";
    }

    //int numThreads = 10; // number of threads for assembly

    ////command line
    //gsCmdLine cmd("Solves Navier-Stokes problem with an isogeometric discretization.");
    //cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    //cmd.addInt("r", "uniformRefine",
    //    "Number of Uniform h-refinement steps to perform before solving", numRefine);
    //cmd.addInt("s", "plotSamples",
    //    "Number of sample points to use for plotting", plot_pts);

    //try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (numRefine<0)
    {
        gsInfo << "Number of refinements must be non-negative, quitting.\n";
        return -1;
    }

    gsAssemblerOptions opt;
    opt.dirStrategy = dirichlet::elimination;
    opt.intStrategy = iFace::glue;

    gsInfo << "Solving turbulent flow for 2d blade profile.\n";


    // ========================================= Define geometry =========================================
    real_t rotation_of_blade;
    real_t flow_rate; //5.73 - optimal, 8.7 - maximal, 1.58 - minimal
    real_t uniformity_param;

    switch (setting_of_blade_rotation) {
      case 1: //maximal
         rotation_of_blade = 10; //in degrees
         flow_rate = 8.7;
         uniformity_param = 0.005;   //parameter for mesh
      break;

      case 2: //minimal
        rotation_of_blade = -15; //in degrees
        flow_rate = 1.58;
        uniformity_param = 0.01;
      break;

      default: //default is optimal
         rotation_of_blade = 0;
         flow_rate = 5.73;
         uniformity_param = 0.01;
      break;
    };

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
    gsVector<> right_angle((num_bladeprofiles));
    right_angle.setConstant(num_bladeprofiles, 90);
    gsVector<> rotate_angle((num_bladeprofiles));
    rotate_angle.setConstant(num_bladeprofiles, rotation_of_blade);

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
    angle = (right_angle - (angle + rotate_angle))*PI/180;
    rotation_center_x << 0.494758396, 0.469497406, 0.444542963, 0.417724545, 0.390108787, 0.361175154, 0.330805204;
    rotation_center_y << 0.060495569, 0.028225794, 0.00125711, -0.006884641, -0.010228889, -0.010435203, -0.00079539;

    gsVector<> length_x1(num_bladeprofiles);
    length_x1 << -0.186559, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
    length_x1 *= 2.;
    real_t length_x2 = 0.433544;


    /*real_t length_x1_fixed = -0.8;
    real_t length_x2_fixed;
    if (geometry_between)
        length_x2_fixed = 1;
    else
        length_x2_fixed = 2;*/
    //  real_t rotate_angle = 0.0 * PI/180;

    //=======================Input velocities==================================================================
        real_t omega = 0.;//2*PI*538.0/60.0;   // value 538 cycles per minute is given
        real_t omega_rot = 2*PI*538.0/60.0;
        //real_t gravity = 9.81;
        //real_t etah = 0.945;
        //real_t H = 9.88929;
        //real_t Hn = H*etah;

        //real_t flow_rate = 5.73; //5.73 - optimal, 8.7 - maximal, 1.58 - minimal
        real_t D_out = 0.5 * 2;
        real_t d_out = 0.122926 * 2;
        real_t vel_mer_out = flow_rate / (PI*(D_out*D_out - d_out*d_out) / 4.0);
        real_t vel_mer_in = vel_mer_out;                //meridial velocity

        gsVector<> velocity_absolute_x(num_bladeprofiles);
        gsVector<> velocity_absolute_y(num_bladeprofiles);
        gsVector<> velocity_absolute_y_data(num_bladeprofiles);
        gsVector<> velocity_relative_x(num_bladeprofiles);
        gsVector<> velocity_relative_y(num_bladeprofiles);
        gsVector<> v_u(num_bladeprofiles);
        gsVector<> v_target(num_bladeprofiles);
        gsVector<> v_t_input(num_bladeprofiles);

        gsVector<> velocity_blade(num_bladeprofiles);    //angular velocity
        for (unsigned i = 0; i<num_bladeprofiles; i++) {
            velocity_blade(i) = -rr(i)*omega;
        }

        for (unsigned i = 0; i<num_bladeprofiles; i++) {
            v_u(i) = rr(i)*omega_rot;
        }

        switch(setting_of_guide_vanes)
        {
            case 1: //56° guide vanes for opt. runner blade, 68° guide vanes for max. runner blade
            {
                switch(setting_of_blade_rotation)
                {
                    case 1: // 68° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 8.04, 9.22, 10.14, 11.03, 12.01, 13.51, 14.29;
                        velocity_absolute_y_data << 5.99, 5.91, 5.29, 4.95, 4.8, 4.72, 4.27;
                        for (unsigned i = 0; i<num_bladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                    default: //56° guide vanes, optimal rotation of runner blade
                    {
                        velocity_relative_x << 6.45, 6.75, 7.3, 7.77, 8.31, 8.31, 8.92;
                        velocity_absolute_y_data << 7.44, 7.21, 6.27, 5.72, 5.4, 5.17, 4.48;
                        for (unsigned i = 0; i<num_bladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                }
            }
            break;

            case 2: //60° guide vanes for opt. runner blade, 72° guide vanes for max. runner blade
            {
                switch(setting_of_blade_rotation)
                {
                    case 1: // 72° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 8.21, 9.47, 10.48, 11.44, 12.52, 14.12, 15.01;
                        velocity_absolute_y_data << 4.88, 4.84, 4.38, 4.12, 4.03, 3.99, 3.64;
                        for (unsigned i = 0; i<num_bladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;

                    default: // 60° guide vanes, maximal rotation of runner blade
                    {
                        velocity_relative_x << 6.54, 6.88, 7.5, 8.03, 8.64, 8.64, 9.33;
                        velocity_absolute_y_data << 6.39, 6.24, 5.47, 5.03, 4.8, 4.63, 4.05;
                        for (unsigned i = 0; i<num_bladeprofiles; i++) {
                           velocity_relative_y(i) =  v_u(i) - velocity_absolute_y_data(i);
                        }
                    }
                    break;
                }
            }
            break;

            default: //compute with formula
            {
                angle_input = angle - leading_angle;
                //      gsInfo << angle*180/PI;
                //      gsInfo<<leading_angle*180/PI;
                //    gsInfo<< "\n angle between input velocity and x-axis (degrees) \n";
                //    gsInfo << angle_input*180/PI;
                for (unsigned i = 0; i<num_bladeprofiles; i++) {
                velocity_relative_x(i) = vel_mer_in;

                // //-------leading angle possibility
                //velocity_relative_y(i) = vel_mer_in*math::tan(angle_input(i));

                // //-------formula possibility
                //v_target(i) =(0.4*rr(i)+0.1)*vel_mer_out;
                //v_t_input(i) = v_target(i) + ((etah*gravity*Hn)/v_u(i));
                //velocity_relative_y(i) =  v_u(i) - v_t_input(i);

                // //-------data with deviation
                //velocity_absolute_y_data << -0.28569, 4.35998, 8.07131, 11.2644, 14.018, 16.5, 18.7736;  //94.98568% of water flow angle with formula
                //velocity_absolute_y_data << -0.27613, 4.18687, 7.64825, 10.5193, 12.9108, 15, 16.8598;  //91.80894% of water flow angle with formula
                //velocity_absolute_y_data << -0.27971, 4.25139, 7.80418, 10.7908, 13.3101, 15.5357, 17.5373; //93
                switch(setting_of_blade_rotation)
                {
                    case 1: //maximal
                    //velocity_absolute_y_data << -0.28273, 4.33645, 8.104, 11.4163, 14.3233, 16.9817, 19.4471; //94
                    velocity_absolute_y_data << -0.74342, 3.67571, 7.26132, 10.38162, 13.08805, 15.53386, 17.77616; //92
                    break;

                    default: //optimal
                    velocity_absolute_y_data << -0.28272, 4.30592, 7.93756, 11.0259, 13.6596, 16.0096, 18.1491; //94


                    break;
                }

                velocity_relative_y(i) =  velocity_absolute_y_data(i);
                }
            }
            break;
        }

        for (unsigned i = 0; i<num_bladeprofiles; i++) {
            velocity_absolute_x(i) = velocity_relative_x(i);
            velocity_absolute_y(i) = velocity_relative_y(i) + velocity_blade(i);
        }

/*
      //gsInfo << "vel_mer_in " << vel_mer_in;

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

    //========================================= Define problem =========================================

    gsInfo << "===================================================\n";
    std::ostringstream strs_profile;
    std::string profile;
    strs_profile << index_of_profile;
    profile = "_blade" + strs_profile.str();
    gsInfo << "profil" + profile;
    gsInfo << "\n";

    gsBoundaryConditions<> bcInfo;
    defineBCs_NS(bcInfo, velocity_absolute_x(index_of_profile), velocity_absolute_y(index_of_profile), velocity_blade(index_of_profile), geometry_between);
    gsFunctionExpr<> f("0", "0", 2);

    gsMultiPatch<real_t> patches;
    if (geometry_between)
    {
        patches = BSplineProfile2DBetweenPatch<real_t>(index_of_profile,
            length_x1[index_of_profile],//length_x1_fixed,
            length_x2, //length_x2_fixed,
            ((2 * PI*rr[index_of_profile]) / num_blades),
            camber_x[index_of_profile],
            camber_y[index_of_profile],
            leading_angle[index_of_profile],
            trailing_angle[index_of_profile],
            thickness_x[index_of_profile],
            thickness_y[index_of_profile],
            ending_offset[index_of_profile],
            output_angle[index_of_profile],
            radius[index_of_profile],
            chord_length[index_of_profile],
            angle[index_of_profile],
            rotation_center_x[index_of_profile],
            rotation_center_y[index_of_profile],
            uniformity_param);

        patches.addInterface(0, boundary::north, (size_t)0, boundary::south);
        patches.addInterface(2, boundary::north, 2, boundary::south);
    }
    else
    {
        patches = BSplineProfile2DPatch<real_t>(index_of_profile,
            length_x1[index_of_profile],//length_x1_fixed,
            length_x2, //length_x2_fixed,
            ((2 * PI*rr[index_of_profile]) / num_blades),
            camber_x[index_of_profile],
            camber_y[index_of_profile],
            leading_angle[index_of_profile],
            trailing_angle[index_of_profile],
            thickness_x[index_of_profile],
            thickness_y[index_of_profile],
            ending_offset[index_of_profile],
            output_angle[index_of_profile],
            radius[index_of_profile],
            chord_length[index_of_profile],
            angle[index_of_profile],
            rotation_center_x[index_of_profile],
            rotation_center_y[index_of_profile]);

        patches.addInterface(0, boundary::north, 1, boundary::south);
        patches.addInterface(4, boundary::north, 5, boundary::south);
    }

    patches.addAutoBoundaries();
    gsInfo << patches << "\n";

    // ========================================= Define basis =========================================
    gsMultiBasis<> tbasis(patches); // basis for RANS equations


    //uniform refinement:
    for (int i = 0; i < numRefine; ++i)
        tbasis.uniformRefine();
    /*gsMatrix<> box_u(2, 2);
    box_u << 0, 1, 0, 0;
    gsMatrix<> box_v(2, 2);
    box_v << 0, 0, 0, 1;
    tbasis.refine(0, box_v);
    tbasis.refine(1, box_u);
    tbasis.refine(1, box_v);
    tbasis.refine(2, box_u);
    tbasis.refine(2, box_v);*/

    //not uniform refinement:

    //refineBasis(tbasis, numRefine, numRefineLocal, geometry_between);
    //uniform refine

    //case without refineBasisUniform
    //refineBasisZones(tbasis, numRefineLocal_v, numKnot_v, numRefineLocal_u0, numRefineLocal_u1_s, numRefineLocal_u1_e, numRefineLocal_u2, numKnot_u0, numKnot_u1_s, numKnot_u1_e, numKnot_u2); //for geometry_between = true

    //refine uniform in v from the knot, then local refinement
    refineBasisUniformZones(tbasis, numRefineUniformLocal_v, numUniformKnot_v);
    /*gsMatrix<> box_u0(2, 2);
    box_u0 << 0, 1, 0, 0;
    tbasis.refine(1, box_u0);*/
    //tbasis.refine(2, box_u0);
    refineBasisZones(tbasis, numRefineLocal_v, math::pow(2, numRefineUniformLocal_v) * numKnot_v, numRefineLocal_u0, numRefineLocal_u1_s, numRefineLocal_u1_e, numRefineLocal_u2, numKnot_u0, numKnot_u1_s, numKnot_u1_e, numKnot_u2); //for geometry_between = true
                                                                                                                                                                                                                                       //one more local refinement in v from 1st knot (in already refined mesh)
    refineBasisZones(tbasis, numRefineLocalFirstKnot_v, 1, 0, 0, 0, 0, 0, 0, 0, 0); //for geometry_between = true
                                                                                    //-------------------------------------------------------------------------------------------------
    bool plotMeshes = true;
    if (plotMeshes)
    {
        gsMesh<> mesh;
        makeMesh(tbasis.at(0), mesh, 10);
        patches.patch(0).evaluateMesh(mesh);
        gsWriteParaview(mesh, "meshPatch0" + profile);
        gsMesh<> mesh1;
        makeMesh(tbasis.at(1), mesh1, 10);
        patches.patch(1).evaluateMesh(mesh1);
        gsWriteParaview(mesh1, "meshPatch1" + profile);
        gsMesh<> mesh2;
        makeMesh(tbasis.at(2), mesh2, 10);
        patches.patch(2).evaluateMesh(mesh2);
        gsWriteParaview(mesh2, "meshPatch2" + profile);
    }

    std::vector< gsMultiBasis<> >  discreteBases;
    discreteBases.push_back(tbasis);//Basis for velocity
    discreteBases.push_back(tbasis);//Basis for pressure
    discreteBases[0].degreeElevate(1); //elevate the velocity space

                                       //    for (int j = 0; j < patches.nPatches(); j++)
                                       //        {
                                       //            const gsBasis<> & basis = (discreteBases.at(0)).piece(j);
                                       //            gsInfo << "basis.maxDegree() = " << basis.maxDegree() << "\n";

                                       //            const gsBasis<> & basisPressure = (discreteBases.at(1)).piece(j);
                                       //            gsInfo << "basisPressure.maxDegree() = " << basisPressure.maxDegree() << "\n";

                                       //        }

    gsInfo << "Velocity basis degree: " << discreteBases[0].degree() << "\n\n";
    file << "Velocity basis degree: " << discreteBases[0].degree() << "\n";

    // ========================================= Compute Steady NS =========================================
    uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);
    uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver

    gsField<> pressureSteady;
    if (!loadIC)
    {
        gsInfo << "Solving Steady case: \n";
        gsInfo << "numDofs: " << navStokes.numDofs() << "\n";
        //getchar();

        navStokes.initialize(); // steady solver
        navStokes.solve(numIterSteadyNS, 1e-5);

        gsField<> velocitySteady = navStokes.constructSolution(0);
        gsWriteParaview<>(velocitySteady, "profile_SteadyVelocity" + profile, plot_pts);//, true);
        pressureSteady = navStokes.constructSolution(1);
        gsWriteParaview<>(pressureSteady, "profile_SteadyPressure" + profile, plot_pts);

        gsFileData<real_t> fd;
        fd << navStokes.getSolution();
        fd.save("NSsteadySolution.xml");
    }

    //-------- read from file--------------
    //gsInfo << "Solving Steady case: \n";
    //gsInfo << "\nSteady solution loaded from xml file.\n\n";
    //uwbINSPde<real_t> NSpde(patches, bcInfo, f, viscositySteady);
    //uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);
    //uwbINSSolverSteady<real_t> navStokes(params); // steady coupled solver
    //gsFileData<> fdRead("NSsteadySolution.xml");
    //gsMatrix<real_t> NSSteadySolution = *(fdRead.getFirst< gsMatrix<real_t> >());
    //navStokes.setSolution(NSSteadySolution);
    //---------------------

    /*gsInfo << "Computing initial velocity by solving Stokes problem...\n";
    uwbINSPde<real_t> NSpde(*patches, bcInfo, f, viscosity);
    uwbINSSolverParams<real_t> params(NSpde, discreteBases, opt);
    uwbINSSolverSteady<real_t> stokesSolver(params);
    gsInfo << "numDofs: " << stokesSolver.numDofs() << "\n";
    stokesSolver.initialize();
    stokesSolver.setStokesSolution();*/

    //================================= wall distance estimation ==================================================
    real_t inletWidth = (((2 * PI*rr[index_of_profile]) / num_blades) / 2.0);
    gsVector<real_t> Uin(2);
    Uin << velocity_absolute_x(index_of_profile), velocity_absolute_y(index_of_profile);
    real_t uFreeStream = Uin.norm();
    real_t Re = uFreeStream * inletWidth / viscosity;

    //vector of the sides of the patches from which the wall distance is computed
    std::vector<boxSide> distanceSides;
    distanceSides.push_back(boundary::north);
    distanceSides.push_back(boundary::south);

    //vector of indexes of the patches corresponding to distanceSides
    //length of the vector distancePatches must be equal to the length of vector distanceSides
    gsVector<int> distancePatches(2);
    distancePatches << 1, 1;

    int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
    real_t maxYplus = 2.5; //maximum dimensionless wall distance which is accepted

                           //table with wall distance information for every chosen side is printed, if the second from the end input parameter is set as true
                           //estimation of the wall distance will be computed, if the last input parameter is set as true
    real_t wallDistance = navStokes.computeDimensionlessWallDistance(distancePatches, distanceSides, viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
    gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";
    file << "\nminimum wallDistance = " << wallDistance << "\n";

    //==================================== compute aspect ratio  ==================================================
    real_t maxAspectRatio = navStokes.getAssembler()->getBlockAssembler().computeAspectRatio();
    gsInfo << "maxAspectRatio = " << maxAspectRatio << "\n";
    file << "maxAspectRatio = " << maxAspectRatio << "\n";
    // ========================================= Define turbulence solver =========================================
    gsBoundaryConditions<> bcInfoTurb;

    real_t kInConst = 1.5 * math::pow(uFreeStream * turbIntensity, 2); // (3/2)*(UI)^2
    gsInfo << "\nkInConst = " << kInConst << "\n";
    file << "\nkInConst = " << kInConst << "\n";
    real_t oInConst = kInConst / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at inlet
    gsInfo << "oInConst = " << oInConst << "\n";
    file << "oInConst = " << oInConst << "\n";

    real_t kBlade = 1.5 * math::pow(velocity_blade(index_of_profile) * turbIntensity, 2);
    real_t oBlade;
    if (omega) {
        oBlade = kBlade / (viscosity * viscosityRatio); // need to satisfy nu_T / nu approximately at the wall
    }
    else {
        real_t beta = 0.0708;
        oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2));
    }
    gsInfo << "kBlade = " << kBlade << "\n";
    gsInfo << "oBlade = " << oBlade << "\n\n";
    file << "kBlade = " << kBlade << "\n";
    file << "oBlade = " << oBlade << "\n\n";

    defineBCs_TM(bcInfoTurb, kInConst, kBlade, oInConst, oBlade, geometry_between);

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
//    gsStopwatch time;
//    real_t Tassembly;
//    real_t Tsolve;
    if (!loadIC)
    {
        gsInfo << "\nComputing steady linearized k-omega...\n";
        gsInfo << "initialization...\n";
        gsField<> velocity = navStokes.constructSolution(0);
        gsInfo << "TMnumDofs = " << koMapper.freeSize() << "\n";
        turbSolver.initialize(velocity);
        //Tassembly = time.stop();
        //gsInfo << "Assembly time:" << Tassembly << "\n";
        //time.restart();
        turbSolver.solve(numIterKOmegaSteady, 1e-5); // solution change norm tol = 10^(-5)
        //Tsolve = time.stop();
        //gsInfo << "Solve time:" << Tsolve << "\n";
        gsField<> kOmegaSteady = turbSolver.constructSolution();
        gsWriteParaview<>(kOmegaSteady, "profile_kOmegaSteady" + profile, plot_pts, true);
    }

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
        koParams.settings().set(constantsINS::TMsupg, true); // set SUPG
        koParams.settings().set(constantsINS::tauStabTypeSUPG, tauSUPGType); //set formula for stabilization parameter tau
    }
    if (AFC)
        koParams.settings().set(constantsINS::TMafc, true);

    uwbTMSolverKOmega<real_t> turbSolver_unsteady(koParams);
    if (loadIC)
        turbSolver_unsteady.setInitialCondition(TMSolution);
    else
        turbSolver_unsteady.setInitialCondition(turbSolver.getSolution());//(koInitial);

    if (tmEvaluator == "koSST")
        solvePoissonEquation(patches, turbSolver_unsteady, profile, geometry_between, plot_pts);

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
    //time.restart();
    if (loadIC)
        ransSolver.setInitialCondition(RANSSolution);
    else
        ransSolver.setInitialCondition(navStokes.getSolution());
    ransSolver.initialize();
    //Tassembly = time.stop();

    gsInfo << "numDofs: " << ransSolver.numDofs() << "\n";
    file << "numDofs: " << ransSolver.numDofs() << "\n";

    //time.restart();

    real_t pitch = ((2 * PI*rr[index_of_profile]) / num_blades);
    gsInfo << "pitch = " << pitch << "\n";
    file << "pitch = " << pitch << "\n";

    //=== solving with solution change relative norm as stopping criterion ============
    //ransSolver.solve(numIterRANS, 1e-5); // solution change norm tol = 10^(-5)

    //=== solving with objective function change relative value as stopping criterion==
    gsVector<real_t> normalStream(2);
    normalStream << 0.0, 1.0;
    normalStream = normalStream / normalStream.norm();

    //sides at which is computed the lift
    std::vector<boxSide> profileSides;
    profileSides.push_back(boundary::north); // sunction side
    profileSides.push_back(boundary::south); // pressure side
                                             //patches at which is computed the lift (must have the same length as the vector of profileSides!!)
    gsVector<int> profilePatchNumbers(2);
    profilePatchNumbers << 1, 1;

    //sides at which is computed the velocity part of the objective function
    std::vector<boxSide> outSides;
    outSides.push_back(boundary::east);
    //patches at which is computed the velocity part of the objective function (must have the same length as the vector of outSides!!)
    gsVector<int> outPatchNumbers(1);
    outPatchNumbers << 2;

    gsField<> uSol = ransSolver.constructSolution(0);
    gsField<> pSol = ransSolver.constructSolution(1);

    real_t angularVelocity = 2 * PI*538.0 / 60.0;
    real_t uCircumferential = -angularVelocity * rr[index_of_profile];
    // target relative tangential velocity
    real_t uTangentialTarget = -(0.4 * rr[index_of_profile] + 0.1) * vel_mer_out - uCircumferential;

    uwbObjectiveFunctionEvaluator<real_t> objFunction(uSol, pSol, discreteBases, patches);
    objFunction.initialize(profilePatchNumbers, profileSides, normalStream, rr[index_of_profile],
        outPatchNumbers, outSides, uTangentialTarget);
    // weights for lift, velocity part at outlet, pressure distribution at profile
    // sum of the weights must equal to 1!!!
    objFunction.setWeights(0.5, 0.0, 0.5, 0.0);
    // set parameters necessary for evaluation of the objective function
    objFunction.setTargetPressureAtProfile(-0.15, 0.15);
    //gsVector<real_t> pTargetProfile(2);
    //pTargetProfile << -10, 10; //must be ordered according to the profileSides
    //objFunction.setTargetPressureAtProfile(pTargetProfile);
    //--- set initial lift if it is not known yet ---
    //objFunction.computeInitialLift();
    //real_t initialLift = objFunction.getInitialLift();
    //--- set initial lift if it is already known ---
    //objFunction.setInitialLift(initialLift);

    //------------ set initial objective parts ----------------
    //objFunction.setInitialObjectiveParts(initialObjectiveParts);

    ransSolver.solveWithObjective(objFunction, tolObjRelVal, tolRelNorm, maxNumIterRANS, minNumIterRANS);
    //------------ get objective parts ----------------
    gsField<> uSolution = ransSolver.constructSolution(0);
    gsField<> pSolution = ransSolver.constructSolution(1);
    objFunction.updateSolution(uSolution, pSolution);
    objFunction.setTargetPressureAtProfile(-0.15, 0.15);
    objFunction.evaluate();
    gsVector<real_t> objectiveParts = objFunction.getObjectives();
    std::ofstream fileObjective;
    fileObjective.open("InitialObjectiveParts.txt");
    fileObjective << objectiveParts << "\n";
    fileObjective.close();
    std::ofstream fileTargetValues;
    fileTargetValues.open("TargetPressure.txt");
    fileTargetValues << objFunction.getPressureTargetProfile() << "\n";
    fileTargetValues.close();

    if (plot)
        ransSolver.plotObjFunValue("profile_ObjectiveFnGraph");
    //---------------------------------------------------------------------------------------------------

    //Tsolve = time.stop();


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
    fd.save("profile" + profile + "_RANS_solution.xml");
    fd_TM.save("profile" + profile + "_TM_solution.xml");
    //=====================================================================================================

    // Optionally plot solution in paraview
    if (plot)
    {
        gsInfo << "Plotting in Paraview...\n";

        gsField<> velocity = ransSolver.constructSolution(0);
        gsField<> pressure = ransSolver.constructSolution(1);
        gsField<> kOmega = turbSolver_unsteady.constructSolution();

        gsWriteParaview<>(velocity, "profile_velocity" + profile, plot_pts);
        gsWriteParaview<>(pressure, "profile_pressure" + profile, plot_pts);
        gsWriteParaview<>(kOmega, "profile_komega" + profile, plot_pts);
        ransSolver.plotTurbulentViscosity("profile_turb_viscosity" + profile, plot_pts);
    }

    return 0;
}

template<class T> void defineBCs_NS(gsBoundaryConditions<T>& bcInfo, real_t velocity_x, real_t velocity_y, real_t velocity_blade, bool geometry)
{

    std::ostringstream strs_vel_x;
    std::ostringstream strs_vel_y;
    std::ostringstream strs_velblade;

    strs_vel_x << velocity_x;
    std::string str_x = strs_vel_x.str();
    strs_vel_y << velocity_y;
    std::string str_y = strs_vel_y.str();
    strs_velblade << velocity_blade;
    std::string str2 = strs_velblade.str();

    gsFunctionExpr<> Uin(str_x, str_y, 2);
    gsFunctionExpr<> Ublade("0", str2, 2);

    if (geometry)
    {
        bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(1, boundary::north, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
    }
    else
    {
        bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(3, boundary::north, condition_type::dirichlet, Ublade, 0);
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uin, 0);
        bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, Uin, 0);
    }

}

template<class T> void defineBCs_TM(gsBoundaryConditions<T>& bcInfoTurb, real_t kIn, real_t kBlade, real_t oIn, real_t oBlade, bool geometry)
{
    // Boundary conditions
    gsFunctionExpr<> Kin(util::to_string(kIn), 2);
    gsFunctionExpr<> Oin(util::to_string(oIn), 2);
    gsFunctionExpr<> Kblade(util::to_string(kBlade), 2);
    gsFunctionExpr<> Oblade(util::to_string(oBlade), 2);

    if (geometry)
    {
        bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Kblade, 0);
        bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Kblade, 0);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Kin, 0);



        bcInfoTurb.addCondition(1, boundary::south, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(1, boundary::north, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Oin, 1);
    }
    else
    {
        bcInfoTurb.addCondition(2, boundary::south, condition_type::dirichlet, Kblade, 0);
        bcInfoTurb.addCondition(3, boundary::north, condition_type::dirichlet, Kblade, 0);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Kin, 0);
        bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Kin, 0);


        bcInfoTurb.addCondition(2, boundary::south, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(3, boundary::north, condition_type::dirichlet, Oblade, 1);
        bcInfoTurb.addCondition(0, boundary::west, condition_type::dirichlet, Oin, 1);
        bcInfoTurb.addCondition(1, boundary::west, condition_type::dirichlet, Oin, 1);
    }
}

    template<class T>
    void solvePoissonEquation(gsMultiPatch<T> patches, uwbTMSolverKOmega<T>& turbSolver_unsteady, std::string profile, bool geometry, int plot_pts)
    {
        bool plotMeshes = true;

        int numRefinePoisson;
        if (geometry)
            numRefinePoisson = 4;
        else
            numRefinePoisson = 3;

        int numRefineUniformLocal_vPoisson = 0;
        int numUniformKnot_vPoisson = 1;

        int numRefineLocal_vPoisson = 0;//4;
        int numKnot_vPoisson = 2;
        int numRefineLocal_u0Poisson = 0;//4;
        int numKnot_u0Poisson = 1;
        int numRefineLocal_u1_sPoisson = 0;//4;
        int numKnot_u1_sPoisson = 1;
        int numRefineLocal_u1_ePoisson = 0;
        int numKnot_u1_ePoisson = 1;
        int numRefineLocal_u2Poisson = 0;
        int numKnot_u2Poisson = 1;

        int numRefineLocalFirstKnot_vPoisson = 0;//5; //after refinement one more time from the first knot (refined net) in v

        gsMultiBasis<> tbasisPoisson(patches); // basis for RANS equations
        for (int i = 0; i < numRefinePoisson; ++i)
            tbasisPoisson.uniformRefine();
        refineBasisUniformZones(tbasisPoisson, numRefineUniformLocal_vPoisson, numUniformKnot_vPoisson);
        gsMatrix<> box_u0Poisson(2, 2);
        box_u0Poisson << 0, 1, 0, 0;
        tbasisPoisson.refine(1, box_u0Poisson);
        refineBasisZones(tbasisPoisson, numRefineLocal_vPoisson, math::pow(2, numRefineUniformLocal_vPoisson) * numKnot_vPoisson, numRefineLocal_u0Poisson, numRefineLocal_u1_sPoisson, numRefineLocal_u1_ePoisson, numRefineLocal_u2Poisson, numKnot_u0Poisson, numKnot_u1_sPoisson, numKnot_u1_ePoisson, numKnot_u2Poisson); //for geometry_between = true
                                                                                                                                                                                                                                           //one more local refinement in v from 1st knot (in already refined mesh)
        refineBasisZones(tbasisPoisson, numRefineLocalFirstKnot_vPoisson, 1, 0, 0, 0, 0, 0, 0, 0, 0);
        if (plotMeshes)
        {
            gsMesh<> mesh;
            makeMesh(tbasisPoisson.at(0), mesh, 10);
            patches.patch(0).evaluateMesh(mesh);
            gsWriteParaview(mesh, "meshPoissonPatch0" + profile);
            gsMesh<> mesh1;
            makeMesh(tbasisPoisson.at(1), mesh1, 10);
            patches.patch(1).evaluateMesh(mesh1);
            gsWriteParaview(mesh1, "meshPoissonPatch1" + profile);
            gsMesh<> mesh2;
            makeMesh(tbasisPoisson.at(2), mesh2, 10);
            patches.patch(2).evaluateMesh(mesh2);
            gsWriteParaview(mesh2, "meshPoissonPatch2" + profile);
        }

        gsFunctionExpr<real_t> fw("1", 2);
        gsFunctionExpr<real_t> gw("0", 2);
        gsFunctionExpr<real_t> wallw("0.0", 2);
        gsBoundaryConditions<real_t> bcInfow;
        bcInfow.addCondition(0, boundary::north, condition_type::neumann, gw);
        bcInfow.addCondition(0, boundary::south, condition_type::neumann, gw);
        bcInfow.addCondition(0, boundary::west, condition_type::neumann, gw);
        bcInfow.addCondition(2, boundary::north, condition_type::neumann, gw);
        bcInfow.addCondition(2, boundary::south, condition_type::neumann, gw);
        bcInfow.addCondition(2, boundary::east, condition_type::neumann, gw);
        bcInfow.addCondition(1, boundary::north, condition_type::dirichlet, wallw);
        bcInfow.addCondition(1, boundary::south, condition_type::dirichlet, wallw);

        gsInfo << "\nSolving Poisson equation.\n";

        turbSolver_unsteady.setPoissonSolution(patches, tbasisPoisson, bcInfow, fw, true, plot_pts);

        gsInfo << "Poisson equation resolved.\n\n";
        turbSolver_unsteady.plotWallDistance("profile_wall_distance" + profile, plot_pts);
    }

template<class T> void refineBasis(gsMultiBasis<T>& tbasis, int numRefineUniform, int numRefineLocal, bool geometry)
{
    if (geometry)
    {


        gsMatrix<> box_v(2, 2);
        gsMatrix<> box_u(2, 2);

        //uniform refine
        for (int i = 0; i < numRefineUniform; i++)
        {
            tbasis.uniformRefine();
        }

        const gsTensorBSplineBasis<2, real_t>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
        real_t uRefineKnot = basis->knot(1, basis->degree(1) + 1); // first non-zero knot in direction v
                                                                   //gsInfo << "uRefineKnot: " << uRefineKnot << "\n";

                                                                   //refine in v, near uper wall
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_v << 0, 0, 1 - (uRefineKnot / math::pow(2, i)), 1;

            tbasis.refine(0, box_v);
            tbasis.refine(1, box_v);
            tbasis.refine(2, box_v);
        }


        //refine in v, near bottom wall
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_v << 0, 0, 0, uRefineKnot / math::pow(2, i);

            tbasis.refine(0, box_v);
            tbasis.refine(1, box_v);
            tbasis.refine(2, box_v);
        }

        //refine in u, before blade
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_u << 1 - uRefineKnot / math::pow(2, i), 1, 0, 0;

            tbasis.refine(0, box_u);
        }

        //refine in u, blade
        for (int i = 0; i < numRefineLocal; i++)
        {
            box_u << 0, uRefineKnot / math::pow(2, i), 0, 0;

            tbasis.refine(1, box_u);
        }
    }
    else
    {
        int numLocRefBeforeBlade = 11;
        int numLocRefBeginningOfBlade = 4;
        int numLocRefEndOfBlade = 6;
        int numLocRefAfterBlade = 7;
        int numLocRefBladePart = 7;
        int numLocRefPeriodicPart = 0;//7;
        int numLocRefInput = 3;
        int numLocRefOutput = 3;

        real_t refKnotBeforeBlade = 0.1;
        real_t refKnotBeginningOfBlade = 0.1;
        real_t refKnotEndOfBlade = 0.1;
        real_t refKnotAfterBlade = 0.3;
        real_t refKnotBladePart = 0.2;
        real_t refKnotPeriodicPart = 0.1;
        real_t refKnotInput = 0.1;
        real_t refKnotOutput = 0.05;

        //================================================

        for (int i = 0; i < numRefineUniform; ++i)
            tbasis.uniformRefine();

        gsMatrix<> box_u0(2, 2);
        gsMatrix<> box_u1(2, 2);
        gsMatrix<> box_u2(2, 2);
        gsMatrix<> box_v0(2, 2);
        gsMatrix<> box_v1(2, 2);
        gsMatrix<> box_v2(2, 2);

        //refine patch 4 and 5
        box_u2 << 0, 1, 0, 0;
        tbasis.refine(4, box_u2);
        tbasis.refine(5, box_u2);

        //refine patch 0, 2 and 4
        /*box_v2 << 0, 0, wallRefineKnot, 1 - wallRefineKnot;
        tbasis.refine(0, box_v2);
        tbasis.refine(2, box_v2);
        tbasis.refine(4, box_v2);*/

        //refine in u before blade
        for (int i = 0; i < numLocRefBeforeBlade; i++)
        {
            box_u0 << 1 - (refKnotBeforeBlade / math::pow(2, i)), 1, 0, 0;

            tbasis.refine(0, box_u0);
            tbasis.refine(1, box_u0);
        }

        //refine in u at the beginning of blade
        for (int i = 0; i < numLocRefBeginningOfBlade; i++)
        {
            box_u0 << 0, (refKnotBeginningOfBlade / math::pow(2, i + 2)), 0, 0;

            tbasis.refine(2, box_u0);
            tbasis.refine(3, box_u0);
        }

        //refine in u at the end of blade
        for (int i = 0; i < numLocRefEndOfBlade; i++)
        {
            box_u0 << 1 - (refKnotEndOfBlade / math::pow(2, i + 2)), 1, 0, 0;

            tbasis.refine(2, box_u0);
            tbasis.refine(3, box_u0);
        }

        //refine in u after blade
        for (int i = 0; i < numLocRefAfterBlade; i++)
        {
            box_u1 << 0, refKnotAfterBlade / math::pow(2, i + 1), 0, 0;

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

        //periodic part
        for (int i = 0; i < numLocRefPeriodicPart; i++)
        {
            box_v0 << 0, 0, 0, refKnotPeriodicPart / math::pow(2, i);
            box_v1 << 0, 0, 1 - (refKnotPeriodicPart / math::pow(2, i)), 1;
            tbasis.refine(0, box_v1);
            tbasis.refine(2, box_v1);
            tbasis.refine(4, box_v1);
            tbasis.refine(1, box_v0);
            tbasis.refine(3, box_v0);
            tbasis.refine(5, box_v0);
        }

        //blade part
        for (int i = 0; i < numLocRefBladePart; i++)
        {
            box_v0 << 0, 0, 0, refKnotBladePart / math::pow(2, i);
            box_v1 << 0, 0, 1 - (refKnotBladePart / math::pow(2, i)), 1;
            tbasis.refine(0, box_v0);
            tbasis.refine(2, box_v0);
            tbasis.refine(4, box_v0);
            tbasis.refine(1, box_v1);
            tbasis.refine(3, box_v1);
            tbasis.refine(5, box_v1);
        }

        //refine input
        for (int i = 0; i < numLocRefInput; i++)
        {
            box_u0 << 0, (refKnotInput / math::pow(2, i)), 0, 0;
            tbasis.refine(0, box_u0);
            tbasis.refine(1, box_u0);
        }

        //refine output
        for (int i = 0; i < numLocRefOutput; i++)
        {
            box_u0 << 1 - (refKnotOutput / math::pow(2, i)), 1, 0, 0;
            tbasis.refine(4, box_u0);
            tbasis.refine(5, box_u0);
        }
    }
}

template<class T> void refineBasisUniformZones(gsMultiBasis<T>& tbasis, int numRefineUniformLocal_v,
    int numUniformKnot_v)
{
    gsMatrix<> box_v_bottom(2, 2);
    gsMatrix<> box_v_up(2, 2);

    for (int i = 0; i < numRefineUniformLocal_v; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_0 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
        gsVector<> uRefineKnot_v_start((i + 1)*(numUniformKnot_v)+1);
        gsVector<> uRefineKnot_v_end((i + 1)*(numUniformKnot_v)+1);
        uRefineKnot_v_start.setZero((i + 1)*(numUniformKnot_v)+1);
        uRefineKnot_v_end.setZero((i + 1)*(numUniformKnot_v)+1);
        int sizeKnots_v = basis_0->knots(1).size() - 1;
        for (int k = 0; k < (i + 1)*(numUniformKnot_v)+1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_v_start(k) = basis_0->knot(1, basis_0->degree(1) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
            uRefineKnot_v_end(k) = basis_0->knot(1, sizeKnots_v - (basis_0->degree(1) + k)); //last numKnot knots and 1 knot in direction v
        }

        for (int j = 0; j < (i + 1)*(numUniformKnot_v); j++)
        {
            box_v_bottom << 0, 0, uRefineKnot_v_start(j), uRefineKnot_v_start(j + 1);
            box_v_up << 0, 0, uRefineKnot_v_end(j + 1), uRefineKnot_v_end(j);

            tbasis.refine(0, box_v_up);
            tbasis.refine(1, box_v_up);
            tbasis.refine(2, box_v_up);
            tbasis.refine(0, box_v_bottom);
            tbasis.refine(1, box_v_bottom);
            tbasis.refine(2, box_v_bottom);
        }
    }
}

template<class T> void refineBasisZones(gsMultiBasis<T>& tbasis,
    int numRefineLocal_v, int numKnot_v,
    int numRefineLocal_u0, int numRefineLocal_u1_s, int numRefineLocal_u1_e, int numRefineLocal_u2, int numKnot_u0, int numKnot_u1_s, int numKnot_u1_e, int numKnot_u2)
{



    gsMatrix<> box_v_bottom(2, 2);
    gsMatrix<> box_v_up(2, 2);
    gsMatrix<> box_u0(2, 2);
    gsMatrix<> box_u1_s(2, 2);
    gsMatrix<> box_u1_e(2, 2);
    gsMatrix<> box_u2(2, 2);






    //knots in v direction assumption !!!SAME FOR ALL PATCHES
    gsVector<> uRefineKnot_v_start(numKnot_v + 1);
    gsVector<> uRefineKnot_v_end(numKnot_v + 1);
    gsVector<> uRefineKnot_u0(numKnot_u0 + 1);
    gsVector<> uRefineKnot_u1_s(numKnot_u1_s + 1);
    gsVector<> uRefineKnot_u1_e(numKnot_u1_e + 1);
    gsVector<> uRefineKnot_u2(numKnot_u2 + 1);


    //refine in v, near bottom and upper wall
    for (int i = 0; i < numRefineLocal_v; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_0 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
        uRefineKnot_v_start.setZero(numKnot_v + 1);
        uRefineKnot_v_end.setZero(numKnot_v + 1);
        int sizeKnots_v = basis_0->knots(1).size() - 1;
        for (int k = 0; k < numKnot_v + 1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_v_start(k) = basis_0->knot(1, basis_0->degree(1) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction
            uRefineKnot_v_end(k) = basis_0->knot(1, sizeKnots_v - (basis_0->degree(1) + k)); //last numKnot knots and 1 knot in direction v
        }

        for (int j = 0; j < numKnot_v; j++)
        {
            box_v_bottom << 0, 0, uRefineKnot_v_start(j), uRefineKnot_v_start(j + 1);
            box_v_up << 0, 0, uRefineKnot_v_end(j + 1), uRefineKnot_v_end(j);

            tbasis.refine(0, box_v_up);
            tbasis.refine(1, box_v_up);
            tbasis.refine(2, box_v_up);
            tbasis.refine(0, box_v_bottom);
            tbasis.refine(1, box_v_bottom);
            tbasis.refine(2, box_v_bottom);
        }
    }



    //first patch, refine in u -> before blade
    for (int i = 0; i < numRefineLocal_u0; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_u0 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(0))); //basis of the first patch
        uRefineKnot_u0.setZero(numKnot_u0 + 1);
        int sizeKnots = basis_u0->knots(0).size() - 1;

        for (int k = 0; k < numKnot_u0 + 1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_u0(k) = basis_u0->knot(0, sizeKnots - (basis_u0->degree(0) + k)); // first numKnot knots and 0 knot in direction u first patch      }


        }


        for (int j = 0; j < numKnot_u0; j++)
        {
            box_u0 << uRefineKnot_u0(j + 1), uRefineKnot_u0(j), 0, 0;

            tbasis.refine(0, box_u0);
        }
    }


    //second patch, refine in u -> leading edge of blade
    for (int i = 0; i < numRefineLocal_u1_s; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_u1_s = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(1))); //basis of the first patch

        uRefineKnot_u1_s.setZero(numKnot_u1_s + 1);
        for (int k = 0; k < numKnot_u1_s + 1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_u1_s(k) = basis_u1_s->knot(0, basis_u1_s->degree(0) + k); // first numKnot knots and 0 knot in direction v - all patches have same knots in v direction        }
        }

        for (int j = 0; j < numKnot_u1_s; j++)
        {
            box_u1_s << uRefineKnot_u1_s(j), uRefineKnot_u1_s(j + 1), 0, 0;

            tbasis.refine(1, box_u1_s);
        }
    }

    //second patch, refine in u -> end
    for (int i = 0; i < numRefineLocal_u1_e; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_u1_e = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(1))); //basis of the first patch
        uRefineKnot_u1_e.setZero(numKnot_u1_e + 1);
        int sizeKnots_1 = basis_u1_e->knots(0).size() - 1;

        for (int k = 0; k < numKnot_u1_e + 1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_u1_e(k) = basis_u1_e->knot(0, sizeKnots_1 - (basis_u1_e->degree(0) + k)); // first numKnot knots and 0 knot in direction u first patch      }


        }


        for (int j = 0; j < numKnot_u1_e; j++)
        {
            box_u1_e << uRefineKnot_u1_e(j + 1), uRefineKnot_u1_e(j), 0, 0;

            tbasis.refine(1, box_u1_e);
        }
    }

    //last patch, refine in u -> start
    for (int i = 0; i < numRefineLocal_u2; i++)
    {
        // in each refinement step take the first numKnot_v non_zero knots
        const gsTensorBSplineBasis<2, real_t>*  basis_u2 = dynamic_cast<const gsTensorBSplineBasis<2, real_t>*>(&(tbasis.basis(2))); //basis of the first patch

        uRefineKnot_u2.setZero(numKnot_u2 + 1);
        for (int k = 0; k < numKnot_u2 + 1; k++) //gsVector of first numKnots knots
        {
            uRefineKnot_u2(k) = basis_u2->knot(0, basis_u2->degree(0) + k); // first numKnot knots and 0 knot in direction u - all patches have same knots in v direction        }
        }

        for (int j = 0; j < numKnot_u2; j++)
        {
            box_u2 << uRefineKnot_u2(j), uRefineKnot_u2(j + 1), 0, 0;

            tbasis.refine(2, box_u2);
        }
    }


}

template<class TT>
gsMultiPatch<TT> BSplineProfile2DPatch(TT const & index, TT const & length_x1, TT const & length_x2, TT const & pitch, TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
    TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY)

{
    //----------------set parameters for blade profile----------------
    int num_samples = 30;
    //int num_samples3d = 30;
    gsVector<TT> vec(2);
    /*vec(0) = rotationCenterX*(chordLength);
    vec(1) = rotationCenterY*(chordLength);*/
    vec(0) = 0.0;
    vec(1) = 0.0;
    gsMatrix<TT> mat(2, 2);
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

    vec(0) = rotationCenterX;
    vec(1) = rotationCenterY;
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



    gsBSpline < TT > bs = pBladeProfile->getPressureSide();
    gsBSpline < TT > bp = pBladeProfile->getSuctionSide();
    gsMatrix < TT > cp_bp(num_cpblade, 2);
    gsMatrix < TT > cp_bs(num_cpblade, 2);

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
    real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ +pitch / 2.0;
    real_t width_y1 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ -pitch / 2.0;
    //    real_t length_x1 = -0.8;
    //    real_t length_x2 = 2.0;
    real_t length_end = 0.0;
    //real_t length = math::abs(length_x1) + math::abs(length_x2);
    real_t point_x = bp.coef(bp.coefsSize() - 2, 0) + (((-bp.coef(bp.coefsSize() - 2, 0) + bp.coef(bp.coefsSize() - 1, 0))*(bp.coef(bp.coefsSize() - 2, 1) - width_y2)) / (bp.coef(bp.coefsSize() - 2, 1) - bp.coef(bp.coefsSize() - 1, 1)));
    real_t point_y = (bp.coef(0, 0) + (((bp.coef(0, 0) - bp.coef(bp.coefsSize() - 1, 0))*(-bp.coef(0, 1) + width_y1)) / (bp.coef(0, 1) - bp.coef(bp.coefsSize() - 1, 1))) + bp.coef(0, 0)) / 2.0;

    //control points of pressure side
    for (unsigned i = 0; i < bp.coefsSize(); i++) {
        cp_bp(i, 0) = bp.coef(i, 0);
        cp_bp(i, 1) = bp.coef(i, 1);
    }
    //control points of suction side
    for (unsigned i = 0; i < bs.coefsSize(); i++) {
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
    gsMultiPatch<TT> mp;



    //initial data for patches without blade profile - patches are linear -> elevate
    gsKnotVector<TT> kv_u(0, 1, 0, 2);
    gsKnotVector<TT> kv_v(0, 1, 0, 2);
    gsTensorBSplineBasis<2, real_t> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, real_t> basis_blade(kvfit, kv_v);



    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(4, 2);
    coef_patch0 << length_x1, (width_y1 + width_y2) / 2.0,
        bp.coef(0, 0), bp.coef(0, 1),
        length_x1, width_y2,
        point_y, width_y2;


    gsTensorBSpline<2, real_t> patch0(basis, coef_patch0);
    patch0.degreeElevate(degree - 1, -1);




    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(4, 2);
    coef_patch1 << length_x1, width_y1,
        point_y, width_y1,
        length_x1, (width_y1 + width_y2) / 2.0,
        bs.coef(0, 0), bs.coef(0, 1);


    gsTensorBSpline<2, real_t> patch1(basis, coef_patch1);
    patch1.degreeElevate(degree - 1, -1);
    //--------------------------------patch 4-------------------------------------------
    gsMatrix<TT> coef_patch4(4, 2);
    coef_patch4 << bp.coef(num_cpblade - 1, 0), bp.coef(num_cpblade - 1, 1),
        // length_x2 + length_end, bp.coef(num_cpblade-1, 1),
        length_x2 + length_end, (width_y1 + width_y2) / 2.0,
        point_x, width_y2,
        length_x2 + length_end, width_y2;

    gsTensorBSpline<2, real_t> patch4(basis, coef_patch4);
    patch4.degreeElevate(degree - 1, -1);


    //--------------------------------patch 5-------------------------------------------
    gsMatrix<TT> coef_patch5(4, 2);
    coef_patch5 << point_x, width_y1,
        length_x2 + length_end, width_y1,
        bs.coef(num_cpblade - 1, 0), bs.coef(num_cpblade - 1, 1),
        //  length_x2 + length_end, bs.coef(num_cpblade - 1, 1),
        length_x2 + length_end, (width_y1 + width_y2) / 2.0;

    gsTensorBSpline<2, real_t> patch5(basis, coef_patch5);
    patch5.degreeElevate(degree - 1, -1);


    real_t insertKnotNum = trunc(length_end) + 1;
    if (insertKnotNum > 1.0) {
        for (real_t i = 1 / (insertKnotNum); i < 1.0; i += 1 / (insertKnotNum)) {
            patch4.insertKnot(i, 0, 1);
            patch5.insertKnot(i, 0, 1);
        }
    }
    gsInfo << "patches 0,1,4,5 was created \n";
    //--------------------------------patch 2-------------------------------------------

    gsMatrix<TT> chordalbs(1, num_cpblade);
    chordalbs = chordalParameterization(cp_bs);

    gsMatrix<TT> coef_patch2(2 * num_cpblade, 2);
    for (unsigned i = 0; i < num_cpblade; i++) {
        coef_patch2(i, 0) = cp_bs(i, 0);
        coef_patch2(i, 1) = cp_bs(i, 1);
    }


    for (unsigned i = num_cpblade; i < 2 * num_cpblade; i++) {
        coef_patch2(i, 0) = point_y + chordalbs(i - num_cpblade) * (point_x - point_y);
        coef_patch2(i, 1) = width_y2;
    }


    gsTensorBSpline<2, real_t> patch2(basis_blade, coef_patch2);
    patch2.degreeElevate(degree - 1, 1);

    if (index == 1) {
        patch0.coef(7, 0) = bs.coef(0, 0) + (((-bs.coef(0, 0) + bs.coef(1, 0))*(bs.coef(0, 1) - patch0.coef(7, 1))) / (bs.coef(0, 1) - bs.coef(1, 1))) - chordLength*0.01;
        patch0.coef(11, 0) = patch0.coef(7, 0);
        patch2.coef(8, 0) = patch0.coef(7, 0);
        patch2.coef(16, 0) = patch0.coef(11, 0);
    };


    //--------------------------------patch 3-------------------------------------------
    gsMatrix<TT> chordalbp = chordalParameterization(cp_bp);
    gsMatrix<TT> coef_patch3(2 * num_cpblade, 2);
    for (unsigned i = num_cpblade; i < 2 * num_cpblade; i++) {
        coef_patch3(i, 0) = bp.coef(i - num_cpblade, 0);
        coef_patch3(i, 1) = bp.coef(i - num_cpblade, 1);
    }
    for (unsigned i = 0; i < num_cpblade; i++) {
        coef_patch3(i, 0) = point_y + chordalbp(i) * (point_x - point_y);
        coef_patch3(i, 1) = width_y1;
    }


    gsTensorBSpline<2, real_t> patch3(basis_blade, coef_patch3);
    patch3.degreeElevate(degree - 1, 1);

    mp.addPatch(patch0);
    mp.addPatch(patch1);
    mp.addPatch(patch2);
    mp.addPatch(patch3);
    mp.addPatch(patch4);
    mp.addPatch(patch5);



    mp.addInterface(0, boundary::south, 1, boundary::north);
    mp.addInterface(0, boundary::east, 2, boundary::west);
    mp.addInterface(1, boundary::east, 3, boundary::west);
    mp.addInterface(2, boundary::east, 4, boundary::west);
    mp.addInterface(3, boundary::east, 5, boundary::west);
    mp.addInterface(4, boundary::south, 5, boundary::north);
    mp.addAutoBoundaries();


    bool plotMeshes = false;
    if (plotMeshes) {
        gsMesh<> mesh0;
        patch0.controlNet(mesh0);

        gsWriteParaview(mesh0, "cpPatch0");

        gsMesh<> mesh1;
        patch1.controlNet(mesh1);

        gsWriteParaview(mesh1, "cpPatch1");

        gsMesh<> mesh2;
        patch2.controlNet(mesh2);

        gsWriteParaview(mesh2, "cpPatch2");

        gsMesh<> mesh3;
        patch3.controlNet(mesh3);

        gsWriteParaview(mesh3, "cpPatch3");

        gsMesh<> mesh4;
        patch4.controlNet(mesh4);
        gsWriteParaview(mesh4, "cpPatch4");

        gsMesh<> mesh5;
        patch5.controlNet(mesh5);
        gsWriteParaview(mesh5, "cpPatch5");



    }
    return mp;
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
    TT const & rotationCenterY,
    TT const & uniformity_param)
{


    //----------------set parameters for blade profile----------------
    bool plot = false;
    int num_samples = 30;
    gsVector<TT> vec(2);
    //gsInfo << pitch << "\n ";
    vec(0) = rotationCenterX;
    vec(1) = rotationCenterY;
    gsMatrix<TT> mat(2, 2);
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
    gsBSpline < TT > bs = pBladeProfile->getPressureSide();
    gsBSpline < TT > bp = pBladeProfile->getSuctionSide();
    gsMatrix < TT > cp_bp(num_cpblade, 2);
    gsMatrix < TT > cp_bp_pom(num_cpblade, 2);
    gsMatrix < TT > cp_bs(num_cpblade, 2);

    //---------------set parameters for boundary of patches-----------------------
    //real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

    //control points of pressure side
    for (unsigned i = 0; i < num_cpblade; i++) {
        cp_bp(i, 0) = bp.coef(i, 0);
        cp_bp(i, 1) = bp.coef(i, 1) + pitch;
        cp_bp_pom(i, 0) = bp.coef(i, 0);
        cp_bp_pom(i, 1) = bp.coef(i, 1);
    }
    //control points of suction side
    for (unsigned i = 0; i < num_cpblade; i++) {
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
    gsMatrix<TT> coef_patchAll(56, 2);
    coef_patchAll.setZero(56, 2);
    gsMatrix<TT> a_cp(14, 2);
    gsMatrix<TT> b_cp(14, 2);
    gsMatrix<TT> c_cp(4, 2);
    gsMatrix<TT> d_cp(4, 2);

    real_t ystart_coor = -((cp_bs(0, 1)*cp_bs(7, 0) - cp_bs(0, 0)*cp_bs(7, 1) - cp_bs(0, 1)*length_x1 + cp_bs(7, 1)*length_x1) / (
        cp_bs(0, 0) - cp_bs(7, 0)));

    gsMatrix<TT> coef_patchStart(4, 2);
    coef_patchStart << length_x1, ystart_coor,
        3.0*length_x1 / 4.0 + cp_bs(0, 0) / 4.0, 3.0*ystart_coor / 4.0 + cp_bs(0, 1) / 4.0,
        length_x1 / 4.0 + 3.0*cp_bs(0, 0) / 4.0, ystart_coor / 4.0 + 3.0*cp_bs(0, 1) / 4.0,
        cp_bs(0, 0), cp_bs(0, 1);

    for (unsigned int i = 0; i<4; i++)
    {
        a_cp(i, 0) = coef_patchStart(i, 0);
        a_cp(i, 1) = coef_patchStart(i, 1);
        b_cp(i, 0) = coef_patchStart(i, 0);
        b_cp(i, 1) = coef_patchStart(i, 1) + pitch;
    }

    for (unsigned int i = 4; i<10; i++)
    {
        a_cp(i, 0) = cp_bs(i - 3, 0);
        a_cp(i, 1) = cp_bs(i - 3, 1);
        b_cp(i, 0) = cp_bp(i - 3, 0);
        b_cp(i, 1) = cp_bp(i - 3, 1);
    }

    real_t yend_coor = -((cp_bs(0, 1)*cp_bs(7, 0) - cp_bs(0, 0)*cp_bs(7, 1) - cp_bs(0, 1)*length_x2 + cp_bs(7, 1)*length_x2) / (
        cp_bs(0, 0) - cp_bs(7, 0)));
    gsMatrix<TT> coef_patchEnd(4, 2);
    coef_patchEnd << cp_bs(7, 0), cp_bs(7, 1),
        length_x2 / 4.0 + 3.0*cp_bs(7, 0) / 4.0, yend_coor / 4.0 + 3.0*cp_bs(7, 1) / 4.0,
        3.0*length_x2 / 4.0 + cp_bs(7, 0) / 4.0, 3.0*yend_coor / 4.0 + cp_bs(7, 1) / 4.0,
        length_x2, yend_coor;
    for (unsigned int i = 10; i<14; i++)
    {
        a_cp(i, 0) = coef_patchEnd(i - 10, 0);
        a_cp(i, 1) = coef_patchEnd(i - 10, 1);
        b_cp(i, 0) = coef_patchEnd(i - 10, 0);
        b_cp(i, 1) = coef_patchEnd(i - 10, 1) + pitch;
    }
    c_cp << length_x1, a_cp(0, 1),
        length_x1, 1.0*b_cp(0, 1) / 4.0 + 3.0*a_cp(0, 1) / 4.0,
        length_x1, 3.0*b_cp(0, 1) / 4.0 + 1.0*a_cp(0, 1) / 4.0,
        length_x1, b_cp(0, 1);
    d_cp << length_x2, a_cp(13, 1),
        length_x2, 1.0*b_cp(13, 1) / 4.0 + 3.0*a_cp(13, 1) / 4.0,
        length_x2, 3.0*b_cp(13, 1) / 4.0 + 1.0*a_cp(13, 1) / 4.0,
        length_x2, b_cp(13, 1);
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
    kv_uu.insert(0.2 / 3.0, 3);
    kv_uu.insert(1 - (0.2 / 3.0), 3);
    gsInfo << kv_uu;
    gsMultiPatch<TT> * boundaries4 = new gsMultiPatch<TT>;
    boundaries4->addPatch(gsBSpline<TT>(kv_vv, c_cp));
    boundaries4->addPatch(gsBSpline<TT>(kv_uu, a_cp));
    boundaries4->addPatch(gsBSpline<TT>(kv_vv, d_cp));
    boundaries4->addPatch(gsBSpline<TT>(kv_uu, b_cp));
    gsCoonsPatch<TT> patchAll = coonsPatch(*boundaries4);
    patchAll.compute();
    mp.addPatch(patchAll.result());

    //=================================optimization===========================================

    gsVector<TT> area_vec(7);
    area_vec.setZero(7);
    area_vec << 0.1, 0.25, 0.5, 0.5, 0.5, 0.75, 1.0;//1,1,1,1,0.75,0.75,0.75;

    real_t orthogonality = 0.0;
    real_t skewness = 0.0;
    real_t eccentricity = 0.0;
    real_t intersection = 0.0;
    real_t uniformity = uniformity_param;//0.01;// 0.25;
    //real_t uniformity = 0.005;// max setting for RB rotation
    real_t area = area_vec(index);
    real_t length = 0;
    real_t epsilon = 1e-7;

    gsQualityMeasure<TT> optimization(mp.patch(0));
    //    real_t opt_val = optimization.functional(orthogonality, skewness,
    //                                             eccentricity, uniformity,
    //                                             length, area,
    //                                             intersection, epsilon);
    optimization.optimize(orthogonality, skewness,
        eccentricity, uniformity,
        length, area,
        intersection, epsilon);


    //    gsInfo << "Value of functional: "
    //           << opt_val
    //           << "\n";

    if (plot)
    {
        gsFileData<> fileData;
        fileData << mp.patch(0);
        std::string out;
        out = "optimize_blade" + util::to_string(index) + ".xml";
        fileData.dump(out);
        gsMesh<TT> mesh;
        makeMesh(mp.patch(0).basis(), mesh, 10);
        mp.patch(0).evaluateMesh(mesh);
        out = "optimize_bladeMesh" + util::to_string(index);
        gsWriteParaview(mesh, out);
        gsMesh<TT> mesh2;
        mp.patch(0).controlNet(mesh2);
        out = "optimize_bladeControlNet" + util::to_string(index);
        gsWriteParaview(mesh2, out);
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
    coef_patch0.setZero(16, 2);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            coef_patch0(i * 4 + j, 0) = coefs(i * 4 + j + i * 10, 0);
            coef_patch0(i * 4 + j, 1) = coefs(i * 4 + j + i * 10, 1);
        }
    }
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    for (real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch0.insertKnot(knot, 0);
        patch0.insertKnot(knot, 1);
    }

    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(num_cpblade * 4, 2);
    coef_patch1.setZero(num_cpblade * 4, 2);

    for (int i = 0; i < 4; i++)
    {
        for (unsigned j = 0; j < num_cpblade; j++)
        {
            coef_patch1(i*num_cpblade + j, 0) = coefs(i * 4 + j + 3 + i * 10, 0);
            coef_patch1(i*num_cpblade + j, 1) = coefs(i * 4 + j + 3 + i * 10, 1);
        }
    }
    //gsInfo << "\n pressure side \n";
    //gsInfo << coef_patch1<< "\n";
    gsTensorBSpline<2, TT> patch1(basis_blade, coef_patch1);
    for (real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch1.insertKnot(knot, 1);
    }

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2(16, 2);
    coef_patch2.setZero(16, 2);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            coef_patch2(i * 4 + j, 0) = coefs(i * 4 + j + 10 + i * 10, 0);
            coef_patch2(i * 4 + j, 1) = coefs(i * 4 + j + 10 + i * 10, 1);
        }
    }
    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    for (real_t knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch2.insertKnot(knot, 0);
        patch2.insertKnot(knot, 1);
    }



    mpFinal.addPatch(patch0);
    mpFinal.addPatch(patch1);
    mpFinal.addPatch(patch2);

    mpFinal.addInterface(0, boundary::east, 1, boundary::west);
    mpFinal.addInterface(1, boundary::east, 2, boundary::west);
    mpFinal.addAutoBoundaries();

    return mpFinal;
}


