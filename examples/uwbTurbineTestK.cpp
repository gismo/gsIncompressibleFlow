/** @file uwbbl ades.cpp
    Author(s): B. Bastl, K. Michalkova
*/
#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif
#include <gsIO/gsIOUtils.h>
#include <iostream>
#include <gismo.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <gsOptimizer/gsQualityMeasure.h>

//#include "gsQualityMeasure.h"
#include "../jku/gsMotorUtils.h"
#include "../jku/gsQualityMeasure2.h"
#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"
#include "uwbDraftTube.h"
#include "uwbProfileOptimization.h"
#include "uwbHydraulicProfile.h"
//#include "uwbKaplanTurbineRunnerBlade.h"
//#include "uwbKaplanTurbineRunnerWheelDomain.h"
//#include "uwbKaplanTurbineGuideVane.h"
#include "gsModeling/gsCoonsPatch.h"
#include <../jku/gsMotorUtils.h>

using namespace gismo;

//=================================================================================
//------------------prototype of functions-----------------------------------------
//compute all permutations with repetitions of the given vector input and save them to matrix permutations,
template<class TT>
void permutate (const gsVector<TT>& input,
                gsMatrix<TT>& permutations,
                gsVector<TT>& index,
                int
                depth,
                int & count, int length);

//save all permutations with repetitions of the given vector input to the matrix
template<class TT>
gsMatrix<TT> permutation(gsVector<TT> input, int length);

template<class TT>
void combine(gsVector<> &input, gsMatrix<> &combinations, gsVector<TT> & data, int start, int end,
                     int index, int & row, int length);
template<class TT>
gsMatrix<TT> combination(gsVector<TT> input, int length);


//save the plus, negative and zero number of points in gsVector
template<class TT>
gsVector<TT> returnJacobianDeterminant(const gsGeometry<TT>& geom, const int points = 100,
                              const bool savePoints = false,
                              const std::string& output = "",
                              const int number = 0);

//return geometry of the space between 2D profiles
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
gsMultiPatch<TT> BSplineProfile2DBetweenPatchStart(TT const & index,
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
                                                   bool & refine);
template<class TT>
gsMultiPatch<TT> DomainBetweenBladeProfiles1(TT const & index,
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
    TT const & uniformity_param,
    std::vector<TT> const & kvfit_knots,
    bool const & coarse,
    gsVector<TT> const & geom_Params, bool & refine);

//optimizeCycle
void optimizeCycle(int iterations,
                   real_t num_weights,
                   real_t num_func,
                   int num_of_func_comb,
                   int jacPts,
                   std::string filename,
                   gsMultiPatch<> geom);

gsMatrix<> readFile(std::string inputFileName);

int main(int argc, char *argv[])
{
    //=======================================================================================
    // File content:
    // 1-GEOMETRY + linear optimization
    // 2-SAVE DATA TO SCV
    // 3-READ DATA TO CSV
    // 4-TEST SEVERAL OPTIONS
    //=======================================================================================

    //=================================GEOMETRY==========================================
    unsigned num_cores = 2; //number of cores for paralelization
    unsigned num_blades = 4;
    unsigned num_bladeprofiles = 7;
    bool readData = false;
    bool saveData = false;

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
     gsVector<> length_x1(num_bladeprofiles);
    real_t PI = 3.14;

    rr << 0.175, 0.229, 0.283, 0.338, 0.392, 0.446, 0.5;
    camber_x << 0.425908805, 0.419, 0.416039882, 0.416576484, 0.416716376, 0.419203131, 0.43280567; // 0.411886743 na originalne druhe pozici
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
    chord_length << 0.36686418, 0.42813715,   0.486298848, 0.545457321, 0.595371262, 0.615928719, 0.584588959;
    angle << 53.696487,       41.265848, 33.007703, 27.603276, 24.437586, 22.893162, 21.162381;
    angle = angle*PI/180;
    rotation_center_x << 0.494758396, 0.469497406, 0.444542963, 0.417724545, 0.390108787, 0.361175154, 0.330805204;
    rotation_center_y << 0.060495569, 0.028225794, 0.00125711, -0.006884641, -0.010228889, -0.010435203, -0.00079539;
    length_x1 << -0.186559, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
    length_x1 *= 2.;
    real_t length_x2_fixed = 0.433544;

    //choose num of weights
    unsigned num_weights = 25;
    //choose num of functions
    unsigned num_func = 6;
    int num_of_func_comb = 2;
    int jacPts=10000;
    real_t intersection = 0.0;
    real_t epsilon = 1e-7;

    //cycle of blade profiles save in csv
    #pragma omp parallel for num_threads(num_cores)
    for(int iop=0;iop<num_bladeprofiles;iop++)
    {
        //int iop = 1;

        unsigned int index_of_profile = iop;
        std::ostringstream strs_profile;
        std::string profile;
        strs_profile << index_of_profile;
        profile = strs_profile.str();
        std::string filename;
        std::string filenamein;
        filename = "dataOptimizeProfil" + profile + ".csv"; //write to filename
        filenamein = MOTOR_DATA_DIR "uwb-pilsen/result_dataOptimizeProfil" + profile + ".csv"; //read from filenamein


        gsMultiPatch<> patchInitial;
        //patches = BSplineProfile2DPatch<real_t>(index_of_profile,length_x1_fixed,length_x2_fixed,((2*PI*rr[index_of_profile])/num_blades), camber_x[index_of_profile], camber_y[index_of_profile],  leading_angle[index_of_profile],  trailing_angle[index_of_profile],  thickness_x[index_of_profile],  thickness_y[index_of_profile],  ending_offset[index_of_profile],  output_angle[index_of_profile],  radius[index_of_profile],  chord_length[index_of_profile],
        //                                                             angle[index_of_profile], rotation_center_x[index_of_profile],  rotation_center_y[index_of_profile]);
        bool refine = true;
        //patchInitial = BSplineProfile2DBetweenPatchStart<real_t>(index_of_profile,length_x1_fixed,length_x2_fixed,((2*PI*rr[index_of_profile])/num_blades), camber_x[index_of_profile], camber_y[index_of_profile],  leading_angle[index_of_profile],  trailing_angle[index_of_profile],  thickness_x[index_of_profile],  thickness_y[index_of_profile],  ending_offset[index_of_profile],  output_angle[index_of_profile],  radius[index_of_profile],  chord_length[index_of_profile],
        //                                 angle[index_of_profile], rotation_center_x[index_of_profile],  rotation_center_y[index_of_profile], refine);


        std::vector<real_t> kvfit_knots = {0.,0.,0.,0.,0.2,0.4,0.6,0.8,1.,1.,1.,1.};
       //  std::vector<real_t> kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};
        gsVector<real_t> geom_Params(1);
        geom_Params << 0.;

        patchInitial = DomainBetweenBladeProfiles1<real_t>(index_of_profile,length_x1[index_of_profile],length_x2_fixed,((2*PI*rr[index_of_profile])/num_blades), camber_x[index_of_profile], camber_y[index_of_profile],  leading_angle[index_of_profile],  trailing_angle[index_of_profile],  thickness_x[index_of_profile],  thickness_y[index_of_profile],  ending_offset[index_of_profile],  output_angle[index_of_profile],  radius[index_of_profile],  chord_length[index_of_profile],
                                         angle[index_of_profile], rotation_center_x[index_of_profile],  rotation_center_y[index_of_profile], 0.,kvfit_knots,false, geom_Params, refine);

        //first step is linear optimization
        int iterations = 5;
        real_t orthogonality = 0.0;
        real_t skewness = 0.0;
        real_t eccentricity = 0.0;
        //real_t intersection = 0.0;
        real_t uniformity = 0.1;
        real_t area = 0.0;
        real_t length = 0.1;
        //real_t epsilon = 1e-7;
        //int jacPts = 10000;
        bool dumped = false;

        gsCmdLine cmd("Optimization");
        //cmd.addString("I", "input", "Input", input);
        //cmd.addString("f", "output", "Output prefix", output);
        cmd.addInt("i", "iterations", "Number of iterations", iterations);
        cmd.addReal("o", "orthogonality",
                    "Weight of quality measure: orthogonality", orthogonality);
        cmd.addReal("s", "skewness",
                     "Weight of quality measure: skewness", skewness);
        cmd.addReal("e", "eccentricity",
                    "Weight of quality measure: eccentricity", eccentricity);
        cmd.addReal("u", "uniformity",
                     "Weight of quality measure: uniformity", uniformity);
        cmd.addReal("L", "length",
                     "Weight of quality measure: length functional", length);
        cmd.addReal("a", "area",
                     "Weight of quality measure: area", area);
        cmd.addReal("n", "intersection",
                     "Weight of quality measure: self-intersection", intersection);
        cmd.addReal("p", "epsilon", "Self intersection variable.", epsilon);
        cmd.addInt("j", "numJacPts",
                    "Number of points where we sample Jacobian determinant",
                    jacPts);
        cmd.addSwitch("dumped", "Use dumped method.", dumped);
        //try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


        gsQualityMeasure<> optimizationIt(patchInitial.patch(0));
        optimizationIt.optimize(orthogonality, skewness,
                              eccentricity, uniformity,
                              length, area,
                              intersection, epsilon);

        //====================SAVE in cvs and optimization cycle===================================
        if(saveData)
        {
            real_t opt_orto = optimizationIt.functional(1, 0,
                                                     0, 0,
                                                     0, 0,
                                                     0, 0);
            real_t opt_skew = optimizationIt.functional(0, 1,
                                                     0, 0,
                                                     0, 0,
                                                     0, 0);
            real_t opt_ecce = optimizationIt.functional(0, 0,
                                                     1, 0,
                                                     0, 0,
                                                     0, 0);
            real_t opt_unif = optimizationIt.functional(0, 0,
                                                     0, 1,
                                                     0, 0,
                                                     0, 0);
            real_t opt_length = optimizationIt.functional(0, 0,
                                                     0, 0,
                                                     1, 0,
                                                     0, 0);
            real_t opt_area = optimizationIt.functional(0, 0,
                                                     0, 0,
                                                     0, 1,
                                                     0, 0);
            real_t opt_All = optimizationIt.functional(orthogonality, skewness,
                                  eccentricity, uniformity,
                                  length, area,
                                  intersection, epsilon);


            gsInfo<< "\n\n------------------------------------------------------------\n"
            "Orthogonality func initial: " << opt_orto << "\n\n"
            "Skewness: func initial " << opt_skew << "\n\n"
            "Eccentricity: func initial " << opt_ecce << "\n\n"
            "Uniformity: func initial " << opt_unif << "\n\n"
            "Length: func initial " << opt_length << "\n\n"
            "Area: func initial " << opt_area << "\n\n"
            "All: func initial " << opt_All << "\n\n"

                      "------------------------------------------------------------\n\n";

            std::string output("results_optimized_initial" + profile);
            gsVector<> JacResult(3);
            JacResult.setZero(3);
            JacResult = returnJacobianDeterminant(patchInitial.patch(0), jacPts, true, output, 0);


            //save initial
            gsFileData<> fileData;
            fileData << patchInitial.patch(0);
            std::string out;
            out = "optimize_initial_blade_"+profile+".xml";
            fileData.dump(out);
            gsMesh<> mesh;
            makeMesh(patchInitial.patch(0).basis(), mesh, 10);
            patchInitial.patch(0).evaluateMesh(mesh);
            out = "optimize_initial_blade_"+profile+"_Mesh";
            gsWriteParaview(mesh, out);
            gsMatrix<> coefs = patchInitial.patch(0).coefs();
            coefs.transposeInPlace();
            out = "optimize_initial_blade_"+profile+"_ControlPoints";
            gsWriteParaviewPoints(coefs, out);



            optimizeCycle(iterations, num_weights, num_func, num_of_func_comb, jacPts, filename, patchInitial);

            // close the output file

        }

        //====================READ data===================================
        if(readData)
        {
           //gsMultiPatch<> patches = patchInitial;
           gsMatrix<> data;
           data = readFile(filenamein);
           gsInfo << data;
           int a;
           std::cin >> a;




           for (int i = 0; i < data.rows(); i++)
           {
                 gsMultiPatch<> patches = patchInitial;
                 gsQualityMeasure<> optimization(patches.patch(0));

                 int count = data(i,0);

                 int iterations = data(i,1);

                 real_t orthogonality = data(i,2);

                 real_t skewness = data(i,3);

                 real_t eccentricity = data(i,4);

                 real_t uniformity = data(i,5);

                 real_t area = data(i,7);

                 real_t length = data(i,6);



                 int it = 1;
                 while (it != iterations + 1)
                 {


                      gsInfo<< "\n\n------------------------------------------------------------\n"
                      "Possibility: " << count << "\n\n"
                               "Iteration of all: " << it << "/" << iterations << "\n\n"
                      "Orthogonality: " << orthogonality << "\n\n"
                      "Skewness: " << skewness << "\n\n"
                      "Eccentricity: " << eccentricity << "\n\n"
                      "Uniformity: " << uniformity << "\n\n"
                      "Length: " << length << "\n\n"
                      "Area: " << area << "\n\n"
                      "Intersection: " << intersection << "\n\n"
                      "Epsilon: " << epsilon << "\n\n"
                                "------------------------------------------------------------\n\n";


                      optimization.optimize(orthogonality, skewness,
                                            eccentricity, uniformity,
                                            length, area,
                                            intersection, epsilon);

                      if(it + 1 == iterations + 1)
                      {
                      std::string output2("optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it));
                      gsVector<> JacResult(3);
                      JacResult=returnJacobianDeterminant(patches.patch(0), jacPts, false, output2, it);

                      gsFileData<> fileData;
                      fileData << patches.patch(0);
                      std::string out;
                      out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+".xml";
                      fileData.dump(out);
                      gsMesh<> mesh;
                      makeMesh(patches.patch(0).basis(), mesh, 10);
                      patches.patch(0).evaluateMesh(mesh);
                      out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"_Mesh";
                      gsWriteParaview(mesh, out);
                      gsMatrix<> coefs = patches.patch(0).coefs();
                      coefs.transposeInPlace();
                      out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"_ControlPoints";
                      gsWriteParaviewPoints(coefs, out);
                      gsMesh<> mesh2;
                      patches.patch(0).controlNet(mesh2);
                      out = "optimize_blade"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"ControlNet";
                      gsWriteParaview(mesh2,out);
                      }

                      it++;

                }
            }
        }


            //=======================================TEST SEVERAL OPTIONS==================================================================


//             iterations = 1;
//             orthogonality = 0.0;
//            skewness = 0.0;
//             eccentricity = 0.0;
//            intersection = 0.0;
//            uniformity = 0.0;
//             area =  1000000000000;
//            length = 10000000000;
//            epsilon = 1e-7;
//            jacPts = 10000;
//             dumped = false;
//            int count = 84351;
//            int it = 1;


//            while (it != iterations + 1)
//            {


//                 gsInfo<< "\n\n------------------------------------------------------------\n"
//                 "Orthogonality: " << orthogonality << "\n\n"
//                 "Skewness: " << skewness << "\n\n"
//                 "Eccentricity: " << eccentricity << "\n\n"
//                 "Uniformity: " << uniformity << "\n\n"
//                 "Length: " << length << "\n\n"
//                 "Area: " << area << "\n\n"
//                 "Intersection: " << intersection << "\n\n"
//                 "Epsilon: " << epsilon << "\n\n"
//                           "------------------------------------------------------------\n\n";

//                 optimizationIt.optimize(orthogonality, skewness,
//                                       eccentricity, uniformity,
//                                       length, area,
//                                       intersection, epsilon);


//                 std::string output2("optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it));
//                 gsVector<> JacResult(3);
//                 JacResult=returnJacobianDeterminant(patchInitial.patch(0), jacPts, false, output2, it);

//                 gsFileData<> fileData;
//                 fileData << patchInitial.patch(0);
//                 std::string out;
//                 out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+".xml";
//                 fileData.dump(out);
//                 gsMesh<> mesh;
//                 makeMesh(patchInitial.patch(0).basis(), mesh, 10);
//                 patchInitial.patch(0).evaluateMesh(mesh);
//                 out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"_Mesh";
//                 gsWriteParaview(mesh, out);
//                 gsMatrix<> coefs = patchInitial.patch(0).coefs();
//                 coefs.transposeInPlace();
//                 out = "optimize_blade_"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"_ControlPoints";
//                 gsWriteParaviewPoints(coefs, out);
//                 gsMesh<> mesh2;
//                 patchInitial.patch(0).controlNet(mesh2);
//                 out = "optimize_blade"+profile+"_possibility_"+ util::to_string(count)+"_it_"+util::to_string(it)+"ControlNet";
//                 gsWriteParaview(mesh2,out);

//                 it++;
//            }

    }

    //========================OLD OPTIMIZATION========================================================

    //                    if(opt_new > 10000.0)
    //                    {
    //                        it = iterations;
    //                        continue;
    //                    }

    //                    std::ostringstream strs_o;
    //                    std::ostringstream strs_e;
    //                    std::ostringstream strs_a;
    //                      std::string output2;

    //                    strs_o << j;
    //                    strs_e << it;
    //                    strs_a << iop;
    //                    output2 = "optimize_blade" + strs_a.str()+ "_" +strs_o.str()+"_" +strs_e.str();
    //                     //saveData2(*geom, output2, 0);


    //    gsWriteParaview(patches,"patches");
    //    system("paraview patches.pvd");
    //    return 0;

    //    std::string filename("Profile2DInRectangular.xml");
    //    gsFileData<> write;
    //    write << patches;
    //    write.dump(filename);
    //    std::string input("Profile2DInRectangular.xml");

    //    gsFileData<> data(input);
    //    gsGeometry<>::uPtr geom;
    //    if (data.has< gsGeometry<> >())
    //    {
    //        geom = data.getFirst< gsGeometry<> >();
    //   }
    //    if (!geom)
    //   {
    //        gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
    //    }



    //    gsInfo << "Value of functional: "
    //           << opt_old0
    //           << "\n";



    //        gsVector<> possibillities(6);
    //        possibillities.setZero(6);
    //        possibillities << 7,2,4,3,2,2;

    //        real_t comb_of_poss_num = 0;
    //        for(int i=0; i<possibillities.rows(); i++)
    //        {
    //            comb_of_poss_num += possibillities(i);
    //        }

    //        gsInfo << "\n comb_of_poss_num: " << comb_of_poss_num << " \n";
    //        gsMatrix<> index_mat(comb_of_poss_num,2);
    //        index_mat.setZero(comb_of_poss_num,2);



    //    //combinations given by brute force
    //    index_mat<< 9,2,
    //                19,1,
    //                39,2,
    //                81,5,
    //                119,1,
    //                123,4,
    //                2626,100,
    //                6,2,
    //                1376,2,
    //                6,2,
    //                276,3,
    //                401,3,
    //                1376,2,
    //                5,2,
    //                525,8,
    //                900,7,
    //                5,2,
    //                900,7,
    //                5,2,
    //                425,13;



    //    real_t comb_of_poss_k = 0;
    //    for(int i=0; i<iop-1; i++)
    //    {
    //        comb_of_poss_k += possibillities(i);
    //    }
    //    for(unsigned k=comb_of_poss_k;k<comb_of_poss_k+possibillities(iop-1);k++)
    //    {
    //        gsInfo << "\n k index: " << k << "\n";
    //        int j = index_mat(k,0);
    //        int iterations = index_mat(k,1);

    //        if (data.has< gsGeometry<> >())
    //        {
    //            geom = data.getFirst< gsGeometry<> >();
    //        }
    //        int it = 0;
    //        //opt_old=opt_old0;
    //        while (it != iterations )
    //        {
    //             gsInfo << "iteration: " << it << " / " << iterations << "\n";
    //             gsInfo << "j: " << j << " / " << "combination" << "\n";

    //             orthogonality =comb(j,0);
    //             skewness =  comb(j,1);
    //             eccentricity = comb(j,2);
    //             uniformity = comb(j,3);
    //                area = comb(j,4);
    //                length = 0.0;
    //                gsInfo<< "\n\n------------------------------------------------------------\n"
    //                "Orthogonality: " << orthogonality << "\n\n"
    //                "Skewness: " << skewness << "\n\n"
    //                "Eccentricity: " << eccentricity << "\n\n"
    //                "Uniformity: " << uniformity << "\n\n"
    //                "Length: " << length << "\n\n"
    //                "Area: " << area << "\n\n"
    //                "Intersection: " << intersection << "\n\n"
    //                "Epsilon: " << epsilon << "\n\n"
    //                          "------------------------------------------------------------\n\n";
    //                gsQualityMeasure<> optimization(*geom);
    //                optimization.optimize(orthogonality, skewness,
    //                                      eccentricity, uniformity,
    //                                      length, area,
    //                                      intersection, epsilon);
    //                real_t opt_new = optimization.functional(orthogonality, skewness,
    //                                                             eccentricity, uniformity, length,
    //                                                             area, intersection, epsilon);


    //                gsInfo << "Value of functional: "
    //                          <<  opt_new
    //                          << "\n";
    //                //std::cin.get();
    //                if(opt_new > 10000.0)
    //                {
    //                    it = iterations;
    //                    continue;
    //                }

    //                std::ostringstream strs_o;
    //                std::ostringstream strs_e;
    //                std::ostringstream strs_a;
    //                  std::string output2;

    //                strs_o << j;
    //                strs_e << it;
    //                strs_a << iop;
    //                output2 = "optimize_blade" + strs_a.str()+ "_" +strs_o.str()+"_" +strs_e.str();
    //                 //saveData2(*geom, output2, 0);



    //                gsInfo<<output2<<"\n";
    //                JacResult=returnJacobianDeterminant(*geom, jacPts, false, output2, it);
    //                gsInfo << JacResult;
    //                //std::cin.get();


    //                        gsFileData<> fileData;
    //                        fileData << *geom;
    //                        std::string out = output2 + ".xml";
    //                        fileData.dump(out);
    //                        gsMesh<> mesh;
    //                        makeMesh(geom->basis(), mesh, 10);
    //                        geom->evaluateMesh(mesh);
    //                        out = output2 + "Mesh";
    //                        gsWriteParaview(mesh, out);
    //                        gsMatrix<> coefs = geom->coefs();
    //                        coefs.transposeInPlace();
    //                        out = output2 + "ControlPoints";
    //                        gsWriteParaviewPoints(coefs, out);
    //                        gsMesh<> mesh2;
    //                        geom->controlNet(mesh2);
    //                        out = output2 + "ControlNet";
    //                        gsWriteParaview(mesh2,out);



    //                it++;


    //                }
    //    }
    //  }
    //    real_t opt_new;
    //    real_t opt_old;

    //    //choose num of weights that change
    //    unsigned num_weights = 5;
    //    gsVector<> weights(num_weights);
    //    weights.setZero(num_weights);

    //    for (int i = 0; i < weights.rows(); i++)
    //    {
    //        weights(i) = i*1.0/(num_weights-1); //(uniform 0..1)
    //    }

    //    gsInfo<< "\n----------------------\n";
    //    gsInfo<< "\n weights \n";
    //    gsInfo<< weights;
    //    gsInfo<< "\n----------------------\n";
    //    //matrix of all permuattions with repetations
    //    gsMatrix<> comb(math::pow(num_weights,num_weights),num_weights);
    //    comb.setZero(math::pow(num_weights,num_weights),num_weights);
    //    comb = combination(weights);

    //    gsVector<> min(4);
    //    min.setZero(4);
    //    real_t minimum = opt_old0;
    //    min(1)=minimum;

    //    //compute geometries for all given combinations
    //    for(unsigned j=1;j<math::pow(5,5);j++)
    //    {
    //        if (data.has< gsGeometry<> >())
    //        {
    //            geom = data.getFirst< gsGeometry<> >();
    //        }
    //        int it = 1;
    //        opt_old=opt_old0;
    //        while (it != iterations )
    //        {
    //            gsInfo << "iteration: " << it << " / " << iterations << "\n";
    //            gsInfo << "j: " << j << " / " << "combination" << "\n";

    //            orthogonality =comb(j,0);
    //            skewness =  comb(j,1);
    //            eccentricity = comb(j,2);
    //            uniformity = comb(j,3);
    //            area = comb(j,4);
    //            length = 0.0;
    //            gsInfo<< "\n\n------------------------------------------------------------\n"
    //            "Orthogonality: " << orthogonality << "\n\n"
    //            "Skewness: " << skewness << "\n\n"
    //            "Eccentricity: " << eccentricity << "\n\n"
    //            "Uniformity: " << uniformity << "\n\n"
    //            "Length: " << length << "\n\n"
    //            "Area: " << area << "\n\n"
    //            "Intersection: " << intersection << "\n\n"
    //            "Epsilon: " << epsilon << "\n\n"
    //                      "------------------------------------------------------------\n\n";
    //            gsQualityMeasure<> optimization(*geom);
    //            optimization.optimize(orthogonality, skewness,
    //                                  eccentricity, uniformity,
    //                                  length, area,
    //                                  intersection, epsilon);
    //            opt_new = optimization.functional(orthogonality, skewness,
    //                                                         eccentricity, uniformity, length,
    //                                                         area, intersection, epsilon);


    //            gsInfo << "Value of functional: "
    //                      <<  opt_new
    //                      << "\n";
    //            //std::cin.get();
    //            if(opt_new > 10000.0)
    //            {
    //                it = iterations;
    //                continue;
    //            }

    //            std::ostringstream strs_o;
    //            std::ostringstream strs_e;
    //            std::ostringstream strs_a;
    //              std::string output2;

    //            strs_o << j;
    //            strs_e << it;
    //            strs_a << iop;
    //            output2 = "optimize_blade" + strs_a.str()+ "_" +strs_o.str()+"_" +strs_e.str();
    //             //saveData2(*geom, output2, 0);



    //            gsInfo<<output2<<"\n";
    //            JacResult=returnJacobianDeterminant(*geom, jacPts, false, output2, it);
    //            gsInfo << JacResult;
    //            //std::cin.get();
    //            if(JacResult(0)!=0 && JacResult(1)!=0)
    //            {
    //                    it++;
    //                    opt_old = opt_new;
    //                    continue;
    //            }

    //            if (opt_new/opt_old < 0.98)
    //            {
    //                   it++;
    //                   opt_old = opt_new;
    //                   continue;
    //            }
    //            else
    //            {
    //                if(opt_new < min(1))
    //                {
    //                    min(0)=j;
    //                    min(1)=opt_new;
    //                    min(2)=JacResult(0);
    //                    min(3)=JacResult(1);

    //                    gsFileData<> fileData;
    //                    fileData << *geom;
    //                    std::string out = output2 + ".xml";
    //                    fileData.dump(out);
    //                    gsMesh<> mesh;
    //                    makeMesh(geom->basis(), mesh, 10);
    //                    geom->evaluateMesh(mesh);
    //                    out = output2 + "Mesh";
    //                    gsWriteParaview(mesh, out);
    //                    gsMatrix<> coefs = geom->coefs();
    //                    coefs.transposeInPlace();
    //                    out = output2 + "ControlPoints";
    //                    gsWriteParaviewPoints(coefs, out);
    //                    gsMesh<> mesh2;
    //                    geom->controlNet(mesh2);
    //                    out = output2 + "ControlNet";
    //                    gsWriteParaview(mesh2,out);


    //                }



    //                it=iterations;

    //            }

    //            }
    //}

    //    final_weights(iop,0)=min(0);
    //    final_weights(iop,1)=min(1);
    //    final_weights(iop,2)=min(2);
    //    final_weights(iop,3)=min(3);
    //    gsInfo<<final_weights;
    //    }//endfor iop

    //    gsInfo<<final_weights;
    //    return 0;
    //}

    return 0;
}

//=================================================================================
//------------------functions-----------------------------------------

template<class TT>
void permutate (const gsVector<TT>& input,
                gsMatrix<TT>& permutations,
                gsVector<TT>& index,
                int depth,
                int & count, int length)
{
    if (depth == length)
    {
        for (int i = 0; i < length; ++i)
        {
            permutations(count,i) = input(index(i));
        }
        ++count;
        return;
    }

    for (int i = 0; i < input.rows(); ++i)
    {
        index(depth) = i;
        permutate(input, permutations, index, depth + 1, count, length);
    }
}

//save k_length permutations with repetitions of the given vector input to the matrix
template<class TT>
gsMatrix<TT> permutation(gsVector<TT> input, int length)
{
    gsMatrix<TT> permutation(math::pow(input.rows(),length),length);
    permutation.setZero(math::pow(input.rows(),length),length);
    gsVector<TT> index(length);
    int count = 0;
    int depth = 0;
    permutate(input, permutation, index, depth, count, length);
    //gsInfo << "\nTotal combinations with repetitions: " << count;
    return permutation;
}

/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Staring and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
template<class TT>
void combine(gsVector<>& input, gsMatrix<>& combinations, gsVector<TT>& data, int start, int end,
                     int index, int & row, int length)
{
    if (index == length)
    {
        for (int j=0; j<length; ++j)
        {
            combinations(row, j) = data(j);
        }
        ++row;
        return;
    }

    for (int i = start; i <= end && end-i+1 >= length-index; ++i)
    {
        data(index) = input(i);
        combine(input, combinations, data, i+1, end, index+1, row, length);
    }
}


//save k_length permutations with repetitions of the given vector input to the matrix
template<class TT>
gsMatrix<TT> combination(gsVector<TT> input, int length)
{
    gsMatrix<TT> combination(factorial(input.rows())/(factorial(length)*factorial(input.rows()-length)),length);
    combination.setZero(factorial(input.rows())/(factorial(length)*factorial(input.rows()-length)),length);
    gsVector<TT> data(length);
    int start = 0;
    int end = input.rows()-1;
    int row = 0;
    int index = 0;
    combine(input, combination, data, start, end, index, row, length);
    //gsInfo << "\nTotal combinations with repetitions: " << count;
    return combination;
}

//save the plus, negative and zero number of points in gsVector
template<class TT>
gsVector<TT> returnJacobianDeterminant(const gsGeometry<TT>& geom, const int points,
                              const bool savePoints,
                              const std::string& output,
                              const int number)
{
    gsInfo << "Checking Jacobian determinant ..." << "\n";
    gsMatrix<TT> para  = geom.support();
    gsVector<TT> c0 = para.col(0);
    gsVector<TT> c1 = para.col(1);
    gsMatrix<TT> pts = uniformPointGrid(c0, c1, points);
    gsMatrix<TT> plus(pts.rows(), pts.cols());
    plus.setZero(pts.rows(), pts.cols());
    gsMatrix<TT> minus(pts.rows(), pts.cols());
    minus.setZero(pts.rows(), pts.cols());
    gsMatrix<TT> zero(pts.rows(), pts.cols());
    zero.setZero(pts.rows(), pts.cols());
    int plusCounter = 0;
    int minusCounter = 0;
    int zeroCounter = 0;

    // calculating the determinant of the Jacobian
    for (int col = 0; col < pts.cols(); col++)
    {
        const real_t determinant = geom.jacobian(pts.col(col)).determinant();
        gsMatrix<TT> tmp = geom.eval(pts.col(col));
        if (determinant < 0)
        {
            minus.col(minusCounter) = tmp.col(0);
            minusCounter++;
        }
        else if (determinant > 0)
        {
            plus.col(plusCounter) = tmp.col(0);
            plusCounter++;
        }
        else
        {
            zero.col(zeroCounter) = tmp.col(0);
            zeroCounter++;
        }
    }
    gsInfo << "Number of points with \n"
              << "  - positive sign: " << plusCounter << "\n"
              << "  - negative sign: " << minusCounter << "\n"
              << "  - zero sign: " << zeroCounter << "\n" << "\n";
    if (savePoints)
    {
        if (0 < plusCounter)
        {
            std::string out = output + "PositivePoints" + util::to_string(number);
            gsWriteParaviewPoints(plus, out);
        }
        if (0 < minusCounter)
        {
            std::string out = output + "NegativePoints" + util::to_string(number);
            gsWriteParaviewPoints(minus, out);
        }
        if (0 < zeroCounter)
        {
            std::string out = output + "ZeroPoints" + util::to_string(number);
            gsWriteParaviewPoints(zero, out);
        }
    }
    gsVector<TT> result(3);
    result << plusCounter,minusCounter,zeroCounter;
    return result;
}

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
                                       TT const & rotationCenterY)
{
    //----------------set parameters for blade profile----------------
    int num_samples = 30;
    gsVector<TT> vec(2);
    vec (0) = rotationCenterX;
    vec (1) = rotationCenterY;
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
     //---------------transform given profile----------------------
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

    //---------------set parameters for boundary of patches-----------------------
    real_t width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

    //control points of pressure side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bp(i, 0) = bp.coef(i, 0) ;
        cp_bp(i, 1) = bp.coef(i, 1) + pitch;
    }
    //control points of suction side
    for(unsigned i = 0; i < num_cpblade; i++){
        cp_bs(i, 0) = bs.coef(i, 0);
        cp_bs(i, 1) = bs.coef(i, 1);
    }

    gsMultiPatch<TT> mp;
    //initial data for patches without blade profile - patches are linear -> elevate
    gsKnotVector<TT> kv_u(0, 1, 0, 2);
    gsKnotVector<TT> kv_v(0, 1, 0, 2);
    gsTensorBSplineBasis<2, real_t> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, real_t> basis_blade(kvfit, kv_v);



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





      for (unsigned int i = 0; i<4;i++){
        a_cp(i,0)=coef_patchStart(i,0);
        a_cp(i,1)=coef_patchStart(i,1);
        b_cp(i,0)=coef_patchStart(i,0);
        b_cp(i,1)=coef_patchStart(i,1)+pitch;
     }
 /*   for (unsigned int i = 0; i<4;i++){
        a_cp(i,0)=patch0.coef(i,0);
        a_cp(i,1)=patch0.coef(i,1);
        b_cp(i,0)=patch1.coef(i+12,0);
        b_cp(i,1)=patch1.coef(i+12,1)+pitch;
    }
*/
    for (unsigned int i = 4; i<10;i++){
        a_cp(i,0)=cp_bs(i-3,0);
        a_cp(i,1)=cp_bs(i-3,1);
        b_cp(i,0)=cp_bp(i-3,0);
        b_cp(i,1)=cp_bp(i-3,1);
    }
     /*
    for (unsigned int i = 10; i<14;i++){
        a_cp(i,0)=patch4.coef(i-10,0);
        a_cp(i,1)=patch4.coef(i-10,1);
        b_cp(i,0)=patch5.coef(i+12-10,0);
        b_cp(i,1)=patch5.coef(i+12-10,1)+pitch;
    }
    */
    real_t yend_coor = -((cp_bs(0,1)*cp_bs(7,0) - cp_bs(0,0)*cp_bs(7,1) - cp_bs(0,1)*length_x2 + cp_bs(7,1)*length_x2)/(
                          cp_bs(0,0) - cp_bs(7,0)));
    gsMatrix<TT> coef_patchEnd (4, 2);
    coef_patchEnd << cp_bs(7,0),cp_bs(7,1),
                     length_x2/4.0 + 3.0*cp_bs(7,0)/4.0, yend_coor/4.0 + 3.0*cp_bs(7,1)/4.0,
                     3.0*length_x2/4.0 + cp_bs(7,0)/4.0, 3.0*yend_coor/4.0 + cp_bs(7,0)/4.0,
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
    boundaries4->addPatch(gsBSpline<>( kv_vv, c_cp));
    boundaries4->addPatch(gsBSpline<>( kv_uu, a_cp));
    boundaries4->addPatch(gsBSpline<>( kv_vv, d_cp));
    boundaries4->addPatch(gsBSpline<>( kv_uu, b_cp));
    gsCoonsPatch<real_t> patchAll = coonsPatch(*boundaries4);
    patchAll.compute();
    mp.addPatch(patchAll.result());
    bool plotMeshes = true;
    if (plotMeshes) {
        gsMesh<> meshAll;
        patchAll.result().controlNet(meshAll);
        gsWriteParaview(meshAll,"ProfilePatchAll");
   }
    return mp;
}


template<class TT>
gsMultiPatch<TT> BSplineProfile2DBetweenPatchStart(TT const & index,
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
                                                   bool & refine)
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



    gsBSpline<TT> boundaryc( kv_vv, c_cp);
    gsBSpline<TT> boundarya( kv_uu, a_cp);
    gsBSpline<TT> boundaryd( kv_vv, d_cp);
    gsBSpline<TT> boundaryb( kv_uu, b_cp);

    if(refine)
    {
        for(real_t knot = 0.2; knot < 1.0; knot += 0.2)
        {
            boundaryc.insertKnot(knot);
            boundaryd.insertKnot(knot);
        }

        for(real_t knot = 0.2; knot < 1.0; knot += 0.2)
        {
            boundarya.insertKnot(knot*(0.2/3.0));
            boundarya.insertKnot(1-(knot*(0.2/3.0)));
            boundaryb.insertKnot(knot*(0.2/3.0));
            boundaryb.insertKnot(1-(knot*(0.2/3.0)));
        }
    }

    boundaries4->addPatch(boundaryc);
    boundaries4->addPatch(boundarya);
    boundaries4->addPatch(boundaryd);
    boundaries4->addPatch(boundaryb);
    gsCoonsPatch<TT> patchAll = coonsPatch(*boundaries4);
    patchAll.compute();


    mp.addPatch(patchAll.result());



//   //=================================optimization===========================================

//    gsVector<TT> area_vec(7);
//    area_vec.setZero(7);
//    area_vec << 0,1,0.25,0.25,0,0,0;

//    real_t orthogonality = 0.0;
//    real_t skewness = 0.0;
//    real_t eccentricity = 0.0;
//    real_t intersection = 0.0;
//    real_t uniformity = 0.25;
//    real_t area = area_vec(index);
//    real_t length = 0;
//    real_t epsilon = 1e-7;

//    gsQualityMeasure<TT> optimization(mp.patch(0));
//    real_t opt_val = optimization.functional(orthogonality, skewness,
//                                             eccentricity, uniformity,
//                                             length, area,
//                                             intersection, epsilon);
//    optimization.optimize(orthogonality, skewness,
//                          eccentricity, uniformity,
//                          length, area,
//                          intersection, epsilon);
//    gsInfo << "Value of functional: "
//           << opt_val
//           << "\n";

//    if(plot)
//    {
//        gsFileData<> fileData;
//        fileData << mp.patch(0);
//        std::string out;
//        out = "optimize_blade"+ util::to_string(index) +".xml";
//        fileData.dump(out);
//        gsMesh<TT> mesh;
//        makeMesh(mp.patch(0).basis(), mesh, 10);
//        mp.patch(0).evaluateMesh(mesh);
//        out = "optimize_bladeMesh" +  util::to_string(index) ;
//        gsWriteParaview(mesh, out);
//        gsMesh<TT> mesh2;
//        mp.patch(0).controlNet(mesh2);
//        out = "optimize_bladeControlNet" +  util::to_string(index);
//        gsWriteParaview(mesh2,out);
//    }

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

    return mp;
}

void optimizeCycle(int iterations, real_t num_weights, real_t num_func, int num_of_func_comb, int jacPts, std::string filename, gsMultiPatch<> geom)
{
    std::ofstream outputFile(filename);

    // write the file headers
    outputFile << "number" << ","
               << "iteration" << ","
               << "w_orthogonality" << ","
               << "w_skewness" << ","
               << "w_eccentricity" << ","
               << "w_uniformity" << ","
               << "w_length" << ","
               << "w_area" << ","
               << "func_orthogonality" << ","
               << "func_skewness" << ","
               << "func_eccentricity" << ","
               << "func_uniformity" << ","
               << "func_length" << ","
               << "func_area" << ","
               << "func_All" << ","
               << "Jacobi plus" << ","
               << "Jacobi minus" << ","
               << "Jacobi zero" <<"\n" ;


    real_t intersection = 0.0;
    real_t epsilon = 1e-7;

    gsVector<> weights(num_weights);
    weights.setZero(num_weights);
    for (int i = 0; i < num_weights; i++)
    {
       weights(i) = math::pow(10,-math::floor(num_weights/2.0)+i); //(uniform 0..1)
    }
    gsInfo<< "\n----------------------\n";
    gsInfo<< "\n weights \n";
    gsInfo<< weights;
    gsInfo<< "\n----------------------\n";

    //matrix of all comb of weights

    int num_weights_comb_all = math::pow(num_weights,num_of_func_comb);
    gsMatrix<> weight_comb(num_weights_comb_all,num_of_func_comb);
    weight_comb.setZero(num_weights_comb_all,num_of_func_comb);

    weight_comb = permutation(weights,num_of_func_comb);

    gsInfo<< "\n----------------------\n";
    gsInfo<< "\n weights_comb \n";
    gsInfo<< weight_comb;
    gsInfo<< "\n----------------------\n";

    gsVector<> index_func(num_func);
    for(int i = 0; i < num_func; i++)
    {
        index_func(i) = i;
    }
    int num_func_comb_all = factorial(num_func)/(factorial(num_func-num_of_func_comb)*factorial(num_of_func_comb));
    gsMatrix<> functionals_comb(num_func_comb_all,num_of_func_comb);
    functionals_comb.setZero(num_func_comb_all,num_of_func_comb);
    gsInfo<< num_func_comb_all << "\n";

    functionals_comb = combination(index_func,num_of_func_comb);

    gsInfo<< "\n----------------------\n";
    gsInfo<< "\n functionals_comb \n";
    gsInfo<< functionals_comb;
    gsInfo<< "\n----------------------\n";

    gsVector<> comb(num_func);
    int number = 1;

    //optimization for given func comb and weights
    for (unsigned f = 0; f < num_func_comb_all; f++)
    {
        for (unsigned w = 0; w < num_weights_comb_all; w++)
        {
           gsQualityMeasure<> optimization(geom.patch(0));
           comb.setZero(num_func);
           for (int k = 0; k < num_of_func_comb; k++)
           {
                comb(functionals_comb(f,k)) = weight_comb(w,k);
           }
           int iteration = 1;
           int it = iteration;
           while (it != iterations)
           {
                number++;
                gsInfo << "iteration: " << it << " / " << iterations << "\n";
                gsInfo << "f: " << f << " / " << num_func_comb_all << "\n";
                gsInfo << "w: " << w << " / " <<  num_weights_comb_all << "\n";
                gsInfo << "comb: " << comb << " / " << "combination" << "\n";

                real_t orthogonality =comb(0);
                real_t skewness =  comb(1);
                real_t eccentricity = comb(2);
                real_t uniformity = comb(3);
                real_t area = comb(4);
                real_t length = comb(5);
                gsInfo<< "\n\n------------------------------------------------------------\n"
                "Orthogonality: " << orthogonality << "\n\n"
                "Skewness: " << skewness << "\n\n"
                "Eccentricity: " << eccentricity << "\n\n"
                "Uniformity: " << uniformity << "\n\n"
                "Area: " << area << "\n\n"
                "Length: " << length << "\n\n"

                          "------------------------------------------------------------\n\n";

                optimization.optimize(orthogonality, skewness,
                                      eccentricity, uniformity,
                                      length, area,
                                      intersection, epsilon);

                real_t opt_orto = optimization.functional(1, 0,
                                                         0, 0,
                                                         0, 0,
                                                         0, 0);
                real_t opt_skew = optimization.functional(0, 1,
                                                         0, 0,
                                                         0, 0,
                                                         0, 0);
                real_t opt_ecce = optimization.functional(0, 0,
                                                         1, 0,
                                                         0, 0,
                                                         0, 0);
                real_t opt_unif = optimization.functional(0, 0,
                                                         0, 1,
                                                         0, 0,
                                                         0, 0);
                real_t opt_length = optimization.functional(0, 0,
                                                         0, 0,
                                                         1, 0,
                                                         0, 0);
                real_t opt_area = optimization.functional(0, 0,
                                                         0, 0,
                                                         0, 1,
                                                         0, 0);
                real_t opt_All = optimization.functional(orthogonality, skewness,
                                      eccentricity, uniformity,
                                      length, area,
                                      intersection, epsilon);

                std::string output2("pom");
                gsVector<> JacResult(3);
                JacResult.setZero(3);
                JacResult=returnJacobianDeterminant(geom.patch(0), jacPts, false, output2, it);


                outputFile << number << ","
                           << it << ","
                           << orthogonality << ","
                           << skewness << ","
                           << eccentricity << ","
                           << uniformity << ","
                           << length << ","
                           << area << ","
                           << opt_orto << ","
                           << opt_skew << ","
                           << opt_ecce << ","
                           << opt_unif << ","
                           << opt_length << ","
                           << opt_area << ","
                           << opt_All << ","
                           << JacResult(0) << ","
                           << JacResult(1) << ","
                           << JacResult(2) <<"\n" ;


                it++;

            }
        }
    }


    return;

}



gsMatrix<> readFile(std::string inputFileName)
{
    std::vector<std::vector<double> > data;
    std::ifstream inputFile(inputFileName);
    int l = 0;
    std::string head;
    getline(inputFile, head);
    while (inputFile)
    {
        l++;
        std::string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#')
        {
           std::istringstream ss(s);
           std::vector<double> record;
           while (ss)
           {
                std::string line;
                if (!getline(ss, line, ','))
                   break;
                try
                {
                   record.push_back(std::stof(line));
                }
                catch (const std::invalid_argument e)
                {
                   std::cout << "NaN found in file " << inputFileName << " line " << l
                        << "\n";
                   e.what();
                }
           }
           data.push_back(record);
       }
   }

   if (!inputFile.eof()) {
       std::cerr << "Could not read file " << inputFileName << "\n";
       //__throw_invalid_argument("File not found.");
   }

   gsMatrix<> dataOut(data.size(),data[0].size());  //all rows has same num of columns
   dataOut.setZero(data.size(),data[0].size());
   for (size_t i=0; i<data.size(); ++i)
   {
       for (size_t j=0; j<data[i].size(); ++j)
       {
           dataOut(i,j) = data[i][j];
       }

   }

    return dataOut;
}


template<class TT>
gsMultiPatch<TT> DomainBetweenBladeProfiles1(TT const & index,
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
    TT const & uniformity_param,
    std::vector<TT> const & kvfit_knots,
    bool const & coarse,
    gsVector<TT> const & geom_Params, bool & refine)
{


    //----------------set parameters for blade profile----------------
    bool plot = true;
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
    //gsKnotVector<TT> kvfit(0, 1, 4, 4);
    //std::vector<TT> kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};
    gsKnotVector<TT> kvfit(kvfit_knots);
    gsKnotVector<TT> kvlinear(0,1,0, 2);
    unsigned num_cpblade = kvfit.size() - kvfit.degree() - 1;
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
    //TT width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

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

    gsInfo << "\n pressure side \n";
    gsInfo << cp_bp;
    gsInfo << "\n pressure side pom \n";
    //gsInfo << cp_bp_pom;
    gsInfo << "\n suction side \n";
    gsInfo << cp_bs;


    gsMultiPatch<TT> mp;
    // compute discrete coons patch to optimize
    gsMatrix<TT> coef_patchAll((kvfit.degree() + 1)*(num_cpblade + kvfit.degree() + kvfit.degree()), 2);
    coef_patchAll.setZero((kvfit.degree() + 1)*(num_cpblade + kvfit.degree() + kvfit.degree()), 2);
    gsMatrix<TT> a_cp(num_cpblade + kvfit.degree() + kvfit.degree(), 2);
    gsMatrix<TT> b_cp(num_cpblade + kvfit.degree() + kvfit.degree(), 2);
    gsMatrix<TT> c_cp(2, 2);
    gsMatrix<TT> d_cp(2, 2);
   a_cp.setZero(num_cpblade + kvfit.degree() + kvfit.degree(), 2);
    b_cp.setZero(num_cpblade + kvfit.degree() + kvfit.degree(), 2);
    c_cp.setZero(2, 2);
    d_cp.setZero(2, 2);

    int cpk = num_cpblade - 1;//last index of cp profile


    TT ystart_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x1 + cp_bs(cpk, 1)*length_x1) / (
        cp_bs(0, 0) - cp_bs(cpk, 0)));

    gsMatrix<TT> coef_patchStart(2, 2);
    coef_patchStart << length_x1, ystart_coor,
                cp_bs(0, 0), cp_bs(0, 1);

    gsBSpline<TT> patchStart(kvlinear, coef_patchStart);
    patchStart.degreeElevate(kvfit.degree()-1);



    for (int i = 0; i<kvfit.degree(); i++)
    {
        a_cp(i, 0) = patchStart.coef(i, 0);
        a_cp(i, 1) = patchStart.coef(i, 1);
        b_cp(i, 0) = patchStart.coef(i, 0);
        b_cp(i, 1) = patchStart.coef(i, 1) + pitch;
    }

    gsInfo << "\n --------- \n";
    gsInfo << a_cp << "\n";
    gsInfo << b_cp << "\n";
    gsInfo << "\n -------- \n";

    for (unsigned int i = kvfit.degree(); i< num_cpblade + kvfit.degree() ; i++)
    {
        a_cp(i, 0) = cp_bs(i - kvfit.degree(), 0);
        a_cp(i, 1) = cp_bs(i - kvfit.degree(), 1);
        b_cp(i, 0) = cp_bp(i - kvfit.degree(), 0);
        b_cp(i, 1) = cp_bp(i - kvfit.degree(), 1);
    }

    gsInfo << "\n --------- \n";
    gsInfo << a_cp << "\n";
    gsInfo << b_cp << "\n";
    gsInfo << "\n -------- \n";
    TT yend_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x2 + cp_bs(cpk, 1)*length_x2) / (
        cp_bs(0, 0) - cp_bs(cpk, 0)));
    gsMatrix<TT> coef_patchEnd(2, 2);
    coef_patchEnd << cp_bs(cpk, 0), cp_bs(cpk, 1),
                length_x2, yend_coor;

     gsInfo << "\n endpatch \n";

    gsBSpline<TT> patchEnd(kvlinear, coef_patchEnd);
    patchEnd.degreeElevate(kvfit.degree()-1);

    gsInfo << patchEnd.coefs() << "\n";
    gsInfo << "\n -------- \n";
    for (unsigned int i = num_cpblade + kvfit.degree(); i< num_cpblade + 2* kvfit.degree(); i++)
    {

        a_cp(i, 0) = patchEnd.coef(i - (num_cpblade + kvfit.degree() )+1, 0);
        a_cp(i, 1) = patchEnd.coef(i - (num_cpblade + kvfit.degree() )+1, 1);
        b_cp(i, 0) = patchEnd.coef(i -(num_cpblade + kvfit.degree() )+1, 0);
        b_cp(i, 1) = patchEnd.coef(i - (num_cpblade + kvfit.degree() )+1, 1) + pitch;
    }

    gsInfo << "\n --------- \n";
    gsInfo << a_cp << "\n";
    gsInfo << b_cp << "\n";
    gsInfo << "\n -------- \n";

    gsKnotVector<TT> kv_uu_start(kvfit_knots);
gsInfo << kv_uu_start;
//        kv_uu.insert(kvfit.at(kvfit.degree() + 1) / 3.0, kvfit.degree());
//        kv_uu.insert(1 - ((1-kvfit.at(kvfit.size() - kvfit.degree() - 1)) / 3.0), kvfit.degree());
    kv_uu_start.transform(1.0/3.0,2.0/3.0);
    gsInfo << kv_uu_start;
    gsKnotVector<TT> kv_uu(0,1,0,kvfit.degree()+1);

    for(unsigned i = 1; i<kv_uu_start.size()-1; i++)
    {
        kv_uu.insert(kv_uu_start.at(i));
        gsInfo << kv_uu.at(i);
    }

    gsInfo << "kv_uu" << kv_uu << "\n";
     gsInfo << "a_cp" << a_cp.rows() << "\n";


    gsBSpline<TT> a_spline(kv_uu, a_cp);
    gsBSpline<TT> b_spline(kv_uu, b_cp);

    gsWriteParaview(a_spline,"a_spline");
    gsWriteParaview(b_spline,"b_spline");

    gsKnotVector<TT> kv_vv(0, 1, 0, kvfit.degree()+1);
    c_cp << length_x1, a_cp(0, 1),
        length_x1, b_cp(0, 1);
    gsBSpline<TT> c_spline(kvlinear, c_cp);
    c_spline.degreeElevate(kvfit.degree()-1);

    gsInfo << "\n --------- \n";
    gsInfo << c_spline.coefs() << "\n";
    gsInfo << "\n -------- \n";
    d_cp << length_x2, a_cp(num_cpblade + kvfit.degree() + kvfit.degree() - 1, 1),
          length_x2, b_cp(num_cpblade + kvfit.degree() + kvfit.degree() - 1, 1);
    gsInfo << "\n --------- \n";
    gsInfo << d_cp << "\n";
    gsInfo << "\n -------- \n";
    gsBSpline<TT> d_spline(kvlinear, d_cp);
    d_spline.degreeElevate(kvfit.degree()-1);

    gsWriteParaview(c_spline,"d_spline");
    gsWriteParaview(d_spline,"c_spline");

    gsInfo << length_x1;
    gsInfo << "\n";
    gsInfo << length_x2;
    gsInfo << "\n";
    //gsInfo << width_y1;
    gsInfo << "\n";
    //gsInfo << width_y2;
    gsInfo << "\n";
    gsInfo<< "\n-------------------------\n";
    gsInfo << a_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << b_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << c_cp;
    gsInfo<< "\n-------------------------\n";
    gsInfo << d_cp;



//        gsMultiPatch<TT> * boundaries4 = new gsMultiPatch<TT>;
//        boundaries4->addPatch(c_spline);
//        boundaries4->addPatch(a_spline);
//        boundaries4->addPatch(d_spline);
//        boundaries4->addPatch(b_spline);
//        gsCoonsPatch<TT> patchAll = coonsPatch(*boundaries4);
//        patchAll.compute();
//        mp.addPatch(patchAll.result());



discreteCoonsPatch(a_spline.coefs(),b_spline.coefs(),c_spline.coefs(),d_spline.coefs(),coef_patchAll,true);
    gsTensorBSplineBasis<2, TT> patchAll_basis(kv_uu, kv_vv);
    gsTensorBSpline<2,TT> patchAll(patchAll_basis,coef_patchAll);
    gsWriteParaview(patchAll,"patchAll",true,true);


    if(refine)
    {
        std::vector<real_t> knot_insert_knots = {0,0.25,0.5,0.75,1};
        gsKnotVector<real_t> knot_insert(knot_insert_knots);
        gsKnotVector<real_t> knot_start = knot_insert;
        gsKnotVector<real_t> knot_end = knot_insert;
        knot_start.transform(0.,1/3.);
        knot_end.transform(2/3.,1);

        for (int i = 1; i < knot_insert_knots.size() - 1; i++)
        {
            patchAll.insertKnot(knot_insert.at(i),1);
            patchAll.insertKnot(knot_start.at(i),0);
            patchAll.insertKnot(knot_end.at(i),0);
        }

    }

    mp.addPatch(patchAll);



    gsWriteParaview(mp,"mp",true,true);

//=================================optimization===========================================

    gsVector<TT> area_vec(7);
    area_vec.setZero(7);
    area_vec << 0.5, 0.25, 0.5, 1, 0.5, 0.75, 1.0;//1,1,1,1,0.75,0.75,0.75;

    TT orthogonality = 0.7;
    TT skewness = 0.0;
    TT eccentricity = 0.0;
    TT intersection = 0.0;
    TT uniformity = 0.01;//uniformity_param+0.5;//0.01;// 0.25;
    //TT uniformity = 0.005;// max seTing for RB rotation
    TT area = area_vec(index);
    TT length = 0.0;
    TT epsilon = 1e-7;

    gsQualityMeasure<TT> optimization(mp.patch(0));
        TT opt_val = optimization.functional(orthogonality, skewness,
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
    gsKnotVector<TT> kv_u(0, 1, 0, kvfit.degree()+1);
    gsKnotVector<TT> kv_v(0, 1, 0, kvfit.degree()+1);
    gsTensorBSplineBasis<2, TT> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, TT> basis_blade(kvfit, kv_v);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0((kvfit.degree()+1)*(kvfit.degree()+1), 2);
    coef_patch0.setZero((kvfit.degree()+1)*(kvfit.degree()+1), 2);

    for (int i = 0; i < kvfit.degree()+1; i++)
    {
        for (int j = 0; j < kvfit.degree()+1; j++)
        {
            coef_patch0(i * (kvfit.degree()+1) + j, 0) = coefs(i * (kvfit.degree()+1) + j + i * (num_cpblade + kvfit.degree() - 1), 0);
            coef_patch0(i * (kvfit.degree()+1) + j, 1) = coefs(i * (kvfit.degree()+1) + j + i * (num_cpblade + kvfit.degree() - 1), 1);
        }
    }
    gsTensorBSpline<2, TT> patch0(basis, coef_patch0);
    for (TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch0.insertKnot(knot, 0);
        patch0.insertKnot(knot, 1);
    }

    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1(num_cpblade * (kvfit.degree()+1), 2);
    coef_patch1.setZero(num_cpblade * (kvfit.degree()+1), 2);

    for (int i = 0; i < (kvfit.degree()+1); i++)
    {
        for (unsigned j = 0; j < num_cpblade; j++)
        {
            coef_patch1(i*num_cpblade + j, 0) = coefs(i * (kvfit.degree()+1) + j + kvfit.degree() + i * (num_cpblade + kvfit.degree() - 1), 0);
            coef_patch1(i*num_cpblade + j, 1) = coefs(i * (kvfit.degree()+1) + j + kvfit.degree() + i * (num_cpblade + kvfit.degree() - 1), 1);
        }
    }
    //gsInfo << "\n pressure side \n";
    //gsInfo << coef_patch1<< "\n";
    gsTensorBSpline<2, TT> patch1(basis_blade, coef_patch1);
    for (TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch1.insertKnot(knot, 1);
    }

    //--------------------------------patch 2-------------------------------------------
    gsMatrix<TT> coef_patch2((kvfit.degree()+1)*(kvfit.degree()+1), 2);
    coef_patch2.setZero((kvfit.degree()+1)*(kvfit.degree()+1), 2);

    for (int i = 0; i < (kvfit.degree()+1); i++)
    {
        for (int j = 0; j < (kvfit.degree()+1); j++)
        {
            coef_patch2(i * (kvfit.degree()+1) + j, 0) = coefs(i * (kvfit.degree()+1) + j + (num_cpblade + kvfit.degree() - 1) + i * (num_cpblade + kvfit.degree() - 1), 0);
            coef_patch2(i * (kvfit.degree()+1) + j, 1) = coefs(i * (kvfit.degree()+1) + j + (num_cpblade + kvfit.degree() - 1) + i * (num_cpblade + kvfit.degree() - 1), 1);
        }
    }
    gsTensorBSpline<2, TT> patch2(basis, coef_patch2);
    for (TT knot = 0.2; knot < 1.0; knot += 0.2)
    {
        patch2.insertKnot(knot, 0);
        patch2.insertKnot(knot, 1);
    }



    mpFinal.addPatch(patch0);
    mpFinal.addPatch(patch1);
    mpFinal.addPatch(patch2);

    mpFinal.addInterface(0, boundary::east, 1, boundary::west);
    mpFinal.addInterface(1, boundary::east, 2, boundary::west);
    //mpFinal.addAutoBoundaries();

    //periodic
   mpFinal.addInterface(0, boundary::south,(size_t) 0, boundary::north);
   mpFinal.addInterface(2, boundary::north, 2, boundary::south);
   mpFinal.addAutoBoundaries();

    return mp;
}
