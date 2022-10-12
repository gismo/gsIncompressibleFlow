#ifndef UWBKAPLANTURBINERUNNERWHEELDOMAIN_H
#define UWBKAPLANTURBINERUNNERWHEELDOMAIN_H

#include <gsIO/gsIOUtils.h>
//#include "gsQualityMeasure.h"

using namespace gismo;

const double PI = 3.14159265358979323846264338327950288419716939937510;

void saveData(const gsGeometry<>& geom,
              const std::string& output,
              const int number)
{
    gsFileData<> fileData;
    fileData << geom;

    std::string out = output + util::to_string(number) + ".xml";
    fileData.dump(out);

    gsMesh<> mesh;
    makeMesh(geom.basis(), mesh, 10);
    geom.evaluateMesh(mesh);
    out = output + "Mesh" + util::to_string(number);
    gsWriteParaview(mesh, out);

//    gsMatrix<TT> coefs = geom.coefs();
//    coefs.transposeInPlace();
//    out = output + "ControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(coefs, out);

//    gsMatrix<index_t> boundary = geom.basis().boundary();
//    gsMatrix<> b(geom.geoDim(), boundary.rows());
//    for (int row = 0; row < boundary.rows(); row++)
//    {
//        b.col(row) = coefs.col((boundary)(row, 0));
//    }

//    out = output + "BoundaryControlPoints" + util::to_string(number);
//    gsWriteParaviewPoints(b, out);

}

void checkJacobianDeterminant(const gsGeometry<>& geom,
                              const int points = 100,
                              const bool savePoints = false,
                              const std::string& output = "",
                              const int number = 0)
{
    gsInfo << "Checking Jacobian determinant ..." << "\n";

    gsMatrix<> para  = geom.support();
    gsVector<> c0 = para.col(0);
    gsVector<> c1 = para.col(1);
    gsMatrix<> pts = uniformPointGrid(c0, c1, points);

    gsMatrix<> plus(pts.rows(), pts.cols());
    plus.setZero();
    gsMatrix<> minus(pts.rows(), pts.cols());
    minus.setZero();
    gsMatrix<> zero(pts.rows(), pts.cols());
    zero.setZero();

    int plusCounter = 0;
    int minusCounter = 0;
    int zeroCounter = 0;

    // calculating the determinant of the Jacobian
    for (int col = 0; col < pts.cols(); col++)
    {
        const real_t determinant = geom.jacobian(pts.col(col)).determinant();
        gsMatrix<> tmp = geom.eval(pts.col(col));

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
}


template <class T>
class KaplanTurbineRunnerWheelDomain
{
public:
    // Constructors and destructors
    KaplanTurbineRunnerWheelDomain() {}
    KaplanTurbineRunnerWheelDomain(HydraulicProfile<T> HydraulicProfile, KaplanTurbineRunnerBlade<T> RunnerBlade, int numBladesInRunnerWheel);
    ~KaplanTurbineRunnerWheelDomain() {}

    // Access function
    //void setSuctionSideSurface(gsTensorBSpline<2,T> bspl) { mySuctionSideSurface = bspl; }
    //gsTensorBSpline<2,T> getSuctionSideSurface() { return mySuctionSideSurface; }
    gsMultiPatch<T> getRunnerDomain() { return m_runnerDomain; }

    // Other functions
    int compute(bool plot, bool print_info = false);
    int computeFrontandBackExtensionSurfaces(bool plot, bool print_info = false);
    int rotateFlowPassageBoundary(bool plot);
    int computeDividerSurfaces(bool plot, bool print_info = false);
    int computeFrontAndBackSurfaces(bool plot, bool print_info = false);
    //int computeTopAndBottomSurfaces(bool plot, bool print_info);
    int computeTopOrBottomSurface(gsTensorBSpline<2,T> surface_front, gsTensorBSpline<2,T> surface_front_rotated, gsTensorBSpline<2,T> surface_frontmiddle, gsTensorBSpline<2,T> surface_frontmost,
                                                                     gsTensorBSpline<2,T> surface_inner, bool side, int domain_index, int num_sample_pars1, int num_sample_pars2, bool plot, bool print_info = false);
    int computeVolumes(bool plot, bool print_info = false);

    gsBSpline<T> computeCircularArc(gsMatrix<T> starting_point, gsMatrix<T> end_point);
    gsMatrix<T> findClosestPointsOnSurface(gsMatrix<T> point, gsTensorBSpline<2, T> surface, gsMatrix<T> initial_pars, T error, int max_iter, bool print_info = false);
    int maxErrorForSurfaceApproximation(gsMatrix<T> points, gsTensorBSpline<2, T> surface, gsMatrix<T> pars, T & max_error, T & avg_error, bool print_info = false);
    int surfaceApproximationIterativeWithBoundary(gsMatrix<T> points, gsKnotVector<T> kv1, gsKnotVector<T> kv2, gsMatrix<T> controlnet_initial, gsMatrix<T> sampling_par_points_surface, gsTensorBSpline<2, T> & surface,
                                         T & max_approx_error, T & avg_approx_error, index_t max_iter, bool print_info = false);
    int surfaceApproximationIterativeFairing(gsMatrix<T> & control_net, int dim1, int dim2, T gamma);
    int selectInitialControlNet(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom, gsTensorBSpline<2, T> surface_parameter_domain,
                            gsTensorBSpline<2, T> surface, gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info = false);
    int selectInitialControlNet(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom,
                                gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info = false);
    int selectInitialControlNet2(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom,
                                gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info = false);
    int controlNetFairing(gsTensorBSpline<2, T> & surface_initial, int num, bool print_info = false);
    int controlNetFairing2(gsTensorBSpline<2, T> & surface_initial, int num, bool print_info = false);

private:

    int m_numBladeprofiles;
    int m_numBladesInRunnerWheel;

    KaplanTurbineRunnerBlade<T> m_runnerBlade;
    HydraulicProfile<T> m_hydraulicProfile;

    gsTensorBSpline<2, T> m_surfacePressureSide;
    gsTensorBSpline<2, T> m_surfaceSuctionSide;
    gsTensorBSpline<2, T> m_surfacePressureSideRotated;
    gsTensorBSpline<2, T> m_surfaceInner;
    gsTensorBSpline<2, T> m_surfaceOuter;
    gsTensorBSpline<2, T> m_surfaceFrontExtension;
    gsTensorBSpline<2, T> m_surfaceBackExtension;
    gsTensorBSpline<2, T> m_surfaceFrontExtensionRotated;
    gsTensorBSpline<2, T> m_surfaceBackExtensionRotated;
    gsTensorBSpline<2, T> m_surfaceFrontMiddleDivider;
    gsTensorBSpline<2, T> m_surfaceMiddleBackDivider;
    gsTensorBSpline<2, T> m_surfaceFront;
    gsTensorBSpline<2, T> m_surfaceBack;
    gsTensorBSpline<2, T> m_surfaceTopFrontDomain;
    gsTensorBSpline<2, T> m_surfaceBottomFrontDomain;
    gsTensorBSpline<2, T> m_surfaceTopMiddleDomain;
    gsTensorBSpline<2, T> m_surfaceBottomMiddleDomain;
    gsTensorBSpline<2, T> m_surfaceTopBackDomain;
    gsTensorBSpline<2, T> m_surfaceBottomBackDomain;

    gsTensorBSpline<3, T> m_volumeFrontDomain;
    gsTensorBSpline<3, T> m_volumeMiddleDomain;
    gsTensorBSpline<3, T> m_volumeBackDomain;

    gsMultiPatch<T> m_runnerDomain;

    gsBSpline<T> m_surfaceFrontExtensionProfiles[20];
    gsBSpline<T> m_surfaceBackExtensionProfiles[20];
};

template<class T>
KaplanTurbineRunnerWheelDomain<T>::KaplanTurbineRunnerWheelDomain(HydraulicProfile<T> HydraulicProfile, KaplanTurbineRunnerBlade<T> RunnerBlade, int numBladesInRunnerWheel) {

    m_hydraulicProfile = HydraulicProfile;
    m_runnerBlade = RunnerBlade;
    m_surfaceInner = m_hydraulicProfile.getHydraulicProfileInner();
    m_surfaceOuter = m_hydraulicProfile.getHydraulicProfileOuter();
    m_numBladesInRunnerWheel = numBladesInRunnerWheel;
    m_numBladeprofiles = m_runnerBlade.getNumBladeProfiles();
    m_surfacePressureSide = m_runnerBlade.getPressureSideSurfaceAfterTrimming();
    m_surfaceSuctionSide = m_runnerBlade.getSuctionSideSurfaceAfterTrimming();
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::compute(bool plot, bool print_info) {

    gsInfo << "\nRunner wheel domain construction\n";
    gsInfo << "================================\n";
    gsInfo << "Computation of front and back extension surfaces ...";
    computeFrontandBackExtensionSurfaces(plot, print_info);
    gsInfo << "Done.\n";
    rotateFlowPassageBoundary(plot);
    gsInfo << "Necessary surfaces rotated ...\n";
    gsInfo << "Computation of divider surfaces between front/middle and middle/back patches ...\n";
    computeDividerSurfaces(plot, print_info);
    gsInfo << "Done.\n";
    gsInfo << "Computation of front most and back most surfaces ...\n";
    computeFrontAndBackSurfaces(plot, print_info);
    gsInfo << "Done.\n";
    gsInfo << "Approximation of up and bottom surfaces (on inner and outer surface of a turbine) ...\n";
    //computeTopAndBottomSurfaces(plot, print_info);
    computeTopOrBottomSurface(m_surfaceFrontExtension, m_surfaceFrontExtensionRotated, m_surfaceFrontMiddleDivider, m_surfaceFront, m_surfaceInner, 0, 0, 30, 20, plot, print_info);
    computeTopOrBottomSurface(m_surfaceFrontExtension, m_surfaceFrontExtensionRotated, m_surfaceFrontMiddleDivider, m_surfaceFront, m_surfaceOuter, 1, 0, 10, 50, plot, print_info);
    computeTopOrBottomSurface(m_surfaceSuctionSide, m_surfacePressureSideRotated, m_surfaceMiddleBackDivider, m_surfaceFrontMiddleDivider, m_surfaceInner, 0, 1, 30, 20, plot, print_info);
    computeTopOrBottomSurface(m_surfaceSuctionSide, m_surfacePressureSideRotated, m_surfaceMiddleBackDivider, m_surfaceFrontMiddleDivider, m_surfaceOuter, 1, 1, 15, 30, plot, print_info);
    computeTopOrBottomSurface(m_surfaceBackExtension, m_surfaceBackExtensionRotated, m_surfaceBack, m_surfaceMiddleBackDivider, m_surfaceInner, 0, 2, 10, 20, plot, print_info);
    computeTopOrBottomSurface(m_surfaceBackExtension, m_surfaceBackExtensionRotated, m_surfaceBack, m_surfaceMiddleBackDivider, m_surfaceOuter, 1, 2, 10, 30, plot, print_info);
    gsInfo << "Done.\n";

    gsInfo << "Computation of final volumes of front, middle and back patches ...\n";
    computeVolumes(plot, print_info);
    gsInfo << "Done.\n";

    return 0;

}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeFrontandBackExtensionSurfaces(bool plot, bool print_info) {

    int num_sample_pars = 30;

    // ********************************
    // Construction of extension curves
    // ********************************
    gsMatrix<T> pars_for_extension_blade(2,4);
    gsMatrix<T> pars_for_extension_HR(2,2);
    gsMatrix<T> suction_surface_corner_points(3,4);
    gsMatrix<T> surface_inner_points(3,2);
    gsMatrix<T> surface_outer_points(3,2);
    gsMatrix<T> suction_surface_corner_ders(6,4);
    gsMatrix<T> pressure_surface_corner_ders(6,4);
    gsMatrix<T> surface_inner_ders(6,2);
    gsMatrix<T> surface_outer_ders(6,2);
    gsTensorBSpline<2, T> suction_side_surface = m_runnerBlade.getSuctionSideSurfaceAfterTrimming();
    gsTensorBSpline<2, T> pressure_side_surface = m_runnerBlade.getPressureSideSurfaceAfterTrimming();
    gsVector<T> vec1(3);
    gsVector<T> vec2(3);
    gsVector<T> direction_vector1_FU(3);
    gsVector<T> direction_vector2_FU(3);
    gsVector<T> direction_vector1_FB(3);
    gsVector<T> direction_vector2_FB(3);
    gsVector<T> direction_vector1_BU(3);
    gsVector<T> direction_vector2_BU(3);
    gsVector<T> direction_vector1_BB(3);
    gsVector<T> direction_vector2_BB(3);
    gsMatrix<T> coefs(4, 3);
    gsKnotVector<T> kv(0, 1, 0, 4);
    gsBSpline<T> curve_front_bottom;
    gsBSpline<T> curve_front_up;
    gsBSpline<T> curve_back_bottom;
    gsBSpline<T> curve_back_up;

    pars_for_extension_blade << 0, 0, 1, 1,
            0, 1, 0, 1;
    suction_side_surface.deriv_into(pars_for_extension_blade, suction_surface_corner_ders);
    suction_side_surface.eval_into(pars_for_extension_blade, suction_surface_corner_points);
    pressure_side_surface.deriv_into(pars_for_extension_blade, pressure_surface_corner_ders);

    // front bottom curve
    vec1(0) = suction_surface_corner_ders(0,0);
    vec1(1) = suction_surface_corner_ders(2,0);
    vec1(2) = suction_surface_corner_ders(4,0);
    vec2(0) = suction_surface_corner_ders(1,0);
    vec2(1) = suction_surface_corner_ders(3,0);
    vec2(2) = suction_surface_corner_ders(5,0);
    cross(vec1, vec2, direction_vector2_FB);

    gsMatrix<T> suction_intersection_point_on_HR = m_runnerBlade.getSuctionIntersectionPointinParameterDomainonInnerHR();
    T a1 = suction_intersection_point_on_HR(1,0);
    T b1 = suction_intersection_point_on_HR(1, suction_intersection_point_on_HR.cols()-1) + (suction_intersection_point_on_HR(1, suction_intersection_point_on_HR.cols()-1) - suction_intersection_point_on_HR(1, 0))/3;
    //gsInfo << a1 << "\n";
    //gsInfo << b1 << "\n";
    pars_for_extension_HR  << 0, 1,
                              a1, b1;
    m_surfaceInner.deriv_into(pars_for_extension_HR, surface_inner_ders);
    m_surfaceInner.eval_into(pars_for_extension_HR, surface_inner_points);
    direction_vector1_FB(0) = surface_inner_ders(0,0);
    direction_vector1_FB(1) = surface_inner_ders(2,0);
    direction_vector1_FB(2) = surface_inner_ders(4,0);

    coefs.row(0) = surface_inner_points.col(0);
    coefs.row(1) = surface_inner_points.col(0) + direction_vector1_FB/6;
    coefs.row(2) = suction_surface_corner_points.col(0) - 2*direction_vector2_FB/3;
    coefs.row(3) = suction_surface_corner_points.col(0);
    curve_front_bottom = gsBSpline<T>( kv, coefs);

    //front up curve
    vec1(0) = suction_surface_corner_ders(0,2);
    vec1(1) = suction_surface_corner_ders(2,2);
    vec1(2) = suction_surface_corner_ders(4,2);
    vec2(0) = suction_surface_corner_ders(1,2);
    vec2(1) = suction_surface_corner_ders(3,2);
    vec2(2) = suction_surface_corner_ders(5,2);
    cross(vec1, vec2, direction_vector2_FU);

    m_surfaceOuter.deriv_into(pars_for_extension_HR, surface_outer_ders);
    m_surfaceOuter.eval_into(pars_for_extension_HR, surface_outer_points);
    direction_vector1_FU(0) = surface_outer_ders(0,0);
    direction_vector1_FU(1) = surface_outer_ders(2,0);
    direction_vector1_FU(2) = surface_outer_ders(4,0);
    coefs.row(0) = surface_outer_points.col(0);
    coefs.row(1) = surface_outer_points.col(0) + direction_vector1_FU/6;
    coefs.row(2) = suction_surface_corner_points.col(2) - 2*direction_vector2_FU/3;
    coefs.row(3) = suction_surface_corner_points.col(2);
    curve_front_up = gsBSpline<T>( kv, coefs);

    // back bottom curve
    direction_vector1_BB(0) = (suction_surface_corner_ders(1,1) + pressure_surface_corner_ders(1,1))/4;
    direction_vector1_BB(1) = (suction_surface_corner_ders(3,1) + pressure_surface_corner_ders(3,1))/4;
    direction_vector1_BB(2) = (suction_surface_corner_ders(5,1) + pressure_surface_corner_ders(5,1))/4;
    direction_vector2_BB(0) = surface_inner_ders(0,1)/2;
    direction_vector2_BB(1) = surface_inner_ders(2,1)/2;
    direction_vector2_BB(2) = surface_inner_ders(4,1)/2;
    coefs.row(0) = suction_surface_corner_points.col(1);
    coefs.row(1) = suction_surface_corner_points.col(1) + direction_vector1_BB/3;
    coefs.row(2) = surface_inner_points.col(1) - direction_vector2_BB/3;
    coefs.row(3) = surface_inner_points.col(1);
    curve_back_bottom = gsBSpline<T>( kv, coefs);

    // back up curve
    direction_vector1_BU(0) = (suction_surface_corner_ders(1,3) + pressure_surface_corner_ders(1,3))/4;
    direction_vector1_BU(1) = (suction_surface_corner_ders(3,3) + pressure_surface_corner_ders(3,3))/4;
    direction_vector1_BU(2) = (suction_surface_corner_ders(5,3) + pressure_surface_corner_ders(5,3))/4;
    direction_vector2_BU(0) = surface_outer_ders(0,1)/4;
    direction_vector2_BU(1) = surface_outer_ders(2,1)/4;
    direction_vector2_BU(2) = surface_outer_ders(4,1)/4;
    coefs.row(0) = suction_surface_corner_points.col(3);
    coefs.row(1) = suction_surface_corner_points.col(3) + direction_vector1_BU/6;
    coefs.row(2) = surface_outer_points.col(1) - direction_vector2_BU/3;
    coefs.row(3) = surface_outer_points.col(1);
    curve_back_up = gsBSpline<T>( kv, coefs);

    if (plot) {
        gsWriteParaview( curve_front_bottom, "extension_curve_front_bottom", 100);
        gsWriteParaview( curve_front_up, "extension_curve_front_up", 100);
        gsWriteParaview( curve_back_bottom, "extension_curve_back_bottom", 100);
        gsWriteParaview( curve_back_up, "extension_curve_back_up", 100);
    }

    // *******************************************************
    // Projecting extension curves on inner and outer surfaces
    // *******************************************************
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);
    gsBSpline<T> curve_FU;
    gsBSpline<T> curve_FB;
    gsBSpline<T> curve_BB;
    gsBSpline<T> curve_BU;
    gsMatrix<T> parameter_points(1, num_sample_pars);
    gsMatrix<T> parameter_points_front_up(2, num_sample_pars);
    gsMatrix<T> parameter_points_front_bottom(2, num_sample_pars);
    gsMatrix<T> parameter_points_back_up(2, num_sample_pars);
    gsMatrix<T> parameter_points_back_bottom(2, num_sample_pars);
    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = (T) i/(num_sample_pars-1);
    }
    gsMatrix<T> points_on_curve_front_up(3, num_sample_pars);
    gsMatrix<T> points_on_curve_front_up_projected(3, num_sample_pars);
    gsMatrix<T> points_on_curve_front_bottom(3, num_sample_pars);
    gsMatrix<T> points_on_curve_front_bottom_projected(3, num_sample_pars);
    gsMatrix<T> points_on_curve_back_up(3, num_sample_pars);
    gsMatrix<T> points_on_curve_back_up_projected(3, num_sample_pars);
    gsMatrix<T> points_on_curve_back_bottom(3, num_sample_pars);
    gsMatrix<T> points_on_curve_back_bottom_projected(3, num_sample_pars);
    curve_front_up.eval_into(parameter_points, points_on_curve_front_up);
    curve_front_bottom.eval_into(parameter_points, points_on_curve_front_bottom);
    curve_back_up.eval_into(parameter_points, points_on_curve_back_up);
    curve_back_bottom.eval_into(parameter_points, points_on_curve_back_bottom);
    gsVector<T> init_sol(2);
    init_sol << 0.0, 0.9;
    gsVector<T> sol(2);
    int num_iter_front_up = 0;
    int max_num_iter_front_up = 0;
    int num_iter_front_bottom = 0;
    int max_num_iter_front_bottom = 0;
    int num_iter_back_up = 0;
    int max_num_iter_back_up = 0;
    int num_iter_back_bottom = 0;
    int max_num_iter_back_bottom = 0;

    for (index_t i = 0; i < num_sample_pars; i++) {
        gsVector<T> point_on_curveFU = points_on_curve_front_up.col(i);
        init_sol << parameter_points(i)/3, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_curveFU, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_front_up, print_info);
        parameter_points_front_up(0, i) = sol(0);
        parameter_points_front_up(1, i) = sol(1);
        if (max_num_iter_front_up < num_iter_front_up) {
            max_num_iter_front_up = num_iter_front_up;
        }
        gsVector<T> point_on_curveFB = points_on_curve_front_bottom.col(i);
        init_sol << parameter_points(i)/4, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_curveFB, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_front_bottom, print_info);
        parameter_points_front_bottom(0, i) = sol(0);
        parameter_points_front_bottom(1, i) = sol(1);
        if (max_num_iter_front_bottom < num_iter_front_bottom) {
            max_num_iter_front_bottom = num_iter_front_bottom;
        }
        gsVector<T> point_on_curveBU = points_on_curve_back_up.col(i);
        init_sol << 0.8, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_curveBU, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_back_up, print_info);
        parameter_points_back_up(0, i) = sol(0);
        parameter_points_back_up(1, i) = sol(1);
        if (max_num_iter_back_up < num_iter_back_up) {
            max_num_iter_back_up = num_iter_back_up;
        }
        gsVector<T> point_on_curveBB = points_on_curve_back_bottom.col(i);
        minimizePointSurfaceDistanceviaNR(point_on_curveBB, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_back_bottom, print_info);
        parameter_points_back_bottom(0, i) = sol(0);
        parameter_points_back_bottom(1, i) = sol(1);
        if (max_num_iter_back_bottom < num_iter_back_bottom) {
            max_num_iter_back_bottom = num_iter_back_bottom;
        }
    }
    m_surfaceOuter.eval_into(parameter_points_front_up, points_on_curve_front_up_projected);
    points_on_curve_front_up_projected.col(0) = points_on_curve_front_up.col(0);
    points_on_curve_front_up_projected.col(num_sample_pars-1) = points_on_curve_front_up.col(num_sample_pars-1);
    m_surfaceOuter.eval_into(parameter_points_back_up, points_on_curve_back_up_projected);
    points_on_curve_back_up_projected.col(0) = points_on_curve_back_up.col(0);
    points_on_curve_back_up_projected.col(num_sample_pars-1) = points_on_curve_back_up.col(num_sample_pars-1);
    m_surfaceInner.eval_into(parameter_points_front_bottom, points_on_curve_front_bottom_projected);
    points_on_curve_front_bottom_projected.col(0) = points_on_curve_front_bottom.col(0);
    points_on_curve_front_bottom_projected.col(num_sample_pars-1) = points_on_curve_front_bottom.col(num_sample_pars-1);
    m_surfaceInner.eval_into(parameter_points_back_bottom, points_on_curve_back_bottom_projected);
    points_on_curve_back_bottom_projected.col(0) = points_on_curve_back_bottom.col(0);
    points_on_curve_back_bottom_projected.col(num_sample_pars-1) = points_on_curve_back_bottom.col(num_sample_pars-1);

    if (plot) {
        gsWriteParaviewPoints<real_t>( points_on_curve_front_up_projected, "extension_curve_front_up_projectedpoints");
        gsWriteParaviewPoints<real_t>( points_on_curve_front_bottom_projected, "extension_curve_front_bottom_projectedpoints");
        gsWriteParaviewPoints<real_t>( points_on_curve_back_up_projected, "extension_curve_back_up_projectedpoints");
        gsWriteParaviewPoints<real_t>( points_on_curve_back_bottom_projected, "extension_curve_back_bottom_projectedpoints");
    }

    gsMatrix<T> points_on_curve_front_up_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_FU(1, num_sample_pars);
    points_on_curve_front_up_projected2 = points_on_curve_front_up_projected.transpose();
    par_points_FU = centripetalParameterization(points_on_curve_front_up_projected2);
    gsMatrix<T> tangent(3,1);
    tangent(0) = direction_vector1_FU(0);
    tangent(1) = direction_vector1_FU(1);
    tangent(2) = direction_vector1_FU(2);
    //gsInfo << tangent << "\n";
    //gsInfo << kvfit << "\n";
    //gsInfo << par_points_FU << "\n";
    //gsInfo << par_points_FU.cols() << "\n";
    //gsInfo << points_on_curve_front_up_projected2 << "\n";
    //gsInfo << points_on_curve_front_up_projected2.rows() << "\n";
    curve_FU = curveFittingWithBoundaryAndInputTangent(points_on_curve_front_up_projected2, par_points_FU, kvfit, tangent);
    gsMatrix<T> curve_FU_control_points = curve_FU.coefs();
    if (print_info) {
        gsInfo << "Front-up extension curve computed.\n";
    }
    //gsInfo << curve_FU << "\n";

    gsMatrix<T> points_on_curve_front_bottom_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_FB(1, num_sample_pars);
    points_on_curve_front_bottom_projected2 = points_on_curve_front_bottom_projected.transpose();
    par_points_FB = centripetalParameterization(points_on_curve_front_bottom_projected2);
    tangent(0) = direction_vector1_FB(0);
    tangent(1) = direction_vector1_FB(1);
    tangent(2) = direction_vector1_FB(2);
    curve_FB = curveFittingWithBoundaryAndInputTangent(points_on_curve_front_bottom_projected2, par_points_FB, kvfit, tangent);
    gsMatrix<T> curve_FB_control_points = curve_FB.coefs();
    if (print_info) {
        gsInfo << "Front-bottom extension curve computed.\n";
    }
    //gsInfo << curve_FB << "\n";

    gsMatrix<T> points_on_curve_back_bottom_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_BB(1, num_sample_pars);
    points_on_curve_back_bottom_projected2 = points_on_curve_back_bottom_projected.transpose();
    par_points_BB = centripetalParameterization(points_on_curve_back_bottom_projected2);
    tangent(0) = direction_vector1_BB(0);
    tangent(1) = direction_vector1_BB(1);
    tangent(2) = direction_vector1_BB(2);
    curve_BB = curveFittingWithBoundaryAndInputTangent(points_on_curve_back_bottom_projected2, par_points_BB, kvfit, tangent);
    gsMatrix<T> curve_BB_control_points = curve_BB.coefs();
    if (print_info) {
        gsInfo << "Back-bottom extension curve computed.\n";
    }
    //gsInfo << curve_BB << "\n";

    gsMatrix<T> points_on_curve_back_up_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_BU(1, num_sample_pars);
    points_on_curve_back_up_projected2 = points_on_curve_back_up_projected.transpose();
    par_points_BU = centripetalParameterization(points_on_curve_back_up_projected2);
    tangent(0) = direction_vector1_BU(0);
    tangent(1) = direction_vector1_BU(1);
    tangent(2) = direction_vector1_BU(2);
    curve_BU = curveFittingWithBoundaryAndInputTangent(points_on_curve_back_up_projected2, par_points_BU, kvfit, tangent);
    gsMatrix<T> curve_BU_control_points = curve_BU.coefs();
    if (print_info) {
        gsInfo << "Back-up extension curve computed.\n";
    }
    //gsInfo << curve_BU << "\n";

    if (plot) {
        gsWriteParaview( curve_FU, "extension_curve_front_up_approximatedonsurface", 100);
        gsWriteParaview( curve_FB, "extension_curve_front_bottom_approximatedonsurface", 100);
        gsWriteParaview( curve_BB, "extension_curve_back_bottom_approximatedonsurface", 100);
        gsWriteParaview( curve_BU, "extension_curve_back_up_approximatedonsurface", 100);
    }

    // ************************************************
    // Construction of front surface (before the blade)
    // ************************************************
    //num_bladeprofiles = m_runnerBlade.getNumBladeProfiles();
    gsKnotVector<T> kvbase_front(0, 1, 0, 2);
    gsMatrix<T> coefs_front(2, 3);
    gsBSpline<T> curve_front;
    coefs_front.row(0) = curve_FB_control_points.row(0);
    coefs_front.row(1) = curve_FU_control_points.row(0);
    curve_front = gsBSpline<>( kvbase_front, coefs_front);
    //gsInfo << curve_front << "\n";
    int deg_dif = kvloft.degree()-kvbase_front.degree();
    for (index_t i = 0; i < deg_dif; i++) {
        curve_front.degreeElevate();
    }
    //gsInfo << curve_front << "\n";
    std::vector<T> res;
    kvloft.difference(curve_front.knots(), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_front.insertKnot(res[i]);
    }
    if (print_info) {
        gsInfo << "The front-most curve of the front extension surface computed.\n";
    }
    //gsInfo << curve_front << "\n";

    gsBSpline<T> curve_front2;
    gsMatrix<T> cp_curve_front2(kvloft.size()-kvloft.degree()-1, 3);
    gsMatrix<T> control_net_suctionside;
    control_net_suctionside = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).coefs();
    //gsInfo << control_net_suctionside << "\n";
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        cp_curve_front2.row(i) = control_net_suctionside.row(i);
    }
    //gsInfo << cp_curve_front2.rows() << ", " << cp_curve_front2.cols() << "\n";
    //gsInfo << cp_curve_front2 << "\n";
    curve_front2 = gsBSpline<T> (kvloft, cp_curve_front2);
    if (print_info) {
        gsInfo << "The back-most curve of the front extension surface computed.\n";
    }

    if (plot) {
        gsWriteParaview( curve_front, "curve_front", 100);
        gsWriteParaview( curve_front2, "curve_front2", 100);
    }

    //gsInfo << curve_FB_control_points << "\n\n";
    //gsInfo << curve_FU_control_points << "\n\n";
    //gsInfo << curve_front.coefs() << "\n\n";
    //gsInfo << curve_front2.coefs() << "\n\n";

    gsVector<T> direction_vector1(3);
    gsVector<T> direction_vector2(3);
    //gsBSpline<T> surface_front_profiles[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_front_profiles;
    gsBSpline<T> curve_aux;
    gsMatrix<T> section_params(1, m_numBladeprofiles);
    section_params = m_runnerBlade.getSectionParameters();
    gsMatrix<T> pars_for_suctionside_profiles(2, m_numBladeprofiles);
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        pars_for_suctionside_profiles(0,i) = section_params(i);
        pars_for_suctionside_profiles(1,i) = 0;
    }
    gsMatrix<T> suction_side_profiles_ders(6, m_numBladeprofiles);
    suction_side_surface.deriv_into(pars_for_suctionside_profiles, suction_side_profiles_ders);
    gsMatrix<T> curve_front_profile_points(3, m_numBladeprofiles);
    gsMatrix<T> curve_front2_profile_points(3, m_numBladeprofiles);
    curve_front.eval_into(section_params, curve_front_profile_points);
    curve_front2.eval_into(section_params, curve_front2_profile_points);
    kvfit.difference(kv, res);
    surface_front_profiles.push_back(curve_FB);
    for (index_t i = 1; i < m_numBladeprofiles-1; i++) {
        direction_vector1(0) = (1-section_params(i)) * direction_vector1_FB(0) + section_params(i) * direction_vector1_FU(0);
        direction_vector1(1) = (1-section_params(i)) * direction_vector1_FB(1) + section_params(i) * direction_vector1_FU(1);
        direction_vector1(2) = (1-section_params(i)) * direction_vector1_FB(2) + section_params(i) * direction_vector1_FU(2);
        vec1(0) = suction_side_profiles_ders(0,i);
        vec1(1) = suction_side_profiles_ders(2,i);
        vec1(2) = suction_side_profiles_ders(4,i);
        vec2(0) = suction_side_profiles_ders(1,i);
        vec2(1) = suction_side_profiles_ders(3,i);
        vec2(2) = suction_side_profiles_ders(5,i);
        cross(vec1, vec2, direction_vector2);
        coefs.row(0) = curve_front_profile_points.col(i);
        coefs.row(1) = curve_front_profile_points.col(i) + direction_vector1/6;
        coefs.row(2) = curve_front2_profile_points.col(i) - direction_vector2_BU/3;
        coefs.row(3) = curve_front2_profile_points.col(i);
        curve_aux = gsBSpline<T>( kv, coefs);
        for (unsigned i = 0; i < res.size(); i++) {
            curve_aux.insertKnot(res[i]);
        }
        //surface_front_profiles[i] = curve_aux;
        surface_front_profiles.push_back(curve_aux);
    }
    //surface_front_profiles[0] = curve_FB;
    //surface_front_profiles[m_numBladeprofiles-1] = curve_FU;
    surface_front_profiles.push_back(curve_FU);
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        m_surfaceFrontExtensionProfiles[i] = surface_front_profiles[i];
    }

    gsTensorBSpline<2, T> surface_front;
    computeLoftSurface(surface_front_profiles, kvfit, m_numBladeprofiles, kvloft, section_params, surface_front);
    m_surfaceFrontExtension = surface_front;

    if (print_info) {
        gsInfo << "Front extension surface computed.\n";
    }

    if (plot) {
        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_front_profiles[i]);
        }
        gsWriteParaview( curves3d, "surface_front_profiles", 100);

        gsWriteParaview( surface_front, "surface_front", 5000);

        gsMesh<> surface_front_mesh;
        surface_front.controlNet(surface_front_mesh);
        gsWriteParaview( surface_front_mesh, "surface_front_mesh");
    }

    // ***********************************************
    // Construction of back surface (behind the blade)
    // ***********************************************
    gsKnotVector<T> kvbase_back(0, 1, 0, 2);
    gsMatrix<T> coefs_back2(2, 3);
    gsBSpline<T> curve_back2;
    int last_index = kvfit.size() - kvfit.degree() - 1 - 1;
    coefs_back2.row(0) = curve_BB_control_points.row(last_index);
    coefs_back2.row(1) = curve_BU_control_points.row(last_index);
    curve_back2 = gsBSpline<T>( kvbase_back, coefs_back2);
    //gsInfo << curve_back2 << "\n";
    deg_dif = kvloft.degree()-kvbase_back.degree();
    for (index_t i = 0; i < deg_dif; i++) {
        curve_back2.degreeElevate();
    }
    //gsInfo << curve_back2 << "\n";
    kvloft.difference(curve_back2.knots(), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_back2.insertKnot(res[i]);
    }
    if (print_info) {
        gsInfo << "The back-most curve of the back extension surface computed.\n";
    }
    //gsInfo << curve_back2 << "\n";

    gsBSpline<T> curve_back;
    gsMatrix<T> cp_curve_back(kvloft.size()-kvloft.degree()-1, 3);
    control_net_suctionside = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).coefs();
    //gsInfo << control_net_suctionside << "\n";
    last_index = control_net_suctionside.rows();
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        cp_curve_back.row(i) = control_net_suctionside.row(last_index - m_numBladeprofiles + i);
    }
    //gsInfo << cp_curve_back.rows() << ", " << cp_curve_back.cols() << "\n";
    //gsInfo << cp_curve_back << "\n";
    curve_back = gsBSpline<T> (kvloft, cp_curve_back);
    if (print_info) {
        gsInfo << "The front-most curve of the back extension surface computed.\n";
    }

    gsMatrix<T> surface_back_control_net(cp_curve_back.rows() * curve_BB_control_points.rows(), 3);
    surface_back_control_net.setZero();
    //gsInfo << curve_BB_control_points << "\n\n";
    //gsInfo << curve_BU_control_points << "\n\n";
    //gsInfo << curve_back.coefs() << "\n\n";
    //gsInfo << curve_back2.coefs() << "\n\n";
    //gsVector<> direction_vector1(3);
    //gsVector<> direction_vector2(3);
    //num_bladeprofiles = runnerBlade.getNumBladeProfiles();
    //gsBSpline<T> surface_back_profiles[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_back_profiles;
    //gsBSpline<> curve_aux;
    //gsMatrix<> section_params(1, num_bladeprofiles);
    //section_params = runnerBlade.getSectionParameters();
    //gsMatrix<> pars_for_suctionside_profiles(2, num_bladeprofiles);
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        pars_for_suctionside_profiles(0,i) = section_params(i);
        pars_for_suctionside_profiles(1,i) = 1;
    }
    gsMatrix<T> pressure_side_profiles_ders(6, m_numBladeprofiles);
    suction_side_surface.deriv_into(pars_for_suctionside_profiles, suction_side_profiles_ders);
    pressure_side_surface.deriv_into(pars_for_suctionside_profiles, pressure_side_profiles_ders);
    gsMatrix<T> curve_back_profile_points(3, m_numBladeprofiles);
    gsMatrix<T> curve_back2_profile_points(3, m_numBladeprofiles);
    curve_back.eval_into(section_params, curve_back_profile_points);
    curve_back2.eval_into(section_params, curve_back2_profile_points);
    kvfit.difference(kv, res);
    surface_back_profiles.push_back(curve_BB);
    for (index_t i = 1; i < m_numBladeprofiles-1; i++) {
        direction_vector2(0) = (1-section_params(i)) * direction_vector2_BB(0) + section_params(i) * direction_vector2_BU(0);
        direction_vector2(1) = (1-section_params(i)) * direction_vector2_BB(1) + section_params(i) * direction_vector2_BU(1);
        direction_vector2(2) = (1-section_params(i)) * direction_vector2_BB(2) + section_params(i) * direction_vector2_BU(2);
        //vec1(0) = suction_side_profiles_ders(0,i);
        //vec1(1) = suction_side_profiles_ders(2,i);
        //vec1(2) = suction_side_profiles_ders(4,i);
        //vec2(0) = suction_side_profiles_ders(1,i);
        //vec2(1) = suction_side_profiles_ders(3,i);
        //vec2(2) = suction_side_profiles_ders(5,i);
        //cross(vec1, vec2, direction_vector1);
        direction_vector1(0) = (suction_side_profiles_ders(1,i) + pressure_side_profiles_ders(1,i))/8;
        direction_vector1(1) = (suction_side_profiles_ders(3,i) + pressure_side_profiles_ders(3,i))/8;
        direction_vector1(2) = (suction_side_profiles_ders(5,i) + pressure_side_profiles_ders(5,i))/8;
        coefs.row(0) = curve_back_profile_points.col(i);
        coefs.row(1) = curve_back_profile_points.col(i) + direction_vector1/3;
        coefs.row(2) = curve_back2_profile_points.col(i) - direction_vector2_BU/6;
        coefs.row(3) = curve_back2_profile_points.col(i);
        curve_aux = gsBSpline<T>( kv, coefs);
        for (unsigned i = 0; i < res.size(); i++) {
            curve_aux.insertKnot(res[i]);
        }
        //surface_back_profiles[i] = curve_aux;
        surface_back_profiles.push_back(curve_aux);
    }
    //surface_back_profiles[0] = curve_BB;
    //surface_back_profiles[m_numBladeprofiles-1] = curve_BU;
    surface_back_profiles.push_back(curve_BU);
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        m_surfaceBackExtensionProfiles[i] = surface_back_profiles[i];
    }

    gsTensorBSpline<2, T> surface_back;
    computeLoftSurface(surface_back_profiles, kvfit, m_numBladeprofiles, kvloft, section_params, surface_back);
    m_surfaceBackExtension = surface_back;

    if (plot) {
        gsWriteParaview( curve_back2, "curve_back2", 100);
        gsWriteParaview( curve_back, "curve_back", 100);

        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_back_profiles[i]);
        }
        gsWriteParaview( curves3d, "surface_back_profiles", 100);

        gsWriteParaview( surface_back, "surface_back", 5000);

        gsMesh<> surface_back_mesh;
        surface_back.controlNet(surface_back_mesh);
        gsWriteParaview( surface_back_mesh, "surface_back_mesh");
    }

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::rotateFlowPassageBoundary(bool plot) {

    //int numBladesInRunnerWheel = 4;
    gsVector<T> turbine_axis(3);
    turbine_axis << 0.0, 0.0, 1.0;
    gsTensorBSpline<2, T> surface_front_rotated = m_surfaceFrontExtension;
    surface_front_rotated.rotate(2*PI/m_numBladesInRunnerWheel, turbine_axis);
    m_surfaceFrontExtensionRotated = surface_front_rotated;

    gsTensorBSpline<2, T> pressure_side_surface_rotated = m_surfacePressureSide;
    pressure_side_surface_rotated.rotate(2*PI/m_numBladesInRunnerWheel, turbine_axis);
    m_surfacePressureSideRotated = pressure_side_surface_rotated;

    gsTensorBSpline<2, T> surface_back_rotated = m_surfaceBackExtension;
    surface_back_rotated.rotate(2*PI/m_numBladesInRunnerWheel, turbine_axis);
    m_surfaceBackExtensionRotated = surface_back_rotated;

    if (plot) {
        gsWriteParaview( surface_front_rotated, "surface_front_rotated", 5000);
        gsWriteParaview( pressure_side_surface_rotated, "pressure_side_surface_rotated", 5000);
        gsWriteParaview( surface_back_rotated, "surface_back_rotated", 5000);
    }

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeDividerSurfaces(bool plot, bool print_info) {

    int num_sample_pars = 30;

    // *********************************************
    // Divider surface between front and midle patch
    // *********************************************
    gsVector<T> vec1(3);
    gsVector<T> vec2(3);
    gsMatrix<T> coefs(4, 3);
    gsMatrix<T> tangent(3,1);
    gsMatrix<T> dersStart(2,3);
    gsMatrix<T> dersEnd(2,3);
    gsVector<T> direction_vector1_FM_B1(3);
    gsVector<T> direction_vector2_FM_B2(3);
    gsVector<T> direction_vector1_FM_U1(3);
    gsVector<T> direction_vector2_FM_U2(3);
    gsVector<T> direction_vector1_MB_B1(3);
    gsVector<T> direction_vector2_MB_B2(3);
    gsVector<T> direction_vector1_MB_U1(3);
    gsVector<T> direction_vector2_MB_U2(3);
    gsMatrix<T> control_net_surface_front((m_surfaceFrontExtension.coefs()).rows(), 3);
    gsMatrix<T> control_net_surface_back((m_surfaceBackExtension.coefs()).rows(), 3);
    gsMatrix<T> control_net_pressureside_rotated;
    gsMatrix<T> control_net_surface_front_rotated;
    gsMatrix<T> control_net_surface_back_rotated;
    gsMatrix<T> control_net_suctionside;
    gsBSpline<T> curve_frontmiddle_bottom;
    gsBSpline<T> curve_frontmiddle_up;
    gsBSpline<T> curve_middleback_bottom;
    gsBSpline<T> curve_middleback_up;

    control_net_surface_front.setZero();
    control_net_surface_front = m_surfaceFrontExtension.coefs();
    control_net_surface_back.setZero();
    control_net_surface_back = m_surfaceBackExtension.coefs();
    control_net_pressureside_rotated.setZero();
    control_net_pressureside_rotated = m_surfacePressureSideRotated.coefs();
    control_net_surface_front_rotated.setZero();
    control_net_surface_front_rotated = m_surfaceFrontExtensionRotated.coefs();
    control_net_surface_back_rotated.setZero();
    control_net_surface_back_rotated = m_surfaceBackExtensionRotated.coefs();
    control_net_suctionside = m_surfaceSuctionSide.coefs();

    gsKnotVector<T> kv(0, 1, 0, 4);
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);

    gsMatrix<T> parameter_points(1, num_sample_pars);
    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = (T) i/(num_sample_pars-1);
    }

    // frontmiddle bottom curve
    int num_cp_frontsurface_u = m_surfaceFrontExtension.knots(0).size() - m_surfaceFrontExtension.degree(0) - 1;
    int num_cp_frontsurface_v = m_surfaceFrontExtension.knots(1).size() - m_surfaceFrontExtension.degree(1) - 1;
    //vec1 = control_net_surface_front.row((curve_FB_control_points.rows()-2)*cp_curve_front2.rows()) - control_net_surface_front.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows());
    //vec2 = control_net_suctionside.row(cp_curve_front2.rows()) - control_net_suctionside.row(0);
    //direction_vector1_FM_B1 = (vec1 + vec2)/2;
    //gsInfo << direction_vector1_FM_B1 << "\n";
    //vec1 = control_net_surface_front_rotated.row((curve_FB_control_points.rows()-2)*cp_curve_front2.rows()) - control_net_surface_front_rotated.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows());;
    //vec2 = control_net_pressureside_rotated.row(cp_curve_front2.rows()) - control_net_pressureside_rotated.row(0);
    //direction_vector2_FM_B2 = (vec1 + vec2)/2;
    //coefs.row(0) = control_net_suctionside.row(0);
    //coefs.row(1) = control_net_suctionside.row(0) + 50*direction_vector1_FM_B1.transpose()/3;
    //coefs.row(2) = control_net_pressureside_rotated.row(0) + 50*direction_vector2_FM_B2.transpose()/3;
    //coefs.row(3) = control_net_pressureside_rotated.row(0);
    vec1 = control_net_surface_front.row((num_cp_frontsurface_v-2)*num_cp_frontsurface_u) - control_net_surface_front.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);
    vec2 = control_net_suctionside.row(num_cp_frontsurface_u) - control_net_suctionside.row(0);
    direction_vector1_FM_B1 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/6;
    //gsInfo << direction_vector1_FM_B1 << "\n";
    vec1 = control_net_surface_front_rotated.row((num_cp_frontsurface_v-2)*num_cp_frontsurface_u) - control_net_surface_front_rotated.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);;
    vec2 = control_net_pressureside_rotated.row(num_cp_frontsurface_u) - control_net_pressureside_rotated.row(0);
    direction_vector2_FM_B2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/6;
    coefs.row(0) = control_net_suctionside.row(0);
    coefs.row(1) = control_net_suctionside.row(0) + direction_vector1_FM_B1.transpose()/3;
    coefs.row(2) = control_net_pressureside_rotated.row(0) + direction_vector2_FM_B2.transpose()/3;
    coefs.row(3) = control_net_pressureside_rotated.row(0);
    //gsInfo << coefs << "\n";
    curve_frontmiddle_bottom = gsBSpline<T>( kv, coefs);
    if (plot) {
        gsWriteParaview( curve_frontmiddle_bottom, "curve_frontmiddle_bottom", 100);
    }

    gsMatrix<T> parameter_points_frontmiddle_bottom(2, num_sample_pars);
    gsMatrix<T> points_on_curve_frontmiddle_bottom(3, num_sample_pars);
    gsMatrix<T> points_on_curve_frontmiddle_bottom_projected(3, num_sample_pars);
    int num_iter_frontmiddle_bottom = 0;
    int max_num_iter_frontmiddle_bottom = 0;
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    curve_frontmiddle_bottom.eval_into(parameter_points, points_on_curve_frontmiddle_bottom);
    for (index_t i = 0; i < num_sample_pars; i++) {
        gsVector<T> point_on_curveFM_B = points_on_curve_frontmiddle_bottom.col(i);
        init_sol << 0.4, 0.45;
        minimizePointSurfaceDistanceviaNR(point_on_curveFM_B, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_frontmiddle_bottom, print_info);
        parameter_points_frontmiddle_bottom(0, i) = sol(0);
        parameter_points_frontmiddle_bottom(1, i) = sol(1);
        if (max_num_iter_frontmiddle_bottom < num_iter_frontmiddle_bottom) {
            max_num_iter_frontmiddle_bottom = num_iter_frontmiddle_bottom;
        }
        if ((num_iter_frontmiddle_bottom > 100) && (print_info)) {
            gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
        }
    }
    if (print_info) {
        gsInfo << "Max number of NR iterations: " << max_num_iter_frontmiddle_bottom << "\n";
    }
    m_surfaceInner.eval_into(parameter_points_frontmiddle_bottom, points_on_curve_frontmiddle_bottom_projected);
    points_on_curve_frontmiddle_bottom_projected.col(0) = points_on_curve_frontmiddle_bottom.col(0);
    points_on_curve_frontmiddle_bottom_projected.col(num_sample_pars-1) = points_on_curve_frontmiddle_bottom.col(num_sample_pars-1);
    if (plot) {
        gsWriteParaviewPoints<real_t>( points_on_curve_frontmiddle_bottom_projected, "curve_frontmiddle_bottom_projectedpoints");
    }

    gsMatrix<T> points_on_curve_frontmiddle_bottom_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_FM_B(1, num_sample_pars);
    gsBSpline<T> curve_FM_B;
    points_on_curve_frontmiddle_bottom_projected2 = points_on_curve_frontmiddle_bottom_projected.transpose();
    par_points_FM_B = centripetalParameterization(points_on_curve_frontmiddle_bottom_projected2);
    tangent(0) = direction_vector1_FM_B1(0);
    tangent(1) = direction_vector1_FM_B1(1);
    tangent(2) = direction_vector1_FM_B1(2);
    curve_FM_B = curveFittingWithBoundaryAndInputTangent(points_on_curve_frontmiddle_bottom_projected2, par_points_FM_B, kvfit, tangent);
    //gsMatrix<> curve_FM_B_control_points = curve_FM_B.coefs();
    //gsInfo << curve_FM_B << "\n";
    if (plot) {
        gsWriteParaview( curve_FM_B, "curve_frontmiddle_bottom_approximatedonsurface", 100);
    }

    // frontmiddle up
    //vec1 = control_net_surface_front.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows()-1) - control_net_surface_front.row(curve_FB_control_points.rows()*cp_curve_front2.rows()-1);
    //vec2 = control_net_suctionside.row(2*cp_curve_front2.rows()-1) - control_net_suctionside.row(cp_curve_front2.rows()-1);
    //direction_vector1_FM_U1 = (vec1 + vec2)/2;
    //gsInfo << direction_vector1_FM_U1 << "\n";
    //vec1 = control_net_surface_front_rotated.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows()-1) - control_net_surface_front_rotated.row(curve_FB_control_points.rows()*cp_curve_front2.rows()-1);;
    //vec2 = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1) - control_net_pressureside_rotated.row(cp_curve_front2.rows()-1);
    //direction_vector2_FM_U2 = (vec1 + vec2)/2;
    //coefs.row(0) = control_net_suctionside.row(cp_curve_front2.rows()-1);
    //coefs.row(1) = control_net_suctionside.row(cp_curve_front2.rows()-1) + 50*direction_vector1_FM_U1.transpose()/3;
    //coefs.row(2) = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1) + 50*direction_vector2_FM_U2.transpose()/3;
    //coefs.row(3) = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1);
    vec1 = control_net_surface_front.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u-1) - control_net_surface_front.row(num_cp_frontsurface_v*num_cp_frontsurface_u-1);
    vec2 = control_net_suctionside.row(2*num_cp_frontsurface_u-1) - control_net_suctionside.row(num_cp_frontsurface_u-1);
    direction_vector1_FM_U1 = 2*(vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/3;
    //gsInfo << direction_vector1_FM_U1 << "\n";
    vec1 = control_net_surface_front_rotated.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u-1) - control_net_surface_front_rotated.row(num_cp_frontsurface_v*num_cp_frontsurface_u-1);;
    vec2 = control_net_pressureside_rotated.row(2*num_cp_frontsurface_u-1) - control_net_pressureside_rotated.row(num_cp_frontsurface_u-1);
    direction_vector2_FM_U2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/21;
    coefs.row(0) = control_net_suctionside.row(num_cp_frontsurface_u-1);
    coefs.row(1) = control_net_suctionside.row(num_cp_frontsurface_u-1) + direction_vector1_FM_U1.transpose()/3;
    coefs.row(2) = control_net_pressureside_rotated.row(num_cp_frontsurface_u-1) + direction_vector2_FM_U2.transpose()/3;
    coefs.row(3) = control_net_pressureside_rotated.row(num_cp_frontsurface_u-1);
    //gsInfo << coefs << "\n";
    curve_frontmiddle_up = gsBSpline<T>( kv, coefs);
    if (plot) {
        gsWriteParaview( curve_frontmiddle_up, "curve_frontmiddle_up", 100);
    }

    gsMatrix<T> parameter_points_frontmiddle_up(2, num_sample_pars);
    gsMatrix<T> points_on_curve_frontmiddle_up(3, num_sample_pars);
    gsMatrix<T> points_on_curve_frontmiddle_up_projected(3, num_sample_pars);
    int num_iter_frontmiddle_up = 0;
    int max_num_iter_frontmiddle_up = 0;
    curve_frontmiddle_up.eval_into(parameter_points, points_on_curve_frontmiddle_up);
    for (index_t i = 0; i < num_sample_pars; i++) {
        gsVector<T> point_on_curveFM_U = points_on_curve_frontmiddle_up.col(i);
        init_sol << 0.4, 0.3+parameter_points(i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_curveFM_U, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-6), 200, num_iter_frontmiddle_up, print_info);
        parameter_points_frontmiddle_up(0, i) = sol(0);
        parameter_points_frontmiddle_up(1, i) = sol(1);
        if (max_num_iter_frontmiddle_up < num_iter_frontmiddle_up) {
            max_num_iter_frontmiddle_up = num_iter_frontmiddle_up;
        }
        if ((num_iter_frontmiddle_up > 100) && (print_info)) {
            gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
        }
    }
    if (print_info) {
        gsInfo << "Max number of NR iterations: " << max_num_iter_frontmiddle_up << "\n";
    }
    m_surfaceOuter.eval_into(parameter_points_frontmiddle_up, points_on_curve_frontmiddle_up_projected);
    points_on_curve_frontmiddle_up_projected.col(0) = points_on_curve_frontmiddle_up.col(0);
    points_on_curve_frontmiddle_up_projected.col(num_sample_pars-1) = points_on_curve_frontmiddle_up.col(num_sample_pars-1);
    if (plot) {
        gsWriteParaviewPoints<real_t>( points_on_curve_frontmiddle_up_projected, "curve_frontmiddle_up_projectedpoints");
    }

    gsMatrix<T> points_on_curve_frontmiddle_up_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_FM_U(1, num_sample_pars);
    gsBSpline<T> curve_FM_U;
    points_on_curve_frontmiddle_up_projected2 = points_on_curve_frontmiddle_up_projected.transpose();
    par_points_FM_U = centripetalParameterization(points_on_curve_frontmiddle_up_projected2);
    tangent(0) = direction_vector1_FM_U1(0);
    tangent(1) = direction_vector1_FM_U1(1);
    tangent(2) = direction_vector1_FM_U1(2);
    //curve_FM_U = curveFittingWithBoundaryAndInputTangent(points_on_curve_frontmiddle_up_projected2, par_points_FM_U, kvfit, tangent);
    dersStart.row(0) = points_on_curve_frontmiddle_up_projected2.row(0);
    dersStart.row(1) = direction_vector1_FM_U1;
    dersEnd.row(0) = points_on_curve_frontmiddle_up_projected2.row(points_on_curve_frontmiddle_up_projected2.rows()-1);
    dersEnd.row(1) = -direction_vector2_FM_U2/2;
    curve_FM_U = curveFittingWithEndDerivatives(points_on_curve_frontmiddle_up_projected2, par_points_FM_U, kvfit, dersStart, dersEnd);
    //gsMatrix<> curve_FM_U_control_points = curve_FM_U.coefs();
    //gsInfo << curve_FM_U << "\n";
    if (plot) {
        gsWriteParaview( curve_FM_U, "curve_frontmiddle_up_approximatedonsurface", 100);
    }
    gsInfo << "Vypocteny okrajove krivky pro delici plochu ...\n";

    // frontmiddle profiles
    //num_bladeprofiles = runnerBlade.getNumBladeProfiles();
    //gsBSpline<T> surface_frontmiddle_profiles[m_numBladeprofiles];
    //gsBSpline<T> surface_frontmiddle_profiles_initial[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_frontmiddle_profiles;
    std::vector<gsBSpline<T> > surface_frontmiddle_profiles_initial;
    gsMatrix<T> parameter_points_aux(2, num_sample_pars);
    gsMatrix<T> points_on_curve_aux(3, num_sample_pars);
    gsMatrix<T> points_on_curve_aux_projected(3, num_sample_pars);
    gsVector<T> cylinder_radii = m_runnerBlade.getCylinderRadius();
    gsTensorBSpline<2, T> cylinder;
    int num_iter_aux = 0;
    int max_num_iter_aux = 0;
    gsMatrix<T> points_on_curve_aux_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_aux(1, num_sample_pars);
    gsBSpline<T> curve_aux;
    gsBSpline<T> curve_aux2;
    gsMatrix<T> pars(2,2);
    gsMatrix<T> pars2(2,2);
    gsMatrix<T> points_on_surfacefrontprofile(3,2);
    gsMatrix<T> ders_on_surfacefrontprofile(6,2);
    gsMatrix<T> points_on_surfacefrontrotatedprofile(3,2);
    gsMatrix<T> ders_on_surfacefrontrotatedprofile(6,2);
    gsMatrix<T> points_on_surfacebladeprofile(3,2);
    gsMatrix<T> ders_on_surfacebladeprofile(6,2);
    gsMatrix<T> points_on_surfacebladerotatedprofile(3,2);
    gsMatrix<T> ders_on_surfacebladerotatedprofile(6,2);
    //BladeProfile<T> * blade_profile;
    //gsBSpline<T> surface_front_profiles[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_front_profiles;
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        //surface_front_profiles[i] = m_surfaceFrontExtensionProfiles[i];
        surface_front_profiles.push_back(m_surfaceFrontExtensionProfiles[i]);
    }
    gsMatrix<T> sec_pars(1,m_numBladeprofiles);
    sec_pars = m_runnerBlade.getSectionParameters();
    //kvfit.difference(kv, res);
    surface_frontmiddle_profiles.push_back(curve_FM_B);
    surface_frontmiddle_profiles_initial.push_back(curve_FM_B);
    for (index_t i = 1; i < m_numBladeprofiles-1; i++) {
        //gsInfo << i << "\n";
        /*pars << 0.0, 1.0;
        surface_front_profiles[i].eval_into(pars, points_on_surfacefrontprofile);
        surface_front_profiles[i].deriv_into(pars, ders_on_surfacefrontprofile);
        blade_profile = m_runnerBlade.getBladeProfile(i);
        (blade_profile->getSuctionSide3d()).eval_into(pars, points_on_surfacebladeprofile);
        (blade_profile->getSuctionSide3d()).deriv_into(pars, ders_on_surfacebladeprofile);
        vec1(0) = - ders_on_surfacefrontprofile(0,1);
        vec1(1) = - ders_on_surfacefrontprofile(1,1);
        vec1(2) = - ders_on_surfacefrontprofile(2,1);
        vec2(0) = ders_on_surfacebladeprofile(0,0);
        vec2(1) = ders_on_surfacebladeprofile(1,0);
        vec2(2) = ders_on_surfacebladeprofile(2,0);*/
        pars << sec_pars(i), sec_pars(i),
                0.0, 1.0;
        //gsInfo << pars << "\n";
        m_surfaceFrontExtension.eval_into(pars, points_on_surfacefrontprofile);
        m_surfaceFrontExtension.deriv_into(pars, ders_on_surfacefrontprofile);
        //gsInfo << points_on_surfacefrontprofile << "\n\n";
        //gsInfo << ders_on_surfacefrontprofile << "\n\n";
        m_surfaceSuctionSide.eval_into(pars, points_on_surfacebladeprofile);
        m_surfaceSuctionSide.deriv_into(pars, ders_on_surfacebladeprofile);
        //gsInfo << points_on_surfacebladeprofile << "\n\n";
        //gsInfo << ders_on_surfacebladeprofile << "\n\n";
        vec1(0) = - ders_on_surfacefrontprofile(1,1);
        vec1(1) = - ders_on_surfacefrontprofile(3,1);
        vec1(2) = - ders_on_surfacefrontprofile(5,1);
        vec2(0) = ders_on_surfacebladeprofile(1,0);
        vec2(1) = ders_on_surfacebladeprofile(3,0);
        vec2(2) = ders_on_surfacebladeprofile(5,0);
        direction_vector1_FM_B1 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/10;
        //gsInfo << direction_vector1_FM_B1 << "\n";

        pars2 << sec_pars(i), sec_pars(i),
                0.0, 1.0;
        //gsInfo << pars2 << "\n";
        m_surfaceFrontExtensionRotated.eval_into(pars2, points_on_surfacefrontrotatedprofile);
        m_surfaceFrontExtensionRotated.deriv_into(pars2, ders_on_surfacefrontrotatedprofile);
        //gsInfo << points_on_surfacefrontrotatedprofile << "\n\n";
        //gsInfo << ders_on_surfacefrontrotatedprofile << "\n\n";
        m_surfacePressureSideRotated.eval_into(pars2, points_on_surfacebladerotatedprofile);
        m_surfacePressureSideRotated.deriv_into(pars2, ders_on_surfacebladerotatedprofile);
        //gsInfo << points_on_surfacebladerotatedprofile << "\n\n";
        //gsInfo << ders_on_surfacebladerotatedprofile << "\n\n";
        vec1(0) = ders_on_surfacefrontrotatedprofile(1,1);
        vec1(1) = ders_on_surfacefrontrotatedprofile(3,1);
        vec1(2) = ders_on_surfacefrontrotatedprofile(5,1);
        vec2(0) = - ders_on_surfacebladerotatedprofile(1,0);
        vec2(1) = - ders_on_surfacebladerotatedprofile(3,0);
        vec2(2) = - ders_on_surfacebladerotatedprofile(5,0);
        direction_vector2_FM_B2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/10;
        //gsInfo << "a\n";

        coefs.row(0) = points_on_surfacefrontprofile.col(1);
        coefs.row(1) = points_on_surfacefrontprofile.col(1) + direction_vector1_FM_B1/3;
        coefs.row(2) = points_on_surfacefrontrotatedprofile.col(1) - direction_vector2_FM_B2/3;
        coefs.row(3) = points_on_surfacefrontrotatedprofile.col(1);
        //gsInfo << coefs << "\n";
        curve_aux = gsBSpline<T>( kv, coefs);
        curve_aux.eval_into(parameter_points, points_on_curve_aux);
        //surface_frontmiddle_profiles_initial[i] = curve_aux;
        surface_frontmiddle_profiles_initial.push_back(curve_aux);
        constructCylinderParametrization(-0.5, 0.3, cylinder_radii(i+1), cylinder);
        if (i==1) {
            gsWriteParaviewPoints<real_t>( points_on_curve_aux, "surface_frontmiddle_profile_inital_points");
        }
        for (index_t i = 0; i < num_sample_pars; i++) {
            gsVector<T> point_on_curve_aux = points_on_curve_aux.col(i);
            init_sol << 0.4, 0.3+parameter_points(i)/4;
            minimizePointSurfaceDistanceviaNR(point_on_curve_aux, cylinder, init_sol, sol, static_cast<real_t>(1e-6), 200, num_iter_aux, print_info);
            parameter_points_aux(0, i) = sol(0);
            parameter_points_aux(1, i) = sol(1);
            if (max_num_iter_aux < num_iter_aux) {
                max_num_iter_aux = num_iter_aux;
            }
            if ((num_iter_aux > 100) && (print_info)) {
                gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
            }
        }
        if (print_info) {
            gsInfo << "Max number of NR iterations: " << max_num_iter_aux << "\n";
        }
        cylinder.eval_into(parameter_points_aux, points_on_curve_aux_projected);
        points_on_curve_aux_projected2 = points_on_curve_aux_projected.transpose();
        par_points_aux = centripetalParameterization(points_on_curve_aux_projected2);
        tangent(0) = direction_vector1_FM_B1(0);
        tangent(1) = direction_vector1_FM_B1(1);
        tangent(2) = direction_vector1_FM_B1(2);
        curve_aux2 = curveFittingWithBoundaryAndInputTangent(points_on_curve_aux_projected2, par_points_aux, kvfit, tangent);
    //    for (unsigned i = 0; i < res.size(); i++) {
    //        curve_aux.insertKnot(res[i]);
    //    }
        //surface_frontmiddle_profiles[i] = curve_aux2;
        surface_frontmiddle_profiles.push_back(curve_aux2);
    }
    //surface_frontmiddle_profiles[0] = curve_FM_B;
    //surface_frontmiddle_profiles[m_numBladeprofiles-1] = curve_FM_U;
    //surface_frontmiddle_profiles_initial[0] = curve_FM_B;
    //surface_frontmiddle_profiles_initial[m_numBladeprofiles-1] = curve_FM_U;
    surface_frontmiddle_profiles.push_back(curve_FM_U);
    surface_frontmiddle_profiles_initial.push_back(curve_FM_U);
    //gsInfo << "a\n";

    gsTensorBSpline<2, T> surface_frontmiddle;
    computeLoftSurface(surface_frontmiddle_profiles, kvfit, m_numBladeprofiles, kvloft, sec_pars, surface_frontmiddle);

    int num1 = kvloft.size() - kvloft.degree() - 1;
    int num2 = kvfit.size() - kvfit.degree() - 1;
    gsMatrix<T> control_net_pom = surface_frontmiddle.coefs();
    for (index_t i = 0; i < num1; i++) {
        control_net_pom.row(i) = m_surfaceSuctionSide.coefs().row(i);
        control_net_pom.row((num2-1)*num1+i) = m_surfacePressureSideRotated.coefs().row(i);
    }
    gsTensorBSplineBasis<2,T> basis(kvloft, kvfit);
    surface_frontmiddle = gsTensorBSpline<2,T> (basis, control_net_pom);
    m_surfaceFrontMiddleDivider = surface_frontmiddle;

    //curves3d = {};
    if (plot) {
        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_frontmiddle_profiles[i]);
        }
        gsWriteParaview( curves3d, "surface_frontmiddle_profiles", 100);
        //curves3d = {};
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_frontmiddle_profiles_initial[i]);
        }
        gsWriteParaview( curves3d, "surface_frontmiddle_profiles_inital", 100);

        gsWriteParaview( surface_frontmiddle, "surface_frontmiddle", 15000);

        gsMesh<> surface_frontmiddle_mesh;
        surface_frontmiddle.controlNet(surface_frontmiddle_mesh);
        gsWriteParaview( surface_frontmiddle_mesh, "surface_frontmiddle_mesh");
    }

    // *********************************************
    // Divider surface between middle and back patch
    // *********************************************
    // middleback bottom curve
    int num_cp_backsurface_u = m_surfaceBackExtension.knots(0).size() - m_surfaceBackExtension.degree(0) - 1;
    int num_cp_backsurface_v = m_surfaceBackExtension.knots(1).size() - m_surfaceBackExtension.degree(1) - 1;
    //vec1 = control_net_surface_front.row((curve_FB_control_points.rows()-2)*cp_curve_front2.rows()) - control_net_surface_front.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows());
    //vec2 = control_net_suctionside.row(cp_curve_front2.rows()) - control_net_suctionside.row(0);
    //direction_vector1_FM_B1 = (vec1 + vec2)/2;
    //gsInfo << direction_vector1_FM_B1 << "\n";
    //vec1 = control_net_surface_front_rotated.row((curve_FB_control_points.rows()-2)*cp_curve_front2.rows()) - control_net_surface_front_rotated.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows());;
    //vec2 = control_net_pressureside_rotated.row(cp_curve_front2.rows()) - control_net_pressureside_rotated.row(0);
    //direction_vector2_FM_B2 = (vec1 + vec2)/2;
    //coefs.row(0) = control_net_suctionside.row(0);
    //coefs.row(1) = control_net_suctionside.row(0) + 50*direction_vector1_FM_B1.transpose()/3;
    //coefs.row(2) = control_net_pressureside_rotated.row(0) + 50*direction_vector2_FM_B2.transpose()/3;
    //coefs.row(3) = control_net_pressureside_rotated.row(0);
    vec1 = control_net_surface_back.row(num_cp_backsurface_u) - control_net_surface_front.row(0);
    vec2 = control_net_suctionside.row((num_cp_backsurface_v-2)*num_cp_backsurface_u) - control_net_suctionside.row((num_cp_backsurface_v-1)*num_cp_backsurface_u);
    direction_vector1_MB_B1 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/4;
    vec1 = control_net_surface_back_rotated.row(num_cp_backsurface_u) - control_net_surface_front_rotated.row(0);;
    vec2 = control_net_pressureside_rotated.row((num_cp_backsurface_v-2)*num_cp_backsurface_u) - control_net_pressureside_rotated.row((num_cp_backsurface_v-1)*num_cp_backsurface_u);
    direction_vector2_MB_B2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/2;
    coefs.row(0) = control_net_suctionside.row((num_cp_backsurface_v-1)*num_cp_backsurface_u);
    coefs.row(1) = control_net_suctionside.row((num_cp_backsurface_v-1)*num_cp_backsurface_u) + direction_vector1_MB_B1.transpose()/3;
    coefs.row(2) = control_net_pressureside_rotated.row((num_cp_backsurface_v-1)*num_cp_backsurface_u) - direction_vector2_MB_B2.transpose()/3;
    coefs.row(3) = control_net_pressureside_rotated.row((num_cp_backsurface_v-1)*num_cp_backsurface_u);
    //gsInfo << coefs << "\n";
    curve_middleback_bottom = gsBSpline<T>( kv, coefs);
    if (plot) {
        gsWriteParaview( curve_middleback_bottom, "curve_middleback_bottom", 100);
    }

    gsMatrix<T> parameter_points_middleback_bottom(2, num_sample_pars);
    gsMatrix<T> points_on_curve_middleback_bottom(3, num_sample_pars);
    gsMatrix<T> points_on_curve_middleback_bottom_projected(3, num_sample_pars);
    int num_iter_middleback_bottom = 0;
    int max_num_iter_middleback_bottom = 0;
    curve_middleback_bottom.eval_into(parameter_points, points_on_curve_middleback_bottom);
    for (index_t i = 0; i < num_sample_pars; i++) {
        gsVector<T> point_on_curveMB_B = points_on_curve_middleback_bottom.col(i);
        init_sol << 0.8, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_curveMB_B, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-10), 200, num_iter_middleback_bottom, print_info);
        parameter_points_middleback_bottom(0, i) = sol(0);
        parameter_points_middleback_bottom(1, i) = sol(1);
        if (max_num_iter_middleback_bottom < num_iter_middleback_bottom) {
            max_num_iter_middleback_bottom = num_iter_middleback_bottom;
        }
        if ((num_iter_middleback_bottom > 100) && (print_info)) {
            gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
        }
    }
    if (print_info) {
        gsInfo << "Max number of NR iterations: " << max_num_iter_middleback_bottom << "\n";
    }
    m_surfaceInner.eval_into(parameter_points_middleback_bottom, points_on_curve_middleback_bottom_projected);
    points_on_curve_middleback_bottom_projected.col(0) = points_on_curve_middleback_bottom.col(0);
    points_on_curve_middleback_bottom_projected.col(num_sample_pars-1) = points_on_curve_middleback_bottom.col(num_sample_pars-1);
    if (plot) {
        gsWriteParaviewPoints<real_t>( points_on_curve_middleback_bottom_projected, "curve_middleback_bottom_projectedpoints");
    }

    gsMatrix<T> points_on_curve_middleback_bottom_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_MB_B(1, num_sample_pars);
    gsBSpline<T> curve_MB_B;
    points_on_curve_middleback_bottom_projected2 = points_on_curve_middleback_bottom_projected.transpose();
    par_points_MB_B = centripetalParameterization(points_on_curve_middleback_bottom_projected2);
    tangent(0) = direction_vector1_MB_B1(0)/3;
    tangent(1) = direction_vector1_MB_B1(1)/3;
    tangent(2) = direction_vector1_MB_B1(2)/3;
    curve_MB_B = curveFittingWithBoundary(points_on_curve_middleback_bottom_projected2, par_points_MB_B, kvfit, print_info);
    //gsMatrix<> curve_FM_B_control_points = curve_FM_B.coefs();
    //gsInfo << curve_FM_B << "\n";
    if (plot) {
        gsWriteParaview( curve_MB_B, "curve_middleback_bottom_approximatedonsurface", 100);
    }

    // middleback up curve
    //vec1 = control_net_surface_front.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows()-1) - control_net_surface_front.row(curve_FB_control_points.rows()*cp_curve_front2.rows()-1);
    //vec2 = control_net_suctionside.row(2*cp_curve_front2.rows()-1) - control_net_suctionside.row(cp_curve_front2.rows()-1);
    //direction_vector1_FM_U1 = (vec1 + vec2)/2;
    //gsInfo << direction_vector1_FM_U1 << "\n";
    //vec1 = control_net_surface_front_rotated.row((curve_FB_control_points.rows()-1)*cp_curve_front2.rows()-1) - control_net_surface_front_rotated.row(curve_FB_control_points.rows()*cp_curve_front2.rows()-1);;
    //vec2 = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1) - control_net_pressureside_rotated.row(cp_curve_front2.rows()-1);
    //direction_vector2_FM_U2 = (vec1 + vec2)/2;
    //coefs.row(0) = control_net_suctionside.row(cp_curve_front2.rows()-1);
    //coefs.row(1) = control_net_suctionside.row(cp_curve_front2.rows()-1) + 50*direction_vector1_FM_U1.transpose()/3;
    //coefs.row(2) = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1) + 50*direction_vector2_FM_U2.transpose()/3;
    //coefs.row(3) = control_net_pressureside_rotated.row(cp_curve_front2.rows()-1);
    vec1 = control_net_surface_back.row(2*num_cp_backsurface_u-1) - control_net_surface_back.row(num_cp_backsurface_u-1);
    vec2 = control_net_suctionside.row((num_cp_backsurface_v-1)*num_cp_backsurface_u-1) - control_net_suctionside.row(num_cp_backsurface_v*num_cp_backsurface_u-1);
    direction_vector1_MB_U1 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/2;
    //gsInfo << direction_vector1_FM_U1 << "\n";
    vec1 = control_net_surface_back_rotated.row(2*num_cp_backsurface_u-1) - control_net_surface_back_rotated.row(num_cp_backsurface_u-1);;
    vec2 = control_net_pressureside_rotated.row((num_cp_backsurface_v-1)*num_cp_backsurface_u-1) - control_net_pressureside_rotated.row(num_cp_backsurface_v*num_cp_backsurface_u-1);
    direction_vector2_MB_U2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/2;
    coefs.row(0) = control_net_suctionside.row(num_cp_backsurface_v*num_cp_backsurface_u-1);
    coefs.row(1) = control_net_suctionside.row(num_cp_backsurface_v*num_cp_backsurface_u-1) + 5*direction_vector1_MB_U1.transpose()/3;
    coefs.row(2) = control_net_pressureside_rotated.row(num_cp_backsurface_v*num_cp_backsurface_u-1) + 5*direction_vector2_MB_U2.transpose()/3;
    coefs.row(3) = control_net_pressureside_rotated.row(num_cp_backsurface_v*num_cp_backsurface_u-1);
    //gsInfo << coefs << "\n";
    curve_middleback_up = gsBSpline<T>( kv, coefs);
    if (plot) {
        gsWriteParaview( curve_middleback_up, "curve_middleback_up", 100);
    }

    gsMatrix<T> parameter_points_middleback_up(2, num_sample_pars);
    gsMatrix<T> points_on_curve_middleback_up(3, num_sample_pars);
    gsMatrix<T> points_on_curve_middleback_up_projected(3, num_sample_pars);
    int num_iter_middleback_up = 0;
    int max_num_iter_middleback_up = 0;
    curve_middleback_up.eval_into(parameter_points, points_on_curve_middleback_up);
    for (index_t i = 0; i < num_sample_pars; i++) {
        gsVector<T> point_on_curveMB_U = points_on_curve_middleback_up.col(i);
        init_sol << 0.8, 0.1+parameter_points(i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_curveMB_U, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-6), 200, num_iter_middleback_up, print_info);
        parameter_points_middleback_up(0, i) = sol(0);
        parameter_points_middleback_up(1, i) = sol(1);
        if (max_num_iter_middleback_up < num_iter_middleback_up) {
            max_num_iter_middleback_up = num_iter_middleback_up;
        }
        if ((num_iter_middleback_up > 100) && (print_info)) {
            gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
        }
    }
    if (print_info) {
        gsInfo << "Max number of NR iterations: " << max_num_iter_middleback_up << "\n";
    }
    m_surfaceOuter.eval_into(parameter_points_middleback_up, points_on_curve_middleback_up_projected);
    points_on_curve_middleback_up_projected.col(0) = points_on_curve_middleback_up.col(0);
    points_on_curve_middleback_up_projected.col(num_sample_pars-1) = points_on_curve_middleback_up.col(num_sample_pars-1);
    if (plot) {
        gsWriteParaviewPoints<real_t>( points_on_curve_middleback_up_projected, "curve_middleback_up_projectedpoints");
    }

    gsMatrix<T> points_on_curve_middleback_up_projected2(num_sample_pars, 3);
    gsMatrix<T> par_points_MB_U(1, num_sample_pars);
    gsBSpline<T> curve_MB_U;
    points_on_curve_middleback_up_projected2 = points_on_curve_middleback_up_projected.transpose();
    par_points_MB_U = centripetalParameterization(points_on_curve_middleback_up_projected2);
    tangent(0) = direction_vector1_MB_U1(0);
    tangent(1) = direction_vector1_MB_U1(1);
    tangent(2) = direction_vector1_MB_U1(2);
    curve_MB_U = curveFittingWithBoundary(points_on_curve_middleback_up_projected2, par_points_MB_U, kvfit, print_info);
    //gsMatrix<> curve_FM_U_control_points = curve_FM_U.coefs();
    //gsInfo << curve_FM_U << "\n";
    if (plot) {
        gsWriteParaview( curve_MB_U, "curve_middleback_up_approximatedonsurface", 100);
    }
    //gsInfo << "Vypocteny okrajove krivky pro delici plochu ...\n";

    // middleback profiles
    //num_bladeprofiles = runnerBlade.getNumBladeProfiles();
    //gsBSpline<T> surface_middleback_profiles[m_numBladeprofiles];
    //gsBSpline<T> surface_middleback_profiles_initial[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_middleback_profiles;
    std::vector<gsBSpline<T> > surface_middleback_profiles_initial;
    //gsMatrix<T> parameter_points_aux(2, num_sample_pars);
    //gsMatrix<T> points_on_curve_aux(3, num_sample_pars);
    //gsMatrix<T> points_on_curve_aux_projected(3, num_sample_pars);
    //gsVector<T> cylinder_radii = m_runnerBlade.getCylinderRadius();
    //gsTensorBSpline<2, T> cylinder;
    num_iter_aux = 0;
    max_num_iter_aux = 0;
    //gsMatrix<T> points_on_curve_aux_projected2(num_sample_pars, 3);
    //gsMatrix<T> par_points_aux(1, num_sample_pars);
    //gsBSpline<T> curve_aux;
    //gsBSpline<T> curve_aux2;
    //gsMatrix<T> pars(2,2);
    //gsMatrix<T> pars2(2,2);
    gsMatrix<T> points_on_surfacebackprofile(3,2);
    gsMatrix<T> ders_on_surfacebackprofile(6,2);
    gsMatrix<T> points_on_surfacebackrotatedprofile(3,2);
    gsMatrix<T> ders_on_surfacebackrotatedprofile(6,2);
    //gsMatrix<T> points_on_surfacebladeprofile(3,2);
    //gsMatrix<T> ders_on_surfacebladeprofile(6,2);
    //gsMatrix<T> points_on_surfacebladerotatedprofile(3,2);
    //gsMatrix<T> ders_on_surfacebladerotatedprofile(6,2);
    //BladeProfile<T> * blade_profile;
    //gsBSpline<T> surface_back_profiles[m_numBladeprofiles];
    std::vector<gsBSpline<T> > surface_back_profiles;
    for (index_t i = 0; i < m_numBladeprofiles; i++) {
        //surface_back_profiles[i] = m_surfaceBackExtensionProfiles[i];
        surface_back_profiles.push_back(m_surfaceBackExtensionProfiles[i]);
    }
    //gsMatrix<T> sec_pars(1,m_numBladeprofiles);
    //sec_pars = m_runnerBlade.getSectionParameters();
    //kvfit.difference(kv, res);
    surface_middleback_profiles.push_back(curve_MB_B);
    surface_middleback_profiles_initial.push_back(curve_MB_B);
    for (index_t i = 1; i < m_numBladeprofiles-1; i++) {
        //gsInfo << i << "\n";
        pars << sec_pars(i), sec_pars(i),
                0.0, 1.0;
        //gsInfo << pars << "\n";
        m_surfaceBackExtension.eval_into(pars, points_on_surfacebackprofile);
        m_surfaceBackExtension.deriv_into(pars, ders_on_surfacebackprofile);
        //gsInfo << points_on_surfacebackprofile << "\n\n";
        //gsInfo << ders_on_surfacebackprofile << "\n\n";
        m_surfaceSuctionSide.eval_into(pars, points_on_surfacebladeprofile);
        m_surfaceSuctionSide.deriv_into(pars, ders_on_surfacebladeprofile);
        //gsInfo << points_on_surfacebladeprofile << "\n\n";
        //gsInfo << ders_on_surfacebladeprofile << "\n\n";
        vec1(0) = ders_on_surfacebackprofile(1,0);
        vec1(1) = ders_on_surfacebackprofile(3,0);
        vec1(2) = ders_on_surfacebackprofile(5,0);
        vec2(0) = - ders_on_surfacebladeprofile(1,1);
        vec2(1) = - ders_on_surfacebladeprofile(3,1);
        vec2(2) = - ders_on_surfacebladeprofile(5,1);
        direction_vector1_MB_B1 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/2;
        //gsInfo << direction_vector1_MB_B1 << "\n";

        pars2 << sec_pars(i), sec_pars(i),
                0.0, 1.0;
        //gsInfo << pars2 << "\n";
        m_surfaceBackExtensionRotated.eval_into(pars2, points_on_surfacebackrotatedprofile);
        m_surfaceBackExtensionRotated.deriv_into(pars2, ders_on_surfacebackrotatedprofile);
        //gsInfo << points_on_surfacebackrotatedprofile << "\n\n";
        //gsInfo << ders_on_surfacebackrotatedprofile << "\n\n";
        m_surfacePressureSideRotated.eval_into(pars2, points_on_surfacebladerotatedprofile);
        m_surfacePressureSideRotated.deriv_into(pars2, ders_on_surfacebladerotatedprofile);
        //gsInfo << points_on_surfacebladerotatedprofile << "\n\n";
        //gsInfo << ders_on_surfacebladerotatedprofile << "\n\n";
        vec1(0) = - ders_on_surfacebackrotatedprofile(1,0);
        vec1(1) = - ders_on_surfacebackrotatedprofile(3,0);
        vec1(2) = - ders_on_surfacebackrotatedprofile(5,0);
        vec2(0) = ders_on_surfacebladerotatedprofile(1,1);
        vec2(1) = ders_on_surfacebladerotatedprofile(3,1);
        vec2(2) = ders_on_surfacebladerotatedprofile(5,1);
        direction_vector2_MB_B2 = (vec1/euclideanNorm(vec1) + vec2/euclideanNorm(vec2))/2;
        //gsInfo << "a\n";

        coefs.row(0) = points_on_surfacebackprofile.col(0);
        coefs.row(1) = points_on_surfacebackprofile.col(0) + direction_vector1_MB_B1/3;
        coefs.row(2) = points_on_surfacebackrotatedprofile.col(0) - direction_vector2_MB_B2/3;
        coefs.row(3) = points_on_surfacebackrotatedprofile.col(0);
        //gsInfo << coefs << "\n";
        curve_aux = gsBSpline<T>( kv, coefs);
        curve_aux.eval_into(parameter_points, points_on_curve_aux);
        //surface_middleback_profiles_initial[i] = curve_aux;
        surface_middleback_profiles_initial.push_back(curve_aux);
        constructCylinderParametrization(-0.5, 0.3, cylinder_radii(i+1), cylinder);
        if (i==1) {
            gsWriteParaviewPoints<real_t>( points_on_curve_aux, "surface_middleback_profile_inital_points");
        }
        for (index_t i = 0; i < num_sample_pars; i++) {
            gsVector<T> point_on_curve_aux = points_on_curve_aux.col(i);
            init_sol << 0.9, 0.25;
            minimizePointSurfaceDistanceviaNR(point_on_curve_aux, cylinder, init_sol, sol, static_cast<real_t>(1e-6), 200, num_iter_aux, print_info);
            parameter_points_aux(0, i) = sol(0);
            parameter_points_aux(1, i) = sol(1);
            if (max_num_iter_aux < num_iter_aux) {
                max_num_iter_aux = num_iter_aux;
            }
            if ((num_iter_aux > 100) && (print_info)) {
                gsInfo << "Point with index " << i << " out of " << num_sample_pars << "does not converged. Solution: " << sol << "\n";
            }
        }
        if (print_info) {
            gsInfo << "Max number of NR iterations: " << max_num_iter_aux << "\n";
        }
        cylinder.eval_into(parameter_points_aux, points_on_curve_aux_projected);
        if (i==1) {
            gsWriteParaviewPoints<real_t>( points_on_curve_aux_projected, "surface_middleback_profile_projected_points");
            gsWriteParaview( cylinder, "cylinder1", 5000);
        }
        points_on_curve_aux_projected2 = points_on_curve_aux_projected.transpose();
        par_points_aux = centripetalParameterization(points_on_curve_aux_projected2);
        tangent(0) = direction_vector1_FM_B1(0);
        tangent(1) = direction_vector1_FM_B1(1);
        tangent(2) = direction_vector1_FM_B1(2);
        curve_aux2 = curveFittingWithBoundary(points_on_curve_aux_projected2, par_points_aux, kvfit, print_info);
    //    for (unsigned i = 0; i < res.size(); i++) {
    //        curve_aux.insertKnot(res[i]);
    //    }
        //surface_middleback_profiles[i] = curve_aux2;
        surface_middleback_profiles.push_back(curve_aux2);
    }
    //surface_middleback_profiles[0] = curve_MB_B;
    //surface_middleback_profiles[m_numBladeprofiles-1] = curve_MB_U;
    //surface_middleback_profiles_initial[0] = curve_MB_B;
    //surface_middleback_profiles_initial[m_numBladeprofiles-1] = curve_MB_U;
    surface_middleback_profiles.push_back(curve_MB_U);
    surface_middleback_profiles_initial.push_back(curve_MB_U);
    //gsInfo << "a\n";

    gsTensorBSpline<2, T> surface_middleback;
    computeLoftSurface(surface_middleback_profiles, kvfit, m_numBladeprofiles, kvloft, sec_pars, surface_middleback);

    //int num1 = kvloft.size() - kvloft.degree() - 1;
    //int num2 = kvfit.size() - kvfit.degree() - 1;
    gsMatrix<T> control_net_pom2 = surface_middleback.coefs();
    for (index_t i = 0; i < num1; i++) {
        control_net_pom2.row(i) = m_surfaceBackExtension.coefs().row(i);
        control_net_pom2.row((num2-1)*num1+i) = m_surfaceBackExtensionRotated.coefs().row(i);
    }
    //gsTensorBSplineBasis<2,T> basis(kvloft, kvfit);
    surface_middleback = gsTensorBSpline<2,T> (basis, control_net_pom2);

    m_surfaceMiddleBackDivider = surface_middleback;

    //curves3d = {};
    if (plot) {
        std::vector<gsGeometry<>*> curves3d;
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_middleback_profiles[i]);
        }
        gsWriteParaview( curves3d, "surface_middleback_profiles", 100);
        //curves3d = {};
        curves3d.clear();
        for (int i = 0; i < m_numBladeprofiles; i++) {
            curves3d.push_back(&surface_middleback_profiles_initial[i]);
        }
        gsWriteParaview( curves3d, "surface_middleback_profiles_inital", 100);

        gsWriteParaview( surface_middleback, "surface_middleback", 5000);

        gsMesh<> surface_middleback_mesh;
        surface_middleback.controlNet(surface_middleback_mesh);
        gsWriteParaview( surface_middleback_mesh, "surface_middleback_mesh");
    }


    return 0;
}

template<class T>
gsBSpline<T> KaplanTurbineRunnerWheelDomain<T>::computeCircularArc(gsMatrix<T> starting_point, gsMatrix<T> end_point) {

    gsMatrix<T> direction1(1,3);
    gsMatrix<T> direction2(1,3);
    gsMatrix<T> coefs(4,3);
    T radius;
    T norm_direction;

    radius = sqrt(pow(starting_point(0),2) + pow(starting_point(1),2));
    direction1(0) = starting_point(1)/radius;
    direction1(1) = -starting_point(0)/radius;
    direction1(2) = 0.0;
    direction2(0) = end_point(1)/radius;
    direction2(1) = -end_point(0)/radius;
    direction2(2) = 0.0;
    norm_direction = - radius * (4.0/3.0)*(tan((2*PI/m_numBladesInRunnerWheel)/4));
    coefs.row(0) = starting_point;
    coefs.row(1) = starting_point + norm_direction * direction1;
    coefs.row(2) = end_point - norm_direction * direction2;
    coefs.row(3) = end_point;

    gsKnotVector<T> kv(0,1,0,4);
    gsBSpline<T> curve;
    curve = gsBSpline<T> ( kv, coefs);

    return curve;
}


template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeFrontAndBackSurfaces(bool plot, bool print_info) {

    gsMatrix<T> surface_front_extension_controlnet(m_surfaceFrontExtension.coefs().rows(), m_surfaceFrontExtension.coefs().cols());
    gsMatrix<T> surface_front_extension_rotated_controlnet(m_surfaceFrontExtension.coefs().rows(), m_surfaceFrontExtension.coefs().cols());
    gsMatrix<T> surface_frontmost_controlnet(m_surfaceFrontExtension.coefs().rows(), m_surfaceFrontExtension.coefs().cols());
    gsMatrix<T> surface_back_extension_controlnet(m_surfaceBackExtension.coefs().rows(), m_surfaceBackExtension.coefs().cols());
    gsMatrix<T> surface_back_extension_rotated_controlnet(m_surfaceBackExtension.coefs().rows(), m_surfaceBackExtension.coefs().cols());
    gsMatrix<T> surface_backmost_controlnet(m_surfaceBackExtension.coefs().rows(), m_surfaceBackExtension.coefs().cols());
    gsMatrix<T> starting_point(1,3);
    gsMatrix<T> end_point(1,3);
    std::vector<T> res;
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    int num_cp_frontsurface_u = m_surfaceFrontExtension.knots(0).size() - m_surfaceFrontExtension.degree(0) - 1;
    int num_cp_frontsurface_v = m_surfaceFrontExtension.knots(1).size() - m_surfaceFrontExtension.degree(1) - 1;
    surface_front_extension_controlnet = m_surfaceFrontExtension.coefs();
    surface_front_extension_rotated_controlnet = m_surfaceFrontExtensionRotated.coefs();
    surface_back_extension_controlnet = m_surfaceBackExtension.coefs();
    surface_back_extension_rotated_controlnet = m_surfaceBackExtensionRotated.coefs();
    gsVector<T> divided_differences(num_cp_frontsurface_u);
    T dist;
    T distall;
    gsMatrix<T> vector(1,3);

    // ********************************************
    // FRONT-MOST SURFACE OF THE FRONT DOMAIN PATCH
    // ********************************************

    // front-most surface bottom curve
    starting_point = surface_front_extension_controlnet.row(0);
    end_point = surface_front_extension_rotated_controlnet.row(0);
    gsBSpline<T> curve_frontmost_bottom;
    curve_frontmost_bottom = computeCircularArc(starting_point, end_point);
    kvfit.difference(curve_frontmost_bottom.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_frontmost_bottom.insertKnot(res[i]);
    }

    // front-most surface upper curve
    starting_point = surface_front_extension_controlnet.row(num_cp_frontsurface_u-1);
    end_point = surface_front_extension_rotated_controlnet.row(num_cp_frontsurface_u-1);
    gsBSpline<T> curve_frontmost_up;
    curve_frontmost_up = computeCircularArc(starting_point, end_point);
    kvfit.difference(curve_frontmost_up.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_frontmost_up.insertKnot(res[i]);
    }

    // front-most surface
    divided_differences.setZero();

    vector = surface_front_extension_controlnet.row(num_cp_frontsurface_u-1) - surface_front_extension_controlnet.row(0);
    distall = euclideanNorm(vector);
    for (index_t i = 1; i < num_cp_frontsurface_u; i++) {
        vector = surface_front_extension_controlnet.row(i) - surface_front_extension_controlnet.row(0);
        dist = euclideanNorm(vector);
        divided_differences(i) = dist/distall;
    }
    //gsInfo << divided_differences << "\n";

    for (index_t i = 0; i < num_cp_frontsurface_v; i++) {
        for (index_t j = 0; j < num_cp_frontsurface_u; j++) {
            if (i==0) {
                surface_frontmost_controlnet.row(j) = surface_front_extension_controlnet.row(j);
            }
            else if (i == num_cp_frontsurface_v-1) {
                surface_frontmost_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+j) = surface_front_extension_rotated_controlnet.row(j);
            }
            else {
                surface_frontmost_controlnet.row(i*num_cp_frontsurface_u+j) = (1-divided_differences(j))*curve_frontmost_bottom.coefs().row(i) + divided_differences(j)*curve_frontmost_up.coefs().row(i);
            }
        }
    }
    //for (index_t i = 0; i < num_cp_frontsurface_u; i++) {
    //    gsInfo << surface_front_extension_controlnet.row(i) << "\n";
    //}
    //gsInfo << "\n";
    //for (index_t i = 0; i < num_cp_frontsurface_u; i++) {
    //    gsInfo << surface_front_extension_rotated_controlnet.row(i) << "\n";
    //}

    gsTensorBSpline<2, T> surface_frontmost;
    gsTensorBSplineBasis<2, T> basis(kvloft, kvfit);
    surface_frontmost = gsTensorBSpline<2, T> (basis, surface_frontmost_controlnet);
    m_surfaceFront = surface_frontmost;

    // ********************************************
    // BACK-MOST SURFACE OF THE BACK DOMAIN PATCH
    // ********************************************

    // back-most surface bottom curve
    starting_point = surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);
    end_point = surface_back_extension_rotated_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);
    gsBSpline<T> curve_backmost_bottom;
    curve_backmost_bottom = computeCircularArc(starting_point, end_point);
    kvfit.difference(curve_backmost_bottom.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_backmost_bottom.insertKnot(res[i]);
    }

    // front-most surface upper curve
    starting_point = surface_back_extension_controlnet.row(num_cp_frontsurface_v*num_cp_frontsurface_u-1);
    end_point = surface_back_extension_rotated_controlnet.row(num_cp_frontsurface_v*num_cp_frontsurface_u-1);
    gsBSpline<T> curve_backmost_up;
    curve_backmost_up = computeCircularArc(starting_point, end_point);
    kvfit.difference(curve_backmost_up.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_backmost_up.insertKnot(res[i]);
    }

    // front-most surface
    divided_differences.setZero();
    //gsInfo << "num_cp_frontsurface_u = " << num_cp_frontsurface_u << "\n";
    //gsInfo << "num_cp_frontsurface_v = " << num_cp_frontsurface_v << "\n";
    vector = surface_back_extension_controlnet.row(num_cp_frontsurface_v*num_cp_frontsurface_u-1) - surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);
    distall = euclideanNorm(vector);
    for (index_t i = 1; i < num_cp_frontsurface_u; i++) {
        vector = surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+i) - surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u);
        dist = euclideanNorm(vector);
        divided_differences(i) = dist/distall;
    }
    //gsInfo << divided_differences << "\n";

    for (index_t i = 0; i < num_cp_frontsurface_v; i++) {
        for (index_t j = 0; j < num_cp_frontsurface_u; j++) {
            if (i==0) {
                surface_backmost_controlnet.row(j) = surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+j);
            }
            else if (i == num_cp_frontsurface_v-1) {
                surface_backmost_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+j) = surface_back_extension_rotated_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+j);
            }
            else {
                surface_backmost_controlnet.row(i*num_cp_frontsurface_u+j) = (1-divided_differences(j))*curve_backmost_bottom.coefs().row(i) + divided_differences(j)*curve_backmost_up.coefs().row(i);
            }
        }
    }
    //for (index_t i = 0; i < num_cp_frontsurface_u; i++) {
    //    gsInfo << surface_back_extension_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+i) << "\n";
    //}
    //gsInfo << "\n";
    //for (index_t i = 0; i < num_cp_frontsurface_u; i++) {
    //    gsInfo << surface_back_extension_rotated_controlnet.row((num_cp_frontsurface_v-1)*num_cp_frontsurface_u+i) << "\n";
    //}

    gsTensorBSpline<2, T> surface_backmost;
    //gsTensorBSplineBasis<2, T> basis(kvloft, kvfit);
    surface_backmost = gsTensorBSpline<2, T> (basis, surface_backmost_controlnet);
    m_surfaceBack = surface_backmost;


    if (plot) {
        gsWriteParaview(curve_frontmost_bottom, "curve_frontmost_bottom", 100);
        gsWriteParaview(curve_frontmost_up, "curve_frontmost_up", 100);
        gsWriteParaview(surface_frontmost, "surface_frontmost", 5000);
        gsWriteParaview(curve_backmost_bottom, "curve_backmost_bottom", 100);
        gsWriteParaview(curve_backmost_up, "curve_backmost_up", 100);
        gsWriteParaview(surface_backmost, "surface_backmost", 5000);
    }


    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeTopOrBottomSurface(gsTensorBSpline<2,T> surface_front, gsTensorBSpline<2,T> surface_front_rotated, gsTensorBSpline<2,T> surface_frontmiddle, gsTensorBSpline<2,T> surface_frontmost,
                                                                 gsTensorBSpline<2,T> surface_inner, bool side, int domain_index, int num_sample_pars1, int num_sample_pars2, bool plot, bool print_info) {

    int num_sample_pars = 50;
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvfit2 = surface_front.knots(1);
    std::vector<T> res;
    gsMatrix<T> parameter_points(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front_rotated(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmiddle(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmost(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_front_rotated2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmiddle2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmost2(num_sample_pars, 2);
    gsMatrix<T> par_points_for_approx(1, num_sample_pars);

    gsMatrix<T> points_on_surface_front(3, num_sample_pars);
    gsMatrix<T> points_on_surface_front_rotated(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmiddle(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmost(3, num_sample_pars);
    gsVector<T> point_on_surface(3);

    gsBSpline<T> curve_left_parameterdomain;
    gsBSpline<T> curve_right_parameterdomain;
    gsBSpline<T> curve_top_parameterdomain;
    gsBSpline<T> curve_bottom_parameterdomain;

    if (domain_index == 1) {
        kvfit2.difference(kvfit, res);
        if (res.size() == 0) {
            surface_front.insertKnot(0.1, 1);
            surface_front.insertKnot(0.3, 1);
            surface_front.insertKnot(0.5, 1);
            surface_front.insertKnot(0.7, 1);
            surface_front.insertKnot(0.9, 1);
            surface_front_rotated.insertKnot(0.1, 1);
            surface_front_rotated.insertKnot(0.3, 1);
            surface_front_rotated.insertKnot(0.5, 1);
            surface_front_rotated.insertKnot(0.7, 1);
            surface_front_rotated.insertKnot(0.9, 1);
            kvfit2 = surface_front.knots(1);
        }
    }

    T k;
    if (side == 0) {
        k = 0.0;
    }
    else {
        k = 1.0;
    }
    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = k;
        parameter_points(1,i) = (T) i/(num_sample_pars-1);
    }

    surface_front.eval_into(parameter_points, points_on_surface_front);
    surface_front_rotated.eval_into(parameter_points, points_on_surface_front_rotated);
    if (domain_index < 2) {
        surface_frontmiddle.eval_into(parameter_points, points_on_surface_frontmiddle);
    }
    if (domain_index > 0) {
        surface_frontmost.eval_into(parameter_points, points_on_surface_frontmost);
    }
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    int num_iter_front = 0;
    int max_num_iter_front = 0;
    int num_iter_front_rotated = 0;
    int max_num_iter_front_rotated = 0;
    int num_iter_frontmiddle = 0;
    int max_num_iter_frontmiddle = 0;
    int num_iter_frontmost = 0;
    int max_num_iter_frontmost = 0;

    for (index_t i = 0; i < num_sample_pars; i++) {
        point_on_surface = points_on_surface_front.col(i);
        switch (domain_index) {
            case 0:
                init_sol << parameter_points(1,i)/3, 0.25;
                break;
            case 1:
                init_sol << 0.6, 0.25;
                break;
            case 2:
                init_sol << 0.95, 0.25;
                break;
            default:
                GISMO_ERROR("Wrong identification of domain patch!");
        }
        minimizePointSurfaceDistanceviaNR(point_on_surface, surface_inner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front, print_info);
        parameter_points_surface_front(0, i) = sol(0);
        parameter_points_surface_front(1, i) = sol(1);
        if (max_num_iter_front < num_iter_front) {
            max_num_iter_front = num_iter_front;
        }
        point_on_surface = points_on_surface_front_rotated.col(i);
        switch (domain_index) {
            case 0:
                init_sol << parameter_points(1,i)/4, 0.5;
                break;
            case 1:
                init_sol << 0.6, 0.5;
                break;
            case 2:
                init_sol << 0.9, 0.35;
                break;
            default:
                GISMO_ERROR("Wrong identification of domain patch!");
        }
        minimizePointSurfaceDistanceviaNR(point_on_surface, surface_inner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front_rotated, print_info);
        parameter_points_surface_front_rotated(0, i) = sol(0);
        parameter_points_surface_front_rotated(1, i) = sol(1);
        if (max_num_iter_front_rotated < num_iter_front_rotated) {
            max_num_iter_front_rotated = num_iter_front_rotated;
        }
        if (domain_index < 2) {
            point_on_surface = points_on_surface_frontmiddle.col(i);
            switch (domain_index) {
                case 0:
                    init_sol << 0.4, 0.3+parameter_points(1,i)/4;
                    break;
                case 1:
                    init_sol << 0.8, 0.25+parameter_points(1,i)/4;
                    break;
                case 2:

                    break;
                default:
                    GISMO_ERROR("Wrong identification of domain patch!");
            }
            minimizePointSurfaceDistanceviaNR(point_on_surface, surface_inner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmiddle, print_info);
            parameter_points_surface_frontmiddle(0, i) = sol(0);
            parameter_points_surface_frontmiddle(1, i) = sol(1);
            if (max_num_iter_frontmiddle < num_iter_frontmiddle) {
                max_num_iter_frontmiddle = num_iter_frontmiddle;
            }
        }
        if (domain_index > 0) {
            point_on_surface = points_on_surface_frontmost.col(i);
            switch (domain_index) {
                case 1:
                    init_sol << 0.4, 0.3+parameter_points(1,i)/4;
                    break;
                case 2:
                    init_sol << 0.8, 0.25+parameter_points(1,i)/4;
                    break;
                default:
                    GISMO_ERROR("Wrong identification of domain patch!");
            }

            minimizePointSurfaceDistanceviaNR(point_on_surface, surface_inner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmost, print_info);
            parameter_points_surface_frontmost(0, i) = sol(0);
            parameter_points_surface_frontmost(1, i) = sol(1);
            if (max_num_iter_frontmost < num_iter_frontmost) {
                max_num_iter_frontmost = num_iter_frontmost;
            }
        }

    }
    if (domain_index > 0) {
        parameter_points_surface_frontmost.col(0) = parameter_points_surface_front.col(0);
        parameter_points_surface_frontmost.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(0);
    }
    if (domain_index < 2) {
        parameter_points_surface_frontmiddle.col(0) = parameter_points_surface_front.col(num_sample_pars-1);
        parameter_points_surface_frontmiddle.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(num_sample_pars-1);
    }

    parameter_points_surface_front2 = parameter_points_surface_front.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front2);
    curve_top_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front2, par_points_for_approx, kvfit2, print_info);
    parameter_points_surface_front_rotated2 = parameter_points_surface_front_rotated.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front_rotated2);
    curve_bottom_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front_rotated2, par_points_for_approx, kvfit2, print_info);

    if (domain_index == 2) {
        gsKnotVector<T> kv(0,1,0,2);
        gsMatrix<T> coefs(2,2);
        coefs.row(0) = parameter_points_surface_front.col(num_sample_pars-1);
        coefs.row(1) = parameter_points_surface_front_rotated.col(num_sample_pars-1);
        curve_right_parameterdomain = gsBSpline<T> ( kv, coefs);
        curve_right_parameterdomain.degreeElevate(2);
        kvfit.difference(curve_right_parameterdomain.knots(0), res);
        for (unsigned i = 0; i < res.size(); i++) {
            curve_right_parameterdomain.insertKnot(res[i]);
        }
    }
    else {
        parameter_points_surface_frontmiddle2 = parameter_points_surface_frontmiddle.transpose();
        par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmiddle2);
        curve_right_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmiddle2, par_points_for_approx, kvfit, print_info);
    }

    if (domain_index == 0) {
        gsKnotVector<T> kv(0,1,0,2);
        gsMatrix<T> coefs(2,2);
        coefs.row(0) = parameter_points_surface_front.col(0);
        coefs.row(1) = parameter_points_surface_front_rotated.col(0);
        curve_left_parameterdomain = gsBSpline<T> ( kv, coefs);
        curve_left_parameterdomain.degreeElevate(2);
        kvfit.difference(curve_left_parameterdomain.knots(0), res);
        for (unsigned i = 0; i < res.size(); i++) {
            curve_left_parameterdomain.insertKnot(res[i]);
        }
    }
    else {
        parameter_points_surface_frontmost2 = parameter_points_surface_frontmost.transpose();
        par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmost2);
        curve_left_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmost2, par_points_for_approx, kvfit, print_info);
    }

    gsMatrix<T> surface_bottom_parameter_domain_controlnet(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    gsMatrix<T> surface_bottom_parameter_domain_controlnet2(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    gsMatrix<T> surface_bottom_parameter_domain_controlnet3(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    if ((side == 0) && (domain_index == 1)) {
        springModelPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet, print_info);
    }
    //else if ((side == 1) && (domain_index == 1)) {
    //    springModelPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet2, print_info);
    //    discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet3, print_info);
    //    surface_bottom_parameter_domain_controlnet = (surface_bottom_parameter_domain_controlnet2 + surface_bottom_parameter_domain_controlnet3)/2;
    //}
    //else if ((side == 1) && (domain_index == 2)) {
    //    springModelPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet2, print_info);
    //    discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet3, print_info);
    //    surface_bottom_parameter_domain_controlnet = (surface_bottom_parameter_domain_controlnet2 + surface_bottom_parameter_domain_controlnet3)/2;
    //}
    else {
        discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet, print_info);
    }
    //gsInfo << surface_bottom_parameter_domain_controlnet << "\n";
    gsTensorBSpline<2, T> surface_bottom_parameter_domain;
    gsTensorBSplineBasis<2, T> basis(kvfit2, kvfit);
    surface_bottom_parameter_domain = gsTensorBSpline<2, T> (basis, surface_bottom_parameter_domain_controlnet);
    /*if ((side == 0) && (domain_index == 1)) {
        springModelPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet, print_info);
        surface_bottom_parameter_domain = gsTensorBSpline<2, T> (basis, surface_bottom_parameter_domain_controlnet);
        gsWriteParaview(surface_bottom_parameter_domain, "MD_bottom_surface_parameterdomain", 5000);
        gsMesh<> surface_bottom_parameter_domain_mesh;
        surface_bottom_parameter_domain.controlNet( surface_bottom_parameter_domain_mesh);
        gsWriteParaview(  surface_bottom_parameter_domain_mesh, "MD_bottom_surface_parameterdomain_mesh");
    }*/
    switch (domain_index) {
        case 0:
            if (side == 0) {
                gsWriteParaview(surface_bottom_parameter_domain, "FD_bottom_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "FD_bottom_surface_parameterdomain_mesh");
            }
            else {
                gsWriteParaview(surface_bottom_parameter_domain, "FD_up_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "FD_up_surface_parameterdomain_mesh");
            }
            break;
        case 1:
            if (side == 0) {
                gsWriteParaview(surface_bottom_parameter_domain, "MD_bottom_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "MD_bottom_surface_parameterdomain_mesh");
            }
            else {
                gsWriteParaview(surface_bottom_parameter_domain, "MD_up_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "MD_up_surface_parameterdomain_mesh");
            }
            break;
        case 2:
            if (side == 0) {
                gsWriteParaview(surface_bottom_parameter_domain, "BD_bottom_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "BD_bottom_surface_parameterdomain_mesh");
            }
            else {
                gsWriteParaview(surface_bottom_parameter_domain, "BD_up_surface_parameterdomain", 5000);
                gsMesh<> surface_bottom_parameter_domain_mesh;
                surface_bottom_parameter_domain.controlNet(surface_bottom_parameter_domain_mesh);
                gsWriteParaview( surface_bottom_parameter_domain_mesh, "BD_up_surface_parameterdomain_mesh");
            }
            break;
        default:
            GISMO_ERROR("Wrong identification of domain patch!");
    }

    // Sampling point on bottom surface for approximation
    //int num_sample_pars1 = 30;
    //int num_sample_pars2 = 20;
    gsMatrix<T> sampling_par_points_surface(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    sampling_par_points_surface.setZero();
    gsMatrix<T> surface_bottom_parameter_domain_points(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_bottom_points(3, (num_sample_pars1-2)*(num_sample_pars2-2));
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface(0, (i-1)*(num_sample_pars2-2)+j-1) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface(1, (i-1)*(num_sample_pars2-2)+j-1) = (T) j/(num_sample_pars2-1);
        }
    }
    //gsInfo << sampling_par_points_surface << "\n";
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface, surface_bottom_parameter_domain_points);
    surface_inner.eval_into(surface_bottom_parameter_domain_points, surface_bottom_points);
    gsMatrix<T> sampling_par_points_surface2(2, (num_sample_pars1)*(num_sample_pars2));
    sampling_par_points_surface2.setZero();
    gsMatrix<T> surface_bottom_parameter_domain_points2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_points2(3, (num_sample_pars1)*(num_sample_pars2));
    for (index_t i = 0; i < num_sample_pars1; i++) {
        for (index_t j = 0; j < num_sample_pars2; j++) {
            sampling_par_points_surface2(0, i*(num_sample_pars2)+j) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface2(1, i*(num_sample_pars2)+j) = (T) j/(num_sample_pars2-1);
        }
    }
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface2, surface_bottom_parameter_domain_points2);
    surface_inner.eval_into(surface_bottom_parameter_domain_points2, surface_bottom_points2);

    // Finding suitable parameters for approximated points with the help of centripetal/chordal parameterization and averaging
    gsMatrix<T> sampling_par_points_surface3(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> pom(3, num_sample_pars2);
    gsMatrix<T> pom2(num_sample_pars2, 3);
    gsMatrix<T> pom_pars(1, num_sample_pars2);
    gsMatrix<T> pom_pars_sum(1, num_sample_pars2);
    pom_pars_sum.setZero();
    for (index_t i = 0; i < num_sample_pars1; i++) {
            for (index_t j = 0; j < num_sample_pars2; j++) {
                pom.col(j) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            gsWriteParaviewPoints(pom, "MD_bottom_surface_pompoints");
            pom2 = pom.transpose();
            switch (domain_index) {
                case 0:
                    pom_pars = chordalParameterization(pom2);
                    break;
                case 1:
                    pom_pars = chordalParameterization(pom2);
                    break;
                case 2:
                    pom_pars = chordalParameterization(pom2);
                    break;
                default:
                    GISMO_ERROR("Wrong identification of domain patch!");
            }
            pom_pars_sum += pom_pars;
            //gsInfo << "Parametrizace pro radu " << i << ": \n" << pom_pars << "\n";
        }
    //gsInfo << pom_pars_sum/num_sample_pars1 << "\n";
    gsMatrix<T> pom3(3, num_sample_pars1);
    gsMatrix<T> pom4(num_sample_pars1, 3);
    gsMatrix<T> pom_pars2(1, num_sample_pars1);
    gsMatrix<T> pom_pars_sum2(1, num_sample_pars1);
    pom_pars_sum2.setZero();
    for (index_t j = 0; j < num_sample_pars2; j++) {
            for (index_t i = 0; i < num_sample_pars1; i++) {
                pom3.col(i) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            gsWriteParaviewPoints(pom3, "MD_bottom_surface_pompoints2");
            pom4 = pom3.transpose();
            switch (domain_index) {
                case 0:
                    pom_pars2 = chordalParameterization(pom4);
                    break;
                case 1:
                    pom_pars2 = chordalParameterization(pom4);
                    break;
                case 2:
                    pom_pars2 = chordalParameterization(pom4);
                    break;
                default:
                    GISMO_ERROR("Wrong identification of domain patch!");
            }
            pom_pars_sum2 += pom_pars2;
            //gsInfo << "Parametrizace pro sloupec " << j << ": \n" << pom_pars2 << "\n";
        }
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface3(0, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum2(i)/num_sample_pars2;
            sampling_par_points_surface3(1, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum(j)/num_sample_pars1;
        }
    }
    //gsInfo << pom_pars_sum2/num_sample_pars2 << "\n";

    // Select initial control net for surface approaximation
    gsTensorBSpline<2, T> surface_bottom_initial;
    //selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_parameter_domain, surface_inner, surface_bottom_initial, kvfit, kvfit, 0, print_info);
    selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_initial, kvfit2, kvfit, side, print_info);
    gsMatrix<T> surface_bottom_initial_controlnet = surface_bottom_initial.coefs();
    //gsMatrix<T> difference(surface_bottom_initial_controlnet.rows(), surface_bottom_initial_controlnet.cols());
    //difference = surface_bottom_initial.coefs() - surface_bottom_initial2.coefs();
    //gsInfo << difference << "\n";
    gsMatrix<index_t> boundaryind = surface_bottom_initial.basis().allBoundary();
    gsMatrix<T> boundarypoints(boundaryind.size(), 3);
    for (index_t i = 0; i < boundaryind.size(); i++) {
        boundarypoints.row(i) = surface_bottom_initial_controlnet.row(boundaryind(i));
    }
    gsMatrix<T> surface_bottom_points3 = surface_bottom_points.transpose();
    switch (domain_index) {
        case 0:
            surface_bottom_initial = surfaceFittingLSQWithBoundary(surface_bottom_points3, sampling_par_points_surface3, kvfit2, kvfit, boundarypoints, print_info);
            break;
        case 1:
            surface_bottom_initial = surfaceFittingLSQWithBoundary(surface_bottom_points3, sampling_par_points_surface, kvfit2, kvfit, boundarypoints, print_info);
            surface_bottom_initial_controlnet = surface_bottom_initial.coefs();
            break;
        case 2:
            surface_bottom_initial = surfaceFittingLSQWithBoundary(surface_bottom_points3, sampling_par_points_surface, kvfit2, kvfit, boundarypoints, print_info);
            //surface_bottom_initial_controlnet = surface_bottom_initial.coefs();
            break;
        default:
            GISMO_ERROR("Wrong identification of domain patch!");
    }
    //surface_bottom_initial_controlnet = surface_bottom_initial.coefs();

    // Surface approximation
    T max_approx_error = 1e-4;
    T avg_approx_error = 1e-5;
    gsTensorBSpline<2, T> surface_bottom;
    switch (domain_index) {
        case 0:
            surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit2, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface, surface_bottom, max_approx_error, avg_approx_error, 20, print_info);
            break;
        case 1:
            surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit2, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface, surface_bottom, max_approx_error, avg_approx_error, 10, print_info);
            break;
        case 2:
            if (side == 0){
                surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit2, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface, surface_bottom, max_approx_error, avg_approx_error, 10, print_info);
            }
            else {
                surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit2, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface, surface_bottom, max_approx_error, avg_approx_error, 1, print_info);
            }
            break;
        default:
            GISMO_ERROR("Wrong identification of domain patch!");
    }
    //gsInfo << surface_bottom << "\n";

    //
    switch (domain_index) {
        case 0:
            if (side == 0) {
                m_surfaceBottomFrontDomain = surface_bottom;
                gsWriteParaview(surface_bottom_initial, "FD_bottom_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "FD_bottom_surface_initial_mesh");
                gsWriteParaview(surface_bottom, "FD_bottom_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "FD_bottom_surface_final_mesh");
            }
            else {
                m_surfaceTopFrontDomain = surface_bottom;
                gsWriteParaview(surface_bottom_initial, "FD_up_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "FD_bottom_up_initial_mesh");
                gsWriteParaview(surface_bottom, "FD_up_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "FD_up_surface_final_mesh");
                gsWriteParaviewPoints(surface_bottom_points, "FD_up_surface_samplepoints");
            }
            break;
        case 1:
            m_surfaceSuctionSide = surface_front;
            m_surfacePressureSideRotated = surface_front_rotated;
            if (side == 0) {
                m_surfaceBottomMiddleDomain = surface_bottom;
                //gsWriteParaviewPoints(surface_bottom_points, "MD_bottom_surface_samplepoints");
                gsWriteParaview(surface_bottom_initial, "MD_bottom_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "MD_bottom_surface_initial_mesh");
                gsWriteParaview(surface_bottom, "MD_bottom_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "MD_bottom_surface_final_mesh");
            }
            else {
                m_surfaceTopMiddleDomain = surface_bottom;
                gsWriteParaviewPoints(surface_bottom_points, "MD_up_surface_samplepoints");
                gsWriteParaview(surface_bottom_initial, "MD_up_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "MD_up_surface_initial_mesh");
                gsWriteParaview(surface_bottom, "MD_up_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "MD_up_surface_final_mesh");
            }
            break;
        case 2:
            if (side == 0) {
                m_surfaceBottomBackDomain = surface_bottom;
                gsWriteParaview(surface_bottom_initial, "BD_bottom_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "BD_bottom_surface_initial_mesh");
                gsWriteParaview(surface_bottom, "BD_bottom_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "BD_bottom_surface_final_mesh");
            }
            else {
                m_surfaceTopBackDomain = surface_bottom;
                gsWriteParaview(surface_bottom_initial, "BD_up_surface_initial", 5000);
                gsMesh<> FD_bottom_surface_initial_mesh;
                surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
                gsWriteParaview( FD_bottom_surface_initial_mesh, "BD_up_surface_initial_mesh");
                gsWriteParaview(surface_bottom, "BD_up_surface_final", 5000);
                gsMesh<> FD_bottom_surface_mesh;
                surface_bottom.controlNet(FD_bottom_surface_mesh);
                gsWriteParaview( FD_bottom_surface_mesh, "BD_up_surface_final_mesh");
                gsWriteParaviewPoints(surface_bottom_points, "BD_up_surface_samplepoints");
            }
            break;
        default:
            GISMO_ERROR("Wrong identification of domain patch!");
    }
    //

    return 0;
}

/*
template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeTopAndBottomSurfaces(bool plot, bool print_info) {

    // *****************************
    // FRONT DOMAIN - BOTTOM SURFACE
    // *****************************
    {
    int num_sample_pars = 50;
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);
    gsTensorBSpline<2,T> surface_front = m_surfaceFrontExtension;
    gsTensorBSpline<2,T> surface_front_rotated = m_surfaceFrontExtensionRotated;
    gsTensorBSpline<2,T> surface_frontmiddle = m_surfaceFrontMiddleDivider;
    gsTensorBSpline<2,T> surface_frontmost = m_surfaceFront;
    gsTensorBSpline<2,T> surface_inner = m_surfaceInner;
    gsMatrix<T> parameter_points(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front_rotated(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmiddle(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_front_rotated2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmiddle2(num_sample_pars, 2);
    gsMatrix<T> par_points_for_approx(1, num_sample_pars);

    gsMatrix<T> points_on_surface_front(3, num_sample_pars);
    gsMatrix<T> points_on_surface_front_rotated(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmiddle(3, num_sample_pars);
    gsVector<T> point_on_surface(3);

    gsBSpline<T> curve_right_parameterdomain;
    gsBSpline<T> curve_top_parameterdomain;
    gsBSpline<T> curve_bottom_parameterdomain;

    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = 0.0;
        parameter_points(1,i) = (T) i/(num_sample_pars-1);
    }

    surface_front.eval_into(parameter_points, points_on_surface_front);
    surface_front_rotated.eval_into(parameter_points, points_on_surface_front_rotated);
    surface_frontmiddle.eval_into(parameter_points, points_on_surface_frontmiddle);
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    int num_iter_front = 0;
    int max_num_iter_front = 0;
    int num_iter_front_rotated = 0;
    int max_num_iter_front_rotated = 0;
    int num_iter_frontmiddle = 0;
    int max_num_iter_frontmiddle = 0;

    for (index_t i = 0; i < num_sample_pars; i++) {
        point_on_surface = points_on_surface_front.col(i);
        init_sol << parameter_points(1,i)/3, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front, print_info);
        parameter_points_surface_front(0, i) = sol(0);
        parameter_points_surface_front(1, i) = sol(1);
        if (max_num_iter_front < num_iter_front) {
            max_num_iter_front = num_iter_front;
        }
        point_on_surface = points_on_surface_front_rotated.col(i);
        init_sol << parameter_points(1,i)/4, 0.5;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front_rotated, print_info);
        parameter_points_surface_front_rotated(0, i) = sol(0);
        parameter_points_surface_front_rotated(1, i) = sol(1);
        if (max_num_iter_front_rotated < num_iter_front_rotated) {
            max_num_iter_front_rotated = num_iter_front_rotated;
        }
        point_on_surface = points_on_surface_frontmiddle.col(i);
        init_sol << 0.4, 0.3+parameter_points(1,i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmiddle, print_info);
        parameter_points_surface_frontmiddle(0, i) = sol(0);
        parameter_points_surface_frontmiddle(1, i) = sol(1);
        if (max_num_iter_frontmiddle < num_iter_frontmiddle) {
            max_num_iter_frontmiddle = num_iter_frontmiddle;
        }

    }
    parameter_points_surface_frontmiddle.col(0) = parameter_points_surface_front.col(num_sample_pars-1);
    parameter_points_surface_frontmiddle.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(num_sample_pars-1);
    parameter_points_surface_front2 = parameter_points_surface_front.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front2);
    curve_top_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front2, par_points_for_approx, kvfit, print_info);
    parameter_points_surface_front_rotated2 = parameter_points_surface_front_rotated.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front_rotated2);
    curve_bottom_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front_rotated2, par_points_for_approx, kvfit, print_info);
    parameter_points_surface_frontmiddle2 = parameter_points_surface_frontmiddle.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmiddle2);
    curve_right_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmiddle2, par_points_for_approx, kvfit, print_info);

    std::vector<T> res;
    gsKnotVector<T> kv(0,1,0,2);
    gsMatrix<T> coefs(2,2);
    gsBSpline<T> curve_left_parameterdomain;
    coefs.row(0) = parameter_points_surface_front.col(0);
    coefs.row(1) = parameter_points_surface_front_rotated.col(0);
    curve_left_parameterdomain = gsBSpline<T> ( kv, coefs);
    curve_left_parameterdomain.degreeElevate(2);
    kvfit.difference(curve_left_parameterdomain.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_left_parameterdomain.insertKnot(res[i]);
    }

    gsMatrix<T> surface_bottom_parameter_domain_controlnet(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet, print_info);
    //gsInfo << surface_bottom_parameter_domain_controlnet << "\n";
    gsTensorBSpline<2, T> surface_bottom_parameter_domain;
    gsTensorBSplineBasis<2, T> basis(kvfit, kvfit);
    surface_bottom_parameter_domain = gsTensorBSpline<2, T> (basis, surface_bottom_parameter_domain_controlnet);

    // Sampling point on bottom surface for approximation
    int num_sample_pars1 = 30;
    int num_sample_pars2 = 20;
    gsMatrix<T> sampling_par_points_surface(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_bottom_parameter_domain_points(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_bottom_points(3, (num_sample_pars1-2)*(num_sample_pars2-2));
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface(0, (i-1)*(num_sample_pars2-2)+j-1) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface(1, (i-1)*(num_sample_pars2-2)+j-1) = (T) j/(num_sample_pars2-1);
        }
    }
    //gsInfo << sampling_par_points_surface << "\n";
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface, surface_bottom_parameter_domain_points);
    surface_inner.eval_into(surface_bottom_parameter_domain_points, surface_bottom_points);
    gsMatrix<T> sampling_par_points_surface2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_parameter_domain_points2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_points2(3, (num_sample_pars1)*(num_sample_pars2));
    for (index_t i = 0; i < num_sample_pars1; i++) {
        for (index_t j = 0; j < num_sample_pars2; j++) {
            sampling_par_points_surface2(0, i*(num_sample_pars2)+j) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface2(1, i*(num_sample_pars2)+j) = (T) j/(num_sample_pars2-1);
        }
    }
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface2, surface_bottom_parameter_domain_points2);
    surface_inner.eval_into(surface_bottom_parameter_domain_points2, surface_bottom_points2);

    // Finding suitable parameters for approximated points with the help of centripetal/chordal parameterization and averaging
    gsMatrix<T> sampling_par_points_surface3(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> pom(3, num_sample_pars2);
    gsMatrix<T> pom2(num_sample_pars2, 3);
    gsMatrix<T> pom_pars(1, num_sample_pars2);
    gsMatrix<T> pom_pars_sum(1, num_sample_pars2);
    pom_pars_sum.setZero();
    for (index_t i = 0; i < num_sample_pars1; i++) {
            for (index_t j = 0; j < num_sample_pars2; j++) {
                pom.col(j) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom, "MD_bottom_surface_pompoints");
            pom2 = pom.transpose();
            pom_pars = chordalParameterization(pom2);
            pom_pars_sum += pom_pars;
            //gsInfo << "Parametrizace pro radu " << i << ": \n" << pom_pars << "\n";
        }
    gsMatrix<T> pom3(3, num_sample_pars1);
    gsMatrix<T> pom4(num_sample_pars1, 3);
    gsMatrix<T> pom_pars2(1, num_sample_pars1);
    gsMatrix<T> pom_pars_sum2(1, num_sample_pars1);
    for (index_t j = 0; j < num_sample_pars2; j++) {
            for (index_t i = 0; i < num_sample_pars1; i++) {
                pom3.col(i) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom3, "MD_bottom_surface_pompoints2");
            pom4 = pom3.transpose();
            pom_pars2 = centripetalParameterization(pom4);
            pom_pars_sum2 += pom_pars2;
            //gsInfo << "Parametrizace pro sloupec " << j << ": \n" << pom_pars2 << "\n";
        }
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface3(0, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum2(i)/num_sample_pars2;
            sampling_par_points_surface3(1, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum(j)/num_sample_pars1;
        }
    }

    // Select initial control net for surface approaximation
    gsTensorBSpline<2, T> surface_bottom_initial;
    //selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_parameter_domain, surface_inner, surface_bottom_initial, kvfit, kvfit, 0, print_info);
    selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_initial, kvfit, kvfit, 0, print_info);
    gsMatrix<T> surface_bottom_initial_controlnet = surface_bottom_initial.coefs();
    //gsMatrix<T> difference(surface_bottom_initial_controlnet.rows(), surface_bottom_initial_controlnet.cols());
    //difference = surface_bottom_initial.coefs() - surface_bottom_initial2.coefs();
    //gsInfo << difference << "\n";
    gsMatrix<index_t> boundaryind = surface_bottom_initial.basis().allBoundary();
    gsMatrix<T> boundarypoints(boundaryind.size(), 3);
    for (index_t i = 0; i < boundaryind.size(); i++) {
        boundarypoints.row(i) = surface_bottom_initial_controlnet.row(boundaryind(i));
    }
    gsMatrix<T> surface_bottom_points3 = surface_bottom_points.transpose();
    surface_bottom_initial = surfaceFittingLSQWithBoundary(surface_bottom_points3, sampling_par_points_surface3, kvfit, kvfit, boundarypoints, print_info);
    surface_bottom_initial_controlnet = surface_bottom_initial.coefs();

    // Surface approximation
    T max_approx_error;
    T avg_approx_error;
    gsTensorBSpline<2, T> surface_bottom;
    surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface3, surface_bottom, max_approx_error, avg_approx_error, 10, print_info);
    m_surfaceBottomFrontDomain = surface_bottom;

    //
    gsWriteParaview(surface_bottom_initial, "FD_bottom_surface_initial", 5000);
    gsMesh<> FD_bottom_surface_initial_mesh;
    surface_bottom_initial.controlNet(FD_bottom_surface_initial_mesh);
    gsWriteParaview( FD_bottom_surface_initial_mesh, "FD_bottom_surface_initial_mesh");

    gsWriteParaview(curve_left_parameterdomain, "FD_curve_left_parameterdomain");
    gsWriteParaview(curve_right_parameterdomain, "FD_curve_right_parameterdomain");
    gsWriteParaview(curve_top_parameterdomain, "FD_curve_top_parameterdomain");
    gsWriteParaview(curve_bottom_parameterdomain, "FD_curve_bottom_parameterdomain");
    gsWriteParaview(surface_bottom_parameter_domain, "FD_surface_bottom_parameterdomain", 5000);
    //

    }

    // *****************************
    // FRONT DOMAIN - UPPER SURFACE
    // *****************************
    {
    int num_sample_pars = 50;
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);
    gsTensorBSpline<2,T> surface_front = m_surfaceFrontExtension;
    gsTensorBSpline<2,T> surface_front_rotated = m_surfaceFrontExtensionRotated;
    gsTensorBSpline<2,T> surface_frontmiddle = m_surfaceFrontMiddleDivider;
    gsTensorBSpline<2,T> surface_frontmost = m_surfaceFront;
    gsTensorBSpline<2,T> surface_outer = m_surfaceOuter;
    gsMatrix<T> parameter_points(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front_rotated(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmiddle(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_front_rotated2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmiddle2(num_sample_pars, 2);
    gsMatrix<T> par_points_for_approx(1, num_sample_pars);

    gsMatrix<T> points_on_surface_front(3, num_sample_pars);
    gsMatrix<T> points_on_surface_front_rotated(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmiddle(3, num_sample_pars);
    gsVector<T> point_on_surface(3);

    gsBSpline<T> curve_right_parameterdomain;
    gsBSpline<T> curve_top_parameterdomain;
    gsBSpline<T> curve_bottom_parameterdomain;

    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = 1.0;
        parameter_points(1,i) = (T) i/(num_sample_pars-1);
    }

    surface_front.eval_into(parameter_points, points_on_surface_front);
    surface_front_rotated.eval_into(parameter_points, points_on_surface_front_rotated);
    surface_frontmiddle.eval_into(parameter_points, points_on_surface_frontmiddle);
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    int num_iter_front = 0;
    int max_num_iter_front = 0;
    int num_iter_front_rotated = 0;
    int max_num_iter_front_rotated = 0;
    int num_iter_frontmiddle = 0;
    int max_num_iter_frontmiddle = 0;

    for (index_t i = 0; i < num_sample_pars; i++) {
        point_on_surface = points_on_surface_front.col(i);
        init_sol << parameter_points(1,i)/3, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front, print_info);
        parameter_points_surface_front(0, i) = sol(0);
        parameter_points_surface_front(1, i) = sol(1);
        if (max_num_iter_front < num_iter_front) {
            max_num_iter_front = num_iter_front;
        }
        point_on_surface = points_on_surface_front_rotated.col(i);
        init_sol << parameter_points(1,i)/4, 0.5;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front_rotated, print_info);
        parameter_points_surface_front_rotated(0, i) = sol(0);
        parameter_points_surface_front_rotated(1, i) = sol(1);
        if (max_num_iter_front_rotated < num_iter_front_rotated) {
            max_num_iter_front_rotated = num_iter_front_rotated;
        }
        point_on_surface = points_on_surface_frontmiddle.col(i);
        init_sol << 0.4, 0.3+parameter_points(1,i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceOuter, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmiddle, print_info);
        parameter_points_surface_frontmiddle(0, i) = sol(0);
        parameter_points_surface_frontmiddle(1, i) = sol(1);
        if (max_num_iter_frontmiddle < num_iter_frontmiddle) {
            max_num_iter_frontmiddle = num_iter_frontmiddle;
        }

    }
    parameter_points_surface_frontmiddle.col(0) = parameter_points_surface_front.col(num_sample_pars-1);
    parameter_points_surface_frontmiddle.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(num_sample_pars-1);
    parameter_points_surface_front2 = parameter_points_surface_front.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front2);
    curve_top_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front2, par_points_for_approx, kvfit, print_info);
    parameter_points_surface_front_rotated2 = parameter_points_surface_front_rotated.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front_rotated2);
    curve_bottom_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front_rotated2, par_points_for_approx, kvfit, print_info);
    parameter_points_surface_frontmiddle2 = parameter_points_surface_frontmiddle.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmiddle2);
    curve_right_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmiddle2, par_points_for_approx, kvfit, print_info);

    std::vector<T> res;
    gsKnotVector<T> kv(0,1,0,2);
    gsMatrix<T> coefs(2,2);
    gsBSpline<T> curve_left_parameterdomain;
    coefs.row(0) = parameter_points_surface_front.col(0);
    coefs.row(1) = parameter_points_surface_front_rotated.col(0);
    curve_left_parameterdomain = gsBSpline<T> ( kv, coefs);
    curve_left_parameterdomain.degreeElevate(2);
    kvfit.difference(curve_left_parameterdomain.knots(0), res);
    for (unsigned i = 0; i < res.size(); i++) {
        curve_left_parameterdomain.insertKnot(res[i]);
    }

    gsMatrix<T> surface_up_parameter_domain_controlnet(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_up_parameter_domain_controlnet, print_info);
    //gsInfo << surface_bottom_parameter_domain_controlnet << "\n";
    gsTensorBSpline<2, T> surface_up_parameter_domain;
    gsTensorBSplineBasis<2, T> basis(kvfit, kvfit);
    surface_up_parameter_domain = gsTensorBSpline<2, T> (basis, surface_up_parameter_domain_controlnet);

    // Sampling point on bottom surface for approximation
    int num_sample_pars1 = 10;
    int num_sample_pars2 = 50;
    gsMatrix<T> sampling_par_points_surface(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_up_parameter_domain_points(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_up_points(3, (num_sample_pars1-2)*(num_sample_pars2-2));
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface(0, (i-1)*(num_sample_pars2-2)+j-1) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface(1, (i-1)*(num_sample_pars2-2)+j-1) = (T) j/(num_sample_pars2-1);
        }
    }
    //gsInfo << sampling_par_points_surface << "\n";
    surface_up_parameter_domain.eval_into(sampling_par_points_surface, surface_up_parameter_domain_points);
    surface_outer.eval_into(surface_up_parameter_domain_points, surface_up_points);
    gsMatrix<T> sampling_par_points_surface2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_parameter_domain_points2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_points2(3, (num_sample_pars1)*(num_sample_pars2));
    for (index_t i = 0; i < num_sample_pars1; i++) {
        for (index_t j = 0; j < num_sample_pars2; j++) {
            sampling_par_points_surface2(0, i*(num_sample_pars2)+j) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface2(1, i*(num_sample_pars2)+j) = (T) j/(num_sample_pars2-1);
        }
    }
    surface_up_parameter_domain.eval_into(sampling_par_points_surface2, surface_bottom_parameter_domain_points2);
    surface_outer.eval_into(surface_bottom_parameter_domain_points2, surface_bottom_points2);

    // Finding suitable parameters for approximated points with the help of centripetal/chordal parameterization and averaging
    gsMatrix<T> sampling_par_points_surface3(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> pom(3, num_sample_pars2);
    gsMatrix<T> pom2(num_sample_pars2, 3);
    gsMatrix<T> pom_pars(1, num_sample_pars2);
    gsMatrix<T> pom_pars_sum(1, num_sample_pars2);
    pom_pars_sum.setZero();
    for (index_t i = 0; i < num_sample_pars1; i++) {
            for (index_t j = 0; j < num_sample_pars2; j++) {
                pom.col(j) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom, "MD_bottom_surface_pompoints");
            pom2 = pom.transpose();
            pom_pars = chordalParameterization(pom2);
            pom_pars_sum += pom_pars;
            //gsInfo << "Parametrizace pro radu " << i << ": \n" << pom_pars << "\n";
        }
    gsMatrix<T> pom3(3, num_sample_pars1);
    gsMatrix<T> pom4(num_sample_pars1, 3);
    gsMatrix<T> pom_pars2(1, num_sample_pars1);
    gsMatrix<T> pom_pars_sum2(1, num_sample_pars1);
    for (index_t j = 0; j < num_sample_pars2; j++) {
            for (index_t i = 0; i < num_sample_pars1; i++) {
                pom3.col(i) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom3, "MD_bottom_surface_pompoints2");
            pom4 = pom3.transpose();
            pom_pars2 = centripetalParameterization(pom4);
            pom_pars_sum2 += pom_pars2;
            //gsInfo << "Parametrizace pro sloupec " << j << ": \n" << pom_pars2 << "\n";
        }
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface3(0, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum2(i)/num_sample_pars2;
            sampling_par_points_surface3(1, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum(j)/num_sample_pars1;
        }
    }

    gsTensorBSpline<2, T> surface_up_initial;
    //selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_up_parameter_domain, surface_outer, surface_up_initial, kvfit, kvfit, 1, print_info);
    selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_up_initial, kvfit, kvfit, 1, print_info);
    //controlNetFairing(surface_up_initial, 3, print_info);
    gsMatrix<T> surface_up_initial_controlnet = surface_up_initial.coefs();
    gsMatrix<index_t> boundaryind = surface_up_initial.basis().allBoundary();
    gsMatrix<T> boundarypoints(boundaryind.size(), 3);
    for (index_t i = 0; i < boundaryind.size(); i++) {
        boundarypoints.row(i) = surface_up_initial_controlnet.row(boundaryind(i));
    }
    gsMatrix<T> surface_up_points3 = surface_up_points.transpose();
    surface_up_initial = surfaceFittingLSQWithBoundary(surface_up_points3, sampling_par_points_surface3, kvfit, kvfit, boundarypoints, print_info);
    surface_up_initial_controlnet = surface_up_initial.coefs();

    T max_approx_error;
    T avg_approx_error;
    gsTensorBSpline<2, T> surface_up;
    surfaceApproximationIterativeWithBoundary(surface_up_points, kvfit, kvfit, surface_up_initial_controlnet, sampling_par_points_surface3, surface_up, max_approx_error, avg_approx_error, 5, print_info);
    m_surfaceTopFrontDomain = surface_up;

    //
    gsWriteParaview(surface_up_initial, "FD_up_surface_initial", 5000);
    gsMesh<> FD_up_surface_initial_mesh;
    surface_up_initial.controlNet(FD_up_surface_initial_mesh);
    gsWriteParaview( FD_up_surface_initial_mesh, "FD_up_surface_initial_mesh");
    //

    }

    // ******************************
    // MIDDLE DOMAIN - BOTTOM SURFACE
    // ******************************
    {
    int num_sample_pars = 50;
    gsKnotVector<T> kvfit = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(1);
    gsKnotVector<T> kvloft = (m_runnerBlade.getSuctionSideSurfaceAfterTrimming()).knots(0);
    gsTensorBSpline<2,T> surface_front = m_surfaceSuctionSide;
    gsTensorBSpline<2,T> surface_front_rotated = m_surfacePressureSideRotated;
    gsTensorBSpline<2,T> surface_frontmiddle = m_surfaceMiddleBackDivider;
    gsTensorBSpline<2,T> surface_frontmost = m_surfaceFrontMiddleDivider;
    gsTensorBSpline<2,T> surface_inner = m_surfaceInner;

    surface_front.insertKnot(0.1, 1);
    surface_front.insertKnot(0.3, 1);
    //surface_front.insertKnot(0.5, 1);
    //surface_front.insertKnot(0.7, 1);
    //surface_front.insertKnot(0.9, 1);
    surface_front_rotated.insertKnot(0.1, 1);
    surface_front_rotated.insertKnot(0.3, 1);
    //surface_front_rotated.insertKnot(0.5, 1);
    //surface_front_rotated.insertKnot(0.7, 1);
    //surface_front_rotated.insertKnot(0.9, 1);
    gsKnotVector<T> kvfit2 = surface_front.knots(1);

    gsMatrix<T> parameter_points(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front_rotated(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmiddle(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_frontmost(2, num_sample_pars);
    gsMatrix<T> parameter_points_surface_front2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_front_rotated2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmiddle2(num_sample_pars, 2);
    gsMatrix<T> parameter_points_surface_frontmost2(num_sample_pars, 2);
    gsMatrix<T> par_points_for_approx(1, num_sample_pars);

    gsMatrix<T> points_on_surface_front(3, num_sample_pars);
    gsMatrix<T> points_on_surface_front_rotated(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmiddle(3, num_sample_pars);
    gsMatrix<T> points_on_surface_frontmost(3, num_sample_pars);
    gsVector<T> point_on_surface(3);

    gsBSpline<T> curve_left_parameterdomain;
    gsBSpline<T> curve_right_parameterdomain;
    gsBSpline<T> curve_top_parameterdomain;
    gsBSpline<T> curve_bottom_parameterdomain;

    for (int i = 0; i < num_sample_pars; i++) {
        parameter_points(0,i) = 0.0;
        parameter_points(1,i) = (T) i/(num_sample_pars-1);
    }

    surface_front.eval_into(parameter_points, points_on_surface_front);
    surface_front_rotated.eval_into(parameter_points, points_on_surface_front_rotated);
    surface_frontmiddle.eval_into(parameter_points, points_on_surface_frontmiddle);
    surface_frontmost.eval_into(parameter_points, points_on_surface_frontmost);
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    int num_iter_front = 0;
    int max_num_iter_front = 0;
    int num_iter_front_rotated = 0;
    int max_num_iter_front_rotated = 0;
    int num_iter_frontmiddle = 0;
    int max_num_iter_frontmiddle = 0;
    int num_iter_frontmost = 0;
    int max_num_iter_frontmost = 0;

    for (index_t i = 0; i < num_sample_pars; i++) {
        point_on_surface = points_on_surface_front.col(i);
        init_sol << 0.5, 0.25;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front, print_info);
        parameter_points_surface_front(0, i) = sol(0);
        parameter_points_surface_front(1, i) = sol(1);
        if (max_num_iter_front < num_iter_front) {
            max_num_iter_front = num_iter_front;
        }
        point_on_surface = points_on_surface_front_rotated.col(i);
        init_sol << 0.5, 0.5;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_front_rotated, print_info);
        parameter_points_surface_front_rotated(0, i) = sol(0);
        parameter_points_surface_front_rotated(1, i) = sol(1);
        if (max_num_iter_front_rotated < num_iter_front_rotated) {
            max_num_iter_front_rotated = num_iter_front_rotated;
        }
        point_on_surface = points_on_surface_frontmiddle.col(i);
        init_sol << 0.9, 0.3+parameter_points(1,i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmiddle, print_info);
        parameter_points_surface_frontmiddle(0, i) = sol(0);
        parameter_points_surface_frontmiddle(1, i) = sol(1);
        if (max_num_iter_frontmiddle < num_iter_frontmiddle) {
            max_num_iter_frontmiddle = num_iter_frontmiddle;
        }
        point_on_surface = points_on_surface_frontmost.col(i);
        init_sol << 0.4, 0.3+parameter_points(1,i)/4;
        minimizePointSurfaceDistanceviaNR(point_on_surface, m_surfaceInner, init_sol, sol, static_cast<real_t>(1e-6), 100, num_iter_frontmost, print_info);
        parameter_points_surface_frontmost(0, i) = sol(0);
        parameter_points_surface_frontmost(1, i) = sol(1);
        if (max_num_iter_frontmost < num_iter_frontmost) {
            max_num_iter_frontmost = num_iter_frontmost;
        }
    }
    parameter_points_surface_frontmiddle.col(0) = parameter_points_surface_front.col(num_sample_pars-1);
    parameter_points_surface_frontmiddle.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(num_sample_pars-1);
    parameter_points_surface_frontmost.col(0) = parameter_points_surface_front.col(0);
    parameter_points_surface_frontmost.col(num_sample_pars-1) = parameter_points_surface_front_rotated.col(0);

    parameter_points_surface_front2 = parameter_points_surface_front.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front2);
    curve_top_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front2, par_points_for_approx, kvfit2, print_info);
    parameter_points_surface_front_rotated2 = parameter_points_surface_front_rotated.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_front_rotated2);
    curve_bottom_parameterdomain = curveFittingWithBoundary(parameter_points_surface_front_rotated2, par_points_for_approx, kvfit2, print_info);
    parameter_points_surface_frontmiddle2 = parameter_points_surface_frontmiddle.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmiddle2);
    curve_right_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmiddle2, par_points_for_approx, kvfit, print_info);
    parameter_points_surface_frontmost2 = parameter_points_surface_frontmost.transpose();
    par_points_for_approx = centripetalParameterization(parameter_points_surface_frontmost2);
    curve_left_parameterdomain = curveFittingWithBoundary(parameter_points_surface_frontmost2, par_points_for_approx, kvfit, print_info);

    gsMatrix<T> surface_bottom_parameter_domain_controlnet(curve_top_parameterdomain.coefs().rows() * curve_left_parameterdomain.coefs().rows(), 2);
    discreteCoonsPatch(curve_top_parameterdomain.coefs(), curve_bottom_parameterdomain.coefs(), curve_left_parameterdomain.coefs(), curve_right_parameterdomain.coefs(), surface_bottom_parameter_domain_controlnet, print_info);
    //gsInfo << surface_bottom_parameter_domain_controlnet << "\n";
    gsTensorBSpline<2, T> surface_bottom_parameter_domain;
    gsTensorBSplineBasis<2, T> basis(kvfit2, kvfit);
    surface_bottom_parameter_domain = gsTensorBSpline<2, T> (basis, surface_bottom_parameter_domain_controlnet);

    // Sampling point on bottom surface for approximation
    int num_sample_pars1 = 30;
    int num_sample_pars2 = 20;
    gsMatrix<T> sampling_par_points_surface(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_bottom_parameter_domain_points(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> surface_bottom_points(3, (num_sample_pars1-2)*(num_sample_pars2-2));
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface(0, (i-1)*(num_sample_pars2-2)+j-1) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface(1, (i-1)*(num_sample_pars2-2)+j-1) = (T) j/(num_sample_pars2-1);
        }
    }
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface, surface_bottom_parameter_domain_points);
    surface_inner.eval_into(surface_bottom_parameter_domain_points, surface_bottom_points);
    //gsInfo << "a\n";

    gsMatrix<T> sampling_par_points_surface2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_parameter_domain_points2(2, (num_sample_pars1)*(num_sample_pars2));
    gsMatrix<T> surface_bottom_points2(3, (num_sample_pars1)*(num_sample_pars2));
    for (index_t i = 0; i < num_sample_pars1; i++) {
        for (index_t j = 0; j < num_sample_pars2; j++) {
            sampling_par_points_surface2(0, i*(num_sample_pars2)+j) = (T) i/(num_sample_pars1-1);
            sampling_par_points_surface2(1, i*(num_sample_pars2)+j) = (T) j/(num_sample_pars2-1);
        }
    }
    surface_bottom_parameter_domain.eval_into(sampling_par_points_surface2, surface_bottom_parameter_domain_points2);
    surface_inner.eval_into(surface_bottom_parameter_domain_points2, surface_bottom_points2);
    //gsInfo << "a\n";

    // Finding suitable parameters for approximated points with the help of centripetal/chordal parameterization and averaging
    gsMatrix<T> sampling_par_points_surface3(2, (num_sample_pars1-2)*(num_sample_pars2-2));
    gsMatrix<T> pom(3, num_sample_pars2);
    gsMatrix<T> pom2(num_sample_pars2, 3);
    gsMatrix<T> pom_pars(1, num_sample_pars2);
    gsMatrix<T> pom_pars_sum(1, num_sample_pars2);
    pom_pars_sum.setZero();
    for (index_t i = 0; i < num_sample_pars1; i++) {
            for (index_t j = 0; j < num_sample_pars2; j++) {
                pom.col(j) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom, "MD_bottom_surface_pompoints");
            pom2 = pom.transpose();
            pom_pars = chordalParameterization(pom2);
            pom_pars_sum += pom_pars;
            //gsInfo << "Parametrizace pro radu " << i << ": \n" << pom_pars << "\n";
        }
    gsMatrix<T> pom3(3, num_sample_pars1);
    gsMatrix<T> pom4(num_sample_pars1, 3);
    gsMatrix<T> pom_pars2(1, num_sample_pars1);
    gsMatrix<T> pom_pars_sum2(1, num_sample_pars1);
    for (index_t j = 0; j < num_sample_pars2; j++) {
            for (index_t i = 0; i < num_sample_pars1; i++) {
                pom3.col(i) = surface_bottom_points2.col(i*(num_sample_pars2)+j);
            }
            //gsWriteParaviewPoints(pom3, "MD_bottom_surface_pompoints2");
            pom4 = pom3.transpose();
            pom_pars2 = centripetalParameterization(pom4);
            pom_pars_sum2 += pom_pars2;
            //gsInfo << "Parametrizace pro sloupec " << j << ": \n" << pom_pars2 << "\n";
        }
    for (index_t i = 1; i < num_sample_pars1-1; i++) {
        for (index_t j = 1; j < num_sample_pars2-1; j++) {
            sampling_par_points_surface3(0, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum2(i)/num_sample_pars2;
            sampling_par_points_surface3(1, (i-1)*(num_sample_pars2-2)+j-1) = pom_pars_sum(j)/num_sample_pars1;
        }
    }
    //gsInfo << sampling_par_points_surface << "\n";
    //gsInfo << sampling_par_points_surface3 << "\n";

    gsTensorBSpline<2, T> surface_bottom_initial;
    //selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_parameter_domain, surface_inner, surface_bottom_initial, kvfit, kvfit, 0, print_info);
    selectInitialControlNet(surface_frontmost, surface_frontmiddle, surface_front, surface_front_rotated, surface_bottom_initial, kvfit2, kvfit, 0, print_info);
    gsMatrix<T> surface_bottom_initial_controlnet = surface_bottom_initial.coefs();
    //controlNetFairing2(surface_bottom_initial, 2, print_info);
    gsMatrix<index_t> boundaryind = surface_bottom_initial.basis().allBoundary();
    gsMatrix<T> boundarypoints(boundaryind.size(), 3);
    for (index_t i = 0; i < boundaryind.size(); i++) {
        boundarypoints.row(i) = surface_bottom_initial_controlnet.row(boundaryind(i));
    }
    gsMatrix<T> surface_bottom_points3 = surface_bottom_points.transpose();
    surface_bottom_initial = surfaceFittingLSQWithBoundary(surface_bottom_points3, sampling_par_points_surface3, kvfit2, kvfit, boundarypoints, print_info);
    surface_bottom_initial_controlnet = surface_bottom_initial.coefs();

    T max_approx_error;
    T avg_approx_error;
    gsTensorBSpline<2, T> surface_bottom;
    surfaceApproximationIterativeWithBoundary(surface_bottom_points, kvfit2, kvfit, surface_bottom_initial_controlnet, sampling_par_points_surface3, surface_bottom, max_approx_error, avg_approx_error, 10, print_info);
    m_surfaceBottomMiddleDomain = surface_bottom;
    m_surfaceSuctionSide = surface_front;
    m_surfacePressureSideRotated = surface_front_rotated;

    //
    //gsWriteParaviewPoints(points_on_surface_front, "MD_points_on_surface_front");
    //gsWriteParaviewPoints(points_on_surface_front_rotated, "MD_points_on_surface_front_rotated");
    //gsWriteParaviewPoints(points_on_surface_frontmiddle, "MD_points_on_surface_frontmiddle");
    //gsWriteParaviewPoints(points_on_surface_frontmost, "MD_points_on_surface_frontmost");
    gsWriteParaviewPoints(parameter_points_surface_front, "MD_parameter_points_surface_front");
    gsWriteParaviewPoints(parameter_points_surface_front_rotated, "MD_parameter_points_surface_frontrotated");
    gsWriteParaviewPoints(parameter_points_surface_frontmiddle, "MD_parameter_points_surface_frontmiddle");
    gsWriteParaviewPoints(parameter_points_surface_frontmost, "MD_parameter_points_surface_frontmost");
    gsWriteParaview(surface_bottom_parameter_domain, "MD_bottom_surface_parameterdomain", 5000);

    gsWriteParaviewPoints(surface_bottom_points, "MD_bottom_sampling points");
    gsWriteParaviewPoints(surface_bottom_parameter_domain_points, "MD_bottom_sampling points_parameterdomain");

    gsWriteParaview(surface_bottom_initial, "MD_bottom_surface_initial");
    gsMesh<> MD_bottom_surface_initial_mesh;
    surface_bottom_initial.controlNet(MD_bottom_surface_initial_mesh);
    gsWriteParaview( MD_bottom_surface_initial_mesh, "MD_bottom_surface_initial_mesh");

    gsWriteParaview(curve_left_parameterdomain, "MD_curve_left_parameterdomain");
    gsWriteParaview(curve_right_parameterdomain, "MD_curve_right_parameterdomain");
    gsWriteParaview(curve_top_parameterdomain, "MD_curve_top_parameterdomain");
    gsWriteParaview(curve_bottom_parameterdomain, "MD_curve_bottom_parameterdomain");
    gsWriteParaview(surface_bottom_parameter_domain, "MD_surface_bottom_parameterdomain", 5000);
    //
    }

    if (plot) {
        gsWriteParaview(m_surfaceBottomFrontDomain, "FD_bottom_surface_final", 5000);
        gsMesh<> FD_bottom_surface_final_mesh;
        m_surfaceBottomFrontDomain.controlNet(FD_bottom_surface_final_mesh);
        gsWriteParaview( FD_bottom_surface_final_mesh, "FD_bottom_surface_final_mesh");

        gsWriteParaview(m_surfaceTopFrontDomain, "FD_up_surface_final", 5000);
        gsMesh<> FD_up_surface_final_mesh;
        m_surfaceTopFrontDomain.controlNet(FD_up_surface_final_mesh);
        gsWriteParaview( FD_up_surface_final_mesh, "FD_up_surface_final_mesh");

        gsWriteParaview(m_surfaceBottomMiddleDomain, "MD_bottom_surface_final", 5000);
        gsMesh<> MD_bottom_surface_final_mesh;
        m_surfaceBottomMiddleDomain.controlNet(MD_bottom_surface_final_mesh);
        gsWriteParaview( MD_bottom_surface_final_mesh, "MD_bottom_surface_final_mesh");
    }

    return 0;
}
*/

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::computeVolumes(bool plot, bool print_info) {

    //T alpha = 0.0;
    gsTensorBSpline<3, T> volume_front_domain;
    gsTensorBSpline<3, T> volume_middle_domain;
    gsTensorBSpline<3, T> volume_back_domain;

    // Front domain
    discreteCoonsPatch3D(m_surfaceFrontExtension, m_surfaceFrontExtensionRotated, m_surfaceFront, m_surfaceFrontMiddleDivider, m_surfaceBottomFrontDomain, m_surfaceTopFrontDomain, volume_front_domain, print_info);
    m_volumeFrontDomain = volume_front_domain;

    // Middle domain
    discreteCoonsPatch3D(m_surfaceSuctionSide, m_surfacePressureSideRotated, m_surfaceFrontMiddleDivider, m_surfaceMiddleBackDivider, m_surfaceBottomMiddleDomain, m_surfaceTopMiddleDomain, volume_middle_domain, print_info);
    m_volumeMiddleDomain = volume_middle_domain;

    // Back domain
    /*springModelPatch3D(m_surfaceBackExtension, m_surfaceBackExtensionRotated, m_surfaceMiddleBackDivider, m_surfaceBack, m_surfaceBottomBackDomain, m_surfaceTopBackDomain, volume_back_domain, print_info);
    gsMatrix<T> volume_back_domain_coefs_spring = volume_back_domain.coefs();
    discreteCoonsPatch3D(m_surfaceBackExtension, m_surfaceBackExtensionRotated, m_surfaceMiddleBackDivider, m_surfaceBack, m_surfaceBottomBackDomain, m_surfaceTopBackDomain, volume_back_domain, print_info);
    gsMatrix<T> volume_back_domain_coefs_coons = volume_back_domain.coefs();

    gsMatrix<T> volume_back_domain_coefs(volume_back_domain_coefs_coons.rows(), volume_back_domain_coefs_coons.cols());
    gsKnotVector<T> kv1 = volume_back_domain.knots(0);
    gsKnotVector<T> kv2 = volume_back_domain.knots(1);
    gsKnotVector<T> kv3 = volume_back_domain.knots(2);
    int num1 = kv1.size() - kv1.degree() - 1;
    int num2 = kv2.size() - kv2.degree() - 1;
    int num3 = kv3.size() - kv3.degree() - 1;
    for (index_t j = 0; j < num2; j++) {
        alpha = (T) j/(num2-1);
        for (index_t k = 0; k < num3; k++) {
            for (index_t i = 0; i < num1; i++) {
                volume_back_domain_coefs.row(k*num1*num2+j*num1+i) = (1-alpha) * volume_back_domain_coefs_coons.row(k*num1*num2+j*num1+i) + alpha * volume_back_domain_coefs_spring.row(k*num1*num2+j*num1+i);
            }
        }
    }
    gsTensorBSplineBasis<3, T> basis = volume_back_domain.basis();
    volume_back_domain = gsTensorBSpline<3, T> (basis, volume_back_domain_coefs);
    m_volumeBackDomain = volume_back_domain;
    */
    discreteCoonsPatch3D(m_surfaceBackExtension, m_surfaceBackExtensionRotated, m_surfaceMiddleBackDivider, m_surfaceBack, m_surfaceBottomBackDomain, m_surfaceTopBackDomain, volume_back_domain, print_info);
    m_volumeBackDomain = volume_back_domain;

    // Rotation of final volume such that the turbine axis concides with x-axis (the construction provides volumes for the turbine axis coincident with z-axis)
    gsVector<real_t> rotation_axis(3);
    rotation_axis << 0, 1, 0;
    //m_volumeFrontDomain.rotate(PI/2, rotation_axis);
    //m_volumeMiddleDomain.rotate(PI/2, rotation_axis);
    //m_volumeBackDomain.rotate(PI/2, rotation_axis);


//    /**************************************************************/
//    /* Optimization of patches */
//    gsTensorBSpline<3> domain = m_volumeBackDomain;
//    std::string output("domain_optimization");
//    int iterations = 5;
//    real_t orthogonality = 0;
//    real_t skewness = 0;
//    real_t eccentricity = 0;
//    real_t intersection = 0;
//    real_t uniformity = 1;
//    real_t area = 0;
//    real_t length = 1;
//    real_t epsilon = 1e-7;
//    int jacPts = 10000;
//    bool dumped = false;

//    gsInfo << "------------------------------------------------------------"
//                 "\nIterations: " << iterations << "\n"
//                 "Orthogonality: " << orthogonality << "\n"
//                 "Skewness: " << skewness << "\n"
//                 "Eccentricity: " << eccentricity << "\n"
//                 "Uniformity: " << uniformity << "\n"
//                 "Length: " << length << "\n"
//                 "Area: " << area << "\n"
//                 "Intersection: " << intersection << "\n"
//                 "Epsilon: " << epsilon << "\n"
//                 "Jacobian Points: " << jacPts << "\n"
//                 "Dumped: " << dumped << "\n"
//                 "------------------------------------------------------------"
//                 "\n\n";

//    gsFileData<> fileData;
//    fileData << domain;

//    std::string out = "domain.xml";
//    fileData.dump(out);

//    gsFileData<> data("domain.xml");
//    gsGeometry<>::uPtr geompatch1;
//    if (data.has< gsGeometry<> >())
//    {
//        geompatch1 = data.getFirst< gsGeometry<> >();
//    }

//    if (!geompatch1)
//    {
//        gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
//    }

//    gsQualityMeasure<real_t> optimization(*geompatch1);
//    gsInfo << "Value of functional: "
//              << optimization.functional(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon)
//              << "\n";
//    gsInfo << "Length: "
//           << optimization.functional(0, 0, 0, 0, 1, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Uniformity: "
//           << optimization.functional(0, 0, 0, 1, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Orthogonality: "
//           << optimization.functional(1, 0, 0, 0, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Skewness: "
//           << optimization.functional(0, 1, 0, 0, 0, 0, 0, epsilon)
//           << "\n";
//    //gsInfo << "Area: " // dim=2 only
//    //       << optimization.functional(orthogonality, skewness, eccentricity, uniformity, length, 1, intersection, epsilon)
//    //       << "\n";
//    gsInfo << "Eccentricity: "
//           << optimization.functional(0, 0, 1, 0, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Intersection: "
//           << optimization.functional(0, 0, 0, 0, 0, 0, 1, epsilon)
//           << "\n\n";

//    saveData(*geompatch1, output, 0);
//    checkJacobianDeterminant(*geompatch1, jacPts, true, output, 0);

//    for (int it = 0; it != iterations; it++)
//    {
//        gsInfo << "Iteration: " << it+1 << " / " << iterations << "\n";

//        optimization.optimize(orthogonality, skewness, eccentricity, uniformity,
//                              length, area,
//                              intersection, epsilon, dumped);

//        gsInfo << "Value of functional: "
//                  << optimization.functional(orthogonality, skewness,
//                                             eccentricity, uniformity, length,
//                                             area, intersection, epsilon)
//                  << "\n";

//        saveData(*geompatch1, output, it + 1);
//        checkJacobianDeterminant(*geompatch1, jacPts, true, output, it + 1);
//    }

//    gsInfo << "Length: "
//           << optimization.functional(0, 0, 0, 0, 1, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Uniformity: "
//           << optimization.functional(0, 0, 0, 1, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Orthogonality: "
//           << optimization.functional(1, 0, 0, 0, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Skewness: "
//           << optimization.functional(0, 1, 0, 0, 0, 0, 0, epsilon)
//           << "\n";
//    //gsInfo << "Area: " // dim=2 only
//    //       << optimization.functional(orthogonality, skewness, eccentricity, uniformity, length, 1, intersection, epsilon)
//    //       << "\n";
//    gsInfo << "Eccentricity: "
//           << optimization.functional(0, 0, 1, 0, 0, 0, 0, epsilon)
//           << "\n";
//    gsInfo << "Intersection: "
//           << optimization.functional(0, 0, 0, 0, 0, 0, 1, epsilon)
//           << "\n\n";
//    /***************************************************************/


    // For testing ONLY
    std::string test_surface_file;
    std::string test_surface_mesh_file;
    /*kv1 = m_volumeFrontDomain.knots(0);
    kv2 = m_volumeFrontDomain.knots(1);
    kv3 = m_volumeFrontDomain.knots(2);
    num1 = kv1.size() - kv1.degree() - 1;
    num2 = kv2.size() - kv2.degree() - 1;
    num3 = kv3.size() - kv3.degree() - 1;
    gsMatrix<T> control_net_test(num1*num2,3);
    control_net_test.setZero();
    for (index_t k = 0; k < num3; k++) {
        for (index_t j = 0; j < num2; j++) {
            for (index_t i = 0; i < num1; i++) {
                control_net_test.row(j*num1+i) = m_volumeFrontDomain.coefs().row(k*num1*num2+j*num1+i);
            }
        }
        test_surface_file = "test_surface_" + util::to_string(k);
        test_surface_mesh_file = "test_surface_mesh_" + util::to_string(k);
        gsTensorBSpline<2,T> test_surface (kv1, kv2, control_net_test);
        gsWriteParaview(test_surface, test_surface_file, 10000);
        gsMesh<> test_surface_mesh;
        test_surface.controlNet(test_surface_mesh);
        gsWriteParaview( test_surface_mesh, test_surface_mesh_file);
    }*/

    int dir = 0;
    int n = 10;
    T par;
    gsTensorBSpline<2,T> volume_slice;
    for (index_t i = 1; i < n+1; i++) {
        par = (T) i/(n+1);
        //gsInfo << par << "\n";
        m_volumeFrontDomain.slice(dir, par, volume_slice);
        test_surface_file = "test_surface_" + util::to_string(i);
        test_surface_mesh_file = "test_surface_mesh_" + util::to_string(i);
        gsWriteParaview(volume_slice, test_surface_file, 10000);
        gsMesh<> volume_slice_mesh;
        volume_slice.controlNet(volume_slice_mesh);
        gsWriteParaview( volume_slice_mesh, test_surface_mesh_file);
    }
    gsInfo << "Front domain:\n";
    gsInfo << m_volumeFrontDomain.knots(0) << "\n";
    gsInfo << m_volumeFrontDomain.knots(1) << "\n";
    gsInfo << m_volumeFrontDomain.knots(2) << "\n";
    gsInfo << volume_slice << "\n";

    //


    // Setting up the final multipatch domain
    gsMultiPatch<T> mp;

    mp.addPatch(m_volumeFrontDomain);
    mp.addPatch(m_volumeMiddleDomain);
    mp.addPatch(m_volumeBackDomain);

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::north, 2, boundary::south);
    mp.addAutoBoundaries();

    m_runnerDomain = mp;

    if (plot) {
        gsWriteParaview(m_volumeFrontDomain, "FD_volume", 10000);
        gsMesh<> FD_volume_mesh;
        m_volumeFrontDomain.controlNet(FD_volume_mesh);
        gsWriteParaview( FD_volume_mesh, "FD_volume_mesh");
        gsWriteParaview(m_volumeMiddleDomain, "MD_volume", 10000);
        gsMesh<> MD_volume_mesh;
        m_volumeMiddleDomain.controlNet(MD_volume_mesh);
        gsWriteParaview( MD_volume_mesh, "MD_volume_mesh");
        gsWriteParaview(m_volumeBackDomain, "BD_volume", 10000);
        gsMesh<> BD_volume_mesh;
        m_volumeBackDomain.controlNet(BD_volume_mesh);
        gsWriteParaview( BD_volume_mesh, "BD_volume_mesh");
    }

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::selectInitialControlNet(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom,
                                                               gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info) {

    int num_cp1 = surface_top.knots(1).size() - surface_top.degree(1) - 1;
    int num_cp2 = surface_left.knots(1).size() - surface_left.degree(1) - 1;
    int num_cp3 = surface_top.knots(0).size() - surface_top.degree(0) - 1;
    //gsInfo << num_cp1 << "\n";
    //gsInfo << num_cp2 << "\n";
    //gsInfo << num_cp3 << "\n";

    gsMatrix<T> cp1(num_cp1, 3);
    gsMatrix<T> cp2(num_cp2, 3);
    gsMatrix<T> cp3(num_cp1, 3);
    gsMatrix<T> cp4(num_cp2, 3);

    gsMatrix<T> surface_initial_controlnet(num_cp1 * num_cp2, 3);

    for (index_t i = 0; i < num_cp1; i++) {
        if (side == 0) {
            cp1.row(i) = surface_top.coefs().row(i * num_cp3);
            cp3.row(i) = surface_bottom.coefs().row(i * num_cp3);
        }
        else {
            cp1.row(i) = surface_top.coefs().row((i+1) * num_cp3 - 1);
            cp3.row(i) = surface_bottom.coefs().row((i+1) * num_cp3 - 1);
        }
    }

    for (index_t i = 0; i < num_cp2; i++) {
        if (side == 0) {
            cp2.row(i) = surface_left.coefs().row(i * num_cp3);
            cp4.row(i) = surface_right.coefs().row(i * num_cp3);
        }
        else {
            cp2.row(i) = surface_left.coefs().row((i+1) * num_cp3 - 1);
            cp4.row(i) = surface_right.coefs().row((i+1) * num_cp3 - 1);
        }
    }

    //gsInfo << "top:\n" << cp1 << "\n";
    //gsInfo << "bottom:\n" << cp3 << "\n";
    //gsInfo << "left:\n" << cp2 << "\n";
    //gsInfo << "right:\n" << cp4 << "\n";

    discreteCoonsPatch(cp1, cp3, cp2, cp4, surface_initial_controlnet, print_info);

    gsTensorBSplineBasis<2, T> basis(kv1, kv2);
    surface_initial = gsTensorBSpline<2, T> (basis, surface_initial_controlnet);

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::selectInitialControlNet2(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom,
                                                               gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info) {

    int num_cp1 = surface_top.knots(1).size() - surface_top.degree(1) - 1;
    int num_cp2 = surface_left.knots(1).size() - surface_left.degree(1) - 1;
    int num_cp3 = surface_top.knots(0).size() - surface_top.degree(0) - 1;
    //gsInfo << num_cp1 << "\n";
    //gsInfo << num_cp2 << "\n";
    //gsInfo << num_cp3 << "\n";

    gsMatrix<T> cp1(num_cp1, 3);
    gsMatrix<T> cp2(num_cp2, 3);
    gsMatrix<T> cp3(num_cp1, 3);
    gsMatrix<T> cp4(num_cp2, 3);

    gsMatrix<T> surface_initial_controlnet(num_cp1 * num_cp2, 3);

    for (index_t i = 0; i < num_cp1; i++) {
        if (side == 0) {
            cp1.row(i) = surface_top.coefs().row(i * num_cp3);
            cp3.row(i) = surface_bottom.coefs().row(i * num_cp3);
        }
        else {
            cp1.row(i) = surface_top.coefs().row((i+1) * num_cp3 - 1);
            cp3.row(i) = surface_bottom.coefs().row((i+1) * num_cp3 - 1);
        }
    }

    for (index_t i = 0; i < num_cp2; i++) {
        if (side == 0) {
            cp2.row(i) = surface_left.coefs().row(i * num_cp3);
            cp4.row(i) = surface_right.coefs().row(i * num_cp3);
        }
        else {
            cp2.row(i) = surface_left.coefs().row((i+1) * num_cp3 - 1);
            cp4.row(i) = surface_right.coefs().row((i+1) * num_cp3 - 1);
        }
    }

    //gsInfo << "top:\n" << cp1 << "\n";
    //gsInfo << "bottom:\n" << cp3 << "\n";
    //gsInfo << "left:\n" << cp2 << "\n";
    //gsInfo << "right:\n" << cp4 << "\n";

    springModelPatch(cp1, cp3, cp2, cp4, surface_initial_controlnet, print_info);

    gsTensorBSplineBasis<2, T> basis(kv1, kv2);
    surface_initial = gsTensorBSpline<2, T> (basis, surface_initial_controlnet);

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::selectInitialControlNet(gsTensorBSpline<2, T> surface_left, gsTensorBSpline<2, T> surface_right, gsTensorBSpline<2, T> surface_top, gsTensorBSpline<2, T> surface_bottom, gsTensorBSpline<2, T> surface_parameter_domain,
                                                               gsTensorBSpline<2, T> surface, gsTensorBSpline<2, T> & surface_initial, gsKnotVector<T> kv1, gsKnotVector<T> kv2, bool side, bool print_info) {

    int num_cp1 = surface_top.knots(1).size() - surface_top.degree(1) - 1;
    int num_cp2 = surface_left.knots(1).size() - surface_left.degree(1) - 1;
    int num_cp3 = surface_top.knots(0).size() - surface_top.degree(0) - 1;
    //gsInfo << num_cp1 << "\n";
    //gsInfo << num_cp2 << "\n";
    //gsInfo << num_cp3 << "\n";

    //int num_cp1 = curve_top_parameterdomain.coefs().rows();
    //int num_cp2 = curve_left_parameterdomain.coefs().rows();
    //int num_cp3 = kvloft.size() - kvloft.degree() - 1;
    gsMatrix<T> divided_differences1_1(1,num_cp1);
    gsMatrix<T> divided_differences1_2(1,num_cp1);
    gsMatrix<T> divided_differences2_1(1,num_cp2);
    gsMatrix<T> divided_differences2_2(1,num_cp2);
    gsMatrix<T> cp1(num_cp1, 3);
    gsMatrix<T> cp2(num_cp2, 3);
    gsMatrix<T> cn_pardomain(2, num_cp1*num_cp2);
    gsMatrix<T> cn_pardomain2(2, num_cp1*num_cp2);
    gsMatrix<T> cn(3, num_cp1*num_cp2);
    //gsInfo << "a\n";
    for (index_t i = 0; i < num_cp1; i++) {
        if (side == 0) {
            cp1.row(i) = surface_top.coefs().row(i * num_cp3);
        }
        else {
            cp1.row(i) = surface_top.coefs().row((i+1) * num_cp3 - 1);
        }

    }
    divided_differences1_1 = centripetalParameterization(cp1);
    for (index_t i = 0; i < num_cp1; i++) {
        if (side == 0) {
            cp1.row(i) = surface_bottom.coefs().row(i * num_cp3);
        }
        else {
            cp1.row(i) = surface_bottom.coefs().row((i+1) * num_cp3 - 1);
        }
    }
    divided_differences1_2 = centripetalParameterization(cp1);
    //gsInfo << divided_differences1_1 << "\n";
    for (index_t i = 0; i < num_cp2; i++) {
        if (side == 0) {
            cp2.row(i) = surface_left.coefs().row(i * num_cp3);
        }
        else {
            cp2.row(i) = surface_left.coefs().row((i+1) * num_cp3 - 1);
        }
    }
    divided_differences2_1 = chordalParameterization(cp2);
    for (index_t i = 0; i < num_cp2; i++) {
        if (side == 0) {
            cp2.row(i) = surface_right.coefs().row(i * num_cp3);
        }
        else {
            cp2.row(i) = surface_right.coefs().row((i+1) * num_cp3 - 1);
        }
    }
    divided_differences2_2 = chordalParameterization(cp2);
    //gsInfo << divided_differences2_1 << "\n";
    for (index_t i = 0; i < num_cp2; i++) {
        for (index_t j = 0; j < num_cp1; j++) {
            //cn_pardomain(0, i * num_cp1 + j) = divided_differences1(j);
            //cn_pardomain(1, i * num_cp1 + j) = divided_differences2(i);
            T alfa = (T) i/(num_cp2-1);
            T beta = (T) j/(num_cp1-1);
            cn_pardomain(0, i * num_cp1 + j) = (1-alfa) * divided_differences1_1(j) + alfa * divided_differences1_2(j);
            cn_pardomain(1, i * num_cp1 + j) = (1-beta) * divided_differences2_1(i) + beta * divided_differences2_2(i);
        }
    }
    //gsInfo << "Parametricke hodnoty pro pocatecni ridici sit:\n" << cn_pardomain << "\n";
    surface_parameter_domain.eval_into(cn_pardomain, cn_pardomain2);
    surface.eval_into(cn_pardomain2, cn);
    //gsInfo << "a\n";
    gsMatrix<T> surface_initial_controlnet(num_cp1 * num_cp2, 3);
    for (index_t i = 0; i < num_cp1; i++) {
        if (side == 0) {
            surface_initial_controlnet.row(i) = surface_top.coefs().row(i * num_cp3); // include control points from FrontExtensionSurface
            surface_initial_controlnet.row((num_cp2-1) * num_cp1 + i) = surface_bottom.coefs().row(i * num_cp3); // include control points from FrontExtensionRotatedSurface
        }
        else {
            surface_initial_controlnet.row(i) = surface_top.coefs().row((i+1) * num_cp3 - 1); // include control points from FrontExtensionSurface
            surface_initial_controlnet.row((num_cp2-1) * num_cp1 + i) = surface_bottom.coefs().row((i+1) * num_cp3 - 1); // include control points from FrontExtensionRotatedSurface
        }
    }
    //gsInfo << "a\n";
    for (index_t i = 1; i < num_cp2-1; i++) {
        for (index_t j = 0; j < num_cp1; j++) {
            if (j == 0) {
                if (side == 0) {
                    surface_initial_controlnet.row(i * num_cp1 + j) = surface_left.coefs().row(i * num_cp3);
                }
                else {
                    surface_initial_controlnet.row(i * num_cp1 + j) = surface_left.coefs().row((i+1) * num_cp3 - 1);
                }
            }
            else if (j == num_cp1-1) {
                if (side == 0) {
                    surface_initial_controlnet.row(i * num_cp1 + j) = surface_right.coefs().row(i * num_cp3);
                }
                else {
                    surface_initial_controlnet.row(i * num_cp1 + j) = surface_right.coefs().row((i+1) * num_cp3 - 1);
                }
            }
            else {
                surface_initial_controlnet(i * num_cp1 + j, 0) = cn(0, i * num_cp1 + j);
                surface_initial_controlnet(i * num_cp1 + j, 1) = cn(1, i * num_cp1 + j);
                surface_initial_controlnet(i * num_cp1 + j, 2) = cn(2, i * num_cp1 + j);
            }
        }
    }
    //gsInfo << "Pocatecni ridici sit:\n" << surface_initial_controlnet << "\n";
    //gsTensorBSpline<2, T> surface_initial;

    //Control net fairing
    /*for (index_t i = 1; i < num_cp2-1; i++) {
        for (index_t j = 1; j < num_cp1-1; j++) {
            surface_initial_controlnet(i * num_cp1 + j, 0) = (surface_initial_controlnet(i * num_cp1 + j - 1, 0) + surface_initial_controlnet(i * num_cp1 + j + 1, 0))/2;
            surface_initial_controlnet(i * num_cp1 + j, 1) = (surface_initial_controlnet(i * num_cp1 + j - 1, 1) + surface_initial_controlnet(i * num_cp1 + j + 1, 1))/2;
            surface_initial_controlnet(i * num_cp1 + j, 2) = (surface_initial_controlnet(i * num_cp1 + j - 1, 2) + surface_initial_controlnet(i * num_cp1 + j + 1, 2))/2;
        }
    }
    for (index_t i = 1; i < num_cp2-1; i++) {
        for (index_t j = 1; j < num_cp1-1; j++) {
            surface_initial_controlnet(i * num_cp1 + j, 0) = (surface_initial_controlnet(i * num_cp1 + j - 1, 0) + surface_initial_controlnet(i * num_cp1 + j + 1, 0))/2;
            surface_initial_controlnet(i * num_cp1 + j, 1) = (surface_initial_controlnet(i * num_cp1 + j - 1, 1) + surface_initial_controlnet(i * num_cp1 + j + 1, 1))/2;
            surface_initial_controlnet(i * num_cp1 + j, 2) = (surface_initial_controlnet(i * num_cp1 + j - 1, 2) + surface_initial_controlnet(i * num_cp1 + j + 1, 2))/2;
        }
    }
    */

    gsTensorBSplineBasis<2, T> basis(kv1, kv2);
    surface_initial = gsTensorBSpline<2, T> (basis, surface_initial_controlnet);

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::controlNetFairing(gsTensorBSpline<2, T> & surface_initial, int num, bool print_info) {

    gsMatrix<T> surface_initial_controlnet = surface_initial.coefs();
    int num_cp1 = surface_initial.knots(0).size() - surface_initial.degree(0) - 1;
    int num_cp2 = surface_initial.knots(1).size() - surface_initial.degree(1) - 1;

    for (index_t k = 0; k < num; k++) {
        for (index_t i = 1; i < num_cp2-1; i++) {
            for (index_t j = 1; j < num_cp1-1; j++) {
                surface_initial_controlnet(i * num_cp1 + j, 0) = (surface_initial_controlnet(i * num_cp1 + j - 1, 0) + surface_initial_controlnet(i * num_cp1 + j + 1, 0))/2;
                surface_initial_controlnet(i * num_cp1 + j, 1) = (surface_initial_controlnet(i * num_cp1 + j - 1, 1) + surface_initial_controlnet(i * num_cp1 + j + 1, 1))/2;
                surface_initial_controlnet(i * num_cp1 + j, 2) = (surface_initial_controlnet(i * num_cp1 + j - 1, 2) + surface_initial_controlnet(i * num_cp1 + j + 1, 2))/2;
            }
        }
    }

    gsTensorBSplineBasis<2, T> basis(surface_initial.knots(0), surface_initial.knots(1));
    surface_initial = gsTensorBSpline<2, T> (basis, surface_initial_controlnet);

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::controlNetFairing2(gsTensorBSpline<2, T> & surface_initial, int num, bool print_info) {

    gsMatrix<T> surface_initial_controlnet = surface_initial.coefs();
    int num_cp1 = surface_initial.knots(0).size() - surface_initial.degree(0) - 1;
    int num_cp2 = surface_initial.knots(1).size() - surface_initial.degree(1) - 1;

    for (index_t k = 0; k < num; k++) {
        for (index_t j = 1; j < num_cp1-1; j++) {
            for (index_t i = 1; i < num_cp2-1; i++) {
                surface_initial_controlnet(i * num_cp1 + j, 0) = (surface_initial_controlnet((i-1) * num_cp1 + j, 0) + surface_initial_controlnet((i+1) * num_cp1 + j, 0))/2;
                surface_initial_controlnet(i * num_cp1 + j, 1) = (surface_initial_controlnet((i-1) * num_cp1 + j, 1) + surface_initial_controlnet((i+1) * num_cp1 + j, 1))/2;
                surface_initial_controlnet(i * num_cp1 + j, 2) = (surface_initial_controlnet((i-1) * num_cp1 + j, 2) + surface_initial_controlnet((i+1) * num_cp1 + j, 2))/2;
            }
        }
    }

    gsTensorBSplineBasis<2, T> basis(surface_initial.knots(0), surface_initial.knots(1));
    surface_initial = gsTensorBSpline<2, T> (basis, surface_initial_controlnet);

    return 0;
}

/*
Surface approximation - method based on the paper:
Y. Kineri, M. Wang, H. Lin, T. Maekawa: B--spline surface fitting by iterative geometric interpolation/approximation algorithms. Computer-Aided Design, Vol. 44, No. 7, pp. 697-708
https://doi.org/10.1016/j.cad.2012.02.011
*/
template<class T>
int KaplanTurbineRunnerWheelDomain<T>::surfaceApproximationIterativeWithBoundary(gsMatrix<T> points, gsKnotVector<T> kv1, gsKnotVector<T> kv2, gsMatrix<T> controlnet_initial, gsMatrix<T> sampling_par_points_surface, gsTensorBSpline<2, T> & surface, T & max_approx_error, T & avg_approx_error, index_t max_iter, bool print_info) {

    if (print_info) {
        gsInfo << "\nSurface approximation with iterative method:\n";
        gsInfo << "============================================\n";
    }
    //gsTensorBSplineBasis<2, T> * surfaceBasis = new gsTensorBSplineBasis<2, T> (kv1, kv2);
    gsTensorBSplineBasis<2, T> basis(kv1, kv2);
    //unsigned int num_basis=surfaceBasis->size();
    unsigned int num_basis=basis.size();
    gsMatrix<T> values;
    gsMatrix<index_t> actives;

    gsMatrix<T> controlnet(num_basis, 3);
    //gsTensorBSpline<2, T> surface_bottom;
    controlnet = controlnet_initial;
    surface = gsTensorBSpline<2, T> (basis, controlnet);

    gsMatrix<index_t> boundaryindices = basis.allBoundary();

    if (print_info) {
        gsInfo << "Initialization\n";
    }
    gsMatrix<T> pars_for_approximation(2, points.cols());
    gsMatrix<T> closestPointsOnSurface(3, points.cols());
    pars_for_approximation = findClosestPointsOnSurface(points, surface, sampling_par_points_surface, static_cast<real_t>(1e-8), 100, print_info);
    surface.eval_into(pars_for_approximation, closestPointsOnSurface);
    gsWriteParaviewPoints(closestPointsOnSurface, "MD_bottom_closestpoints");

    bool ok;
    T max_approx_error2;
    T avg_approx_error2;
    maxErrorForSurfaceApproximation(points, surface, pars_for_approximation, max_approx_error2, avg_approx_error2, print_info);
    if (print_info) {
        gsInfo << "Initial errors: max. error = " << max_approx_error2 << ", " << "avg. error = " << avg_approx_error2 << "\n";
    }
    //gsInfo << "Max. approximation error: " << max_approx_error << "\n";
    //gsInfo << "Avg. approximation error: " << avg_approx_error << "\n";

    gsMatrix<T> NUM(num_basis, 3);
    gsMatrix<T> DEN(num_basis, 1);
    NUM.setZero();
    DEN.setZero();
    if (print_info) {
        gsInfo << "Approximation\n";
    }
    for (index_t i = 0; i < max_iter; i++) {
        if ((max_approx_error2 < max_approx_error) || (avg_approx_error2 < avg_approx_error)) {
            if (print_info) {
                gsInfo << "Solution converged to prescribed accuracy (max. error " << max_approx_error << ", avg. error " << avg_approx_error << "):\n max. error = " << max_approx_error2 <<
                          ", avg. error = " << avg_approx_error2 << "\n";
            }
            break;
        }
        if (print_info) {
            gsInfo << "--------------\n";
            gsInfo << "Iteration: " << i+1 << "\n";
        }
        basis.eval_into(pars_for_approximation, values);
        basis.active_into(pars_for_approximation, actives);

        for (index_t j = 0; j < pars_for_approximation.cols(); j++) {
            for (index_t k = 0; k < actives.rows(); k++) {
                NUM(actives(k,j), 0) += values(k,j) * (points(0, j) - closestPointsOnSurface(0,j));
                NUM(actives(k,j), 1) += values(k,j) * (points(1, j) - closestPointsOnSurface(1,j));
                NUM(actives(k,j), 2) += values(k,j) * (points(2, j) - closestPointsOnSurface(2,j));
                DEN(actives(k,j)) += values(k,j);
            }
        }
        for (unsigned j = 0; j < num_basis; j++) {
            ok = true;
            for (index_t l = 0; l < boundaryindices.rows(); l++) {
                if (j == boundaryindices(l)) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                controlnet(j, 0) += NUM(j, 0)/DEN(j);
                controlnet(j, 1) += NUM(j, 1)/DEN(j);
                controlnet(j, 2) += NUM(j, 2)/DEN(j);
            }

        }
        //if (i % 3 == 0) {
        //    surfaceApproximationIterativeFairing(controlnet, dim1, dim2, 2);
        //}
        surface = gsTensorBSpline<2, T> (basis, controlnet);

        maxErrorForSurfaceApproximation(points, surface, pars_for_approximation, max_approx_error2, avg_approx_error2, print_info);
        if (print_info) {
            gsInfo << "Errors: max. error = " << max_approx_error2 << ", " << "avg. error = " << avg_approx_error2 << "\n";
        }
        //gsInfo << "Max. approximation error: " << max_approx_error << "\n";
        //gsInfo << "Avg. approximation error: " << avg_approx_error << "\n";

        if (i < 20) {
            pars_for_approximation = findClosestPointsOnSurface(points, surface, pars_for_approximation, static_cast<real_t>(1e-8), 100, print_info);
            surface.eval_into(pars_for_approximation, closestPointsOnSurface);
        }

        //maxErrorForSurfaceApproximation(surface_bottom_points, surface_bottom, parameter_points_for_approximation, max_approx_error, avg_approx_error, print_info);
        //gsInfo << "Max. approximatio error: " << max_approx_error << "\n";
        //gsInfo << "Avg. approximatio error: " << avg_approx_error << "\n";

    }

    //delete surfaceBasis;

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::surfaceApproximationIterativeFairing(gsMatrix<T> & control_net, int dim1, int dim2, T gamma) {

    T e1 = exp(-gamma);
    T e2 = exp(-2 * gamma);

    for (index_t i = 1; i < dim2-1; i++) {
        for (index_t j = 1; j < dim1-1; j++) {
            control_net.row(i * dim1 + j) = (control_net.row(i * dim1 + j) + e1 * (control_net.row((i-1) * dim1 + j) + control_net.row(i * dim1 + j - 1) + control_net.row((i+1) * dim1 + j) + control_net.row(i * dim1 + j + 1))
                                             + e2 * (control_net.row((i-1) * dim1 + j - 1) + control_net.row((i+1) * dim1 + j - 1) + control_net.row((i-1) * dim1 + j + 1) + control_net.row((i+1) * dim1 + j + 1)))/(1 + 4 * e1 + 4 * e2);
        }
    }

    return 0;
}

template<class T>
int KaplanTurbineRunnerWheelDomain<T>::maxErrorForSurfaceApproximation(gsMatrix<T> points, gsTensorBSpline<2, T> surface, gsMatrix<T> pars, T & max_error, T & avg_error, bool print_info) {

    gsMatrix<T> pointsonsurface(3, points.cols());
    gsVector<T> vector(3);
    T norm;

    max_error = -1;
    avg_error = 0;
    surface.eval_into(pars, pointsonsurface);
    for (index_t i = 0; i < points.cols(); i++) {
        vector(0) = points(0, i) - pointsonsurface(0, i);
        vector(1) = points(1, i) - pointsonsurface(1, i);
        vector(2) = points(2, i) - pointsonsurface(2, i);
        norm = euclideanNorm(vector);
        if (max_error < norm) {
            max_error = norm;
        }
        avg_error += norm;
    }

    avg_error = avg_error/points.cols();

    return 0;

}


template<class T>
gsMatrix<T> KaplanTurbineRunnerWheelDomain<T>::findClosestPointsOnSurface(gsMatrix<T> points, gsTensorBSpline<2, T> surface, gsMatrix<T> initial_pars, T error, int max_iter, bool print_info) {

    int num_iter = 0;
    int max_num_iter = 0;
    gsVector<T> init_sol(2);
    gsVector<T> sol(2);
    gsVector<T> point;
    gsMatrix<T> parameter_points_final(2, points.cols());

    for (index_t i = 0; i < points.cols(); i++) {
        point = points.col(i);
        init_sol(0) = initial_pars(0, i);
        init_sol(1) = initial_pars(1, i);
        minimizePointSurfaceDistanceviaNR(point, surface, init_sol, sol, static_cast<real_t>(error), max_iter, num_iter, print_info);
        parameter_points_final(0, i) = sol(0);
        parameter_points_final(1, i) = sol(1);
        if (max_num_iter < num_iter) {
            max_num_iter = num_iter;
        }
        if (num_iter > max_iter/2) {
            //gsInfo << "Point with index " << i << " out of " << points.cols() << " does not converged. Solution: " << sol << "\n";
        }
    }
    //gsInfo << "Max number of NR iterations: " << max_num_iter << "\n";
    if (print_info) {
        gsInfo << "Max number of NR iterations: " << max_num_iter << "\n";
    }

    return parameter_points_final;

}



#endif // UWBKAPLANTURBINERUNNERWHEELDOMAIN_H
