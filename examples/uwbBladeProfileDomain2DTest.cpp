
#include <iostream>
#include <gismo.h>
#include <math.h>

//#include "../jku/gsMotorUtils.h"
//#include "uwbTurbineUtils.h"
//#include "uwbBladeProfile.h"
//#include "uwbDraftTube.h"
//#include "uwbProfileOptimization.h"
//#include "uwbHydraulicProfile.h"
//#include "uwbKaplanTurbineRunnerBlade.h"
//#include "uwbKaplanTurbineRunnerWheelDomain.h"
//#include "uwbKaplanTurbineGuideVane.h"
#include "uwbRANSExamplesSetting.h"
#include "uwbGeometryCreators.h"
#include "uwbDataFileReader.h"

//#include <gsIO/gsIOUtils.h>
//#include <gsOptimizer/gsQualityMeasure.h>
//#include "../jku/gsMotorUtils.h"
//#include "gsModeling/gsCoonsPatch.h"

using namespace gismo;

/*
gsMatrix<real_t> computeCircleArc2D(gsMatrix<real_t> start, real_t centreX, real_t centreY, real_t alpha)
{
    //constant for circle section

   real_t radius = math::sqrt(math::pow(start(0)-centreX,2)+math::pow(start(1)-centreY,2));
   real_t norm_koef =   radius * (4.0/3.0)*(math::tan(alpha/4));

   gsMatrix<real_t> end(2,1);
   end << centreX + (-centreX + start(0))*( math::cos(alpha)) - (-centreY + start(1)) * (math:: sin(alpha)),
          centreY + (-centreY + start(1))*(math::cos(alpha)) + (-centreX + start(0)) *(math::sin(alpha));

   gsMatrix<real_t> coefsclp(4, 2);
   coefsclp << start(0), start(1),
               start(0) - (norm_koef*(start(1)-centreY))/(math::sqrt(math::pow(start(0)-centreX,2)+math::pow(start(1)-centreY,2))), start(1) + (norm_koef*(start(0)-centreX))/(math::sqrt(math::pow(-start(0)+centreX,2)+math::pow(start(1)-centreY,2))),
               end(0) + (norm_koef*(end(1)-centreY))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))), end(1) - (norm_koef*(end(0)-centreX))/(math::sqrt(math::pow(-end(1)+centreY,2)+math::pow(end(0)-centreX,2))),
               end(0), end(1);

   return coefsclp;
}

gsMatrix<real_t> LSQFergusonSmooth(gsMatrix<real_t> P0, gsMatrix<real_t> T0, gsMatrix<real_t> T1, gsMatrix<real_t> P3) {

    real_t alfa = (-3*P0(0,0)*(T0(0,0) - 4*T0(0,1)*T1(0,0)*T1(0,1) + 4*T0(0,0)*pow(T1(0,1),2)) + 3*P3(0,0)*(T0(0,0) - 4*T0(0,1)*T1(0,0)*T1(0,1) + 4*T0(0,0)*pow(T1(0,1),2)) + 3*(P0(0,1) - P3(0,1))*(4*T0(0,0)*T1(0,0)*T1(0,1) + T0(0,1)*(-5 + 4*pow(T1(0,1),2))))/
            (9 + 16*(-2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2))));
    real_t beta = (3*(-(P0(0,0)*(T1(0,0) + 4*pow(T0(0,1),2)*T1(0,0) - 4*T0(0,0)*T0(0,1)*T1(0,1))) + P3(0,0)*(T1(0,0) + 4*pow(T0(0,1),2)*T1(0,0) - 4*T0(0,0)*T0(0,1)*T1(0,1)) + (P0(0,1) - P3(0,1))*(-5*T1(0,1) + 4*T0(0,1)*(T0(0,0)*T1(0,0) + T0(0,1)*T1(0,1)))))/
            (-9 + 32*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) - 16*pow(T1(0,1),2) + 16*pow(T0(0,1),2)*(-1 + 2*pow(T1(0,1),2)));
    gsInfo << alfa << "\n";
    gsInfo << beta << "\n";

    gsMatrix<real_t> res(4,2);
    res << P0(0,0), P0(0,1),
           P0(0,0) + alfa * T0(0,0), P0(0,1) + alfa * T0(0,1),
           P3(0,0) + beta * T1(0,0), P3(0,1) + beta * T1(0,1),
           P3(0,0), P3(0,1);

    return res;
}

gsMatrix<real_t> LSQFergusonShort(gsMatrix<real_t> P0, gsMatrix<real_t> T0, gsMatrix<real_t> T1, gsMatrix<real_t> P3) {

    real_t alfa = (-(P0(0,0)*(T0(0,0) - T0(0,1)*T1(0,0)*T1(0,1) + T0(0,0)*pow(T1(0,1),2))) + P3(0,0)*(T0(0,0) - T0(0,1)*T1(0,0)*T1(0,1) + T0(0,0)*pow(T1(0,1),2)) + (P0(0,1) - P3(0,1))*(T0(0,0)*T1(0,0)*T1(0,1) + T0(0,1)*(-2 + pow(T1(0,1),2))))/
            (3 - 2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2)));
    real_t beta = -((-(P0(0,0)*(1 + pow(T0(0,1),2))*T1(0,0)) + P0(0,0)*T0(0,0)*T0(0,1)*T1(0,1) + P3(0,0)*(T1(0,0) + pow(T0(0,1),2)*T1(0,0) - T0(0,0)*T0(0,1)*T1(0,1)) + (P0(0,1) - P3(0,1))*(T0(0,0)*T0(0,1)*T1(0,0) + (-2 + pow(T0(0,1),2))*T1(0,1)))/
                    (3 - 2*T0(0,0)*T0(0,1)*T1(0,0)*T1(0,1) + pow(T1(0,1),2) + pow(T0(0,1),2)*(1 - 2*pow(T1(0,1),2))));
    gsInfo << alfa << "\n";
    gsInfo << beta << "\n";
    alfa = abs(alfa);
    beta = abs(beta);

    gsMatrix<real_t> res(4,2);
    res << P0(0,0), P0(0,1),
           P0(0,0) + alfa * T0(0,0), P0(0,1) + alfa * T0(0,1),
           P3(0,0) + beta * T1(0,0), P3(0,1) + beta * T1(0,1),
           P3(0,0), P3(0,1);

    return res;
}

gsMatrix<real_t> AxisDirectionNormed(gsMatrix<real_t> T0, gsMatrix<real_t> T1) {

    gsMatrix<real_t> aux = (T0/T0.norm() + T1/T1.norm());

    return aux/aux.norm();
}

// Change the knot vector of curve1 to be the same as the knot vector of curve2
void UnifyKnotVectors(gsBSpline<real_t> &curve1, gsBSpline<real_t> curve2) {

    std::vector<real_t> knotsdiff;

    gsKnotVector<real_t> kv1 = curve1.knots();
    gsKnotVector<real_t> kv2 = curve2.knots();
    kv2.difference(kv1, knotsdiff);
    for (unsigned i = 0; i < knotsdiff.size(); i++) {
        curve1.insertKnot(knotsdiff[i]);
    }

}

gsMultiPatch<real_t> DomainBetweenBladeProfiles5(real_t const & index, real_t const & length_x1, real_t const & length_x2, real_t const & pitch, real_t const & camberX, real_t const & camberY, real_t const & leadingAngle, real_t const & trailingAngle, real_t const & thicknessX,
                                         real_t const & thicknessY, real_t const & endingOffset, real_t const & outputAngle, real_t const & radius, real_t const & chordLength, real_t const & Angle, real_t const & rotationCenterX, real_t const & rotationCenterY,
                                         real_t const & uniformity_param) {

    //----------------set parameters for blade profile----------------
    real_t offset_distance = 0.03;
    real_t wake_width = 3.0;
    //gsKnotVector<real_t> kvfit(0, 1, 4, 4);
    std::vector<real_t> kvfit_knots = {0.0, 0.0, 0.0, 0.0, 0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617, 1.0, 1.0, 1.0, 1.0};
    gsKnotVector<real_t> kvfit = gsKnotVector<real_t> (kvfit_knots);
    gsKnotVector<real_t> kvcub(0, 1, 0, 4);
    gsKnotVector<real_t> kvlin(0, 1, 0, 2);

    bool plot = true;
    bool plotMeshes = true;
    int num_samples = 100;
    gsVector<real_t> vec(2);
    //gsInfo << pitch << "\n ";
    vec(0) = rotationCenterX;
    vec(1) = rotationCenterY;
    gsMatrix<real_t> mat(2, 2);
    mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
           chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
    gsBSpline<real_t> suction_side_curve;
    gsBSpline<real_t> pressure_side_curve;
    gsBSpline<real_t> suction_side_offset_curve;
    gsBSpline<real_t> pressure_side_offset_curve;
    BladeProfile<real_t> * pBladeProfile = 0;

    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<real_t>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

    //---------------transform given profile----------------------
    suction_side_curve.translate(-vec);
    pressure_side_curve.translate(-vec);
    pressure_side_curve.linearTransform(mat);
    suction_side_curve.linearTransform(mat);
    pBladeProfile->setPressureSide(pressure_side_curve);
    pBladeProfile->setSuctionSide(suction_side_curve);

    pBladeProfile->computeOffset(offset_distance, suction_side_offset_curve, pressure_side_offset_curve, kvfit, num_samples);

    vec(0) = 0.0;
    vec(1) = pitch;
    suction_side_curve.translate(vec);
    suction_side_offset_curve.translate(vec);
    pBladeProfile->setSuctionSide(suction_side_curve);
    pBladeProfile->setSuctionSideOffset(suction_side_offset_curve);

    gsMatrix<real_t> suction_side_cp = suction_side_curve.coefs();
    gsMatrix<real_t> pressure_side_cp = pressure_side_curve.coefs();
    gsMatrix<real_t> suction_side_offset_cp = suction_side_offset_curve.coefs();
    gsMatrix<real_t> pressure_side_offset_cp = pressure_side_offset_curve.coefs();

    // --------------- cross-section curves ---------------------------------------------------------------------------------------------------
    gsMatrix<real_t> aux_cp(2, 2);
    aux_cp << suction_side_offset_cp.row(0),
              suction_side_cp.row(0);
    gsBSpline<real_t> cs_curve1 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve1.degreeElevate(2);
    gsMatrix<real_t> cs_curve1_cp = cs_curve1.coefs();
    aux_cp << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1),
              suction_side_cp.row(suction_side_curve.coefsSize()-1);
    gsBSpline<real_t> cs_curve2 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve2.degreeElevate(2);
    gsMatrix<real_t> cs_curve2_cp = cs_curve2.coefs();
    aux_cp << pressure_side_offset_cp.row(0),
              pressure_side_cp.row(0);
    gsBSpline<real_t> cs_curve3 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve3.degreeElevate(2);
    gsMatrix<real_t> cs_curve3_cp = cs_curve3.coefs();
    aux_cp << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1),
              pressure_side_cp.row(pressure_side_curve.coefsSize()-1);
    gsBSpline<real_t> cs_curve4 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve4.degreeElevate(2);
    gsMatrix<real_t> cs_curve4_cp = cs_curve4.coefs();

    // ---------------- outer boundary -------------------------------------------------------------------------------------------------------
    gsMatrix<real_t> cp_bs = suction_side_curve.coefs();
    real_t ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
    real_t yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

    aux_cp << length_x1, ystart_coor - pitch,
              length_x1, ystart_coor;
    gsBSpline<real_t> left_boundary_curve = gsBSpline<real_t> ( kvlin, aux_cp );
    left_boundary_curve.degreeElevate(2);
    gsMatrix<real_t> left_boundary_curve_cp = left_boundary_curve.coefs();

    //aux_cp << length_x2, yend_coor + ft*pitch,
    //          length_x2, yend_coor - fb*pitch;
    //gsBSpline<real_t> right_boundary_curve = gsBSpline<real_t> ( kvlin, aux_cp );
    //right_boundary_curve.degreeElevate(2);
    gsMatrix<real_t> rightSplitPoint1(1,2), rightSplitPoint2(1,2);
    real_t inserted_knot_right = wake_width*offset_distance/pitch;
    rightSplitPoint1(0,0) = length_x2;
    rightSplitPoint1(0,1) = (1 - inserted_knot_right) * (yend_coor - pitch) + inserted_knot_right * yend_coor;
    rightSplitPoint2(0,0) = length_x2;
    rightSplitPoint2(0,1) = inserted_knot_right * (yend_coor - pitch) + (1 - inserted_knot_right) * yend_coor;
    aux_cp << rightSplitPoint1(0,0), rightSplitPoint1(0,1),
              length_x2, yend_coor - pitch;
    gsBSpline<real_t> right_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    right_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> right_boundary_curve_A_cp = right_boundary_curve_A.coefs();
    aux_cp << rightSplitPoint1(0,0), rightSplitPoint1(0,1),
              rightSplitPoint2(0,0), rightSplitPoint2(0,1);
    gsBSpline<real_t> right_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    right_boundary_curve_B.degreeElevate(2);
    gsMatrix<real_t> right_boundary_curve_B_cp = right_boundary_curve_B.coefs();
    aux_cp << rightSplitPoint2(0,0), rightSplitPoint2(0,1),
              length_x2, yend_coor;
    gsBSpline<real_t> right_boundary_curve_C = gsBSpline<real_t> ( kvlin, aux_cp );
    right_boundary_curve_C.degreeElevate(2);
    gsMatrix<real_t> right_boundary_curve_C_cp = right_boundary_curve_C.coefs();

    aux_cp << length_x1, ystart_coor,
              suction_side_offset_cp(0,0), suction_side_offset_cp(0,1);
    gsBSpline<real_t> top_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    top_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
    aux_cp << suction_side_cp(suction_side_curve.coefsSize()-1,0), suction_side_cp(suction_side_curve.coefsSize()-1,1),
              length_x2, yend_coor;
    gsBSpline<real_t> top_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    top_boundary_curve_B.degreeElevate(2);
    gsMatrix<real_t> top_boundary_curve_B_cp = top_boundary_curve_B.coefs();

    aux_cp << length_x1, ystart_coor - pitch,
              pressure_side_offset_cp(0,0), pressure_side_offset_cp(0,1);
    gsBSpline<real_t> bottom_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    bottom_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
    aux_cp << pressure_side_cp(pressure_side_curve.coefsSize()-1,0), pressure_side_cp(pressure_side_curve.coefsSize()-1,1),
              length_x2, yend_coor - pitch;
    gsBSpline<real_t> bottom_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    bottom_boundary_curve_B.degreeElevate(2);
    gsMatrix<real_t> bottom_boundary_curve_B_cp = bottom_boundary_curve_B.coefs();

    // ---------------- cross-section curves 2 -------------------------------------------------------------------------------------------------
    gsMatrix<real_t> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
    //p0 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
    //t0 = suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-2); t0 = t0/t0.norm();
    //t1 << - length_x2 + length_x1, - yend_coor + ystart_coor; t1 = t1/t1.norm();
    //p3 << rightSplitPoint2(0,0), rightSplitPoint2(0,1);
    //gsMatrix<real_t> cs_curve5_cp(4,2);
    //cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
    //gsInfo << cs_curve5_cp << "\n";
    //gsBSpline<real_t> cs_curve5 = gsBSpline<real_t> ( kvcub, cs_curve5_cp );
    aux_cp << suction_side_offset_cp(suction_side_offset_curve.coefsSize()-1,0), suction_side_offset_cp(suction_side_offset_curve.coefsSize()-1,1),
              rightSplitPoint2(0,0), rightSplitPoint2(0,1);
    gsBSpline<real_t> cs_curve5 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve5.degreeElevate(2);
    gsMatrix<real_t> cs_curve5_cp = cs_curve5.coefs();

    //p0 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
    //t0 = pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-2); t0 = t0/t0.norm();
    //t1 << - length_x2 + length_x1, - yend_coor + ystart_coor; t1 = t1/t1.norm();
    //p3 << rightSplitPoint1(0,0), rightSplitPoint1(0,1);
    //gsMatrix<real_t> cs_curve6_cp(4,2);
    //cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
    //gsInfo << cs_curve6_cp << "\n";
    //gsBSpline<real_t> cs_curve6 = gsBSpline<real_t> ( kvcub, cs_curve6_cp );
    aux_cp << pressure_side_offset_cp(pressure_side_offset_curve.coefsSize()-1,0), pressure_side_offset_cp(pressure_side_offset_curve.coefsSize()-1,1),
              rightSplitPoint1(0,0), rightSplitPoint1(0,1);
    gsBSpline<real_t> cs_curve6 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve6.degreeElevate(2);
    gsMatrix<real_t> cs_curve6_cp = cs_curve6.coefs();

    p0 << pressure_side_offset_cp.row(0);
    dir1 << length_x1 - length_x2, ystart_coor - yend_coor,
    dir2 << ystart_coor - yend_coor, length_x2 - length_x1;
    t0 = AxisDirectionNormed(dir1, dir2);
    //dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
    //dir2 << yend_coor - ystart_coor, length_x1 - length_x2;
    //t1 = AxisDirectionNormed(dir1, dir2);
    t1 << length_x1 - suction_side_offset_cp(0,0), ystart_coor - pitch - suction_side_offset_cp(0,1); t1 = t1/t1.norm();
    p3 << suction_side_offset_cp.row(0);
    gsMatrix<real_t> cs_curve7_cp(4,2);
    cs_curve7_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve7_cp << "\n";
    gsBSpline<real_t> cs_curve7 = gsBSpline<real_t> ( kvcub, cs_curve7_cp );

    //p0 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
    ////t0 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
    //dir1 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - pressure_side_cp.row(pressure_side_curve.coefsSize()-1);
    //dir2 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
    //t0 = AxisDirectionNormed(dir1, dir2);
    ////t1 << - yend_coor + ystart_coor, length_x2 - length_x1; t1 = t1/t1.norm();
    //dir1 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - suction_side_cp.row(suction_side_curve.coefsSize()-1);
    //dir2 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
    //t1 = AxisDirectionNormed(dir1, dir2);
    //p3 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
    //gsMatrix<real_t> cs_curve8_cp(4,2);
    //cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
    //gsInfo << cs_curve8_cp << "\n";
    //gsBSpline<real_t> cs_curve8 = gsBSpline<real_t> ( kvcub, cs_curve8_cp );
    aux_cp << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1),
              suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
    gsBSpline<real_t> cs_curve8 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve8.degreeElevate(2);
    gsMatrix<real_t> cs_curve8_cp = cs_curve8.coefs();


    // ---------------- construction of patches ------------------------------------------------------------------------------------------------
    gsMatrix<real_t> patch1_cp(suction_side_cp.rows() * cs_curve1_cp.rows(), 2);
    discreteCoonsPatch(suction_side_offset_cp, suction_side_cp, cs_curve1_cp, cs_curve2_cp, patch1_cp, true);
    gsTensorBSpline<2, real_t> patch1 = gsTensorBSpline<2, real_t> (kvfit, kvcub, patch1_cp);

    gsMatrix<real_t> patch2_cp(pressure_side_cp.rows() * cs_curve3_cp.rows(), 2);
    discreteCoonsPatch(pressure_side_offset_cp, pressure_side_cp, cs_curve3_cp, cs_curve4_cp, patch2_cp, true);
    gsTensorBSpline<2, real_t> patch2 = gsTensorBSpline<2, real_t> (kvfit, kvcub, patch2_cp);

    gsMatrix<real_t> patch3_cp(top_boundary_curve_A_cp.rows() * left_boundary_curve_cp.rows(), 2);
    discreteCoonsPatch(bottom_boundary_curve_A_cp, top_boundary_curve_A_cp, left_boundary_curve_cp, cs_curve7_cp, patch3_cp, true);
    gsTensorBSpline<2, real_t> patch3 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch3_cp);

    gsMatrix<real_t> patch4_cp(suction_side_offset_cp.rows() * cs_curve7_cp.rows(), 2);
    discreteCoonsPatch(pressure_side_offset_cp, suction_side_offset_cp, cs_curve7_cp, cs_curve8_cp, patch4_cp, true);
    gsTensorBSpline<2, real_t> patch4 = gsTensorBSpline<2, real_t> (kvfit, kvcub, patch4_cp);

    gsMatrix<real_t> patch5_cp(cs_curve5_cp.rows() * cs_curve8_cp.rows(), 2);
    discreteCoonsPatch(cs_curve6_cp, cs_curve5_cp, cs_curve8_cp, right_boundary_curve_B_cp, patch5_cp, true);
    gsTensorBSpline<2, real_t> patch5 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch5_cp);

    gsMatrix<real_t> patch6_cp(cs_curve5_cp.rows() * cs_curve2_cp.rows(), 2);
    discreteCoonsPatch(cs_curve5_cp, top_boundary_curve_B_cp, cs_curve2_cp, right_boundary_curve_C_cp, patch6_cp, true);
    gsTensorBSpline<2, real_t> patch6 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch6_cp);

    gsMatrix<real_t> patch7_cp(cs_curve6_cp.rows() * cs_curve4_cp.rows(), 2);
    discreteCoonsPatch(cs_curve6_cp, bottom_boundary_curve_B_cp, cs_curve4_cp, right_boundary_curve_A_cp, patch7_cp, true);
    gsTensorBSpline<2, real_t> patch7 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch7_cp);

    // ---------------- construction of multipatch ---------------------------------------------------------------------------------------------
    gsMultiPatch<real_t> mpFinal;
    mpFinal.addPatch(patch1);
    mpFinal.addPatch(patch2);
    mpFinal.addPatch(patch3);
    mpFinal.addPatch(patch4);
    mpFinal.addPatch(patch5);
    mpFinal.addPatch(patch6);
    mpFinal.addPatch(patch7);

    mpFinal.addInterface(2, boundary::east, 3, boundary::west);
    mpFinal.addInterface(3, boundary::south, 0, boundary::north);
    mpFinal.addInterface(3, boundary::north, 1, boundary::north);
    mpFinal.addInterface(3, boundary::east, 4, boundary::west);
    mpFinal.addInterface(4, boundary::south, 5, boundary::north);
    mpFinal.addInterface(4, boundary::north, 6, boundary::north);
    // periodic interfaces
    mpFinal.addInterface(2, boundary::south, 6, boundary::north);
    mpFinal.addInterface(0, boundary::west, 1, boundary::west);
    mpFinal.addInterface(5, boundary::south, 6, boundary::south);
    mpFinal.addAutoBoundaries();

    // ---------------- plotting ---------------------------------------------------------------------------------------------------------------
    if (plot) {
        std::vector<gsGeometry<>*> curves;
        curves.clear();
        curves.push_back(&suction_side_curve);
        curves.push_back(&pressure_side_curve);
        curves.push_back(&suction_side_offset_curve);
        curves.push_back(&pressure_side_offset_curve);
        curves.push_back(&cs_curve1);
        curves.push_back(&cs_curve2);
        curves.push_back(&cs_curve3);
        curves.push_back(&cs_curve4);
        curves.push_back(&left_boundary_curve);
        curves.push_back(&right_boundary_curve_A);
        curves.push_back(&right_boundary_curve_B);
        curves.push_back(&right_boundary_curve_C);
        curves.push_back(&top_boundary_curve_A);
        curves.push_back(&top_boundary_curve_B);
        curves.push_back(&bottom_boundary_curve_A);
        curves.push_back(&bottom_boundary_curve_B);
        curves.push_back(&cs_curve5);
        curves.push_back(&cs_curve6);
        curves.push_back(&cs_curve7);
        curves.push_back(&cs_curve8);

        gsWriteParaview( curves, "section_curves", 1000);

        mpFinal.uniformRefine(); mpFinal.uniformRefine();
        gsWriteParaview( mpFinal, "patches", 50000, true);
    }

    return mpFinal;


}

gsMultiPatch<real_t> DomainAroundBladeProfile2(real_t const & index, real_t const & length_x1, real_t const & length_x2, real_t const & pitch, real_t const & camberX, real_t const & camberY, real_t const & leadingAngle, real_t const & trailingAngle, real_t const & thicknessX,
                                         real_t const & thicknessY, real_t const & endingOffset, real_t const & outputAngle, real_t const & radius, real_t const & chordLength, real_t const & Angle, real_t const & rotationCenterX, real_t const & rotationCenterY,
                                         real_t const & uniformity_param) {

    //----------------set parameters for blade profile----------------
    real_t offset_distance = 0.03;
    real_t inserted_knot = 0.3;
    real_t inserted_knot_left_1 = 0.5;
    real_t inserted_knot_left_2 = 0.9;
    real_t inserted_knot_top_1 = 0.25;
    real_t inserted_knot_top_2 = 0.65;
    real_t inserted_knot_right = 0.3;
    real_t fb = 0.5 + camberY + thicknessY;
    real_t ft = 1 - fb;
    //gsKnotVector<real_t> kvfit(0, 1, 4, 4);
    std::vector<real_t> kvfit_knots = {0.0, 0.0, 0.0, 0.0, 0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617, 1.0, 1.0, 1.0, 1.0};
    gsKnotVector<real_t> kvfit = gsKnotVector<real_t> (kvfit_knots);
    gsKnotVector<real_t> kvcub(0, 1, 0, 4);
    gsKnotVector<real_t> kvlin(0, 1, 0, 2);

    bool plot = true;
    bool plotMeshes = true;
    int num_samples = 100;
    gsVector<real_t> vec(2);
    //gsInfo << pitch << "\n ";
    vec(0) = rotationCenterX;
    vec(1) = rotationCenterY;
    gsMatrix<real_t> mat(2, 2);
    mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
           chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
    gsBSpline<real_t> suction_side_curve;
    gsBSpline<real_t> pressure_side_curve;
    gsBSpline<real_t> suction_side_offset_curve;
    gsBSpline<real_t> pressure_side_offset_curve;
    unsigned num_cpblade = 8;
    BladeProfile<real_t> * pBladeProfile = 0;

    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<real_t>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

    //---------------transform given profile----------------------
    suction_side_curve.translate(-vec);
    pressure_side_curve.translate(-vec);
    pressure_side_curve.linearTransform(mat);
    suction_side_curve.linearTransform(mat);
    pBladeProfile->setPressureSide(pressure_side_curve);
    pBladeProfile->setSuctionSide(suction_side_curve);

    pBladeProfile->computeOffset(offset_distance, suction_side_offset_curve, pressure_side_offset_curve, kvfit, num_samples);

    suction_side_curve.insertKnot(inserted_knot, suction_side_curve.degree() + 1);
    pressure_side_curve.insertKnot(inserted_knot, pressure_side_curve.degree() + 1);
    suction_side_offset_curve.insertKnot(inserted_knot, suction_side_offset_curve.degree() + 1);
    pressure_side_offset_curve.insertKnot(inserted_knot, pressure_side_offset_curve.degree() + 1);
    gsInfo << suction_side_curve.knots() << "\n";

    real_t nPointsBefore = 0;
    while ((suction_side_curve.knots())[nPointsBefore] < inserted_knot) { nPointsBefore++; }
    gsInfo << nPointsBefore << "\n";

    // -------------- knot vector for trimmed suction and pressure side curves and their trimmed offset curves --------------------------------
    gsKnotVector<real_t> kvfit2 = suction_side_curve.knots();
    kvfit2.trimLeft(nPointsBefore);
    kvfit2.affineTransformTo(0.0, 1.0);
    gsInfo << kvfit2 << "\n";

    // -------------- knot vector for leading curve --------------------------------------------------------------------------------------------
    gsKnotVector<real_t> kvfit3 = suction_side_curve.knots();
    //kvfit3.trimRight(nPointsBefore + kvfit3.multiplicity(inserted_knot));
    kvfit3.trimRight(kvfit3.size() - kvfit3.multiplicity(inserted_knot) - nPointsBefore);
    kvfit3.affineTransformTo(0.0, 1.0);
    gsInfo << kvfit3 << "\n";
    std::vector<real_t> knots_for_kvfit3(kvfit3.size() + kvfit3.degree() + (kvfit3.size() - 2 * kvfit3.degree() - 2));
    int j = 0;
    for (index_t i = kvfit3.size()-1; i > kvfit3.degree(); i--) { knots_for_kvfit3[j] = math::abs(1 - kvfit3[i])/2; j++; gsInfo << math::abs(1 - kvfit3[i])/2 << "\n"; }
    for (index_t i = 1; i < kvfit3.size(); i++) { knots_for_kvfit3[j] = 0.5 + kvfit3[i]/2; j++; gsInfo << 0.5 + kvfit3[i]/2 << "\n"; }
    kvfit3 = gsKnotVector<real_t> (knots_for_kvfit3);
    gsInfo << kvfit3 << "\n";

    // -------------- trimmed suction and pressure side curves and their trimmed offsets -------------------------------------------------------
    gsInfo << suction_side_curve.coefs() << "\n\n";
    gsInfo << pressure_side_curve.coefs() << "\n\n";
    gsInfo << suction_side_curve.coefsSize() << "\n\n";
    gsMatrix<real_t> suction_side_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2); suction_side_trimmed_cp.setZero();
    gsMatrix<real_t> pressure_side_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
    gsMatrix<real_t> suction_side_offset_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
    gsMatrix<real_t> pressure_side_offset_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
    for (index_t i = nPointsBefore; i < suction_side_curve.coefsSize(); i++) {
        gsInfo << i << "\n";
        suction_side_trimmed_cp.row(i-nPointsBefore) = suction_side_curve.coef(i);
        //gsInfo << suction_side_trimmed_cp << "\n";
        pressure_side_trimmed_cp.row(i-nPointsBefore) = pressure_side_curve.coef(i);
        suction_side_offset_trimmed_cp.row(i-nPointsBefore) = suction_side_offset_curve.coef(i);
        pressure_side_offset_trimmed_cp.row(i-nPointsBefore) = pressure_side_offset_curve.coef(i);
    }
    gsInfo << suction_side_trimmed_cp << "\n\n";
    gsInfo << pressure_side_trimmed_cp << "\n";
    gsBSpline<real_t> suction_side_trimmed = gsBSpline<real_t>( kvfit2, suction_side_trimmed_cp);
    gsBSpline<real_t> pressure_side_trimmed = gsBSpline<real_t>( kvfit2, pressure_side_trimmed_cp);
    gsBSpline<real_t> suction_side_offset_trimmed = gsBSpline<real_t>( kvfit2, suction_side_offset_trimmed_cp);
    gsBSpline<real_t> pressure_side_offset_trimmed = gsBSpline<real_t>( kvfit2, pressure_side_offset_trimmed_cp);
    gsInfo << suction_side_trimmed << "\n";

    // -------------- leading curve -------------------------------------------------------------------------------------------------------------
    gsMatrix<real_t> leading_curve_cp(kvfit3.size()-kvfit3.degree()-1, 2);
    gsMatrix<real_t> leading_offset_curve_cp(kvfit3.size()-kvfit3.degree()-1, 2);
    j = 0;
    for (index_t i = nPointsBefore-1; i > 0; i--) {
        leading_curve_cp.row(j) = suction_side_curve.coef(i);
        leading_offset_curve_cp.row(j) = suction_side_offset_curve.coef(i);
        j++;
    }
    for (index_t i = 0; i < nPointsBefore; i++) {
        leading_curve_cp.row(j) = pressure_side_curve.coef(i);
        leading_offset_curve_cp.row(j) = pressure_side_offset_curve.coef(i);
        j++;
    }
    gsInfo << leading_curve_cp << "\n";
    gsBSpline<real_t> leading_curve = gsBSpline<real_t>( kvfit3, leading_curve_cp);
    gsBSpline<real_t> leading_offset_curve = gsBSpline<real_t>( kvfit3, leading_offset_curve_cp);

    // --------------- cross-section curves ---------------------------------------------------------------------------------------------------
    gsMatrix<real_t> aux_cp(2, 2);
    aux_cp << suction_side_trimmed_cp.row(0),
              suction_side_offset_trimmed_cp.row(0);
    gsBSpline<real_t> cs_curve1 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve1.degreeElevate(2);
    gsMatrix<real_t> cs_curve1_cp = cs_curve1.coefs();
    aux_cp << pressure_side_trimmed_cp.row(0),
              pressure_side_offset_trimmed_cp.row(0);
    gsBSpline<real_t> cs_curve2 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve2.degreeElevate(2);
    gsMatrix<real_t> cs_curve2_cp = cs_curve2.coefs();
    aux_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
              pressure_side_offset_trimmed_cp.row(pressure_side_offset_trimmed.coefsSize()-1);
    gsBSpline<real_t> cs_curve3 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve3.degreeElevate(2);
    gsMatrix<real_t> cs_curve3_cp = cs_curve3.coefs();
    aux_cp << suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1),
              suction_side_offset_trimmed_cp.row(suction_side_offset_trimmed.coefsSize()-1);
    gsBSpline<real_t> cs_curve4 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve4.degreeElevate(2);
    gsMatrix<real_t> cs_curve4_cp = cs_curve4.coefs();

    // ---------------- outer boundary -------------------------------------------------------------------------------------------------------
    gsMatrix<real_t> cp_bs = suction_side_curve.coefs();
    real_t ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
    real_t yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

    gsMatrix<real_t> leftSplitPoint1(1,2), leftSplitPoint2(1,2);
    leftSplitPoint1(0,0) = length_x1;
    leftSplitPoint1(0,1) = (1 - inserted_knot_left_1) * (ystart_coor - fb*pitch) + inserted_knot_left_1 * (ystart_coor + ft*pitch);
    leftSplitPoint2(0,0) = length_x1;
    leftSplitPoint2(0,1) = (1 - inserted_knot_left_2) * (ystart_coor - fb*pitch) + inserted_knot_left_2 * (ystart_coor + ft*pitch);
    aux_cp << leftSplitPoint1(0,0), leftSplitPoint1(0,1),
              length_x1, ystart_coor - fb*pitch;
    gsBSpline<real_t> left_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    left_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> left_boundary_curve_A_cp = left_boundary_curve_A.coefs();
    aux_cp << leftSplitPoint1(0,0), leftSplitPoint1(0,1),
              leftSplitPoint2(0,0), leftSplitPoint2(0,1);
    gsBSpline<real_t> left_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    left_boundary_curve_B.degreeElevate(2);
    UnifyKnotVectors(left_boundary_curve_B, leading_offset_curve);
    gsMatrix<real_t> left_boundary_curve_B_cp = left_boundary_curve_B.coefs();
    aux_cp << leftSplitPoint2(0,0), leftSplitPoint2(0,1),
              length_x1, ystart_coor + ft*pitch;
    gsBSpline<real_t> left_boundary_curve_C = gsBSpline<real_t> ( kvlin, aux_cp );
    left_boundary_curve_C.degreeElevate(2);
    left_boundary_curve_C.insertKnot(0.333); left_boundary_curve_C.insertKnot(0.666);
    gsMatrix<real_t> left_boundary_curve_C_cp = left_boundary_curve_C.coefs();

    //aux_cp << length_x2, yend_coor + ft*pitch,
    //          length_x2, yend_coor - fb*pitch;
    //gsBSpline<real_t> right_boundary_curve = gsBSpline<real_t> ( kvlin, aux_cp );
    //right_boundary_curve.degreeElevate(2);
    gsMatrix<real_t> rightSplitPoint(1,2);
    rightSplitPoint(0,0) = length_x2;
    rightSplitPoint(0,1) = (1 - inserted_knot_right) * (yend_coor - fb*pitch) + inserted_knot_right * (yend_coor + ft*pitch);
    aux_cp << rightSplitPoint(0,0), rightSplitPoint(0,1),
            length_x2, yend_coor - fb*pitch;
    gsBSpline<real_t> right_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    right_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> right_boundary_curve_A_cp = right_boundary_curve_A.coefs();
    aux_cp << length_x2, yend_coor + ft*pitch,
              rightSplitPoint(0,0), rightSplitPoint(0,1);
    gsBSpline<real_t> right_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    right_boundary_curve_B.degreeElevate(2);
    gsMatrix<real_t> right_boundary_curve_B_cp = right_boundary_curve_B.coefs();

    gsMatrix<real_t> topSplitPoint1(1,2), topSplitPoint2(1,2);
    topSplitPoint1(0,0) = (1 - inserted_knot_top_1) * length_x1 + inserted_knot_top_1 * length_x2;
    topSplitPoint1(0,1) = (1 - inserted_knot_top_1) * (ystart_coor + ft*pitch) + inserted_knot_top_1 * (yend_coor + ft*pitch);
    topSplitPoint2(0,0) = (1 - inserted_knot_top_2) * length_x1 + inserted_knot_top_2 * length_x2;
    topSplitPoint2(0,1) = (1 - inserted_knot_top_2) * (ystart_coor + ft*pitch) + inserted_knot_top_2 * (yend_coor + ft*pitch);
    aux_cp << topSplitPoint1(0,0), topSplitPoint1(0,1),
              length_x1, ystart_coor + ft*pitch;
    gsBSpline<real_t> top_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    top_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
    aux_cp << topSplitPoint1(0,0), topSplitPoint1(0,1),
              topSplitPoint2(0,0), topSplitPoint2(0,1);
    gsBSpline<real_t> top_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    top_boundary_curve_B.degreeElevate(2);
    UnifyKnotVectors(top_boundary_curve_B, pressure_side_offset_trimmed);
    gsMatrix<real_t> top_boundary_curve_B_cp = top_boundary_curve_B.coefs();
    aux_cp << topSplitPoint2(0,0), topSplitPoint2(0,1),
              length_x2, yend_coor + ft*pitch;
    gsBSpline<real_t> top_boundary_curve_C = gsBSpline<real_t> ( kvlin, aux_cp );
    top_boundary_curve_C.degreeElevate(2);
    gsMatrix<real_t> top_boundary_curve_C_cp = top_boundary_curve_C.coefs();

    gsMatrix<real_t> bottomSplitPoint1(1,2), bottomSplitPoint2(1,2);
    bottomSplitPoint1(0,0) = (1 - inserted_knot_top_1) * length_x1 + inserted_knot_top_1 * length_x2;
    bottomSplitPoint1(0,1) = (1 - inserted_knot_top_1) * (ystart_coor - fb*pitch) + inserted_knot_top_1 * (yend_coor - fb*pitch);
    bottomSplitPoint2(0,0) = (1 - inserted_knot_top_2) * length_x1 + inserted_knot_top_2 * length_x2;
    bottomSplitPoint2(0,1) = (1 - inserted_knot_top_2) * (ystart_coor - fb*pitch) + inserted_knot_top_2 * (yend_coor - fb*pitch);
    aux_cp << bottomSplitPoint1(0,0), bottomSplitPoint1(0,1),
              length_x1, ystart_coor - fb*pitch;
    gsBSpline<real_t> bottom_boundary_curve_A = gsBSpline<real_t> ( kvlin, aux_cp );
    bottom_boundary_curve_A.degreeElevate(2);
    gsMatrix<real_t> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
    aux_cp << bottomSplitPoint1(0,0), bottomSplitPoint1(0,1),
              bottomSplitPoint2(0,0), bottomSplitPoint2(0,1);
    gsBSpline<real_t> bottom_boundary_curve_B = gsBSpline<real_t> ( kvlin, aux_cp );
    bottom_boundary_curve_B.degreeElevate(2);
    UnifyKnotVectors(bottom_boundary_curve_B, suction_side_offset_trimmed);
    gsMatrix<real_t> bottom_boundary_curve_B_cp = bottom_boundary_curve_B.coefs();
    aux_cp << bottomSplitPoint2(0,0), bottomSplitPoint2(0,1),
              length_x2, yend_coor - fb*pitch;
    gsBSpline<real_t> bottom_boundary_curve_C = gsBSpline<real_t> ( kvlin, aux_cp );
    bottom_boundary_curve_C.degreeElevate(2);
    bottom_boundary_curve_C.insertKnot(0.333); bottom_boundary_curve_C.insertKnot(0.666);
    gsMatrix<real_t> bottom_boundary_curve_C_cp = bottom_boundary_curve_C.coefs();

    // ---------------- cross-section curves 2 -------------------------------------------------------------------------------------------------
    gsMatrix<real_t> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
    p0 << suction_side_offset_trimmed_cp.row(0);
    //t0 << (suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0))/(suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0)).norm();
    dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
    dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
    t0 = AxisDirectionNormed(dir1, dir2);
    //t1 << - yend_coor + ystart_coor, length_x2 - length_x1; t1 = t1/t1.norm();
    dir1 << - yend_coor + ystart_coor, length_x2 - length_x1;
    dir2 << suction_side_offset_trimmed_cp.row(0) - bottomSplitPoint1;
    t1 = AxisDirectionNormed(dir1, dir2);
    p3 << bottomSplitPoint1(0,0), bottomSplitPoint1(0,1);
    gsMatrix<real_t> cs_curve5_cp(4,2);
    cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve5_cp << "\n";
    gsBSpline<real_t> cs_curve5 = gsBSpline<real_t> ( kvcub, cs_curve5_cp );

    p0 << pressure_side_offset_trimmed_cp.row(0);
    dir1 << pressure_side_offset_trimmed_cp.row(0) - pressure_side_trimmed_cp.row(0);
    dir2 << topSplitPoint1 - pressure_side_offset_trimmed_cp.row(0);
    t0 = AxisDirectionNormed(dir1, dir2);
    dir1 << yend_coor - ystart_coor, - length_x2 + length_x1;
    dir2 << suction_side_offset_trimmed_cp.row(0) - topSplitPoint1;
    t1 = AxisDirectionNormed(dir1, dir2);
    p3 << topSplitPoint1(0,0), topSplitPoint1(0,1);
    gsMatrix<real_t> cs_curve6_cp(4,2);
    cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve6_cp << "\n";
    gsBSpline<real_t> cs_curve6 = gsBSpline<real_t> ( kvcub, cs_curve6_cp );
    cs_curve6.insertKnot(0.333); cs_curve6.insertKnot(0.666);
    cs_curve6_cp = cs_curve6.coefs();

    p0 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
    //t0 << (pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1))/(pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1)).norm();
    dir1 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
    dir2 << topSplitPoint2 - pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
    t0 = AxisDirectionNormed(dir1, dir2);
    dir1 << yend_coor - ystart_coor, - length_x2 + length_x1;
    dir2 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - topSplitPoint2;
    t1 = AxisDirectionNormed(dir1, dir2);
    p3 << topSplitPoint2;
    gsMatrix<real_t> cs_curve7_cp(4,2);
    cs_curve7_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve7_cp << "\n";
    gsBSpline<real_t> cs_curve7 = gsBSpline<real_t> ( kvcub, cs_curve7_cp );
    cs_curve7.insertKnot(0.333); cs_curve7.insertKnot(0.666);
    cs_curve7_cp = cs_curve7.coefs();

    p0 << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1);
    //t0 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
    dir1 << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1);
    dir2 << bottomSplitPoint2 - suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1);
    t0 = AxisDirectionNormed(dir1, dir2);
    //t1 << - yend_coor + ystart_coor, length_x2 - length_x1; t1 = t1/t1.norm();
    dir1 << - yend_coor + ystart_coor, length_x2 - length_x1;
    dir2 << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - bottomSplitPoint2;
    t1 = AxisDirectionNormed(dir1, dir2);
    p3 << bottomSplitPoint2;
    gsMatrix<real_t> cs_curve8_cp(4,2);
    cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve8_cp << "\n";
    gsBSpline<real_t> cs_curve8 = gsBSpline<real_t> ( kvcub, cs_curve8_cp );

    p0 << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1);
    //dir1 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
    //dir2 << length_x2 - length_x1, yend_coor - ystart_coor;
    //t0 = AxisDirectionNormed(dir1, dir2);
    t0 << length_x2 - suction_side_offset_trimmed_cp(suction_side_trimmed.coefsSize()-1, 0), yend_coor - fb*pitch - suction_side_offset_trimmed_cp(suction_side_trimmed.coefsSize()-1, 1);
    t1 << length_x1 - length_x2, ystart_coor - yend_coor; t1 = t1/t1.norm();
    p3 << rightSplitPoint;
    gsMatrix<real_t> cs_curve9_cp(4,2);
    cs_curve9_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve9_cp << "\n";
    gsBSpline<real_t> cs_curve9 = gsBSpline<real_t> ( kvcub, cs_curve9_cp );
    cs_curve9.insertKnot(0.333); cs_curve9.insertKnot(0.666);
    cs_curve9_cp = cs_curve9.coefs();

    p3 << leftSplitPoint1;
    //dir1 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
    //dir2 << length_x2 - length_x1, yend_coor - ystart_coor;
    //t0 = AxisDirectionNormed(dir1, dir2);
    t1 << length_x2 - length_x1, yend_coor - ystart_coor; t1 = t1/t1.norm();
    gsMatrix<real_t> par(1,1); par << 0.0;
    dir1 = leading_offset_curve.deriv(par).transpose();
    //dir2 << (suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0))/(suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0)).norm();
    dir2 = cs_curve5.deriv(par).transpose();
    t0 = AxisDirectionNormed(dir1, dir2);
    p0 << suction_side_offset_trimmed_cp.row(0);
    gsMatrix<real_t> cs_curve10_cp(4,2);
    cs_curve10_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve10_cp << "\n";
    gsBSpline<real_t> cs_curve10 = gsBSpline<real_t> ( kvcub, cs_curve10_cp );

    p0 << pressure_side_offset_trimmed_cp.row(0);
    //dir1 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
    //dir2 << length_x2 - length_x1, yend_coor - ystart_coor;
    //t0 = AxisDirectionNormed(dir1, dir2);
    t1 << length_x2 - length_x1, yend_coor - ystart_coor; t1 = t1/t1.norm();
    par << 1.0;
    dir1 = - leading_offset_curve.deriv(par).transpose();
    dir2 << (pressure_side_offset_trimmed_cp.row(0) - pressure_side_trimmed_cp.row(0))/(pressure_side_offset_trimmed_cp.row(0) - pressure_side_trimmed_cp.row(0)).norm();
    t0 = AxisDirectionNormed(dir1, dir2);
    p3 << leftSplitPoint2;
    gsMatrix<real_t> cs_curve11_cp(4,2);
    cs_curve11_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve11_cp << "\n";
    gsBSpline<real_t> cs_curve11 = gsBSpline<real_t> ( kvcub, cs_curve11_cp );

    // ---------------- circular patches behind the blade ----------------------------------------------------------------------------------------
    par << 1.0;
    dir1 << pressure_side_trimmed.deriv(par).transpose();
    dir2 << suction_side_trimmed.deriv(par).transpose();
    dir2 = AxisDirectionNormed(dir1, dir2);
    dir1 = (pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1))/(pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1)).norm();
    real_t alpha = acos(dir1(0,0) * dir2(0,0) + dir1(0,1) * dir2(0,1));
    gsInfo << "alpha = " << alpha << "\n";
    gsInfo << pressure_side_offset_trimmed_cp << "\n\n";
    gsInfo << pressure_side_trimmed_cp << "\n\n";
    gsInfo << suction_side_offset_trimmed_cp << "\n\n";
    gsMatrix<real_t> circularArcPressureSide_cp = computeCircleArc2D(pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1), pressure_side_trimmed_cp(pressure_side_trimmed.coefsSize()-1,0), pressure_side_trimmed_cp(pressure_side_trimmed.coefsSize()-1,1), -alpha);
    gsInfo << "ctvrtkruh1 = " << circularArcPressureSide_cp << "\n";
    gsBSpline<real_t> circularArcPressureSide = gsBSpline<real_t> ( kvcub, circularArcPressureSide_cp );
    gsMatrix<real_t> circularArcSuctionSide_cp = computeCircleArc2D(circularArcPressureSide_cp.row(circularArcPressureSide.coefsSize()-1), suction_side_trimmed_cp(suction_side_trimmed.coefsSize()-1,0), suction_side_trimmed_cp(suction_side_trimmed.coefsSize()-1,1), -alpha);
    gsInfo << "ctvrtkruh2 = " << circularArcSuctionSide_cp << "\n";
    gsBSpline<real_t> circularArcSuctionSide = gsBSpline<real_t> ( kvcub, circularArcSuctionSide_cp );

    aux_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
              circularArcSuctionSide_cp.row(0);
    gsBSpline<real_t> cs_curve12 = gsBSpline<real_t> ( kvlin, aux_cp );
    cs_curve12.degreeElevate(2);
    gsMatrix<real_t> cs_curve12_cp = cs_curve12.coefs();
    gsInfo << cs_curve12_cp << "\n";

    gsMatrix<real_t> endpoint_curve_cp(4,2);
    endpoint_curve_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                         pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                         pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                         pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
    gsInfo << endpoint_curve_cp << "\n";

    gsMatrix<real_t> patch10_cp(endpoint_curve_cp.rows() * cs_curve3_cp.rows(), 2);
    discreteCoonsPatch(endpoint_curve_cp, circularArcPressureSide_cp, cs_curve3_cp, cs_curve12_cp, patch10_cp, true);
    gsTensorBSpline<2, real_t> patch10 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch10_cp);

    gsMatrix<real_t> patch11_cp(endpoint_curve_cp.rows() * cs_curve12_cp.rows(), 2);
    discreteCoonsPatch(endpoint_curve_cp, circularArcSuctionSide_cp, cs_curve12_cp, cs_curve4_cp, patch11_cp, true);
    gsTensorBSpline<2, real_t> patch11 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch11_cp);

    // ---------------- cross-section curves 3 -------------------------------------------------------------------------------------------------
    p0 << circularArcPressureSide_cp.row(circularArcPressureSide.coefsSize()-1);
    t0 << circularArcSuctionSide_cp.row(0) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1); t0 = t0/t0.norm();
    dir1 << length_x1-length_x2, ystart_coor - yend_coor;
    dir2 << 0, -1;
    t1 = AxisDirectionNormed(dir1, dir2);
    p3 << length_x2, yend_coor + ft*pitch;
    gsMatrix<real_t> cs_curve13_cp(4,2);
    cs_curve13_cp = LSQFergusonShort(p0, t0, t1, p3);
    gsInfo << cs_curve13_cp << "\n";
    gsBSpline<real_t> cs_curve13 = gsBSpline<real_t> ( kvcub, cs_curve13_cp );
    cs_curve13.insertKnot(0.333); cs_curve13.insertKnot(0.666);
    cs_curve13_cp = cs_curve13.coefs();

    // ---------------- construction of patches ------------------------------------------------------------------------------------------------
    gsMultiPatch<real_t> *boundaries = new gsMultiPatch<real_t>;

    gsMatrix<real_t> patch1_cp(suction_side_trimmed_cp.rows() * cs_curve1_cp.rows(), 2);
    discreteCoonsPatch(suction_side_trimmed_cp, suction_side_offset_trimmed_cp, cs_curve1_cp, cs_curve4_cp, patch1_cp, true);
    gsTensorBSpline<2, real_t> patch1 = gsTensorBSpline<2, real_t> (kvfit2, kvcub, patch1_cp);

    gsMatrix<real_t> patch2_cp(leading_curve_cp.rows() * cs_curve1_cp.rows(), 2);
    discreteCoonsPatch(leading_curve_cp, leading_offset_curve_cp, cs_curve1_cp, cs_curve2_cp, patch2_cp, true);
    gsTensorBSpline<2, real_t> patch2 = gsTensorBSpline<2, real_t> (kvfit3, kvcub, patch2_cp);

    gsMatrix<real_t> patch3_cp(pressure_side_trimmed_cp.rows() * cs_curve2_cp.rows(), 2);
    discreteCoonsPatch(pressure_side_trimmed_cp, pressure_side_offset_trimmed_cp, cs_curve2_cp, cs_curve3_cp, patch3_cp, true);
    gsTensorBSpline<2, real_t> patch3 = gsTensorBSpline<2, real_t> (kvfit2, kvcub, patch3_cp);

    gsMatrix<real_t> patch4_cp(suction_side_offset_trimmed_cp.rows() * cs_curve5_cp.rows(), 2);
    springModelPatch(suction_side_offset_trimmed_cp, bottom_boundary_curve_B_cp, cs_curve5_cp, cs_curve8_cp, patch4_cp, true);
    gsTensorBSpline<2, real_t> patch4 = gsTensorBSpline<2, real_t> (kvfit2, kvcub, patch4_cp);

    gsMatrix<real_t> patch5_cp(cs_curve5_cp.rows() * cs_curve10_cp.rows(), 2);
    springModelPatch(cs_curve5_cp, left_boundary_curve_A_cp, cs_curve10_cp, bottom_boundary_curve_A_cp, patch5_cp, true);
    gsTensorBSpline<2, real_t> patch5 = gsTensorBSpline<2, real_t> (kvcub, kvcub, patch5_cp);

    gsMatrix<real_t> patch6_cp(leading_offset_curve_cp.rows() * cs_curve10_cp.rows(), 2);
    discreteCoonsPatch(leading_offset_curve_cp, left_boundary_curve_B_cp, cs_curve10_cp, cs_curve11_cp, patch6_cp, true);
    gsTensorBSpline<2, real_t> patch6 = gsTensorBSpline<2, real_t> (kvfit3, kvcub, patch6_cp);

    gsMatrix<real_t> patch7_cp(cs_curve6_cp.rows() * cs_curve11_cp.rows(), 2);
    discreteCoonsPatch(cs_curve6_cp, left_boundary_curve_C_cp, cs_curve11_cp, top_boundary_curve_A_cp, patch7_cp, true);
    gsTensorBSpline<2, real_t> patch7 = gsTensorBSpline<2, real_t> (cs_curve6.knots(), cs_curve11.knots(), patch7_cp);

    gsMatrix<real_t> patch8_cp(pressure_side_offset_trimmed_cp.rows() * cs_curve6_cp.rows(), 2);
    discreteCoonsPatch(pressure_side_offset_trimmed_cp, top_boundary_curve_B_cp, cs_curve6_cp, cs_curve7_cp, patch8_cp, true);
    gsTensorBSpline<2, real_t> patch8 = gsTensorBSpline<2, real_t> (pressure_side_offset_trimmed.knots(), cs_curve6.knots(), patch8_cp);

    gsMatrix<real_t> patch9_cp(cs_curve8_cp.rows() * cs_curve9_cp.rows(), 2);
    discreteCoonsPatch(cs_curve8_cp, right_boundary_curve_A_cp, cs_curve9_cp, bottom_boundary_curve_C_cp, patch9_cp, true);
    gsTensorBSpline<2, real_t> patch9 = gsTensorBSpline<2, real_t> (cs_curve8.knots(), cs_curve9.knots(), patch9_cp);

    gsMatrix<real_t> patch12_cp(circularArcPressureSide_cp.rows() * cs_curve7_cp.rows(), 2);
    discreteCoonsPatch(circularArcPressureSide_cp, top_boundary_curve_C_cp, cs_curve7_cp, cs_curve13_cp, patch12_cp, true);
    gsTensorBSpline<2, real_t> patch12 = gsTensorBSpline<2, real_t> (circularArcPressureSide.knots(), cs_curve7.knots(), patch12_cp);

    gsMatrix<real_t> patch13_cp(circularArcSuctionSide_cp.rows() * cs_curve13_cp.rows(), 2);
    springModelPatch(circularArcSuctionSide_cp, right_boundary_curve_B_cp, cs_curve13_cp, cs_curve9_cp, patch13_cp, true);
    gsTensorBSpline<2, real_t> patch13 = gsTensorBSpline<2, real_t> (circularArcSuctionSide.knots(), cs_curve13.knots(), patch13_cp);

    // ---------------- construction of multipatch ---------------------------------------------------------------------------------------------
    gsMultiPatch<real_t> mpFinal;
    //mpFinal.addPatch(patch1.result());
    mpFinal.addPatch(patch1);
    //mpFinal.addPatch(patch2.result());
    mpFinal.addPatch(patch2);
    //mpFinal.addPatch(patch3.result());
    mpFinal.addPatch(patch3);
    //mpFinal.addPatch(patch4.result());
    mpFinal.addPatch(patch4);
    //mpFinal.addPatch(patch5.result());
    mpFinal.addPatch(patch5);
    //mpFinal.addPatch(patch6.result());
    mpFinal.addPatch(patch6);
    //mpFinal.addPatch(patch7.result());
    mpFinal.addPatch(patch7);
    mpFinal.addPatch(patch8);
    mpFinal.addPatch(patch9);
    mpFinal.addPatch(patch10);
    mpFinal.addPatch(patch11);
    mpFinal.addPatch(patch12);
    mpFinal.addPatch(patch13);

    mpFinal.addInterface(5, boundary::east, 6, boundary::west);
    mpFinal.addInterface(5, boundary::west, 4, boundary::west);
    mpFinal.addInterface(6, boundary::north, 7, boundary::west);
    mpFinal.addInterface(4, boundary::north, 3, boundary::west);
    mpFinal.addInterface(3, boundary::east, 8, boundary::north);
    mpFinal.addInterface(3, boundary::north, 0, boundary::south);
    mpFinal.addInterface(7, boundary::east, 11, boundary::west);
    mpFinal.addInterface(7, boundary::north, 2, boundary::south);
    mpFinal.addInterface(11, boundary::east, 12, boundary::west);
    mpFinal.addInterface(12, boundary::east, 8, boundary::west);
    mpFinal.addInterface(5, boundary::north, 1, boundary::south);
    mpFinal.addInterface(2, boundary::east, 9, boundary::west);
    mpFinal.addInterface(3, boundary::east, 10, boundary::east);
    mpFinal.addInterface(9, boundary::east, 10, boundary::west);
    mpFinal.addInterface(9, boundary::south, 11, boundary::north);
    mpFinal.addInterface(10, boundary::south, 12, boundary::north);
    mpFinal.addInterface(0, boundary::west, 1, boundary::west);
    mpFinal.addInterface(1, boundary::east, 2, boundary::west);
    // periodic interfaces
    mpFinal.addInterface(6, boundary::east, 4, boundary::east);
    mpFinal.addInterface(7, boundary::south, 3, boundary::south);
    mpFinal.addInterface(11, boundary::south, 8, boundary::east);
    mpFinal.addAutoBoundaries();

    // ---------------- plotting ---------------------------------------------------------------------------------------------------------------
    if (plot) {
        std::vector<gsGeometry<>*> curves;
        curves.clear();
        curves.push_back(&suction_side_trimmed);
        curves.push_back(&pressure_side_trimmed);
        curves.push_back(&suction_side_offset_trimmed);
        curves.push_back(&pressure_side_offset_trimmed);
        curves.push_back(&leading_curve);
        curves.push_back(&leading_offset_curve);
        curves.push_back(&cs_curve1);
        curves.push_back(&cs_curve2);
        curves.push_back(&cs_curve3);
        curves.push_back(&cs_curve4);
        curves.push_back(&left_boundary_curve_A);
        curves.push_back(&left_boundary_curve_B);
        curves.push_back(&left_boundary_curve_C);
        curves.push_back(&right_boundary_curve_A);
        curves.push_back(&right_boundary_curve_B);
        curves.push_back(&top_boundary_curve_A);
        curves.push_back(&top_boundary_curve_B);
        curves.push_back(&top_boundary_curve_C);
        curves.push_back(&bottom_boundary_curve_A);
        curves.push_back(&bottom_boundary_curve_B);
        curves.push_back(&bottom_boundary_curve_C);
        curves.push_back(&cs_curve5);
        curves.push_back(&cs_curve6);
        curves.push_back(&cs_curve7);
        curves.push_back(&cs_curve8);
        curves.push_back(&cs_curve9);
        curves.push_back(&cs_curve10);
        curves.push_back(&cs_curve11);
        curves.push_back(&circularArcPressureSide);
        curves.push_back(&circularArcSuctionSide);
        curves.push_back(&cs_curve12);
        curves.push_back(&cs_curve13);

        gsWriteParaview( curves, "section_curves", 100);

        mpFinal.uniformRefine(); mpFinal.uniformRefine();
        gsWriteParaview( mpFinal, "patches", 15000, true);
    }

    if (plotMeshes)
    {
        std::ostringstream strs_patch;
        std::string strpatch;
        gsMultiBasis<real_t> tbasis(mpFinal);
        tbasis.uniformRefine(); tbasis.uniformRefine();
        gsMesh<real_t> mesh;
        for (index_t i = 0; i < tbasis.nPieces(); i++ ) {
            mesh.cleanMesh();
            strs_patch.str().clear();
            strs_patch << i;
            strpatch =strs_patch.str();
            makeMesh(tbasis.at(i), mesh, 10);
            mpFinal.patch(i).evaluateMesh(mesh);
            gsWriteParaview(mesh, "patch_mesh_" + strpatch);
        }
    }


    return mpFinal;

}
*/

int main(int argc, char *argv[])
{
    //bool plot = true; // If set to true, paraview file is generated and launched on exit
    //bool print_info = true; // If set to true, some debug informations will be generated to console

    //int num_sample_pars = 30;
    //unsigned index_of_profile = 0; // i_blade for all profiles; 0-6, do not use 0 for geometry_between = false


    unsigned num_blades = 4;
    unsigned num_bladeprofiles = 7;

    gsVector<real_t> camber_x(num_bladeprofiles);
    gsVector<real_t> camber_y(num_bladeprofiles);
    gsVector<real_t> leading_angle(num_bladeprofiles);
    gsVector<real_t> trailing_angle(num_bladeprofiles);
    gsVector<real_t> thickness_x(num_bladeprofiles);
    gsVector<real_t> thickness_y(num_bladeprofiles);
    gsVector<real_t> ending_offset(num_bladeprofiles);
    gsVector<real_t> output_angle(num_bladeprofiles);
    gsVector<real_t> radius(num_bladeprofiles);
    gsVector<real_t> chord_length(num_bladeprofiles);
    gsVector<real_t> angle(num_bladeprofiles);
    gsVector<real_t> rotation_center_x(num_bladeprofiles);
    gsVector<real_t> rotation_center_y(num_bladeprofiles);
    gsVector<real_t> rr(num_bladeprofiles);
    gsVector<real_t> angle_input(num_bladeprofiles);
    gsVector<real_t> right_angle(num_bladeprofiles);
    right_angle.setConstant(num_bladeprofiles, 90);
    gsVector<real_t> rotate_angle((num_bladeprofiles));
    rotate_angle.setConstant(num_bladeprofiles, 0.0);

    // SETTING PARAMETERS
    rr << 0.175, 0.229, 0.283, 0.338, 0.392, 0.446, 0.5;
    camber_x << 0.425908805, 0.419, 0.416039882, 0.416576484, 0.416716376, 0.419203131, 0.43280567; // 0.411886743 na originalne druhe pozici
    //gsInfo << camber_x << "\n";
    camber_y << 0.095896144, 0.064997436, 0.037083707, 0.02724709, 0.024356984, 0.023262639, 0.019802704;
    leading_angle << 38.47692, 21.399859, 12.641776, 9.28275, 8.38282, 8.338553, 8.446091;
    leading_angle = leading_angle*EIGEN_PI/180;
    trailing_angle << 20.717182, 10.721469, 6.371868, 4.702573, 4.217487, 4.164196, 4.233241;
    trailing_angle = trailing_angle*pi/180;
    thickness_x << 0.281408739, 0.275195139, 0.273161229, 0.272459059, 0.272310247, 0.272053329, 0.271637628;
    thickness_y << 0.064950942, 0.044591162, 0.030435219, 0.020799775, 0.01464564, 0.011033004, 0.009989223;
    //ending_offset << 0.001184, 0.000812, 0.000555, 0.000379, 0.000267, 0.000201, 0.000182;
    ending_offset << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    output_angle << 8.381088, 5.891258, 4.042351, 2.765094, 1.950139, 1.466739, 1.329705;
    output_angle = output_angle*EIGEN_PI/180;
    radius << 0.020637, 0.014575, 0.007271, 0.003687, 0.001955, 0.001012, 0.000828;
    chord_length << 0.36686418, 0.42813715,	0.486298848, 0.545457321, 0.595371262, 0.615928719, 0.584588959;
    angle << 53.696487,	41.265848, 33.007703, 27.603276, 24.437586, 22.893162, 21.162381;
    angle = (right_angle - (angle + rotate_angle))*EIGEN_PI/180;
    rotation_center_x << 0.494758396, 0.469497406, 0.444542963, 0.417724545, 0.390108787, 0.361175154, 0.330805204;
    rotation_center_y << 0.060495569, 0.028225794, 0.00125711, -0.006884641, -0.010228889, -0.010435203, -0.00079539;

    gsVector<real_t> length_x1(num_bladeprofiles);
    length_x1 << -0.186559, -0.177974, -0.169389, -0.160645, -0.152061, -0.143476, -0.134891;
    length_x1 *= 2.;
    real_t length_x2 = 0.433544;

    real_t uniformity_param = 0.01;

    gsInfo << "Reading file geometry. \n";

    gsVector<int> inGeoInt(6);
    gsMatrix<int> refineUniformSettings(0,0);
    gsMatrix<int> refineLocalSettings(0,0);
    gsVector<real_t> geomParams(0);
    gsVector<bool> inGeoBool(2);
    std::vector<real_t> kvfit_knots;


    readInitialGeometry(MOTOR_DATA_DIR "uwb-pilsen/initialGeometry_domain3.txt", inGeoInt, refineUniformSettings, refineLocalSettings, geomParams, inGeoBool,kvfit_knots);

    int geomChoice =  inGeoInt(0);
    int index_of_profile = inGeoInt(1);
    int uniformRefine = inGeoInt(2);
    bool uniform_knot = inGeoBool(0);
    bool coarse = inGeoBool(1);

    if (uniform_knot)
        kvfit_knots = {0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1};
    else
        kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};

    uwbGeometryCreators<real_t> domain;
    gsMultiPatch<real_t> patches;
    gsFileData<real_t> fdpatches;
    std::ostringstream strs_domain;
    strs_domain << geomChoice;

    switch (geomChoice) {
        case 2:
            gsInfo << "DomainAroundBladeProfile2 computation ...\n";
            gsInfo << "=========================================\n";
            patches = domain.DomainAroundBladeProfile2(index_of_profile, length_x1[index_of_profile], length_x2, ((2 * EIGEN_PI*rr[index_of_profile]) / num_blades),
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
                    kvfit_knots,
                    coarse,
                    geomParams);
            break;

        case 3:
            gsInfo << "DomainBetweenBladeProfiles3 computation ...\n";
            gsInfo << "=========================================\n";
            patches = domain.DomainBetweenBladeProfiles3b(index_of_profile, length_x1[index_of_profile], length_x2, ((2 * EIGEN_PI*rr[index_of_profile]) / num_blades),
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
                    uniformity_param,
                    kvfit_knots,
                    coarse,
                    geomParams);
            break;

        case 5:
            gsInfo << "DomainBetweenBladeProfiles5 computation ...\n";
            gsInfo << "===========================================\n";
            patches = domain.DomainBetweenBladeProfiles5(index_of_profile, length_x1[index_of_profile], length_x2, ((2 * EIGEN_PI*rr[index_of_profile]) / num_blades),
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
                    uniformity_param,
                    kvfit_knots,
                    coarse,
                    geomParams);
            break;
        default:
            GISMO_ASSERT((geomChoice != 2) && (geomChoice != 5), "Not supported index of geometry selection!");
    }

    for (int i = 0; i < uniformRefine; i++)
        patches.uniformRefine();

    fdpatches << patches;
    fdpatches.save("domain" + strs_domain.str() + "_multipatch.xml");

    gsInfo << "Konec\n";

    //patches.addInterface(0, boundary::north, 0, boundary::south);

    return 0;

}
