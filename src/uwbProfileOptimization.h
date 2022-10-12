#ifndef UWBPROFILEOPTIMIZATION_H
#define UWBPROFILEOPTIMIZATION_H


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <gismo.h>

#include <sstream>
#include <vector>
#include <math.h>


#include "uwbBladeProfile.h"
#include "uwbTurbineUtils.h"

namespace gismo
{

template<class TT>
gsMultiPatch<TT>  BSplineProfile2D(TT const & camberX, TT const & camberY, TT const & leadingAngle, TT const & trailingAngle, TT const & thicknessX, TT const & thicknessY, TT const & endingOffset, TT const & outputAngle, TT const & radius, TT const & chordLength,
                                                     TT const & Angle, TT const & rotationCenterX, TT const & rotationCenterY)

{
    //----------------set parameters for blade profile----------------
    int num_samples = 30;
    gsBSpline<TT> suction_side_curve;
    gsBSpline<TT> pressure_side_curve;
    gsKnotVector<TT> kvfit(0, 1, 4, 4);
    unsigned num_cpblade = 8;
    BladeProfile<TT> * pBladeProfile = 0;

    //---------------set parameters for boundary of patches-----------------------
    //point [length_x1, width_y1] is the lower left point of boundary rectangle
    //point [length_x2, width_y2] is the upper right point of boundary rectangle
    real_t width_y2 = 0.91412/2 + 0.3;
    real_t width_y1 = -0.91412/2 - 0.5;
    real_t length_x1 = -0.8;
    real_t length_x2 = chordLength + 1.0;
    real_t length = math::abs(length_x1) + math::abs(length_x2);
    unsigned degree = 3;

    //---------------compute blade profile for given parameters----------------------
    pBladeProfile = new BladeProfile<TT>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength,
                                    Angle, rotationCenterX, rotationCenterY);
    pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

    suction_side_curve.rotate(-Angle);
    pressure_side_curve.rotate(-Angle);
    pBladeProfile->setSuctionSide(suction_side_curve);
    pBladeProfile->setPressureSide(pressure_side_curve);

    gsBSpline<TT> bp = pBladeProfile->getPressureSide();
    gsBSpline<TT> bs = pBladeProfile->getSuctionSide();
    gsMatrix<TT> cp_bp (num_cpblade,2);
    gsMatrix<TT> cp_bs (num_cpblade,2);

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

    gsMultiPatch<TT> mp;



    //initial data for patches without blade profile - patches are linear -> elevate
    gsKnotVector<TT> kv_u(0, 1, 0, 2);
    gsKnotVector<TT> kv_v(0, 1, 0, 2);
    gsTensorBSplineBasis<2, TT> basis(kv_u, kv_v);
    gsTensorBSplineBasis<2, TT> basis_blade(kvfit, kv_v);

    //--------------------------------patch 0-------------------------------------------
    gsMatrix<TT> coef_patch0(4, 2);
    coef_patch0 <<  length_x1, 0,
                    bp.coef(0, 0), bp.coef(0, 1),
                    length_x1, width_y2,
                    length_x1 + length/4, width_y2;


    gsTensorBSpline<2, real_t> patch0(basis, coef_patch0);
    patch0.degreeElevate(degree - 2, -1);

    //--------------------------------patch 1-------------------------------------------
    gsMatrix<TT> coef_patch1 (4, 2);
    coef_patch1 <<  length_x1, width_y1,
                    length_x1 + length/4, width_y1,
                    length_x1, 0,
                    bs.coef(0, 0), bs.coef(0, 1);


    gsTensorBSpline<2, TT> patch1(basis, coef_patch1);
    patch1.degreeElevate(degree - 2, -1);
    //--------------------------------patch 4-------------------------------------------
    gsMatrix<TT> coef_patch4(4, 2);
    coef_patch4 <<  bp.coef(num_cpblade-1, 0), bp.coef(num_cpblade-1, 1),
                    length_x2, bp.coef(num_cpblade-1, 1),
                    length_x2 - length/4, width_y2,
                    length_x2, width_y2;

    gsTensorBSpline<2, TT> patch4(basis, coef_patch4);
    patch4.degreeElevate(degree - 2, -1);
    //--------------------------------patch 5-------------------------------------------
    gsMatrix<TT> coef_patch5 (4, 2);
    coef_patch5 <<  length_x2 - length/4, width_y1,
                    length_x2, width_y1,
                    bs.coef(num_cpblade - 1, 0), bs.coef(num_cpblade - 1, 1),
                    length_x2, bs.coef(num_cpblade - 1, 1);
gsInfo << coef_patch5 << "\n";
    gsTensorBSpline<2, real_t> patch5(basis, coef_patch5);
    patch5.degreeElevate(degree - 2, -1);

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
            coef_patch2(i, 0)= length_x1 + (length/4.0) +   chordalbs(i-num_cpblade) * (length/2);
            coef_patch2(i, 1)= width_y2 ;
        }


    gsTensorBSpline<2, TT> patch2(basis_blade, coef_patch2);
    patch2.degreeElevate(degree - 2, 1);

    //--------------------------------patch 3-------------------------------------------
    gsMatrix<TT> chordalbp = chordalParameterization(cp_bp);
    gsMatrix<TT> coef_patch3 (2 * num_cpblade, 2);
    for(unsigned i = num_cpblade; i < 2 * num_cpblade; i++){
            coef_patch3(i, 0)= bp.coef(i -num_cpblade, 0);
            coef_patch3(i, 1)= bp.coef(i -num_cpblade, 1) ;
        }
    for(unsigned i = 0; i < num_cpblade; i++){
            coef_patch3(i, 0)= length_x1 + length/4 +   chordalbp(i) * (length/2);
            coef_patch3(i, 1)= width_y1 ;
        }


    gsTensorBSpline<2, TT> patch3(basis_blade, coef_patch3);
    patch3.degreeElevate(degree - 2, 1);

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


    return mp;
}

//---------------------Example for parameters------------------------
/*
// promenne, ktere lze menit
real_t camberX = 0.41604;               //x_ova souradnice max prohnuti strednice
real_t camberY =  0.0370837;            //y_ova souradnice max prohnuti strednice
real_t leadingAngle = 0.220641;         //vstupni uhel strednice
real_t trailingAngle = 0.11121;         //vystupni uhel strednice
real_t thicknessX = 0.273161;           //x_ova souradnice max tloustky
real_t thicknessY = 0.0304352;          //y_ova souradnice max tloustky
real_t outputAngle = 0.11121;           //vystupni uhel tlousky
real_t radius = 0.007271;               //polomer osk. kr. v nabehu
real_t chordLength = 1.0;               //delka tetivy profilu


// promenne, ktere "nelze" menit, tzn. jen v pripade dalsich zmen v souboru uwbProfileOptimization.h
real_t endingOffset = 0.0;                //odsazeni tloustky
real_t Angle = leadingAngle; //natoceni profilu, nastaveno tak, aby vstupni rychlost nabihala na profil ve smeru osy x
real_t rotationCenterX = 0.0;           //stred otaceni profilu
real_t rotationCenterY = 0.0;           //stred otaceni profilu



 gsMultiPatch<> patches;
 patches = BSplineProfile2D(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength,
                             Angle, rotationCenterX, rotationCenterY);
*/

} // end namespace gismo


#endif // UWBPROFILEOPTIMIZATION

