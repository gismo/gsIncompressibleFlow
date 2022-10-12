/** @file uwbGeometryCreators.h

Author(s): K. Michalkova, E. Turnerova
*/

#pragma once

#include "uwbTurbineUtils.h"
#include "uwbBladeProfile.h"

//#include <gsIO/gsIOUtils.h>
#include <gsOptimizer/gsQualityMeasure.h>

#include "../jku/gsMotorUtils.h"
//#include "../jku/gsQualityMeasure2.h"

//#include "uwbProfileOptimization.h"
#include "gsModeling/gsCoonsPatch.h"
//#include <../../motor/jku/gsMotorUtils.h>

namespace gismo
{

template<class T>
class uwbGeometryCreators
{
public:
    uwbGeometryCreators() { }

    ~uwbGeometryCreators() { }

public:
    gsMultiPatch<T> BSplineBackwardStep2D(T const & a, T const & b, T const & a_in)
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

    gsMultiPatch<T> BSplineBackwardStep2D_4patches(int const & deg, T const & stepHeight, T const & preStepHeight, T const & preStepLength, T const & prePatchLength, T const & length)
    {
        //T domainHeight = stepHeight + preStepHeight;

        gsMultiPatch<T> mp;

        /*mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, 0.0, length, stepHeight));
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0, stepHeight, length, domainHeight));
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-preStepLength, stepHeight, 0.0, domainHeight));
        mp.addPatch(gsNurbsCreator<>::BSplineRectangle(-preStepLength - prePatchLength, stepHeight, -preStepLength, domainHeight));
        */
        //int deg = 1;
        mp.addPatch(BSplineRect(deg, 0.0, 0.0, length, stepHeight));
        mp.addPatch(BSplineRect(deg, 0.0, stepHeight, length, preStepHeight));

        mp.addPatch(BSplineRect(deg, -preStepLength, stepHeight, preStepLength, preStepHeight));


        //mp.addPatch(BSplineRect(deg, -preStepLength, stepHeight, preStepLength, preStepHeight));
        //mp.addPatch(BSplineRect(deg, -preStepLength - prePatchLength, stepHeight, prePatchLength, preStepHeight));

        mp.addInterface(0, boundary::north, 1, boundary::south);
        mp.addInterface(1, boundary::west, 2, boundary::east);
        //mp.addInterface(2, boundary::west, 3, boundary::east);
        mp.addAutoBoundaries();

        return mp;
    }

    gsTensorBSpline<2, T> BSplineRect(int deg, const T llx, const T lly, const T a, const T b) // llx - lower left x, lly - lower left y
    {
        gsKnotVector<T> kv(0, 1, 0, deg + 1); // first, last, inter, mult_end

        int n = kv.size() - deg - 1;
        gsMatrix<T> coef(n*n, 2);

        switch (deg)
        {
        case 1:
        {
            coef << llx + 0, lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + b,
                llx + a, lly + b;
            break;
        }
        case 2:
        {
            coef << llx + 0, lly + 0,
                llx + (a / 2), lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + (b / 2),
                llx + (a / 2), lly + (b / 2),
                llx + a, lly + (b / 2),
                llx + 0, lly + b,
                llx + (a / 2), lly + b,
                llx + a, lly + b;
            break;
        }
        case 3:
        {
            coef << llx + 0, lly + 0,
                llx + (a / 3), lly + 0,
                llx + (2. / 3) * a, lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + (b / 3),
                llx + (a / 3), lly + (b / 3),
                llx + (2. / 3) * a, lly + (b / 3),
                llx + a, lly + (b / 3),
                llx + 0, lly + (2. / 3) * b,
                llx + (a / 3), lly + (2. / 3) * b,
                llx + (2. / 3) * a, lly + (2. / 3) * b,
                llx + a, lly + (2. / 3) * b,
                llx + 0, lly + b,
                llx + (a / 3), lly + b,
                llx + (2. / 3) * a, lly + b,
                llx + a, lly + b;
            break;
        }
        case 4:
        {
            coef << llx + 0, lly + 0,
                llx + (a / 4), lly + 0,
                llx + (2. / 4) * a, lly + 0,
                llx + (3. / 4) * a, lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + (b / 4),
                llx + (a / 4), lly + (b / 4),
                llx + (2. / 4) * a, lly + (b / 4),
                llx + (3. / 4) * a, lly + (b / 4),
                llx + a, lly + (b / 4),
                llx + 0, lly + (2. / 4) * b,
                llx + (a / 4), lly + (2. / 4) * b,
                llx + (2. / 4) * a, lly + (2. / 4) * b,
                llx + (3. / 4) * a, lly + (2. / 4) * b,
                llx + a, lly + (2. / 4) * b,
                llx + 0, lly + (3. / 4) * b,
                llx + (a / 4), lly + (3. / 4) * b,
                llx + (2. / 4) * a, lly + (3. / 4) * b,
                llx + (3. / 4) * a, lly + (3. / 4) * b,
                llx + a, lly + (3. / 4) * b,
                llx + 0, lly + b,
                llx + (a / 4), lly + b,
                llx + (2. / 4) * a, lly + b,
                llx + (3. / 4) * a, lly + b,
                llx + a, lly + b;
            break;
        }
        case 5:
        {
            coef << llx + 0, lly + 0,
                llx + (a / 5), lly + 0,
                llx + (2. / 5) * a, lly + 0,
                llx + (3. / 5) * a, lly + 0,
                llx + (4. / 5) * a, lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + (b / 5),
                llx + (a / 5), lly + (b / 5),
                llx + (2. / 5) * a, lly + (b / 5),
                llx + (3. / 5) * a, lly + (b / 5),
                llx + (4. / 5) * a, lly + (b / 5),
                llx + a, lly + (b / 5),
                llx + 0, lly + (2. / 5) * b,
                llx + (a / 5), lly + (2. / 5) * b,
                llx + (2. / 5) * a, lly + (2. / 5) * b,
                llx + (3. / 5) * a, lly + (2. / 5) * b,
                llx + (4. / 5) * a, lly + (2. / 5) * b,
                llx + a, lly + (2. / 5) * b,
                llx + 0, lly + (3. / 5) * b,
                llx + (a / 5), lly + (3. / 5) * b,
                llx + (2. / 5) * a, lly + (3. / 5) * b,
                llx + (3. / 5) * a, lly + (3. / 5) * b,
                llx + (4. / 5) * a, lly + (3. / 5) * b,
                llx + a, lly + (3. / 5) * b,
                llx + 0, lly + (4. / 5) * b,
                llx + (a / 5), lly + (4. / 5) * b,
                llx + (2. / 5) * a, lly + (4. / 5) * b,
                llx + (3. / 5) * a, lly + (4. / 5) * b,
                llx + (4. / 5) * a, lly + (4. / 5) * b,
                llx + a, lly + (4. / 5) * b,
                llx + 0, lly + b,
                llx + (a / 5), lly + b,
                llx + (2. / 5) * a, lly + b,
                llx + (3. / 5) * a, lly + b,
                llx + (4. / 5) * a, lly + b,
                llx + a, lly + b;
            break;
        }
        case 6:
        {
            coef << llx + 0, lly + 0,
                llx + (a / 6), lly + 0,
                llx + (2. / 6) * a, lly + 0,
                llx + (3. / 6) * a, lly + 0,
                llx + (4. / 6) * a, lly + 0,
                llx + (5. / 6) * a, lly + 0,
                llx + a, lly + 0,
                llx + 0, lly + (b / 6),
                llx + (a / 6), lly + (b / 6),
                llx + (2. / 6) * a, lly + (b / 6),
                llx + (3. / 6) * a, lly + (b / 6),
                llx + (4. / 6) * a, lly + (b / 6),
                llx + (5. / 6) * a, lly + (b / 6),
                llx + a, lly + (b / 6),
                llx + 0, lly + (2. / 6) * b,
                llx + (a / 6), lly + (2. / 6) * b,
                llx + (2. / 6) * a, lly + (2. / 6) * b,
                llx + (3. / 6) * a, lly + (2. / 6) * b,
                llx + (4. / 6) * a, lly + (2. / 6) * b,
                llx + (5. / 6) * a, lly + (2. / 6) * b,
                llx + a, lly + (2. / 6) * b,
                llx + 0, lly + (3. / 6) * b,
                llx + (a / 6), lly + (3. / 6) * b,
                llx + (2. / 6) * a, lly + (3. / 6) * b,
                llx + (3. / 6) * a, lly + (3. / 6) * b,
                llx + (4. / 6) * a, lly + (3. / 6) * b,
                llx + (5. / 6) * a, lly + (3. / 6) * b,
                llx + a, lly + (3. / 6) * b,
                llx + 0, lly + (4. / 6) * b,
                llx + (a / 6), lly + (4. / 6) * b,
                llx + (2. / 6) * a, lly + (4. / 6) * b,
                llx + (3. / 6) * a, lly + (4. / 6) * b,
                llx + (4. / 6) * a, lly + (4. / 6) * b,
                llx + (5. / 6) * a, lly + (4. / 6) * b,
                llx + a, lly + (4. / 6) * b,
                llx + 0, lly + (5. / 6) * b,
                llx + (a / 6), lly + (5. / 6) * b,
                llx + (2. / 6) * a, lly + (5. / 6) * b,
                llx + (3. / 6) * a, lly + (5. / 6) * b,
                llx + (4. / 6) * a, lly + (5. / 6) * b,
                llx + (5. / 6) * a, lly + (5. / 6) * b,
                llx + a, lly + (5. / 6) * b,
                llx + 0, lly + b,
                llx + (a / 6), lly + b,
                llx + (2. / 6) * a, lly + b,
                llx + (3. / 6) * a, lly + b,
                llx + (4. / 6) * a, lly + b,
                llx + (5. / 6) * a, lly + b,
                llx + a, lly + b;
            break;
        }
        default:
            GISMO_ERROR("Degree not implemented.");
            break;
        }

        return gsTensorBSpline<2, T>(kv, kv, give(coef));
    }

    gsMultiPatch<T> BSplineBackwardStep2D1Patch()
    {
        // create knot vector [0,0,0, 0.25,0.5,0.75 1,1,1]
        gsKnotVector<T> tK1(0,1,3,3);
        // create knot vector [0,0,0, 0.5, 1,1,1]
       std::vector<T> knots = {0., 0., 0., 0.0833333, 0.166667, 0.25, 0.333333, 0.333333, 0.416667, 0.5, 0.583333, 0.666667, 0.666667, 0.75, 0.833333, 0.916667, 1., 1.,
                  1.};
        gsKnotVector<T> tK2(knots);

        gsMatrix<T> C(96,2);

        C << -1., 1.,
                -0.875, 1.,
                -0.625, 1.,
                -0.375, 1.,
                -0.125, 1.,
                0., 1.,
                0.3125, 1.,
                0.9375, 1.,
                1.5625, 1.,
                2.1875, 1.,
                2.5, 1.,
                2.8125, 1.,
                3.4375, 1.,
                4.0625, 1.,
                4.6875, 1.,
                5., 1.,
                -1., 0.875,
                -0.875, 0.875,
                -0.625, 0.875,
                -0.375, 0.875,
                -0.125, 0.875,
                0., 0.875,
                0.2734375, 0.859375,
                0.8203125, 0.828125,
                1.3671875, 0.796875,
                1.9140625, 0.765625,
                2.1875, 0.75,
                2.5390625, 0.75,
                3.2421875, 0.75,
                3.9453125, 0.75,
                4.6484375, 0.75,
                5., 0.75,
                -1., 0.625,
                -0.875, 0.625,
                -0.625, 0.625,
                -0.375, 0.625,
                -0.125, 0.625,
                0., 0.625,
                0.1953125, 0.578125,
                0.5859375, 0.484375,
                0.9765625, 0.390625,
                1.3671875, 0.296875,
                1.5625, 0.25,
                1.9921875, 0.25,
                2.8515625, 0.25,
                3.7109375, 0.25,
                4.5703125, 0.25,
                5., 0.25,
                -1., 0.375,
                -0.875, 0.375,
                -0.625, 0.375,
                -0.375, 0.375,
                -0.125, 0.375,
                0., 0.375,
                0.1171875, 0.296875,
                0.3515625, 0.140625,
                0.5859375, -0.015625,
                0.8203125, -0.171875,
                0.9375, -0.25,
                1.4453125, -0.25,
                2.4609375, -0.25,
                3.4765625, -0.25,
                4.4921875, -0.25,
                5., -0.25,
                -1., 0.125,
                -0.875, 0.125,
                -0.625, 0.125,
                -0.375, 0.125,
                -0.125, 0.125,
                0., 0.125,
                0.0390625, 0.015625,
                0.1171875, -0.203125,
                0.1953125, -0.421875,
                0.2734375, -0.640625,
                0.3125, -0.75,
                0.8984375, -0.75,
                2.0703125, -0.75,
                3.2421875, -0.75,
                4.4140625, -0.75,
                5., -0.75,
                -1., 0.,
                -0.875, 0.,
                -0.625, 0.,
                -0.375, 0.,
                -0.125, 0.,
                0., 0.,
                0., -0.125,
                0., -0.375,
                0., -0.625,
                0., -0.875,
                0., -1.,
                0.625, -1.,
                1.875, -1.,
                3.125, -1.,
                4.375, -1.,
                5., -1.;

        gsVector<T> vConstant;
        vConstant.setConstant(C.rows(), 0.25);
        C = C/4;
        C.middleCols(1, 1) = C.col(1) + vConstant;
        for (int i = 0; i < C.rows(); i++)
        {
            if (C(i, 0) > 0)
            {
                C(i, 0) = C(i, 0) * 1.6;
            }
        }

        gsMultiPatch<T> mp;
        gsTensorBSpline<2, T> patch = gsTensorBSpline<2,T>(tK2,tK1, give(C));
        mp.addPatch(patch);
        mp.addAutoBoundaries();

        return mp;
    }

    gsMultiPatch<T>  BSplineCircle2D4Patches(const T & length, T const & width, T const & widthExpand, T const & radius,
                                             T const & centreX, T const & centreY, bool const &prepatch, T const & prepatchWidth, bool periodic = true)
    {
        const real_t PI = 3.141592653589793238463;
        unsigned degree = 3;
        real_t alpha = PI/2.0;
        real_t const K = (4.0/3.0)*(tan(alpha/4.0)); //constant for circle
        bool plotMeshes = false;


        gsMultiPatch<T> mp;




        //initial data for patches
        gsKnotVector<T> kv_cub(0, 1, 0, 4);
        gsKnotVector<T> kv_lin(0, 1, 0, 2);
        gsTensorBSplineBasis<2, T> basis(kv_lin, kv_cub);

        //--------------------------------patch 0-------------------------------------------
        gsMatrix<T> coef_patch0(8, 2);
        coef_patch0 <<  - width/2, - width/2,
                        - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                         - width/2, - (width + widthExpand)/4,
                        - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius * sin(alpha/2)),
                          - width/2, (width + widthExpand)/4,
                        - radius*cos(alpha/2) + centreX +  K * ( - radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius *  sin(alpha/2)),
                        - width/2,  width/2 + widthExpand,
                        - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
        gsTensorBSpline<2, T> patch0(basis, coef_patch0);
        patch0.degreeElevate(degree - 1, 0);


        //--------------------------------patch 1-------------------------------------------
        gsMatrix<T> coef_patch1(8, 2);
        coef_patch1 <<  - width/2, width/2 + widthExpand,
                - radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY,
                        - width/4, width/2 + widthExpand,
                - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                         width/4, width/2 + widthExpand,
                            radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +  K * ( radius * sin(alpha/2)),
                         + width/2, width/2 + widthExpand,

                                    radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
        gsTensorBSpline<2, T> patch1(basis, coef_patch1);
        patch1.degreeElevate(degree - 1, 0);

        //--------------------------------patch 2-------------------------------------------
        gsMatrix<T> coef_patch2(8, 2);
        coef_patch2 <<  - width/2,- width/2,
                       - radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                 - width/4, - width/2,
                - radius*cos(alpha/2) + centreX + K * ( radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
                 width/4, - width/2,
                radius*cos(alpha/2) + centreX + K * (- radius * cos(alpha/2)),  - radius*sin(alpha/2) + centreY +  K * ( - radius * sin(alpha/2)),
                  + width/2,- width/2,
                radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY;

        gsTensorBSpline<2, T> patch2(basis, coef_patch2);
        patch2.degreeElevate(degree - 1, 0);

        //--------------------------------patch 3-------------------------------------------
        gsMatrix<T> coef_patch3(8, 2);
        coef_patch3 << width/2, - width/2,
                       radius*cos(alpha/2) + centreX, - radius*sin(alpha/2) + centreY,
                          width/2, -(width + widthExpand)/4,
                        radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)), - radius*sin(alpha/2) + centreY +   K * (radius*sin(alpha/2)),
                       width/2,  (width + widthExpand)/4,
                        radius*cos(alpha/2) + centreX +  K * ( radius * cos(alpha/2)),  radius*sin(alpha/2) + centreY +   K * (- radius * sin(alpha/2)),
                        width/2,  (width/2 + widthExpand),
                         radius*cos(alpha/2) + centreX, radius*sin(alpha/2) + centreY;
        gsTensorBSpline<2, T> patch3(basis, coef_patch3);
        patch3.degreeElevate(degree - 1, 0);

        //--------------------------------patch 4-------------------------------------------
        gsMatrix<T> coef_patch4(8, 2);
        coef_patch4 << length - width/2, - width/2,
                       width/2, - width/2,
                          length - width/2, - (width + widthExpand)/4,
                        width/2, - (width + widthExpand)/4,
                       length - width/2,  (width + widthExpand)/4,
                        width/2, + (width + widthExpand)/4,
                        length - width/2,  (width/2 + widthExpand),
                       width/2,  (width/2 + widthExpand);
        gsTensorBSpline<2, T> patch4(basis, coef_patch4);
        patch4.degreeElevate(degree - 1, 0);

        //refining last patch according to length/width
        real_t insertKnotNum = trunc(length/width)-1;
        if(insertKnotNum>1.0){
            for(real_t i = 1/(2*insertKnotNum); i < 1.0; i+=1/(2*insertKnotNum)){
            patch4.insertKnot(i,0,1);
            }
        }


        mp.addPatch(patch0);
        mp.addPatch(patch1);
        mp.addPatch(patch2);
        mp.addPatch(patch3);
        mp.addPatch(patch4);

        mp.addInterface(0, boundary::north, 1, boundary::south);
        mp.addInterface(0, boundary::south, 2, boundary::south);
        mp.addInterface(3, boundary::north, 1, boundary::north);
        mp.addInterface(3, boundary::south, 2, boundary::north);
        mp.addInterface(4, boundary::east, 3, boundary::west);
        if (periodic)
        {
            mp.addInterface(1, boundary::west, 2, boundary::west);
            mp.addInterface(4, boundary::north, 4, boundary::south);
        }

        //add start patch
        if(prepatch){
            //--------------------------------patch 5-------------------------------------------
            gsMatrix<T> coef_patch5(8, 2);
            coef_patch5 << -prepatchWidth - width/2, - width/2,
                    - width/2, - width/2,
                     - prepatchWidth - width/2, - (width + widthExpand)/4,
                    - width/2 , - (width + widthExpand)/4,
                      - prepatchWidth - width/2, (width + widthExpand)/4,
                    - width/2 , (width + widthExpand)/4,
                    - prepatchWidth - width/2,  width/2 + widthExpand,
                    - width/2,  width/2 + widthExpand;
            gsTensorBSpline<2, T> patch5(basis, coef_patch5);
            patch5.degreeElevate(degree - 1, 0);

            real_t insertKnotNumPrepatch = trunc(prepatchWidth/width);

            if(insertKnotNumPrepatch > 0.00001){
                for(real_t i = 1/(2*insertKnotNumPrepatch); i < 1.0; i+=1/(2*insertKnotNumPrepatch)){

                patch5.insertKnot(i,0,1);
                }
            }


            mp.addPatch(patch5);
            mp.addInterface(0, boundary::west, 5, boundary::east);
            if (periodic)
                mp.addInterface(5, boundary::north, 5, boundary::south);

            if (plotMeshes){
                gsMesh<> mesh5;
                patch5.controlNet(mesh5);
                gsWriteParaview(mesh5,"cpCirclePatch5");
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

    //8 patches around circle
    gsMultiPatch<T>  BSplineCircle2D(const T & length, T const & width, T const & widthExpand, T const & radius, T const & centreX, T const & centreY, bool const &prepatch, T const & prepatchWidth)
    {
        const real_t PI = 3.141592653589793238463;
        unsigned degree = 3;
        real_t alpha = PI/2.0;
        real_t const K = (4.0/3.0)*(tan(alpha/4.0)); //constant for circle
        bool plotMeshes = true;


        gsMultiPatch<T> mp;

        //initial data for patches
        gsKnotVector<T> kv_cub(0, 1, 0, 4);
        gsKnotVector<T> kv_lin(0, 1, 0, 2);
        gsTensorBSplineBasis<2, T> basis_lc(kv_lin, kv_cub);
        gsTensorBSplineBasis<2, T> basis_lin(kv_lin, kv_lin);

        gsMatrix<T> boundary_02(4,2);
        gsMatrix<T> boundary_24(4,2);
        gsMatrix<T> boundary_46(4,2);
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
        gsMatrix<T> boundary_13(4,2);
        gsMatrix<T> boundary_35(4,2);
        gsMatrix<T> boundary_57(4,2);
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


            gsMatrix<T> coef_patch0(8, 2);
            coef_patch0 <<  - width/2, widthExpand/2,
                            boundary_02(0,0), boundary_02(0,1),
                            - width/2, (4*widthExpand + width)/6,
                            boundary_02(1,0), boundary_02(1,1),
                            - width/2, (5*widthExpand + 2*width)/6,
                            boundary_02(2,0), boundary_02(2,1),
                            - width/2, width/2 + widthExpand,
                            boundary_02(3,0), boundary_02(3,1);
            gsTensorBSpline<2, T> patch0(basis_lc, coef_patch0);
            patch0.degreeElevate(degree - 1, 0);


        //--------------------------------patch 1-------------------------------------------
        gsMatrix<T> coef_patch1(8, 2);
        coef_patch1 <<  - width/2, - width/2,
                         boundary_13(0,0), boundary_13(0,1),
                        - width/2, (widthExpand - 2*width)/6,
                boundary_13(1,0), boundary_13(1,1),
                        -width/2, (2*widthExpand - width)/6,
                boundary_13(2,0), boundary_13(2,1),
                        - width/2, widthExpand/2,
                boundary_13(3,0), boundary_13(3,1);

        gsTensorBSpline<2, T> patch1(basis_lc, coef_patch1);
        patch1.degreeElevate(degree - 1, 0);

        //--------------------------------patch 2-------------------------------------------

        gsMultiPatch<T> boundaries2;
        gsMatrix<T> boundary_2s(4,2);
        gsMatrix<T> boundary_2n(4,2);

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

        gsMultiPatch<T> boundaries3;
        gsMatrix<T> boundary_3s(4,2);
        gsMatrix<T> boundary_3n(4,2);

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
        gsMultiPatch<T> boundaries4;
        gsMatrix<T> boundary_4s(4,2);
        gsMatrix<T> boundary_4n(4,2);

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
        gsMultiPatch<T> boundaries5;
        gsMatrix<T> boundary_5s(4,2);
        gsMatrix<T> boundary_5n(4,2);

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
        gsMatrix<T> coef_patch6(8, 2);
        coef_patch6 << boundary_46(0,0), boundary_46(0,1),
                        width/2, centreY + widthExpand/2,
                         boundary_46(1,0), boundary_46(1,1),
                        width/2, (4*widthExpand + width)/6,
                       boundary_46(2,0), boundary_46(2,1),
                        width/2, (5*widthExpand + 2*width)/6,
                        boundary_46(3,0), boundary_46(3,1),
                        width/2, width/2 + widthExpand;

        gsTensorBSpline<2, T> patch6(basis_lc, coef_patch6);
        patch6.degreeElevate(degree - 1, 0);


        //--------------------------------patch 7-------------------------------------------
        gsMatrix<T> coef_patch7(8, 2);
        coef_patch7 <<  boundary_57(0,0), boundary_57(0,1),
                        + width/2, - width/2,
                        boundary_57(1,0), boundary_57(1,1),
                        + width/2, (widthExpand - 2*width)/6,
                        boundary_57(2,0), boundary_57(2,1),
                        + width/2, (2*widthExpand - width)/6,
                        boundary_57(3,0), boundary_57(3,1),
                        + width/2, widthExpand/2;

        gsTensorBSpline<2, T> patch7(basis_lc, coef_patch7);
        patch7.degreeElevate(degree - 1, 0);


        //--------------------------------patch 8-------------------------------------------
        gsMatrix<T> coef_patch8(4, 2);
        coef_patch8 <<  width/2, widthExpand/2,
                        length - width/2, widthExpand/2,
                        width/2, width/2 + widthExpand,
                        length - width/2, width/2 + widthExpand;
        gsTensorBSpline<2, T> patch8(basis_lin, coef_patch8);
        patch8.degreeElevate(degree - 1, -1);



        //--------------------------------patch 9-------------------------------------------
        gsMatrix<T> coef_patch9(4, 2);
        coef_patch9 <<  width/2, -width/2,
                        length - width/2, -width/2,
                        width/2, widthExpand/2,
                        length - width/2, widthExpand/2;
        gsTensorBSpline<2, T> patch9(basis_lin, coef_patch9);
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
            gsMatrix<T> coef_patch10(4, 2);
            coef_patch10 << - prepatchWidth - width/2, widthExpand/2,
                    - width/2, widthExpand/2,
                    - prepatchWidth - width/2, width/2 + widthExpand,
                     - width/2, width/2 + widthExpand;

            gsTensorBSpline<2, T> patch10(basis_lin, coef_patch10);
            patch10.degreeElevate(degree - 1, -1);



            //--------------------------------patch 11-------------------------------------------
            gsMatrix<T> coef_patch11(4, 2);
            coef_patch11 << - prepatchWidth - width/2, -width/2,
                    - width/2,-width/2,
                    - prepatchWidth - width/2, widthExpand/2,
                     - width/2, widthExpand/2;

            gsTensorBSpline<2, T> patch11(basis_lin, coef_patch11);
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

public:
    gsMultiPatch<T> DomainBetweenBladeProfilesUniformKnotVector(
        T const & index,
        T const & length_x1,
        T const & length_x2,
        T const & pitch,
        T const & camberX,
        T const & camberY,
        T const & leadingAngle,
        T const & trailingAngle,
        T const & thicknessX,
        T const & thicknessY,
        T const & endingOffset,
        T const & outputAngle,
        T const & radius,
        T const & chordLength,
        T const & Angle,
        T const & rotationCenterX,
        T const & rotationCenterY,
        T const & uniformity_param)
    {


        //----------------set parameters for blade profile----------------
        bool plot = false;
        int num_samples = 30;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
            chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        //gsInfo << vec << "\n \n";
        //gsInfo << mat << "\n \n";
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        gsBSpline<T> suction_side_curve_transf;
        gsBSpline<T> pressure_side_curve_transf;
        std::vector<T> kvfit_knots{0,0,0,0,0.2,0.4,0.6,0.8,1,1,1,1};
        gsKnotVector<T> kvfit(kvfit_knots);
        unsigned num_cpblade = 8;
        BladeProfile<T> * pBladeProfile = 0;
        //unsigned degree = 3;
        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
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
        gsBSpline < T > bs = pBladeProfile->getPressureSide();
        gsBSpline < T > bp = pBladeProfile->getSuctionSide();
        gsMatrix < T > cp_bp(num_cpblade, 2);
        gsMatrix < T > cp_bp_pom(num_cpblade, 2);
        gsMatrix < T > cp_bs(num_cpblade, 2);

        //---------------set parameters for boundary of patches-----------------------
        //T width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

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

        gsMultiPatch<T> mp;
        // compute discrete coons patch to optimize
        gsMatrix<T> coef_patchAll(56, 2);
        coef_patchAll.setZero(56, 2);
        gsMatrix<T> a_cp(14, 2);
        gsMatrix<T> b_cp(14, 2);
        gsMatrix<T> c_cp(4, 2);
        gsMatrix<T> d_cp(4, 2);

        T ystart_coor = -((cp_bs(0, 1)*cp_bs(7, 0) - cp_bs(0, 0)*cp_bs(7, 1) - cp_bs(0, 1)*length_x1 + cp_bs(7, 1)*length_x1) / (
            cp_bs(0, 0) - cp_bs(7, 0)));

        gsMatrix<T> coef_patchStart(4, 2);
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

        T yend_coor = -((cp_bs(0, 1)*cp_bs(7, 0) - cp_bs(0, 0)*cp_bs(7, 1) - cp_bs(0, 1)*length_x2 + cp_bs(7, 1)*length_x2) / (
            cp_bs(0, 0) - cp_bs(7, 0)));
        gsMatrix<T> coef_patchEnd(4, 2);
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
        gsKnotVector<T> kv_uu(0, 1, 4, 4);
        gsKnotVector<T> kv_vv(0, 1, 0, 4);
        kv_uu.insert(0.2 / 3.0, 3);
        kv_uu.insert(1 - (0.2 / 3.0), 3);
        gsInfo << kv_uu;
        gsMultiPatch<T> * boundaries4 = new gsMultiPatch<T>;
        boundaries4->addPatch(gsBSpline<T>(kv_vv, c_cp));
        boundaries4->addPatch(gsBSpline<T>(kv_uu, a_cp));
        boundaries4->addPatch(gsBSpline<T>(kv_vv, d_cp));
        boundaries4->addPatch(gsBSpline<T>(kv_uu, b_cp));
        gsCoonsPatch<T> patchAll = coonsPatch(*boundaries4);
        patchAll.compute();
        mp.addPatch(patchAll.result());

        //=================================optimization===========================================

        gsVector<T> area_vec(7);
        area_vec.setZero(7);
        area_vec << 0.1, 0.25, 0.5, 0.5, 0.5, 0.75, 1;//1,1,1,1,0.75,0.75,0.75;

        T orthogonality = 0.0;
        T skewness = 0.0;
        T eccentricity = 0.0;
        T intersection = 0.0;
        T uniformity = uniformity_param;//0.01;// 0.25;
        //T uniformity = 0.005;// max seTing for RB rotation
        T area = area_vec(index);
        T length = 0;
        T epsilon = 1e-7;

        gsQualityMeasure<T> optimization(mp.patch(0));
        //    T opt_val = optimization.functional(orthogonality, skewness,
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
            gsMesh<T> mesh;
            makeMesh(mp.patch(0).basis(), mesh, 10);
            mp.patch(0).evaluateMesh(mesh);
            out = "optimize_bladeMesh" + util::to_string(index);
            gsWriteParaview(mesh, out);
            gsMesh<T> mesh2;
            mp.patch(0).controlNet(mesh2);
            out = "optimize_bladeControlNet" + util::to_string(index);
            gsWriteParaview(mesh2, out);
        }

        //=================================divide gsGeometry into three patches===========================================

        //initial data for patches
        gsMatrix<T> coefs = mp.patch(0).coefs();
        gsMultiPatch<T> mpFinal;
        gsKnotVector<T> kv_u(0, 1, 0, 4);
        gsKnotVector<T> kv_v(0, 1, 0, 4);
        gsTensorBSplineBasis<2, T> basis(kv_u, kv_v);
        gsTensorBSplineBasis<2, T> basis_blade(kvfit, kv_v);

        //--------------------------------patch 0-------------------------------------------
        gsMatrix<T> coef_patch0(16, 2);
        coef_patch0.setZero(16, 2);

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                coef_patch0(i * 4 + j, 0) = coefs(i * 4 + j + i * 10, 0);
                coef_patch0(i * 4 + j, 1) = coefs(i * 4 + j + i * 10, 1);
            }
        }
        gsTensorBSpline<2, T> patch0(basis, coef_patch0);
        for (T knot = 0.2; knot < 1.0; knot += 0.2)
        {
            patch0.insertKnot(knot, 0);
            patch0.insertKnot(knot, 1);
        }

        //--------------------------------patch 1-------------------------------------------
        gsMatrix<T> coef_patch1(num_cpblade * 4, 2);
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
        gsTensorBSpline<2, T> patch1(basis_blade, coef_patch1);
        for (T knot = 0.2; knot < 1.0; knot += 0.2)
        {
            patch1.insertKnot(knot, 1);
        }

        //--------------------------------patch 2-------------------------------------------
        gsMatrix<T> coef_patch2(16, 2);
        coef_patch2.setZero(16, 2);

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                coef_patch2(i * 4 + j, 0) = coefs(i * 4 + j + 10 + i * 10, 0);
                coef_patch2(i * 4 + j, 1) = coefs(i * 4 + j + 10 + i * 10, 1);
            }
        }
        gsTensorBSpline<2, T> patch2(basis, coef_patch2);
        for (T knot = 0.2; knot < 1.0; knot += 0.2)
        {
            patch2.insertKnot(knot, 0);
            patch2.insertKnot(knot, 1);
        }



        mpFinal.addPatch(patch0);
        mpFinal.addPatch(patch1);
        mpFinal.addPatch(patch2);

        gsInfo << "patches computed.\n";
          gsWriteParaview( mpFinal, "patches", 50000, true);

          mpFinal.patch(0);

        mpFinal.addInterface(0, boundary::south,(size_t) 0, boundary::north);
        mpFinal.addInterface(0, boundary::east, 1, boundary::west);

        mpFinal.addInterface(1, boundary::east, 2, boundary::west);
        mpFinal.addInterface(2, boundary::north, 2, boundary::south);

        mpFinal.addAutoBoundaries();

        return mpFinal;
    }


    gsMultiPatch<T> DomainBetweenBladeProfiles1(T const & index,
         T const & length_x1,
         T const & length_x2,
         T const & pitch,
         T const & camberX,
         T const & camberY,
         T const & leadingAngle,
         T const & trailingAngle,
         T const & thicknessX,
         T const & thicknessY,
         T const & endingOffset,
         T const & outputAngle,
         T const & radius,
         T const & chordLength,
         T const & Angle,
         T const & rotationCenterX,
         T const & rotationCenterY,
         T const & uniformity_param,
         std::vector<T> const & kvfit_knots,
         bool const & coarse,
         gsVector<T> const & geom_Params)
     {


         //----------------set parameters for blade profile----------------
         //bool plot = true;
         int num_samples = 30;
         gsVector<T> vec(2);
         //gsInfo << pitch << "\n ";
         vec(0) = rotationCenterX;
         vec(1) = rotationCenterY;
         gsMatrix<T> mat(2, 2);
         mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
             chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
         //gsInfo << vec << "\n \n";
         //gsInfo << mat << "\n \n";
         gsBSpline<T> suction_side_curve;
         gsBSpline<T> pressure_side_curve;
         gsBSpline<T> suction_side_curve_transf;
         gsBSpline<T> pressure_side_curve_transf;
         //gsKnotVector<T> kvfit(0, 1, 4, 4);
         //std::vector<T> kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};
         gsKnotVector<T> kvfit(kvfit_knots);

         gsKnotVector<T> kvlinear(0,1,0, 2);
         gsKnotVector<T> kvcubic(0,1,0, 4);
           gsKnotVector<T> kvuniform(0,1,0, kvfit.degree()+1);
         unsigned num_cpblade = kvfit.size() - kvfit.degree() - 1;
         BladeProfile<T> * pBladeProfile = 0;
         //unsigned degree = 3;
         //---------------compute blade profile for given parameters----------------------
         pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
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
         gsBSpline < T > bs = pBladeProfile->getPressureSide();
         gsBSpline < T > bp = pBladeProfile->getSuctionSide();
         gsMatrix < T > cp_bp(num_cpblade, 2);
         gsMatrix < T > cp_bp_pom(num_cpblade, 2);
         gsMatrix < T > cp_bs(num_cpblade, 2);

         //---------------set parameters for boundary of patches-----------------------
         //T width_y2 = /*bp.coef(bp.coefsSize()-1, 1)/2.0*/ + pitch/2.0;

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

         gsBSpline < T > bs_spline(kvfit,cp_bs);
         gsBSpline < T > bp_spline(kvfit,cp_bp);


         int cpk = num_cpblade - 1;//last index of cp profile


         T ystart_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x1 + cp_bs(cpk, 1)*length_x1) / (
             cp_bs(0, 0) - cp_bs(cpk, 0)));

         gsMatrix<T> coef_patchStart(2, 2);
         coef_patchStart << length_x1, ystart_coor,
                     cp_bs(0, 0), cp_bs(0, 1);

         gsBSpline<T> patchStart(kvlinear, coef_patchStart);
         patchStart.degreeElevate(kvfit.degree()-1);

         gsMatrix<T> coef_patchStart_pitch(2, 2);
         coef_patchStart_pitch << length_x1, ystart_coor + pitch,
                     cp_bp(0, 0), cp_bp(0, 1);

         gsBSpline<T> patchStart_pitch(kvlinear, coef_patchStart_pitch);
         patchStart_pitch.degreeElevate(kvfit.degree()-1);

         T yend_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x2 + cp_bs(cpk, 1)*length_x2) / (
             cp_bs(0, 0) - cp_bs(cpk, 0)));
         gsMatrix<T> coef_patchEnd(2, 2);
         coef_patchEnd << cp_bs(cpk, 0), cp_bs(cpk, 1),
                     length_x2, yend_coor;

          gsInfo << "\n endpatch \n";

         gsBSpline<T> patchEnd(kvlinear, coef_patchEnd);
         patchEnd.degreeElevate(kvfit.degree()-1);

         gsMatrix<T> coef_patchEnd_pitch(2, 2);
         coef_patchEnd_pitch << cp_bp(cpk, 0), cp_bp(cpk, 1),
                     length_x2, yend_coor+pitch;

          gsInfo << "\n endpatch \n";

         gsBSpline<T> patchEnd_pitch(kvlinear, coef_patchEnd_pitch);
         patchEnd_pitch.degreeElevate(kvfit.degree()-1);



 gsMatrix<T> c_cp(2, 2);
         c_cp << length_x1, ystart_coor,
             length_x1, ystart_coor + pitch;
         gsBSpline<T> c_spline(kvlinear, c_cp);
         c_spline.degreeElevate(kvfit.degree()-1);

 gsMatrix<T> d_cp(2, 2);
         d_cp << length_x2, yend_coor,
               length_x2, yend_coor + pitch;
         gsBSpline<T> d_spline(kvlinear, d_cp);
         d_spline.degreeElevate(kvfit.degree()-1);


         gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
         gsMatrix<T> par(1,1);

         gsMatrix<T> leading_coef(2, 2);

         p0 << cp_bs(0,0),cp_bs(0,1);
         par << 0.0;
         dir1 << bs_spline.deriv(par).transpose();
        par << 1.0;
         dir2 = -1*(patchStart.deriv(par).transpose());
         t0 = AxisDirectionNormed(dir1, dir2);
         p3 << cp_bp(0,0),cp_bp(0,1);
         par << 0.0;
         dir1 << bp_spline.deriv(par).transpose();
         dir2 = -1*(patchStart_pitch.deriv(par).transpose());

         t1 = AxisDirectionNormed(dir1, dir2);

         leading_coef =  LSQFergusonShort(p0, t0, t1, p3);

          gsBSpline<T> leading_spline_cubic(kvcubic, leading_coef);
          gsBSpline<T> leading_spline;

          gsMatrix<T> trailing_coef(2, 2);

          p0 << cp_bs(cpk,0),cp_bs(cpk,1);
          par << 1.0;
          dir1 << -1*(bs_spline.deriv(par).transpose());
         par << 0.0;
          dir2 = patchEnd.deriv(par).transpose();
          t0 = -1*AxisDirectionNormed(dir1, dir2);
          p3 << cp_bp(cpk,0),cp_bp(cpk,1);
          par << 1.0;
          dir1 <<-1*(bp_spline.deriv(par).transpose());
         par << 0.0;
          dir2 = patchEnd_pitch.deriv(par).transpose();
          t1 = AxisDirectionNormed(dir1, dir2);

          trailing_coef =  LSQFergusonShort(p0, t0, t1, p3);

           gsBSpline<T> trailing_spline_cubic(kvcubic, trailing_coef);
           gsBSpline<T> trailing_spline;

           gsKnotVector<T> kvaprox;

        if (kvfit.degree() == 3)
        {
            leading_spline = leading_spline_cubic;
            trailing_spline = trailing_spline_cubic;
            kvaprox = kvuniform;
        }
        else
        {
            if (kvfit.degree() > 3)
            {
                leading_spline = leading_spline_cubic;
                trailing_spline = trailing_spline_cubic;
                leading_spline.degreeElevate(kvfit.degree()-3);
                trailing_spline.degreeElevate(kvfit.degree()-3);
                kvaprox = kvuniform;
            }
            else
            {
                kvaprox = kvuniform;
                kvaprox.insert(0.25);
                kvaprox.insert(0.5);
                kvaprox.insert(0.75);

                int num_samples = 100;
                gsMatrix<T> parameter_points(1, num_samples+1);
                gsMatrix<T> leading_spline_points(2, num_samples+1);
                gsMatrix<T> trailing_spline_points(2, num_samples+1);
                gsMatrix<T> leading_spline_parameter_points(num_samples+1,1);
                gsMatrix<T> trailing_spline_parameter_points(num_samples+1,1);

                for (int i = 0; i < num_samples; i++) {
                    parameter_points(0,i) = (1-cos(pi/2*i/num_samples));
                }
                parameter_points(0,num_samples) = 1;

                leading_spline_points = leading_spline_cubic.eval(parameter_points);
                trailing_spline_points = trailing_spline_cubic.eval(parameter_points);

                leading_spline_parameter_points = centripetalParameterization(leading_spline_points);
                        trailing_spline_parameter_points = centripetalParameterization(trailing_spline_points);
               leading_spline = curveFittingWithBoundary(leading_spline_points,leading_spline_parameter_points,kvaprox,false);
               trailing_spline = curveFittingWithBoundary(trailing_spline_points,trailing_spline_parameter_points,kvaprox,false);

            }


        }


        gsInfo << "leading_spline" << leading_spline.coefs() << "\n";
         gsInfo << "trailing_spline" << trailing_spline.coefs() << "\n";



         gsTensorBSplineBasis<2, T> basis(kvuniform, kvaprox);
          gsTensorBSplineBasis<2, T> basis_blade(kvfit, kvaprox);


         //--------------------------------patch 0-------------------------------------------
         gsMatrix<T> coef_patch0((kvuniform.degree()+1)*(kvaprox.size() - kvaprox.degree() - 1), 2);
         coef_patch0.setZero((kvuniform.degree()+1)*(kvaprox.size() - kvaprox.degree() - 1), 2);

         discreteCoonsPatch(patchStart.coefs(),patchStart_pitch.coefs(),c_spline.coefs(),leading_spline.coefs(),coef_patch0,true);
         gsTensorBSpline<2, T> patch0(basis, coef_patch0);

         //--------------------------------patch 2-------------------------------------------
         gsMatrix<T> coef_patch2((kvuniform.degree()+1)*(kvaprox.size() - kvaprox.degree() - 1), 2);
         coef_patch2.setZero((kvuniform.degree()+1)*(kvaprox.size() - kvaprox.degree() - 1), 2);

         discreteCoonsPatch(patchEnd.coefs(),patchEnd_pitch.coefs(),trailing_spline.coefs(),d_spline.coefs(),coef_patch2,true);
         gsTensorBSpline<2, T> patch2(basis, coef_patch2);

         //--------------------------------patch 1-------------------------------------------
         gsMatrix<T> coef_patch1(num_cpblade*(kvaprox.size() - kvaprox.degree() - 1), 2);
         coef_patch0.setZero((num_cpblade)*(kvaprox.size() - kvaprox.degree() - 1), 2);

         discreteCoonsPatch(bs_spline.coefs(),bp_spline.coefs(),leading_spline.coefs(),trailing_spline.coefs(),coef_patch1,true);
         gsTensorBSpline<2, T> patch1(basis_blade, coef_patch1);




         if(!coarse)
        {
            gsInfo << "Special refinement. \n";

            gsInfo << patch1.knots(0) << "\n";

            std::vector<real_t> patch1_kv0_unique = patch1.knots(0).unique();


            //std::vector<real_t> patch1_inserted_knots = smartKnotIdentification2(patch1_kv0_unique);
            std::vector<real_t> patch1_inserted_knots = smartKnotIdentification(patch1_kv0_unique);

            for (size_t i = 0; i < patch1_inserted_knots.size(); i++) {
                if ((patch1_inserted_knots[i] > 0.0) && (patch1_inserted_knots[i] < 1.0))
                patch1.insertKnot(patch1_inserted_knots[i], 0);
            }

            for (T knot = 0.1; knot < 0.9; knot += 0.1)
            {
                gsInfo << knot << "\n";
                patch0.insertKnot(knot, 0);
                patch2.insertKnot(knot, 0);
            }

            for (T knot = 0.2; knot < 1.0; knot += 0.2)

            {
                gsInfo << knot << "\n";
                patch0.insertKnot(knot, 1);
                patch2.insertKnot(knot, 1);
                patch1.insertKnot(knot, 1);
            }
/*



            patch0.insertKnot(0.95, 0);
            patch2.insertKnot(0.05, 0);

            patch0.insertKnot(0.975, 0);
            patch2.insertKnot(0.025, 0);





            patch0.insertKnot(0.9, 1);
            patch2.insertKnot(0.9, 1);
            patch1.insertKnot(0.9, 1);
            patch0.insertKnot(0.1, 1);
            patch2.insertKnot(0.1, 1);
            patch1.insertKnot(0.1, 1);
            patch0.insertKnot(0.95, 1);
            patch2.insertKnot(0.95, 1);
            patch1.insertKnot(0.95, 1);
            patch0.insertKnot(0.05, 1);
            patch2.insertKnot(0.05, 1);
            patch1.insertKnot(0.05, 1);
            patch0.insertKnot(0.975, 1);
            patch2.insertKnot(0.975, 1);
            patch1.insertKnot(0.975, 1);
            patch0.insertKnot(0.025, 1);
            patch2.insertKnot(0.025, 1);
            patch1.insertKnot(0.025, 1);*/


            patch0.insertKnot(0.9, 1);
            patch2.insertKnot(0.9, 1);
            patch1.insertKnot(0.9, 1);
            patch0.insertKnot(0.1, 1);
            patch2.insertKnot(0.1, 1);
            patch1.insertKnot(0.1, 1);
            patch0.insertKnot(0.95, 1);
            patch2.insertKnot(0.95, 1);
            patch1.insertKnot(0.95, 1);
            patch0.insertKnot(0.05, 1);
            patch2.insertKnot(0.05, 1);
            patch1.insertKnot(0.05, 1);
            patch0.insertKnot(0.975, 1);
            patch2.insertKnot(0.975, 1);
            patch1.insertKnot(0.975, 1);
            patch0.insertKnot(0.025, 1);
            patch2.insertKnot(0.025, 1);
            patch1.insertKnot(0.025, 1);
        }

           gsMultiPatch<T> mpFinal;

         mpFinal.addPatch(patch0);
         mpFinal.addPatch(patch1);
         mpFinal.addPatch(patch2);

         mpFinal.addInterface(0, boundary::east, 1, boundary::west);
         mpFinal.addInterface(1, boundary::east, 2, boundary::west);
             //periodic
        mpFinal.addInterface(0, boundary::south,(size_t) 0, boundary::north);
        mpFinal.addInterface(2, boundary::north, 2, boundary::south);
        mpFinal.addAutoBoundaries();

         return mpFinal;
     }




    gsMultiPatch<T> DomainAroundBladeProfile2(T const & index,
                                              T const & length_x1,
                                              T const & length_x2,
                                              T const & pitch,
                                              T const & camberX,
                                              T const & camberY,
                                              T const & leadingAngle,
                                              T const & trailingAngle,
                                              T const & thicknessX,
                                              T const & thicknessY,
                                              T const & endingOffset,
                                              T const & outputAngle,
                                              T const & radius,
                                              T const & chordLength,
                                              T const & Angle,
                                              T const & rotationCenterX,
                                              T const & rotationCenterY,
                                              std::vector<T> const & kvfit_knots,
                                              bool const & coarse,
                                              gsVector<T> const & geom_Params) {

        //----------------set parameters for blade profile----------------
        T offset_distance = geom_Params(0);
        T inserted_knot = geom_Params(1);
        T inserted_knot_left_1 = geom_Params(2);
        T inserted_knot_left_2 = geom_Params(3);
        T inserted_knot_top_1 = geom_Params(4);
        T inserted_knot_top_2 = geom_Params(5);
        T inserted_knot_right = geom_Params(6);
        T fb = 0.5 + camberY + thicknessY;
        T ft = 1 - fb;
        //gsKnotVector<T> kvfit(0, 1, 4, 4);
        gsKnotVector<T> kvfit = gsKnotVector<T> (kvfit_knots);
        gsInfo << "kvfit:" << kvfit << "\n";
        gsKnotVector<T> kvcub(0, 1, 0, 4);
        gsKnotVector<T> kvlin(0, 1, 0, 2);

        bool plot = false;
        //bool plotMeshes = true;
        int num_samples = 100;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
               chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        gsBSpline<T> suction_side_offset_curve;
        gsBSpline<T> pressure_side_offset_curve;
        //unsigned num_cpblade = 8;
        BladeProfile<T> * pBladeProfile = 0;

        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
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

        T nPointsBefore = 0;
        while ((suction_side_curve.knots())[nPointsBefore] < inserted_knot) { nPointsBefore++; }
        gsInfo << nPointsBefore << "\n";

        // -------------- knot vector for trimmed suction and pressure side curves and their trimmed offset curves --------------------------------
        gsKnotVector<T> kvfit2 = suction_side_curve.knots();
        kvfit2.trimLeft(nPointsBefore);
        kvfit2.affineTransformTo(0.0, 1.0);
        gsInfo << kvfit2 << "\n";

        // -------------- knot vector for leading curve --------------------------------------------------------------------------------------------
        gsKnotVector<T> kvfit3 = suction_side_curve.knots();
        //kvfit3.trimRight(nPointsBefore + kvfit3.multiplicity(inserted_knot));
        kvfit3.trimRight(kvfit3.size() - kvfit3.multiplicity(inserted_knot) - nPointsBefore);
        kvfit3.affineTransformTo(0.0, 1.0);
        gsInfo << kvfit3 << "\n";
        std::vector<T> knots_for_kvfit3(kvfit3.size() + kvfit3.degree() + (kvfit3.size() - 2 * kvfit3.degree() - 2));
        int j = 0;
        for (index_t i = kvfit3.size()-1; i > kvfit3.degree(); i--) { knots_for_kvfit3[j] = math::abs(1 - kvfit3[i])/2; j++; gsInfo << math::abs(1 - kvfit3[i])/2 << "\n"; }
        for (unsigned i = 1; i < kvfit3.size(); i++) { knots_for_kvfit3[j] = 0.5 + kvfit3[i]/2; j++; gsInfo << 0.5 + kvfit3[i]/2 << "\n"; }
        kvfit3 = gsKnotVector<T> (knots_for_kvfit3);
        gsInfo << kvfit3 << "\n";

        // -------------- trimmed suction and pressure side curves and their trimmed offsets -------------------------------------------------------
        gsInfo << suction_side_curve.coefs() << "\n\n";
        gsInfo << pressure_side_curve.coefs() << "\n\n";
        gsInfo << suction_side_curve.coefsSize() << "\n\n";
        gsMatrix<T> suction_side_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2); suction_side_trimmed_cp.setZero();
        gsMatrix<T> pressure_side_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
        gsMatrix<T> suction_side_offset_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
        gsMatrix<T> pressure_side_offset_trimmed_cp(kvfit2.size()-kvfit2.degree()-1, 2);
        for (unsigned i = nPointsBefore; i < suction_side_curve.coefsSize(); i++) {
            gsInfo << i << "\n";
            suction_side_trimmed_cp.row(i-nPointsBefore) = suction_side_curve.coef(i);
            //gsInfo << suction_side_trimmed_cp << "\n";
            pressure_side_trimmed_cp.row(i-nPointsBefore) = pressure_side_curve.coef(i);
            suction_side_offset_trimmed_cp.row(i-nPointsBefore) = suction_side_offset_curve.coef(i);
            pressure_side_offset_trimmed_cp.row(i-nPointsBefore) = pressure_side_offset_curve.coef(i);
        }
        gsInfo << suction_side_trimmed_cp << "\n\n";
        gsInfo << pressure_side_trimmed_cp << "\n";
        gsBSpline<T> suction_side_trimmed = gsBSpline<T>( kvfit2, suction_side_trimmed_cp);
        gsBSpline<T> pressure_side_trimmed = gsBSpline<T>( kvfit2, pressure_side_trimmed_cp);
        gsBSpline<T> suction_side_offset_trimmed = gsBSpline<T>( kvfit2, suction_side_offset_trimmed_cp);
        gsBSpline<T> pressure_side_offset_trimmed = gsBSpline<T>( kvfit2, pressure_side_offset_trimmed_cp);
        gsInfo << suction_side_trimmed << "\n";

        // -------------- leading curve -------------------------------------------------------------------------------------------------------------
        gsMatrix<T> leading_curve_cp(kvfit3.size()-kvfit3.degree()-1, 2);
        gsMatrix<T> leading_offset_curve_cp(kvfit3.size()-kvfit3.degree()-1, 2);
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
        gsBSpline<T> leading_curve = gsBSpline<T>( kvfit3, leading_curve_cp);
        gsBSpline<T> leading_offset_curve = gsBSpline<T>( kvfit3, leading_offset_curve_cp);

        // --------------- cross-section curves ---------------------------------------------------------------------------------------------------
        gsMatrix<T> aux_cp(2, 2);
        aux_cp << suction_side_trimmed_cp.row(0),
                  suction_side_offset_trimmed_cp.row(0);
        gsBSpline<T> cs_curve1 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve1.degreeElevate(2);
        gsMatrix<T> cs_curve1_cp = cs_curve1.coefs();
        aux_cp << pressure_side_trimmed_cp.row(0),
                  pressure_side_offset_trimmed_cp.row(0);
        gsBSpline<T> cs_curve2 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve2.degreeElevate(2);
        gsMatrix<T> cs_curve2_cp = cs_curve2.coefs();
        aux_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                  pressure_side_offset_trimmed_cp.row(pressure_side_offset_trimmed.coefsSize()-1);
        gsBSpline<T> cs_curve3 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve3.degreeElevate(2);
        gsMatrix<T> cs_curve3_cp = cs_curve3.coefs();
        aux_cp << suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1),
                  suction_side_offset_trimmed_cp.row(suction_side_offset_trimmed.coefsSize()-1);
        gsBSpline<T> cs_curve4 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve4.degreeElevate(2);
        gsMatrix<T> cs_curve4_cp = cs_curve4.coefs();

        // ---------------- outer boundary -------------------------------------------------------------------------------------------------------
        gsMatrix<T> cp_bs = suction_side_curve.coefs();
        T ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
        T yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

        gsMatrix<T> leftSplitPoint1(1,2), leftSplitPoint2(1,2);
        leftSplitPoint1(0,0) = length_x1;
        leftSplitPoint1(0,1) = (1 - inserted_knot_left_1) * (ystart_coor - fb*pitch) + inserted_knot_left_1 * (ystart_coor + ft*pitch);
        leftSplitPoint2(0,0) = length_x1;
        leftSplitPoint2(0,1) = (1 - inserted_knot_left_2) * (ystart_coor - fb*pitch) + inserted_knot_left_2 * (ystart_coor + ft*pitch);
        aux_cp << leftSplitPoint1(0,0), leftSplitPoint1(0,1),
                  length_x1, ystart_coor - fb*pitch;
        gsBSpline<T> left_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> left_boundary_curve_A_cp = left_boundary_curve_A.coefs();
        /*aux_cp << leftSplitPoint1(0,0), leftSplitPoint1(0,1),
                  leftSplitPoint2(0,0), leftSplitPoint2(0,1);
        gsBSpline<T> left_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve_B.degreeElevate(2);
        */
        gsMatrix<T> aux_cp2(4,2);
        aux_cp2 << leftSplitPoint1(0,0), leftSplitPoint1(0,1),
                   0.9 * leftSplitPoint1(0,0) + 0.1 * leftSplitPoint2(0,0), 0.9 * leftSplitPoint1(0,1) + 0.1 * leftSplitPoint2(0,1),
                   0.1 * leftSplitPoint1(0,0) + 0.9 * leftSplitPoint2(0,0), 0.1 * leftSplitPoint1(0,1) + 0.9 * leftSplitPoint2(0,1),
                   leftSplitPoint2(0,0), leftSplitPoint2(0,1);
        gsBSpline<T> left_boundary_curve_B = gsBSpline<T> ( kvcub, aux_cp2 );
        UnifyKnotVectors(left_boundary_curve_B, leading_offset_curve);
        gsMatrix<T> left_boundary_curve_B_cp = left_boundary_curve_B.coefs();
        aux_cp << leftSplitPoint2(0,0), leftSplitPoint2(0,1),
                  length_x1, ystart_coor + ft*pitch;
        gsBSpline<T> left_boundary_curve_C = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve_C.degreeElevate(2);
        //left_boundary_curve_C.insertKnot(0.333); left_boundary_curve_C.insertKnot(0.666);
        gsMatrix<T> left_boundary_curve_C_cp = left_boundary_curve_C.coefs();
        gsInfo << "left_boundary_curve_B_cp: " << left_boundary_curve_B_cp.rows() << ", " << left_boundary_curve_B_cp.cols() << "\n";

        /*aux_cp << length_x2, yend_coor + ft*pitch,
                  length_x2, yend_coor - fb*pitch;
        gsBSpline<T> right_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve.degreeElevate(2);*/
        gsMatrix<T> rightSplitPoint(1,2);
        rightSplitPoint(0,0) = length_x2;
        rightSplitPoint(0,1) = (1 - inserted_knot_right) * (yend_coor - fb*pitch) + inserted_knot_right * (yend_coor + ft*pitch);
        aux_cp << rightSplitPoint(0,0), rightSplitPoint(0,1),
                length_x2, yend_coor - fb*pitch;
        gsBSpline<T> right_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_A_cp = right_boundary_curve_A.coefs();
        aux_cp << length_x2, yend_coor + ft*pitch,
                  rightSplitPoint(0,0), rightSplitPoint(0,1);
        gsBSpline<T> right_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve_B.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_B_cp = right_boundary_curve_B.coefs();

        gsMatrix<T> topSplitPoint1(1,2), topSplitPoint2(1,2);
        topSplitPoint1(0,0) = (1 - inserted_knot_top_1) * length_x1 + inserted_knot_top_1 * length_x2;
        topSplitPoint1(0,1) = (1 - inserted_knot_top_1) * (ystart_coor + ft*pitch) + inserted_knot_top_1 * (yend_coor + ft*pitch);
        topSplitPoint2(0,0) = (1 - inserted_knot_top_2) * length_x1 + inserted_knot_top_2 * length_x2;
        topSplitPoint2(0,1) = (1 - inserted_knot_top_2) * (ystart_coor + ft*pitch) + inserted_knot_top_2 * (yend_coor + ft*pitch);
        aux_cp << topSplitPoint1(0,0), topSplitPoint1(0,1),
                  length_x1, ystart_coor + ft*pitch;
        gsBSpline<T> top_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
        aux_cp << topSplitPoint1(0,0), topSplitPoint1(0,1),
                  topSplitPoint2(0,0), topSplitPoint2(0,1);
        gsBSpline<T> top_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_B.degreeElevate(2);
        UnifyKnotVectors(top_boundary_curve_B, pressure_side_offset_trimmed);
        gsMatrix<T> top_boundary_curve_B_cp = top_boundary_curve_B.coefs();
        aux_cp << topSplitPoint2(0,0), topSplitPoint2(0,1),
                  length_x2, yend_coor + ft*pitch;
        gsBSpline<T> top_boundary_curve_C = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_C.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_C_cp = top_boundary_curve_C.coefs();

        gsMatrix<T> bottomSplitPoint1(1,2), bottomSplitPoint2(1,2);
        bottomSplitPoint1(0,0) = (1 - inserted_knot_top_1) * length_x1 + inserted_knot_top_1 * length_x2;
        bottomSplitPoint1(0,1) = (1 - inserted_knot_top_1) * (ystart_coor - fb*pitch) + inserted_knot_top_1 * (yend_coor - fb*pitch);
        bottomSplitPoint2(0,0) = (1 - inserted_knot_top_2) * length_x1 + inserted_knot_top_2 * length_x2;
        bottomSplitPoint2(0,1) = (1 - inserted_knot_top_2) * (ystart_coor - fb*pitch) + inserted_knot_top_2 * (yend_coor - fb*pitch);
        aux_cp << bottomSplitPoint1(0,0), bottomSplitPoint1(0,1),
                  length_x1, ystart_coor - fb*pitch;
        gsBSpline<T> bottom_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
        aux_cp << bottomSplitPoint1(0,0), bottomSplitPoint1(0,1),
                  bottomSplitPoint2(0,0), bottomSplitPoint2(0,1);
        gsBSpline<T> bottom_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_B.degreeElevate(2);
        UnifyKnotVectors(bottom_boundary_curve_B, suction_side_offset_trimmed);
        gsMatrix<T> bottom_boundary_curve_B_cp = bottom_boundary_curve_B.coefs();
        aux_cp << bottomSplitPoint2(0,0), bottomSplitPoint2(0,1),
                  length_x2, yend_coor - fb*pitch;
        gsBSpline<T> bottom_boundary_curve_C = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_C.degreeElevate(2);
        //bottom_boundary_curve_C.insertKnot(0.333); bottom_boundary_curve_C.insertKnot(0.666);
        gsMatrix<T> bottom_boundary_curve_C_cp = bottom_boundary_curve_C.coefs();

        // ---------------- cross-section curves 2 -------------------------------------------------------------------------------------------------
        gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
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
        gsMatrix<T> cs_curve5_cp(4,2);
        cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve5_cp << "\n";
        gsBSpline<T> cs_curve5 = gsBSpline<T> ( kvcub, cs_curve5_cp );

        p0 << pressure_side_offset_trimmed_cp.row(0);
        dir1 << pressure_side_offset_trimmed_cp.row(0) - pressure_side_trimmed_cp.row(0);
        dir2 << topSplitPoint1 - pressure_side_offset_trimmed_cp.row(0);
        t0 = AxisDirectionNormed(dir1, dir2);
        dir1 << yend_coor - ystart_coor, - length_x2 + length_x1;
        dir2 << suction_side_offset_trimmed_cp.row(0) - topSplitPoint1;
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << topSplitPoint1(0,0), topSplitPoint1(0,1);
        gsMatrix<T> cs_curve6_cp(4,2);
        cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve6_cp << "\n";
        gsBSpline<T> cs_curve6 = gsBSpline<T> ( kvcub, cs_curve6_cp );
        //cs_curve6.insertKnot(0.333); cs_curve6.insertKnot(0.666);
        cs_curve6_cp = cs_curve6.coefs();

        /*p0 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
        //t0 << (pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1))/(pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1)).norm();
        dir1 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
        dir2 << topSplitPoint2 - pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
        t0 = AxisDirectionNormed(dir1, dir2);
        dir1 << yend_coor - ystart_coor, - length_x2 + length_x1;
        dir2 << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - topSplitPoint2;
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << topSplitPoint2;
        gsMatrix<T> cs_curve7_cp(4,2);
        cs_curve7_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve7_cp << "\n";
        gsBSpline<T> cs_curve7 = gsBSpline<T> ( kvcub, cs_curve7_cp );
        //cs_curve7.insertKnot(0.333); cs_curve7.insertKnot(0.666);
        cs_curve7_cp = cs_curve7.coefs();*/
        aux_cp << pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                  topSplitPoint2;
        gsBSpline<T> cs_curve7 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve7.degreeElevate(2);
        gsMatrix<T> cs_curve7_cp = cs_curve7.coefs();

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
        gsMatrix<T> cs_curve8_cp(4,2);
        cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve8_cp << "\n";
        gsBSpline<T> cs_curve8 = gsBSpline<T> ( kvcub, cs_curve8_cp );

        /*p0 << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1);
        //dir1 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
        //dir2 << length_x2 - length_x1, yend_coor - ystart_coor;
        //t0 = AxisDirectionNormed(dir1, dir2);
        t0 << length_x2 - suction_side_offset_trimmed_cp(suction_side_trimmed.coefsSize()-1, 0), yend_coor - fb*pitch - suction_side_offset_trimmed_cp(suction_side_trimmed.coefsSize()-1, 1);
        t1 << length_x1 - length_x2, ystart_coor - yend_coor; t1 = t1/t1.norm();
        p3 << rightSplitPoint;
        gsMatrix<T> cs_curve9_cp(4,2);
        cs_curve9_cp = LSQFergusonSmooth(p0, t0, t1, p3);
        gsInfo << cs_curve9_cp << "\n";
        gsBSpline<T> cs_curve9 = gsBSpline<T> ( kvcub, cs_curve9_cp );
        //cs_curve9.insertKnot(0.333); cs_curve9.insertKnot(0.666);
        cs_curve9_cp = cs_curve9.coefs();
        */
        aux_cp << suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1),
                  rightSplitPoint;
        gsBSpline<T> cs_curve9 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve9.degreeElevate(2);
        gsMatrix<T> cs_curve9_cp = cs_curve9.coefs();

        p3 << leftSplitPoint1;
        //dir1 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
        //dir2 << length_x2 - length_x1, yend_coor - ystart_coor;
        //t0 = AxisDirectionNormed(dir1, dir2);
        t1 << length_x2 - length_x1, yend_coor - ystart_coor; t1 = t1/t1.norm();
        gsMatrix<T> par(1,1); par << 0.0;
        dir1 = leading_offset_curve.deriv(par).transpose();
        //dir2 << (suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0))/(suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0)).norm();
        dir2 = cs_curve5.deriv(par).transpose();
        t0 = AxisDirectionNormed(dir1, dir2);
        p0 << suction_side_offset_trimmed_cp.row(0);
        gsMatrix<T> cs_curve10_cp(4,2);
        cs_curve10_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve10_cp << "\n";
        gsBSpline<T> cs_curve10 = gsBSpline<T> ( kvcub, cs_curve10_cp );

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
        gsMatrix<T> cs_curve11_cp(4,2);
        cs_curve11_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve11_cp << "\n";
        gsBSpline<T> cs_curve11 = gsBSpline<T> ( kvcub, cs_curve11_cp );

        // ---------------- circular patches behind the blade ----------------------------------------------------------------------------------------
        par << 1.0;
        dir1 << pressure_side_trimmed.deriv(par).transpose();
        dir2 << suction_side_trimmed.deriv(par).transpose();
        dir2 = AxisDirectionNormed(dir1, dir2);
        dir1 = (pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1))/(pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1)).norm();
        T alpha = acos(dir1(0,0) * dir2(0,0) + dir1(0,1) * dir2(0,1));
        gsInfo << "alpha = " << alpha << "\n";
        gsInfo << pressure_side_offset_trimmed_cp << "\n\n";
        gsInfo << pressure_side_trimmed_cp << "\n\n";
        gsInfo << suction_side_offset_trimmed_cp << "\n\n";
        gsMatrix<T> pressure_side_start(2,1);
        pressure_side_start = pressure_side_offset_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
        gsMatrix<T> circularArcPressureSide_cp = computeCircleArc2D(pressure_side_start, pressure_side_trimmed_cp(pressure_side_trimmed.coefsSize()-1,0), pressure_side_trimmed_cp(pressure_side_trimmed.coefsSize()-1,1), -alpha);
        gsInfo << "ctvrtkruh1 = " << circularArcPressureSide_cp << "\n";
        gsBSpline<T> circularArcPressureSide = gsBSpline<T> ( kvcub, circularArcPressureSide_cp );
        gsMatrix<T> circularArcPressureSide_start(2,1);
        circularArcPressureSide_start = circularArcPressureSide_cp.row(circularArcPressureSide.coefsSize()-1);
        gsMatrix<T> circularArcSuctionSide_cp = computeCircleArc2D(circularArcPressureSide_start, suction_side_trimmed_cp(suction_side_trimmed.coefsSize()-1,0), suction_side_trimmed_cp(suction_side_trimmed.coefsSize()-1,1), -alpha);
        gsInfo << "ctvrtkruh2 = " << circularArcSuctionSide_cp << "\n";
        gsBSpline<T> circularArcSuctionSide = gsBSpline<T> ( kvcub, circularArcSuctionSide_cp );

        aux_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                  circularArcSuctionSide_cp.row(0);
        gsBSpline<T> cs_curve12 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve12.degreeElevate(2);
        gsMatrix<T> cs_curve12_cp = cs_curve12.coefs();
        gsInfo << cs_curve12_cp << "\n";

        gsMatrix<T> endpoint_curve_cp(4,2);
        endpoint_curve_cp << pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                             pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                             pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1),
                             pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1);
        gsInfo << endpoint_curve_cp << "\n";

        gsMatrix<T> patch10_cp(endpoint_curve_cp.rows() * cs_curve3_cp.rows(), 2);
        discreteCoonsPatch(endpoint_curve_cp, circularArcPressureSide_cp, cs_curve3_cp, cs_curve12_cp, patch10_cp, true);
        gsTensorBSpline<2, T> patch10 = gsTensorBSpline<2, T> (kvcub, kvcub, patch10_cp);

        gsMatrix<T> patch11_cp(endpoint_curve_cp.rows() * cs_curve12_cp.rows(), 2);
        discreteCoonsPatch(endpoint_curve_cp, circularArcSuctionSide_cp, cs_curve12_cp, cs_curve4_cp, patch11_cp, true);
        gsTensorBSpline<2, T> patch11 = gsTensorBSpline<2, T> (kvcub, kvcub, patch11_cp);

        // ---------------- cross-section curves 3 -------------------------------------------------------------------------------------------------
        p0 << circularArcPressureSide_cp.row(circularArcPressureSide.coefsSize()-1);
        //t0 << circularArcSuctionSide_cp.row(0) - pressure_side_trimmed_cp.row(pressure_side_trimmed.coefsSize()-1); t0 = t0/t0.norm();
        t0 << length_x2 - length_x1, yend_coor - ystart_coor; t0 = t0/t0.norm();
        dir1 << length_x1-length_x2, ystart_coor - yend_coor;
        dir2 << 0, -1;
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << length_x2, yend_coor + ft*pitch;
        gsMatrix<T> cs_curve13_cp(4,2);
        cs_curve13_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve13_cp << "\n";
        gsBSpline<T> cs_curve13 = gsBSpline<T> ( kvcub, cs_curve13_cp );
        //cs_curve13.insertKnot(0.333); cs_curve13.insertKnot(0.666);
        cs_curve13_cp = cs_curve13.coefs();

        // ---------------- construction of patches ------------------------------------------------------------------------------------------------
        //gsMultiPatch<T> *boundaries = new gsMultiPatch<T>;
        /*boundaries->addPatch(cs_curve1);
        boundaries->addPatch(suction_side_trimmed);
        boundaries->addPatch(cs_curve4);
        boundaries->addPatch(suction_side_offset_trimmed);
        gsCoonsPatch<T> patch1 = coonsPatch(*boundaries);
        patch1.compute();
        delete boundaries;
        */

        gsMatrix<T> patch1_cp(suction_side_trimmed_cp.rows() * cs_curve1_cp.rows(), 2);
        //discreteCoonsPatch(suction_side_trimmed_cp, suction_side_offset_trimmed_cp, cs_curve1_cp, cs_curve4_cp, patch1_cp, true);
        discreteCoonsPatch(cs_curve1_cp, cs_curve4_cp, suction_side_trimmed_cp, suction_side_offset_trimmed_cp, patch1_cp, true);
        //gsTensorBSpline<2, T> patch1 = gsTensorBSpline<2, T> (kvfit2, kvcub, patch1_cp);
        gsTensorBSpline<2, T> patch1 = gsTensorBSpline<2, T> (kvcub, kvfit2, patch1_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve1);
        boundaries->addPatch(leading_curve);
        boundaries->addPatch(cs_curve2);
        boundaries->addPatch(leading_offset_curve);
        gsCoonsPatch<T> patch2 = coonsPatch(*boundaries);
        patch2.compute();
        delete boundaries;*/

        gsMatrix<T> patch2_cp(leading_curve_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(leading_curve_cp, leading_offset_curve_cp, cs_curve1_cp, cs_curve2_cp, patch2_cp, true);
        gsTensorBSpline<2, T> patch2 = gsTensorBSpline<2, T> (kvfit3, kvcub, patch2_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve2);
        boundaries->addPatch(pressure_side_trimmed);
        boundaries->addPatch(cs_curve3);
        boundaries->addPatch(pressure_side_offset_trimmed);
        gsCoonsPatch<T> patch3 = coonsPatch(*boundaries);
        patch3.compute();
        delete boundaries;*/

        gsMatrix<T> patch3_cp(pressure_side_trimmed_cp.rows() * cs_curve2_cp.rows(), 2);
        discreteCoonsPatch(pressure_side_trimmed_cp, pressure_side_offset_trimmed_cp, cs_curve2_cp, cs_curve3_cp, patch3_cp, true);
        gsTensorBSpline<2, T> patch3 = gsTensorBSpline<2, T> (kvfit2, kvcub, patch3_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve5);
        boundaries->addPatch(suction_side_offset_trimmed);
        boundaries->addPatch(cs_curve8);
        boundaries->addPatch(bottom_boundary_curve_A);
        gsCoonsPatch<T> patch4 = coonsPatch(*boundaries);
        patch4.compute();
        delete boundaries;*/

        gsMatrix<T> patch4_cp(suction_side_offset_trimmed_cp.rows() * cs_curve5_cp.rows(), 2);
        gsMatrix<T> patch4_cp1(suction_side_offset_trimmed_cp.rows() * cs_curve5_cp.rows(), 2);
        gsMatrix<T> patch4_cp2(suction_side_offset_trimmed_cp.rows() * cs_curve5_cp.rows(), 2);
        springModelPatch(cs_curve5_cp, cs_curve8_cp, suction_side_offset_trimmed_cp, bottom_boundary_curve_B_cp, patch4_cp1, true);
        discreteCoonsPatch(cs_curve5_cp, cs_curve8_cp, suction_side_offset_trimmed_cp, bottom_boundary_curve_B_cp, patch4_cp2, true);
        patch4_cp = (patch4_cp1 + patch4_cp2) / 2;
        //gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvfit2, kvcub, patch4_cp);
        gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvcub, kvfit2, patch4_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve5);
        boundaries->addPatch(leading_offset_curve);
        boundaries->addPatch(cs_curve6);
        boundaries->addPatch(left_boundary_curve);
        gsCoonsPatch<T> patch5 = coonsPatch(*boundaries);
        patch5.compute();
        delete boundaries;*/

        gsMatrix<T> patch5_cp(cs_curve5_cp.rows() * cs_curve10_cp.rows(), 2);
        //springModelPatch(cs_curve5_cp, left_boundary_curve_A_cp, cs_curve10_cp, bottom_boundary_curve_A_cp, patch5_cp, true);
        springModelPatch(cs_curve10_cp, bottom_boundary_curve_A_cp, cs_curve5_cp, left_boundary_curve_A_cp, patch5_cp, true);
        gsTensorBSpline<2, T> patch5 = gsTensorBSpline<2, T> (kvcub, kvcub, patch5_cp);

        gsMatrix<T> patch6_cp(leading_offset_curve_cp.rows() * cs_curve10_cp.rows(), 2);
        discreteCoonsPatch(leading_offset_curve_cp, left_boundary_curve_B_cp, cs_curve10_cp, cs_curve11_cp, patch6_cp, true);
        gsTensorBSpline<2, T> patch6 = gsTensorBSpline<2, T> (kvfit3, kvcub, patch6_cp);

        gsMatrix<T> patch7_cp(cs_curve6_cp.rows() * cs_curve11_cp.rows(), 2);
        discreteCoonsPatch(cs_curve6_cp, left_boundary_curve_C_cp, cs_curve11_cp, top_boundary_curve_A_cp, patch7_cp, true);
        gsTensorBSpline<2, T> patch7 = gsTensorBSpline<2, T> (cs_curve6.knots(), cs_curve11.knots(), patch7_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve6);
        boundaries->addPatch(pressure_side_offset_trimmed);
        boundaries->addPatch(cs_curve7);
        boundaries->addPatch(top_boundary_curve_A);
        gsCoonsPatch<T> patch6 = coonsPatch(*boundaries);
        patch6.compute();
        delete boundaries;*/

        gsMatrix<T> patch8_cp(pressure_side_offset_trimmed_cp.rows() * cs_curve6_cp.rows(), 2);
        discreteCoonsPatch(pressure_side_offset_trimmed_cp, top_boundary_curve_B_cp, cs_curve6_cp, cs_curve7_cp, patch8_cp, true);
        gsTensorBSpline<2, T> patch8 = gsTensorBSpline<2, T> (pressure_side_offset_trimmed.knots(), cs_curve6.knots(), patch8_cp);

        /*boundaries = new gsMultiPatch<T>;
        boundaries->addPatch(cs_curve9);
        boundaries->addPatch(right_boundary_curve_A);
        boundaries->addPatch(bottom_boundary_curve_B);
        boundaries->addPatch(cs_curve8);
        gsCoonsPatch<T> patch7 = coonsPatch(*boundaries);
        patch7.compute();
        delete boundaries;*/

        gsMatrix<T> patch9_cp(cs_curve8_cp.rows() * cs_curve9_cp.rows(), 2);
        discreteCoonsPatch(cs_curve8_cp, right_boundary_curve_A_cp, cs_curve9_cp, bottom_boundary_curve_C_cp, patch9_cp, true);
        gsTensorBSpline<2, T> patch9 = gsTensorBSpline<2, T> (cs_curve8.knots(), cs_curve9.knots(), patch9_cp);

        gsMatrix<T> patch12_cp(circularArcPressureSide_cp.rows() * cs_curve7_cp.rows(), 2);
        discreteCoonsPatch(circularArcPressureSide_cp, top_boundary_curve_C_cp, cs_curve7_cp, cs_curve13_cp, patch12_cp, true);
        gsTensorBSpline<2, T> patch12 = gsTensorBSpline<2, T> (circularArcPressureSide.knots(), cs_curve7.knots(), patch12_cp);

        gsMatrix<T> patch13_cp(circularArcSuctionSide_cp.rows() * cs_curve13_cp.rows(), 2);
        springModelPatch(circularArcSuctionSide_cp, right_boundary_curve_B_cp, cs_curve13_cp, cs_curve9_cp, patch13_cp, true);
        gsTensorBSpline<2, T> patch13 = gsTensorBSpline<2, T> (circularArcSuctionSide.knots(), cs_curve13.knots(), patch13_cp);

        // ---------------- construction of multipatch ---------------------------------------------------------------------------------------------




         if(!coarse)
        {
             gsInfo << "Special refinement. \n";
             gsInfo << patch1.knots(0) << "\n";
             gsInfo << patch1.knots(1) << "\n";
             gsInfo << patch2.knots(0) << "\n";
             gsInfo << patch2.knots(1) << "\n";
             std::vector<real_t> patch1_kv0_unique = patch1.knots(1).unique();
             std::vector<real_t> patch2_kv0_unique_aux = patch2.knots(0).unique();
             std::vector<real_t> patch2_kv0_unique(patch2_kv0_unique_aux.size() / 2 + 1);
             for (size_t i = 0; i < patch2_kv0_unique_aux.size() / 2 + 1; i++)
                 patch2_kv0_unique[i] = patch2_kv0_unique_aux[patch2_kv0_unique_aux.size() / 2 + i];
             for (size_t i = 0; i < patch1_kv0_unique.size(); i++)
                 gsInfo << patch1_kv0_unique[i] << ",";
             gsInfo << "\n";
             for (size_t i = 0; i < patch2_kv0_unique.size(); i++)
                 gsInfo << patch2_kv0_unique[i] << ",";
             gsInfo << "\n";

             std::vector<real_t> patch1_inserted_knots = smartKnotIdentification(patch1_kv0_unique);
             std::vector<real_t> patch2_inserted_knots = smartKnotIdentification(patch2_kv0_unique);

             /*
             std::vector<real_t> patch1_kv0_spans(patch1_kv0_unique.size()-1), patch2_kv0_spans(patch2_kv0_unique.size()-1);
             for (index_t i = 0; i < patch1_kv0_unique.size()-1; i++)
                 patch1_kv0_spans[i] = patch1_kv0_unique[i+1] - patch1_kv0_unique[i];
             for (index_t i = 0; i < patch2_kv0_unique.size()-1; i++)
                 patch2_kv0_spans[i] = patch2_kv0_unique[i+1] - patch2_kv0_unique[i];
             for (index_t i = 0; i < patch1_kv0_spans.size(); i++)
                 gsInfo << patch1_kv0_spans[i] << ",";
             gsInfo << "\n";
             for (index_t i = 0; i < patch2_kv0_spans.size(); i++)
                 gsInfo << patch2_kv0_spans[i] << ",";
             gsInfo << "\n";

             bool changed = false;
             gsVector<real_t> fr1(2), fr2(3);
             fr1 << 0.35, 0.65;
             fr2 << 0.25, 0.33, 0.417;
             std::vector<real_t> patch1_inserted_knots;
             index_t iter = 0;
             auto it1 = patch1_kv0_spans.begin();
             do {
                 changed = false;
                 for (index_t i = 0; i < patch1_kv0_spans.size()-1; i++) {
                     gsInfo << i << "\n";
                     if (patch1_kv0_spans[i+1] <= fr1[1] * patch1_kv0_spans[i]) {
                         gsInfo << "a\n";
                         patch1_inserted_knots.push_back(patch1_kv0_unique[i] + fr2[0] * patch1_kv0_spans[i]);
                         patch1_inserted_knots.push_back(patch1_kv0_unique[i] + (fr2[0]+fr2[1]) * patch1_kv0_spans[i]);
                         it1 = patch1_kv0_spans.insert(it1 + i + 1, fr2[2] * patch1_kv0_spans[i]);
                         patch1_kv0_spans.insert(it1, fr2[1] * patch1_kv0_spans[i]);
                         patch1_kv0_spans[i] = fr2[0] * patch1_kv0_spans[i];
                         it1 = patch1_kv0_spans.begin();
                         for (index_t i = 0; i < patch1_kv0_spans.size(); i++)
                             gsInfo << patch1_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                     else if (patch1_kv0_spans[i+1] < (0.99 * patch1_kv0_spans[i])) {
                         gsInfo << "b\n";
                         patch1_inserted_knots.push_back(patch1_kv0_unique[i] + fr1[0] * patch1_kv0_spans[i]);
                         it1 = patch1_kv0_spans.insert(it1 + i + 1, fr1[1] * patch1_kv0_spans[i]);
                         patch1_kv0_spans[i] = fr1[1] * patch1_kv0_spans[i];
                         it1 = patch1_kv0_spans.begin();
                         for (index_t i = 0; i < patch1_kv0_spans.size(); i++)
                             gsInfo << patch1_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                 }
                 iter++;
                 gsInfo << "iter = " << iter << "\n";
             } while (changed);
             gsInfo << "patch1_inserted_knots: ";
             for (index_t i = 0; i < patch1_inserted_knots.size(); i++)
                 gsInfo << patch1_inserted_knots[i] << ",";
             gsInfo << "\n";
             do {
                 changed = false;
                 for (index_t i = 0; i < patch1_kv0_spans.size()-1; i++) {
                     gsInfo << i << "\n";
                     if ((3 * patch1_kv0_spans[i]) < patch1_kv0_spans[i+1]) {
                         gsInfo << "c\n";
                         patch1_inserted_knots.push_back(patch1_kv0_unique[i+1] + fr1[0] * patch1_kv0_spans[i+1]);
                         T span_old = patch1_kv0_spans[i+1];
                         patch1_kv0_spans[i+1] = fr1[1] * span_old;
                         it1 = patch1_kv0_spans.insert(it1 + i + 1, fr1[0] * span_old);
                         it1 = patch1_kv0_spans.begin();
                         for (index_t i = 0; i < patch1_kv0_spans.size(); i++)
                             gsInfo << patch1_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                 }
                 iter++;
                 gsInfo << "iter = " << iter << "\n";
             } while (changed);
             gsInfo << "patch1_inserted_knots: ";
             for (index_t i = 0; i < patch1_inserted_knots.size(); i++)
                 gsInfo << patch1_inserted_knots[i] << ",";
             gsInfo << "\n";

             std::vector<real_t> patch2_inserted_knots;
             iter = 0;
             auto it2 = patch2_kv0_spans.begin();
             do {
                 changed = false;
                 for (index_t i = 0; i < patch2_kv0_spans.size()-1; i++) {
                     gsInfo << i << "\n";
                     if (patch2_kv0_spans[i+1] <= fr1[1] * patch2_kv0_spans[i]) {
                         gsInfo << "a\n";
                         patch2_inserted_knots.push_back(patch2_kv0_unique[i] + fr2[0] * patch2_kv0_spans[i]);
                         patch2_inserted_knots.push_back(patch2_kv0_unique[i] + (fr2[0]+fr2[1]) * patch2_kv0_spans[i]);
                         it2 = patch2_kv0_spans.insert(it2 + i + 1, fr2[2] * patch2_kv0_spans[i]);
                         patch2_kv0_spans.insert(it2, fr2[1] * patch2_kv0_spans[i]);
                         patch2_kv0_spans[i] = fr2[0] * patch2_kv0_spans[i];
                         it2 = patch2_kv0_spans.begin();
                         for (index_t i = 0; i < patch2_kv0_spans.size(); i++)
                             gsInfo << patch2_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                     else if (patch2_kv0_spans[i+1] < (0.99 * patch2_kv0_spans[i])) {
                         gsInfo << "b\n";
                         patch2_inserted_knots.push_back(patch2_kv0_unique[i] + fr1[0] * patch2_kv0_spans[i]);
                         it2 = patch2_kv0_spans.insert(it2 + i + 1, fr1[1] * patch2_kv0_spans[i]);
                         patch2_kv0_spans[i] = fr1[1] * patch2_kv0_spans[i];
                         it2 = patch2_kv0_spans.begin();
                         for (index_t i = 0; i < patch2_kv0_spans.size(); i++)
                             gsInfo << patch2_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                 }
                 iter++;
                 gsInfo << "iter = " << iter << "\n";
             } while (changed);
             gsInfo << "patch2_inserted_knots: ";
             for (index_t i = 0; i < patch2_inserted_knots.size(); i++)
                 gsInfo << patch2_inserted_knots[i] << ",";
             gsInfo << "\n";
             do {
                 changed = false;
                 for (index_t i = 0; i < patch2_kv0_spans.size()-1; i++) {
                     gsInfo << i << "\n";
                     if ((3 * patch2_kv0_spans[i]) < patch2_kv0_spans[i+1]) {
                         gsInfo << "c\n";
                         patch2_inserted_knots.push_back(patch2_kv0_unique[i+1] + fr1[0] * patch2_kv0_spans[i+1]);
                         T span_old = patch2_kv0_spans[i+1];
                         patch2_kv0_spans[i+1] = fr1[1] * span_old;
                         it2 = patch2_kv0_spans.insert(it2 + i + 1, fr1[0] * span_old);
                         it2 = patch2_kv0_spans.begin();
                         for (index_t i = 0; i < patch2_kv0_spans.size(); i++)
                             gsInfo << patch2_kv0_spans[i] << ",";
                         gsInfo << "\n";
                         changed = true;
                     }
                 }
                 iter++;
                 gsInfo << "iter = " << iter << "\n";
             } while (changed);
             gsInfo << "patch2_inserted_knots: ";
             for (index_t i = 0; i < patch2_inserted_knots.size(); i++)
                 gsInfo << patch2_inserted_knots[i] << ",";
             gsInfo << "\n";
             */

             // patch 2 + patch 6
             /*if (inserted_knot == 0.15) {
                 patch2.insertKnot(0.3, 0);
                 patch2.insertKnot(0.7, 0);
                 patch6.insertKnot(0.3, 0);
                 patch6.insertKnot(0.7, 0);
             }
             else if (inserted_knot == 0.3) {
                 patch2.insertKnot(0.4, 0);
                 patch2.insertKnot(0.6, 0);
                 patch6.insertKnot(0.4, 0);
                 patch6.insertKnot(0.6, 0);
             }*/
             for (size_t i = 0; i < patch2_inserted_knots.size(); i++) {
                 patch2.insertKnot(patch2_inserted_knots[i], 0);
                 patch2.insertKnot(1 - patch2_inserted_knots[i], 0);
                 patch6.insertKnot(patch2_inserted_knots[i], 0);
                 patch6.insertKnot(1 - patch2_inserted_knots[i], 0);
             }

             // patch 1 + patch 4 + patch 3 + patch 8
             /*if (inserted_knot == 0.15) {
                 patch1.insertKnot(0.15, 0);
                 patch4.insertKnot(0.15, 0);
                 patch3.insertKnot(0.15, 0);
                 patch8.insertKnot(0.15, 0);
             }*/
             for (size_t i = 0; i < patch1_inserted_knots.size(); i++) {
                 patch1.insertKnot(patch1_inserted_knots[i], 1);
                 patch4.insertKnot(patch1_inserted_knots[i], 1);
                 patch3.insertKnot(patch1_inserted_knots[i], 0);
                 patch8.insertKnot(patch1_inserted_knots[i], 0);
             }

             // patch 6 + patch 5 + patch 7
             patch6.insertKnot(0.2, 1);
             patch6.insertKnot(0.45, 1);
             patch6.insertKnot(0.7, 1);
             patch5.insertKnot(0.2, 0);
             patch5.insertKnot(0.45, 0);
             patch5.insertKnot(0.7, 0);
             patch7.insertKnot(0.2, 1);
             patch7.insertKnot(0.45, 1);
             patch7.insertKnot(0.7, 1);

             // patch 5 + patch 4 + patch 9
             patch5.insertKnot(0.5, 1);
             patch4.insertKnot(0.5, 0);
             patch9.insertKnot(0.5, 0);

             // patch 9
             patch9.insertKnot(0.4, 1);
             patch13.insertKnot(0.4, 1);
             patch12.insertKnot(0.4, 1);
             patch8.insertKnot(0.4, 1);
             patch7.insertKnot(0.4, 0);

             patch12.insertKnot(0.4,0); //periodic with patch9
             patch10.insertKnot(0.4,0);
             patch11.insertKnot(0.4,0);
             patch13.insertKnot(0.4,0);
        }

        gsMultiPatch<T> mpFinal;
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

        //=================================optimization===========================================

        gsInfo << "Optimalizace 4. platu (mezi offsety) ...\n";
        gsInfo << "----------------------------------------\n";

        T orthogonality = 0.9;
        T skewness = 0.0;
        T eccentricity = 0.0;
        T intersection = 0.0;
        T uniformity = 0.01;
        T area = 1-orthogonality-uniformity;
        T length = 0;
        T epsilon = 1e-7;

        gsQualityMeasure<T> optimization(mpFinal.patch(3));
        //    T opt_val = optimization.functional(orthogonality, skewness,
        //                                             eccentricity, uniformity,
        //                                             length, area,
        //                                             intersection, epsilon);
        //optimization.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);

        gsQualityMeasure<T> optimization2(mpFinal.patch(11));
        optimization2.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);

        gsQualityMeasure<T> optimization3(mpFinal.patch(12));
        optimization3.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);

        //gsQualityMeasure<T> optimization4(mpFinal.patch(5));
        //optimization4.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);

        //    gsInfo << "Value of functional: "
        //           << opt_val
        //           << "\n";
        gsInfo << "Konec optimalizace.\n";


        mpFinal.addInterface(5, boundary::east, 6, boundary::west);
        mpFinal.addInterface(5, boundary::west, 4, boundary::south);
        mpFinal.addInterface(6, boundary::south, 7, boundary::west);
        mpFinal.addInterface(4, boundary::west, 3, boundary::south);
        mpFinal.addInterface(3, boundary::north, 8, boundary::south);
        mpFinal.addInterface(3, boundary::west, 0, boundary::east);
        mpFinal.addInterface(7, boundary::east, 11, boundary::west);
        mpFinal.addInterface(7, boundary::south, 2, boundary::north);
        mpFinal.addInterface(11, boundary::east, 12, boundary::west);
        mpFinal.addInterface(12, boundary::east, 8, boundary::west);
        mpFinal.addInterface(5, boundary::south, 1, boundary::north);
        mpFinal.addInterface(2, boundary::east, 9, boundary::west);
        mpFinal.addInterface(0, boundary::north, 10, boundary::east);
        mpFinal.addInterface(9, boundary::east, 10, boundary::west);
        mpFinal.addInterface(9, boundary::north, 11, boundary::south);
        mpFinal.addInterface(10, boundary::north, 12, boundary::south);
        mpFinal.addInterface(0, boundary::south, 1, boundary::west);
        mpFinal.addInterface(1, boundary::east, 2, boundary::west);
        // periodic interfaces
        mpFinal.addInterface(6, boundary::east, 4, boundary::north);
        mpFinal.addInterface(7, boundary::north, 3, boundary::east);
        mpFinal.addInterface(11, boundary::north, 8, boundary::east);
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

            //mpFinal.uniformRefine(); mpFinal.uniformRefine();
            gsWriteParaview( mpFinal, "patches", 15000, true);
        }

        /*if (plotMeshes)
        {
            std::ostringstream strs_patch;
            std::string strpatch;
            gsMultiBasis<T> tbasis(mpFinal);
            tbasis.uniformRefine(); tbasis.uniformRefine();
            gsMesh<T> mesh;
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
        */

        return mpFinal;

    }

    gsMultiPatch<T> DomainBetweenBladeProfiles3(T const & index,
        T const & length_x1,
        T const & length_x2,
        T const & pitch,
        T const & camberX,
        T const & camberY,
        T const & leadingAngle,
        T const & trailingAngle,
        T const & thicknessX,
        T const & thicknessY,
        T const & endingOffset,
        T const & outputAngle,
        T const & radius,
        T const & chordLength,
        T const & Angle,
        T const & rotationCenterX,
        T const & rotationCenterY,
        T const & uniformity_param,
        std::vector<T> const & kvfit_knots,
        bool const & coarse,
        gsVector<T> const & geom_Params)
    {
        T blade_param1 = 0.3;
        T blade_param2 = 0.6;

         //----------------knot vector seTing----------------
        //std::vector<T> kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};
        //std::vector<T> kvfit_knots = {0.,0.,0.,0.,0.2,0.4,0.6,0.8,1.,1.,1.,1.};
        gsKnotVector<T> kvfit(kvfit_knots);
        gsInfo << "degree of knot: " << kvfit.degree() << "\n";

         //----------------set parameters for blade profile----------------
        //bool plot = false;
        int num_samples = 30;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
            chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        //gsInfo << vec << "\n \n";
        //gsInfo << mat << "\n \n";
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        gsBSpline<T> suction_side_curve_transf;
        gsBSpline<T> pressure_side_curve_transf;


        //gsKnotVector<T> kvfit(0,1,4,4);

        BladeProfile<T> * pBladeProfile = 0;
        //unsigned degree = 3;
        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
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
        pBladeProfile->setSuctionSide(suction_side_curve_transf);
        pBladeProfile->setPressureSide(pressure_side_curve_transf);

        gsBSpline < T > bs = pBladeProfile->getPressureSide();
        gsBSpline < T > bp = pBladeProfile->getSuctionSide();




        GISMO_ASSERT(kvfit.inDomain(blade_param1), "blade_param1 is not in the parametric domain");
        GISMO_ASSERT(kvfit.inDomain(blade_param2), "blade_param2 is not in the parametric domain");

        std::vector<real_t>::const_iterator span1;
        span1 = kvfit.iFind(blade_param1);
        std::vector<real_t>::const_iterator span2;
        span2 = kvfit.iFind(blade_param2);

        //insert blade_param as the boundary of span or in the middle
        T middle_span1 = (span1[0]+span1[1])/2.;
        std::vector<T> diff_span1(3);
        T middle_span2 = (span2[0]+span2[1])/2.;
        std::vector<T> diff_span2(3);
        diff_span1[0] = math::abs(span1[0]-blade_param1);
        diff_span1[1] = math::abs(span1[1]-blade_param1);
        diff_span1[2] = math::abs(middle_span1-blade_param1);
        diff_span2[0] = math::abs(span2[0]-blade_param2);
        diff_span2[1] = math::abs(span2[1]-blade_param2);
        diff_span2[2] = math::abs(middle_span2-blade_param2);

        for(unsigned int i = 0; i < diff_span1.size(); i++)
        {
          gsInfo << "diff_span1" << diff_span1[i] << "\n " ;

          }
        for(unsigned int i = 0; i < diff_span2.size(); i++)
        {
          gsInfo << "diff_span2" << diff_span2[i] << "\n " ;

          }

        T blade_param_new1;
        T blade_param_new2;

        auto min_diff_span1 = std::min_element( diff_span1.begin(), diff_span1.end());

        if (*min_diff_span1 == diff_span1[2] || (span1[0]==0. || span1[1]==1.))
        {
            blade_param_new1 = middle_span1;
            bs.insertKnot(middle_span1,kvfit.degree()+1);
            bp.insertKnot(middle_span1,kvfit.degree()+1);
            kvfit.insert(middle_span1,kvfit.degree()+1);
            gsInfo << "a" << middle_span1 << "\n";
        }
        else
        {
            blade_param_new1 = span1[std::distance(diff_span1.begin(), min_diff_span1)];
            bs.insertKnot(blade_param_new1,kvfit.degree());
            bp.insertKnot(blade_param_new1,kvfit.degree());
            kvfit.insert(blade_param_new1,kvfit.degree());
        }

        auto min_diff_span2 = std::min_element( diff_span2.begin(), diff_span2.end());
        if (*min_diff_span2 == diff_span2[2] || (span2[0]==0. || span2[1]==1.))
        {
            blade_param_new2 = middle_span2;
            bs.insertKnot(middle_span2,kvfit.degree()+1);
            bp.insertKnot(middle_span2,kvfit.degree()+1);
            kvfit.insert(middle_span2,kvfit.degree()+1);
        }
        else
        {
            blade_param_new2 = span2[std::distance(diff_span2.begin(), min_diff_span2)];
            bs.insertKnot(blade_param_new2,kvfit.degree());
            bp.insertKnot(blade_param_new2,kvfit.degree());
            kvfit.insert(blade_param_new2,kvfit.degree());
        }




        gsInfo << "min_diff_span1" << *min_diff_span1 << "\n";
        gsInfo << "min_diff_span2" << *min_diff_span2 << "\n";
        gsInfo << "blade_param_new1" << blade_param_new1 << "\n";
        gsInfo << "blade_param_new2" << blade_param_new2 << "\n";


        //---------------knot vectors of patches-----------------------

   //     gsKnotVector<T> kvuniform(0,1,0,kvfit.degree()+1);
        gsKnotVector<T> kvlinear(0,1,0,2);
        gsKnotVector<T> kvleadingpart=kvfit;
        gsKnotVector<T> kvbladepart=kvfit;
        gsKnotVector<T> kvtrailingpart=kvfit;

        unsigned num_knots_leading = kvfit.degree()+1;
        unsigned num_knots_blade = kvfit.degree()+1;
         unsigned num_knots_trailing;

        unsigned k = 0;

        while (blade_param_new1 != kvfit.at(k))
        {
            num_knots_leading++;
            k++;
        }

        while (blade_param_new2 != kvfit.at(k))
        {
            num_knots_blade++;
            k++;
        }



        num_knots_trailing = kvfit.degree()+1 + kvfit.size() - num_knots_leading -  num_knots_blade + kvfit.degree()+1;
        gsInfo << "num_knots_leading" << num_knots_leading << "\n";
        gsInfo << "num_knots_blade" <<num_knots_blade << "\n";
        gsInfo << "num_knots_trailing" << num_knots_trailing << "\n";


        kvleadingpart.trimRight(kvfit.size() - num_knots_leading);

        kvleadingpart.transform(0,1);

        kvtrailingpart.trimLeft(kvfit.size() - num_knots_trailing);
        kvtrailingpart.transform(0,1);

//        gsInfo << "kvtrailingpart" << kvtrailingpart << "\n";
//        gsInfo << "kvleadingpart" << kvleadingpart << "\n";
        kvbladepart.trimLeft( (num_knots_leading - (kvfit.degree()+1)));
 //        gsInfo << "kvbladepart" << kvbladepart << "\n";
        kvbladepart.trimRight((num_knots_trailing - (kvfit.degree()+1)));
        kvbladepart.transform(0,1);


        std::vector<T> knots_diff1;
        std::vector<T> knots_diff2;

        kvbladepart.difference(kvleadingpart,knots_diff1);

        gsKnotVector<T> kvlbpart = kvleadingpart;

        for(unsigned i = 0; i < knots_diff1.size(); i++)
        {
            kvlbpart.insert(knots_diff1[i]);
        }


        kvtrailingpart.difference(kvbladepart,knots_diff2);

        gsKnotVector<T> kvfinal = kvlbpart;

        for(unsigned i = 0; i < knots_diff2.size(); i++)
        {
            kvfinal.insert(knots_diff2[i]);
        }

       gsInfo << "kvfinal" << kvfinal << "\n";

        std::vector<T> knots_diff_leadp;
         std::vector<T> knots_diff_bladep;
         std::vector<T> knots_diff_trailp;

          kvfinal.difference(kvleadingpart,knots_diff_leadp);
          kvfinal.difference(kvbladepart,knots_diff_bladep);
          kvfinal.difference(kvtrailingpart,knots_diff_trailp);


        //---------------some boundary coefs of patches-----------------------
        gsMultiPatch<T> mpFinal;

        unsigned num_cpblade = bs.coefsSize();
        unsigned num_cpleadingpart = num_knots_leading - kvfit.degree() - 1;
        unsigned num_cpbladepart = num_knots_blade - kvfit.degree() - 1;
        unsigned num_cptrailingpart = num_knots_trailing - kvfit.degree() - 1;
        unsigned num_cpfinal = kvfinal.size() -  kvfinal.degree() - 1;

        gsInfo << "num_cpblade " << num_cpblade << "\n";
        gsInfo << "num_cpleadingpart " << num_cpleadingpart << "\n";
        gsInfo << "num_cpbladepart " << num_cpbladepart << "\n";
        gsInfo << "num_cptrailingpart  " << num_cptrailingpart << "\n";
        gsInfo << "num_cpfinal " << num_cpfinal << "\n";


        gsMatrix < T > cp_bp(num_cpblade, 2);
        gsMatrix < T > cp_bs(num_cpblade, 2);
        //control points of suction,pressure side
        for (unsigned i = 0; i < num_cpblade; i++) {
            cp_bp(i, 0) = bp.coef(i, 0);
            cp_bp(i, 1) = bp.coef(i, 1) + pitch;
            cp_bs(i, 0) = bs.coef(i, 0);
            cp_bs(i, 1) = bs.coef(i, 1);
        }

        gsInfo << "cp_bp" << cp_bp << "\n";
        gsInfo << "cp_bs" << cp_bs << "\n";

        int cpk = num_cpblade-1;
        int cpk_leadingpart = num_cpleadingpart-1;
        int cpk_bladepart = num_cpleadingpart-1+num_cpbladepart;


        T ystart_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x1 + cp_bs(cpk, 1)*length_x1) / (
            cp_bs(0, 0) - cp_bs(cpk, 0)));
        T yend_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x2 + cp_bs(cpk, 1)*length_x2) / (
            cp_bs(0, 0) - cp_bs(cpk, 0)));

        T xlead_coor = (cp_bp(cpk_leadingpart,0) + cp_bs(0,0) + 2.0*length_x1)/4.0;
        T ylead_coor = (cp_bp(cpk_leadingpart,1) + cp_bs(0,1) + 2.0*ystart_coor + pitch)/4.0;

        T xblade_coor = (cp_bp(cpk_bladepart,0) + cp_bs(cpk_leadingpart,0))/2.0;
        T yblade_coor = (cp_bp(cpk_bladepart,1) + cp_bs(cpk_leadingpart,1))/2.0;

        T xtrail_coor = (cp_bp(cpk,0) + cp_bs(cpk_bladepart,0))/2.0;
        T ytrail_coor = (cp_bp(cpk,1) + cp_bs(cpk_bladepart,1))/2.0;

        T xoutlet_coor = (cp_bp(cpk,0) + cp_bs(cpk,0) + 2.0*length_x2)/4.0;
        T youtlet_coor = (cp_bp(cpk,1) + cp_bs(cpk,1) + 2.0*yend_coor + pitch)/4.0;

        //---------------leading_pressure-------------------------------
        gsMatrix<T> leading_pressure_north_coef(num_cpleadingpart, 2);
        gsMatrix<T> leading_pressure_west_coef(2, 2);
        gsMatrix<T> leading_pressure_east_coef(2, 2);
        gsMatrix<T> leading_pressure_south_coef(2, 2);
        leading_pressure_north_coef.setZero(num_cpleadingpart, 2);
        leading_pressure_west_coef.setZero(2, 2);
        leading_pressure_south_coef.setZero(2, 2);
        leading_pressure_east_coef.setZero(2, 2);

        for(unsigned i = 0; i < num_cpleadingpart; i++)
        {
            leading_pressure_north_coef(i, 0) =  cp_bp(i, 0);
            leading_pressure_north_coef(i, 1) =  cp_bp(i, 1);
        }
        leading_pressure_west_coef <<  length_x1, ystart_coor + pitch,
                                      cp_bp(0,0),cp_bp(0,1);
        leading_pressure_east_coef <<  xlead_coor, ylead_coor,
                                      cp_bp(cpk_leadingpart,0),cp_bp(cpk_leadingpart,1);
        leading_pressure_south_coef << length_x1, ystart_coor + pitch,
                                    xlead_coor, ylead_coor;

        gsBSpline<T>  leading_pressure_south(kvlinear,leading_pressure_south_coef);
        leading_pressure_south.degreeElevate(kvfinal.degree()-1);
        for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
        {
             leading_pressure_south.insertKnot(kvfinal.at(i));
        }
         gsBSpline<T>  leading_pressure_west(kvlinear,leading_pressure_west_coef);
         leading_pressure_west.degreeElevate(kvfinal.degree()-1);
         for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
         {
              leading_pressure_west.insertKnot(kvfinal.at(i));
         }
         gsBSpline<T>  leading_pressure_east(kvlinear,leading_pressure_east_coef);
         leading_pressure_east.degreeElevate(kvfinal.degree()-1);
         for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
         {
              leading_pressure_east.insertKnot(kvfinal.at(i));
         }
          gsBSpline<T>  leading_pressure_north(kvleadingpart,leading_pressure_north_coef);

        for(unsigned i = 0; i < knots_diff_leadp.size(); i++)
        {
             leading_pressure_north.insertKnot(knots_diff_leadp[i]);
        }

        gsInfo << "knots diff" << knots_diff_leadp.size();
        for(unsigned i = 0; i < knots_diff_leadp.size(); i++)
        {
            gsInfo << "knots"  << knots_diff_leadp[i] << "\n";
        }

        gsInfo << "input: " << leading_pressure_south.coefs() << "\n       " << leading_pressure_north.coefs()<< "\n      " <<leading_pressure_west.coefs()<< "\n     " <<leading_pressure_east.coefs()<< "\n";
        gsMatrix<T> leading_pressure_mat((num_cpfinal )*(num_cpfinal),2);
        leading_pressure_mat.setZero((num_cpfinal)*(num_cpfinal),2);
        discreteCoonsPatch(leading_pressure_south.coefs(),leading_pressure_north.coefs(),leading_pressure_west.coefs(),leading_pressure_east.coefs(),leading_pressure_mat,true);
        gsTensorBSplineBasis<2, T> leading_pressure_basis(kvfinal, kvfinal);
        gsTensorBSpline<2,T> leading_pressure(leading_pressure_basis,leading_pressure_mat);

        mpFinal.addPatch(leading_pressure);



        //---------------inlet_suction-------------------------------

        gsMatrix<T> inlet_suction_north_coef(num_cpfinal, 2);
        gsMatrix<T> inlet_suction_west_coef(2, 2);
        gsMatrix<T> inlet_suction_east_coef(2, 2);
        gsMatrix<T> inlet_suction_south_coef(2, 2);
        inlet_suction_north_coef.setZero(num_cpfinal, 2);
        inlet_suction_west_coef.setZero(2, 2);
        inlet_suction_south_coef.setZero(2, 2);
        inlet_suction_east_coef.setZero(2, 2);

        inlet_suction_north_coef = leading_pressure_south.coefs();

        inlet_suction_west_coef <<  length_x1,ystart_coor ,
                                    length_x1, ystart_coor + pitch;
        inlet_suction_east_coef <<  cp_bs (0,0), cp_bs (0,1),
                                    xlead_coor, ylead_coor;
        inlet_suction_south_coef << length_x1, ystart_coor,
                                    cp_bs(0,0), cp_bs(0,1);

        gsBSpline<T> inlet_suction_north(kvfinal,inlet_suction_north_coef);
        gsBSpline<T> inlet_suction_west(kvlinear,inlet_suction_west_coef);
        inlet_suction_west.degreeElevate(kvfinal.degree()-1);
        for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
        {
            inlet_suction_west.insertKnot(kvfinal.at(i));
        }

        gsBSpline<T> inlet_suction_east(kvlinear,inlet_suction_east_coef);
        inlet_suction_east.degreeElevate(kvfinal.degree()-1);
        for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
        {
            inlet_suction_east.insertKnot(kvfinal.at(i));
        }

        gsBSpline<T> inlet_suction_south(kvlinear,inlet_suction_south_coef);
        inlet_suction_south.degreeElevate(kvfinal.degree()-1);
        for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
        {
            inlet_suction_south.insertKnot(kvfinal.at(i));
        }


        gsMatrix<T> inlet_suction_mat((num_cpfinal)*(num_cpfinal),2);
        inlet_suction_mat.setZero((num_cpfinal)*(num_cpfinal),2);
        discreteCoonsPatch(inlet_suction_south.coefs(),inlet_suction_north.coefs(),inlet_suction_west.coefs(),inlet_suction_east.coefs(),inlet_suction_mat,true);
        gsTensorBSplineBasis<2, T> inlet_suction_basis(kvfinal, kvfinal);
        gsTensorBSpline<2,T> inlet_suction(inlet_suction_basis,inlet_suction_mat);



        mpFinal.addPatch(inlet_suction);

      gsInfo << "inlet_pressure computed \n";

        //---------------bladepart_pressure-------------------------------

      gsMatrix<T> bladepart_pressure_north_coef(num_cpbladepart, 2);
      gsMatrix<T> bladepart_pressure_west_coef(num_cpfinal, 2);
      gsMatrix<T> bladepart_pressure_east_coef(2, 2);
      gsMatrix<T> bladepart_pressure_south_coef(2, 2);
      bladepart_pressure_north_coef.setZero(num_cpbladepart, 2);
      bladepart_pressure_west_coef.setZero(kvfinal.degree()+1, 2);
      bladepart_pressure_south_coef.setZero(2, 2);
      bladepart_pressure_east_coef.setZero(2, 2);

      for(unsigned i = 0; i < num_cpbladepart ; i++)
      {
          bladepart_pressure_north_coef(i, 0) =  cp_bp(i + num_cpleadingpart , 0);
          bladepart_pressure_north_coef(i, 1) =  cp_bp(i + num_cpleadingpart , 1);
      }
      bladepart_pressure_west_coef = leading_pressure_east.coefs();
      bladepart_pressure_east_coef << xblade_coor, yblade_coor,
                                    cp_bp(cpk_bladepart,0),cp_bp(cpk_bladepart,1);
      bladepart_pressure_south_coef << xlead_coor,ylead_coor,
                                  xblade_coor, yblade_coor;

      gsBSpline<T>  bladepart_pressure_south(kvlinear,bladepart_pressure_south_coef);
      bladepart_pressure_south.degreeElevate(kvfinal.degree()-1);
      for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
      {
           bladepart_pressure_south.insertKnot(kvfinal.at(i));
      }

       gsBSpline<T>  bladepart_pressure_west(kvfinal,bladepart_pressure_west_coef);
       gsBSpline<T>  bladepart_pressure_east(kvlinear,bladepart_pressure_east_coef);
       bladepart_pressure_east.degreeElevate(kvfinal.degree()-1);
       for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
       {
            bladepart_pressure_east.insertKnot(kvfinal.at(i));
       }
        gsBSpline<T>  bladepart_pressure_north(kvbladepart,bladepart_pressure_north_coef);

      for(unsigned i = 0; i < knots_diff_bladep.size(); i++)
      {
           bladepart_pressure_north.insertKnot(knots_diff_bladep[i]);
      }

 gsInfo << "input: bladepart_pressure" << bladepart_pressure_south.coefs() << "\n       " << bladepart_pressure_north.coefs()<< "\n      " <<bladepart_pressure_west.coefs()<< "\n     " <<bladepart_pressure_east.coefs()<< "\n";
      gsMatrix<T> bladepart_pressure_mat((num_cpfinal )*(num_cpfinal),2);
      bladepart_pressure_mat.setZero((num_cpfinal)*(num_cpfinal),2);
      discreteCoonsPatch(bladepart_pressure_south.coefs(),bladepart_pressure_north.coefs(),bladepart_pressure_west.coefs(),bladepart_pressure_east.coefs(),bladepart_pressure_mat,true);
      gsTensorBSplineBasis<2, T> bladepart_pressure_basis(kvfinal, kvfinal);
      gsTensorBSpline<2,T> bladepart_pressure(bladepart_pressure_basis,bladepart_pressure_mat);

      mpFinal.addPatch(bladepart_pressure);


      //---------------leading_suction-------------------------------

    gsMatrix<T> leading_suction_north_coef(num_cpfinal, 2);
    gsMatrix<T> leading_suction_west_coef(num_cpfinal, 2);
    gsMatrix<T> leading_suction_east_coef(2, 2);
    gsMatrix<T> leading_suction_south_coef(num_cpleadingpart, 2);
  leading_suction_north_coef.setZero(num_cpfinal, 2);
  leading_suction_west_coef.setZero(num_cpfinal, 2);
     leading_suction_east_coef.setZero(2, 2);
    leading_suction_south_coef.setZero(num_cpleadingpart, 2);

    for(unsigned i = 0; i < num_cpleadingpart; i++)
    {
        leading_suction_south_coef(i, 0) =  cp_bs(i, 0);
        leading_suction_south_coef(i, 1) =  cp_bs(i, 1);
    }
    leading_suction_west_coef = inlet_suction_east.coefs();
    leading_suction_east_coef <<  cp_bs(cpk_leadingpart,0),cp_bs(cpk_leadingpart,1),
                                   xblade_coor,yblade_coor;
    leading_suction_north_coef = bladepart_pressure_south.coefs();

    gsBSpline<T>  leading_suction_south(kvleadingpart,leading_suction_south_coef);
    for(unsigned i = 0; i < knots_diff_leadp.size(); i++)
    {
         leading_suction_south.insertKnot(knots_diff_leadp[i]);
    }

     gsBSpline<T>  leading_suction_west(kvfinal,leading_suction_west_coef);
           gsBSpline<T>  leading_suction_east(kvlinear,leading_suction_east_coef);
     leading_suction_east.degreeElevate(kvfinal.degree()-1);
     for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
     {
          leading_suction_east.insertKnot(kvfinal.at(i));
     }
      gsBSpline<T>  leading_suction_north(kvfinal,leading_suction_north_coef);

      gsInfo << "input: leading_suction" << leading_suction_south.coefs() << "\n       " << leading_suction_north.coefs()<< "\n      " <<leading_suction_west.coefs()<< "\n     " <<leading_suction_east.coefs()<< "\n";


    gsMatrix<T> leading_suction_mat((num_cpfinal )*(num_cpfinal),2);
    leading_suction_mat.setZero((num_cpfinal)*(num_cpfinal),2);
    discreteCoonsPatch(leading_suction_south.coefs(),leading_suction_north.coefs(),leading_suction_west.coefs(),leading_suction_east.coefs(),leading_suction_mat,true);
    gsTensorBSplineBasis<2, T> leading_suction_basis(kvfinal, kvfinal);
    gsTensorBSpline<2,T> leading_suction(leading_suction_basis,leading_suction_mat);

    mpFinal.addPatch(leading_suction);

        //---------------trailing_pressure-------------------------------

      gsMatrix<T> trailing_pressure_north_coef(num_cptrailingpart, 2);
      gsMatrix<T> trailing_pressure_west_coef(num_cpfinal, 2);
      gsMatrix<T> trailing_pressure_east_coef(2, 2);
      gsMatrix<T> trailing_pressure_south_coef(2, 2);
      trailing_pressure_north_coef.setZero(num_cptrailingpart, 2);
      trailing_pressure_west_coef.setZero(num_cpfinal, 2);
      trailing_pressure_south_coef.setZero(2, 2);
      trailing_pressure_east_coef.setZero(2, 2);

      for(unsigned i = 0; i < num_cptrailingpart; i++)
      {
          trailing_pressure_north_coef(i, 0) =  cp_bp(i + num_cpleadingpart + num_cpbladepart , 0);
          trailing_pressure_north_coef(i, 1) =  cp_bp(i + num_cpleadingpart + num_cpbladepart , 1);
      }
      trailing_pressure_west_coef = bladepart_pressure_east.coefs();
      trailing_pressure_east_coef << xtrail_coor, ytrail_coor,
                                    cp_bp(cpk,0),cp_bp(cpk,1);
      trailing_pressure_south_coef << xblade_coor,yblade_coor,
                                  xtrail_coor, ytrail_coor;

      gsBSpline<T>  trailing_pressure_south(kvlinear,trailing_pressure_south_coef);
      trailing_pressure_south.degreeElevate(kvfinal.degree()-1);
      for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
      {
           trailing_pressure_south.insertKnot(kvfinal.at(i));
      }

       gsBSpline<T>  trailing_pressure_west(kvfinal,trailing_pressure_west_coef);
       gsBSpline<T>  trailing_pressure_east(kvlinear,trailing_pressure_east_coef);
       trailing_pressure_east.degreeElevate(kvfinal.degree()-1);
       for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
       {
            trailing_pressure_east.insertKnot(kvfinal.at(i));
       }
        gsBSpline<T>  trailing_pressure_north(kvtrailingpart,trailing_pressure_north_coef);

      for(unsigned i = 0; i < knots_diff_trailp.size(); i++)
      {
           trailing_pressure_north.insertKnot(knots_diff_trailp[i]);
      }



      gsMatrix<T> trailing_pressure_mat((num_cpfinal )*(num_cpfinal),2);
      trailing_pressure_mat.setZero((num_cpfinal)*(num_cpfinal),2);
      discreteCoonsPatch(trailing_pressure_south.coefs(),trailing_pressure_north.coefs(),trailing_pressure_west.coefs(),trailing_pressure_east.coefs(),trailing_pressure_mat,true);
      gsTensorBSplineBasis<2, T> trailing_pressure_basis(kvfinal, kvfinal);
      gsTensorBSpline<2,T> trailing_pressure(trailing_pressure_basis,trailing_pressure_mat);

      mpFinal.addPatch(trailing_pressure);

      gsInfo << "trailing_suction computed";

      //---------------bladepart_suction-------------------------------

    gsMatrix<T> bladepart_suction_north_coef(num_cpfinal, 2);
    gsMatrix<T> bladepart_suction_west_coef(num_cpfinal, 2);
    gsMatrix<T> bladepart_suction_east_coef(2, 2);
    gsMatrix<T> bladepart_suction_south_coef(num_cpbladepart, 2);
  bladepart_suction_north_coef.setZero(num_cpfinal, 2);
  bladepart_suction_west_coef.setZero(num_cpfinal, 2);
     bladepart_suction_east_coef.setZero(2, 2);
    bladepart_suction_south_coef.setZero(num_cpbladepart, 2);

    for(unsigned i = 0; i < num_cpbladepart; i++)
    {
        bladepart_suction_south_coef(i, 0) =  cp_bs(i + num_cpleadingpart , 0);
        bladepart_suction_south_coef(i, 1) =  cp_bs(i  + num_cpleadingpart, 1);
    }
    bladepart_suction_west_coef = leading_suction_east.coefs();
    bladepart_suction_east_coef <<  cp_bs(cpk_bladepart,0),cp_bs(cpk_bladepart,1),
                                   xtrail_coor,ytrail_coor;
    bladepart_suction_north_coef = trailing_pressure_south.coefs();

    gsBSpline<T>  bladepart_suction_south(kvbladepart,bladepart_suction_south_coef);
    for(unsigned i = 0; i < knots_diff_bladep.size(); i++)
    {
         bladepart_suction_south.insertKnot(knots_diff_bladep[i]);
    }

     gsBSpline<T>  bladepart_suction_west(kvfinal,bladepart_suction_west_coef);
           gsBSpline<T>  bladepart_suction_east(kvlinear,bladepart_suction_east_coef);
     bladepart_suction_east.degreeElevate(kvfinal.degree()-1);
     for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
     {
          bladepart_suction_east.insertKnot(kvfinal.at(i));
     }
      gsBSpline<T>  bladepart_suction_north(kvfinal,bladepart_suction_north_coef);


    gsMatrix<T> bladepart_suction_mat((num_cpfinal )*(num_cpfinal),2);
    bladepart_suction_mat.setZero((num_cpfinal)*(num_cpfinal),2);
    discreteCoonsPatch(bladepart_suction_south.coefs(),bladepart_suction_north.coefs(),bladepart_suction_west.coefs(),bladepart_suction_east.coefs(),bladepart_suction_mat,true);
    gsTensorBSplineBasis<2, T> bladepart_suction_basis(kvfinal, kvfinal);
    gsTensorBSpline<2,T> bladepart_suction(bladepart_suction_basis,bladepart_suction_mat);

    mpFinal.addPatch(bladepart_suction);

 gsInfo << "bladepart_pressure computed \n";
    //---------------outlet_pressure-------------------------------

  gsMatrix<T> outlet_pressure_north_coef(2, 2);
  gsMatrix<T> outlet_pressure_west_coef(num_cpfinal, 2);
  gsMatrix<T> outlet_pressure_east_coef(2, 2);
  gsMatrix<T> outlet_pressure_south_coef(2, 2);
  outlet_pressure_north_coef.setZero(2, 2);
  outlet_pressure_west_coef.setZero(num_cpfinal, 2);
  outlet_pressure_south_coef.setZero(2, 2);
  outlet_pressure_east_coef.setZero(2, 2);


  outlet_pressure_west_coef = trailing_pressure_east.coefs();
  outlet_pressure_north_coef << cp_bp(cpk,0),cp_bp(cpk,1),
                                length_x2, yend_coor + pitch;
  outlet_pressure_east_coef << xoutlet_coor, youtlet_coor,
                                length_x2, yend_coor + pitch;
  outlet_pressure_south_coef << xtrail_coor, ytrail_coor,
                              xoutlet_coor, youtlet_coor;

  gsBSpline<T>  outlet_pressure_south(kvlinear,outlet_pressure_south_coef);
  outlet_pressure_south.degreeElevate(kvfinal.degree()-1);
  for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
  {
       outlet_pressure_south.insertKnot(kvfinal.at(i));
  }

  gsBSpline<T>  outlet_pressure_north(kvlinear,outlet_pressure_north_coef);
  outlet_pressure_north.degreeElevate(kvfinal.degree()-1);
  for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
  {
       outlet_pressure_north.insertKnot(kvfinal.at(i));
  }

   gsBSpline<T>  outlet_pressure_west(kvfinal,outlet_pressure_west_coef);
   gsBSpline<T>  outlet_pressure_east(kvlinear,outlet_pressure_east_coef);
   outlet_pressure_east.degreeElevate(kvfinal.degree()-1);
   for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
   {
        outlet_pressure_east.insertKnot(kvfinal.at(i));
   }


  gsMatrix<T> outlet_pressure_mat((num_cpfinal )*(num_cpfinal),2);
  outlet_pressure_mat.setZero((num_cpfinal)*(num_cpfinal),2);
  discreteCoonsPatch(outlet_pressure_south.coefs(),outlet_pressure_north.coefs(),outlet_pressure_west.coefs(),outlet_pressure_east.coefs(),outlet_pressure_mat,true);
  gsTensorBSplineBasis<2, T> outlet_pressure_basis(kvfinal, kvfinal);
  gsTensorBSpline<2,T> outlet_pressure(outlet_pressure_basis,outlet_pressure_mat);

  mpFinal.addPatch(outlet_pressure);

  //---------------trailing_suction-------------------------------

gsMatrix<T> trailing_suction_north_coef(num_cpfinal, 2);
gsMatrix<T> trailing_suction_west_coef(num_cpfinal, 2);
gsMatrix<T> trailing_suction_east_coef(2, 2);
gsMatrix<T> trailing_suction_south_coef(num_cptrailingpart, 2);
trailing_suction_north_coef.setZero(num_cpfinal, 2);
trailing_suction_west_coef.setZero(num_cpfinal, 2);
 trailing_suction_east_coef.setZero(2, 2);
trailing_suction_south_coef.setZero(num_cptrailingpart, 2);

for(unsigned i = 0; i < num_cptrailingpart; i++)
{
    trailing_suction_south_coef(i, 0) =  cp_bs(i + num_cpleadingpart + num_cpbladepart, 0);
    trailing_suction_south_coef(i, 1) =  cp_bs(i  + num_cpleadingpart + num_cpbladepart, 1);
}
trailing_suction_west_coef = bladepart_suction_east.coefs();
trailing_suction_east_coef <<  cp_bs(cpk,0),cp_bs(cpk,1),
                               xoutlet_coor,youtlet_coor;
trailing_suction_north_coef = outlet_pressure_south.coefs();

gsBSpline<T>  trailing_suction_south(kvtrailingpart,trailing_suction_south_coef);
for(unsigned i = 0; i < knots_diff_trailp.size(); i++)
{
     trailing_suction_south.insertKnot(knots_diff_trailp[i]);
}

 gsBSpline<T>  trailing_suction_west(kvfinal,trailing_suction_west_coef);
       gsBSpline<T>  trailing_suction_east(kvlinear,trailing_suction_east_coef);
 trailing_suction_east.degreeElevate(kvfinal.degree()-1);
 for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
 {
      trailing_suction_east.insertKnot(kvfinal.at(i));
 }
  gsBSpline<T>  trailing_suction_north(kvfinal,trailing_suction_north_coef);


gsMatrix<T> trailing_suction_mat((num_cpfinal )*(num_cpfinal),2);
trailing_suction_mat.setZero((num_cpfinal)*(num_cpfinal),2);
discreteCoonsPatch(trailing_suction_south.coefs(),trailing_suction_north.coefs(),trailing_suction_west.coefs(),trailing_suction_east.coefs(),trailing_suction_mat,true);
gsTensorBSplineBasis<2, T> trailing_suction_basis(kvfinal, kvfinal);
gsTensorBSpline<2,T> trailing_suction(trailing_suction_basis,trailing_suction_mat);

mpFinal.addPatch(trailing_suction);

//---------------outlet-------------------------------

gsMatrix<T> outlet_north_coef(num_cpfinal, 2);
gsMatrix<T> outlet_west_coef(num_cpfinal, 2);
gsMatrix<T> outlet_east_coef(2, 2);
gsMatrix<T> outlet_south_coef(2, 2);
outlet_north_coef.setZero(num_cpfinal, 2);
outlet_west_coef.setZero(num_cpfinal, 2);
outlet_south_coef.setZero(2, 2);
outlet_east_coef.setZero(2, 2);


outlet_west_coef = trailing_suction_east.coefs();
outlet_north_coef = outlet_pressure_east.coefs();
outlet_east_coef << length_x2, yend_coor,
                            length_x2, yend_coor + pitch;
outlet_south_coef << cp_bs(cpk,0), cp_bs(cpk,1),
                         length_x2, yend_coor;

gsBSpline<T>  outlet_south(kvlinear,outlet_south_coef);
outlet_south.degreeElevate(kvfinal.degree()-1);
for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
{
    outlet_south.insertKnot(kvfinal.at(i));
}

gsBSpline<T>  outlet_north(kvfinal,outlet_north_coef);


gsBSpline<T>  outlet_west(kvfinal,outlet_west_coef);
gsBSpline<T>  outlet_east(kvlinear,outlet_east_coef);
outlet_east.degreeElevate(kvfinal.degree()-1);
for(unsigned i = kvfinal.degree()+1; i < kvfinal.size()- (kvfinal.degree()+1); i++)
{
     outlet_east.insertKnot(kvfinal.at(i));
}
gsInfo << "input: bladepart_suction" << outlet_south.coefs() << "\n       " << outlet_north.coefs()<< "\n      " <<outlet_west.coefs()<< "\n     " <<outlet_east.coefs()<< "\n";

gsMatrix<T> outlet_mat((num_cpfinal)*(num_cpfinal),2);
outlet_mat.setZero((num_cpfinal)*(num_cpfinal),2);
discreteCoonsPatch(outlet_south.coefs(),outlet_north.coefs(),outlet_west.coefs(),outlet_east.coefs(),outlet_mat,true);
gsTensorBSplineBasis<2, T> outlet_basis(kvfinal, kvfinal);
gsTensorBSpline<2,T> outlet(outlet_basis,outlet_mat);

mpFinal.addPatch(outlet);



        mpFinal.addInterface(0, boundary::south, 1, boundary::north);
       mpFinal.addInterface(2, boundary::south, 3, boundary::north);
       mpFinal.addInterface(4, boundary::south, 5, boundary::north);
       mpFinal.addInterface(6, boundary::south, 7, boundary::north);
       mpFinal.addInterface(0, boundary::east, 2, boundary::west);
       mpFinal.addInterface(2, boundary::east, 4, boundary::west);
       mpFinal.addInterface(4, boundary::east, 6, boundary::west);
       mpFinal.addInterface(1, boundary::east, 3, boundary::west);
       mpFinal.addInterface(3, boundary::east, 5, boundary::west);

       mpFinal.addInterface(5, boundary::east, 7, boundary::west);
        mpFinal.addInterface(6, boundary::east, 8, boundary::north);
         mpFinal.addInterface(7, boundary::east, 8, boundary::west);



           //periodic
           mpFinal.addInterface(0, boundary::west, 1, boundary::south);
           mpFinal.addInterface(6, boundary::north, 8, boundary::south);

        mpFinal.addAutoBoundaries();


        return mpFinal;
    }

    gsMultiPatch<T> DomainBetweenBladeProfiles3b(T const & index,
        T const & length_x1,
        T const & length_x2,
        T const & pitch,
        T const & camberX,
        T const & camberY,
        T const & leadingAngle,
        T const & trailingAngle,
        T const & thicknessX,
        T const & thicknessY,
        T const & endingOffset,
        T const & outputAngle,
        T const & radius,
        T const & chordLength,
        T const & Angle,
        T const & rotationCenterX,
        T const & rotationCenterY,
        T const & uniformity_param,
        std::vector<T> const & kvfit_knots,
        bool const & coarse,
        gsVector<T> const & geom_Params)
    {
        T blade_param1 = geom_Params(0);
        T blade_param2 = geom_Params(1);

        gsKnotVector<T> kvfit = gsKnotVector<T> (kvfit_knots);
        gsKnotVector<T> kvcub(0, 1, 0, 4);
        gsKnotVector<T> kvlin(0, 1, 0, 2);

        bool plot = true;
        //bool plotMeshes = true;
        int num_samples = 100;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
               chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        BladeProfile<T> * pBladeProfile = 0;

        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
        pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

        //---------------transform given profile-----------------------------------------
        suction_side_curve.translate(-vec);
        pressure_side_curve.translate(-vec);
        pressure_side_curve.linearTransform(mat);
        suction_side_curve.linearTransform(mat);
        pBladeProfile->setPressureSide(pressure_side_curve);
        pBladeProfile->setSuctionSide(suction_side_curve);

        vec(0) = 0.0;
        vec(1) = pitch;
        suction_side_curve.translate(vec);
        pBladeProfile->setSuctionSide(suction_side_curve);

        gsMatrix<T> suction_side_cp = suction_side_curve.coefs();
        gsMatrix<T> pressure_side_cp = pressure_side_curve.coefs();

        // --------------outer boundary - part 1 -----------------------------------------
        gsMatrix<T> aux_cp(2, 2);
        gsMatrix<T> cp_bs = suction_side_curve.coefs();
        T ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
        T yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

        aux_cp << length_x1, ystart_coor - pitch,
                  length_x1, ystart_coor;
        gsBSpline<T> left_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve.degreeElevate(2);
        gsMatrix<T> left_boundary_curve_cp = left_boundary_curve.coefs();

        aux_cp << length_x2, yend_coor - pitch,
                  length_x2, yend_coor;
        gsBSpline<T> right_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_cp = right_boundary_curve.coefs();

        aux_cp << length_x1, ystart_coor,
                  suction_side_cp(0,0), suction_side_cp(0,1);
        gsBSpline<T> top_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();

        aux_cp << suction_side_cp(suction_side_curve.coefsSize()-1,0), suction_side_cp(suction_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor;
        gsBSpline<T> top_boundary_curve_E = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_E.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_E_cp = top_boundary_curve_E.coefs();

        aux_cp << length_x1, ystart_coor - pitch,
                  pressure_side_cp(0,0), pressure_side_cp(0,1);
        gsBSpline<T> bottom_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();

        aux_cp << pressure_side_cp(pressure_side_curve.coefsSize()-1,0), pressure_side_cp(pressure_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor - pitch;
        gsBSpline<T> bottom_boundary_curve_E = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_E.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_E_cp = bottom_boundary_curve_E.coefs();

        // --------------outer boundary - part 2 -----------------------------------------
        // separation of suction side and pressure side curves
        suction_side_curve.insertKnot(blade_param1, suction_side_curve.degree() + 1);
        pressure_side_curve.insertKnot(blade_param1, pressure_side_curve.degree() + 1);
        suction_side_curve.insertKnot(blade_param2, suction_side_curve.degree() + 1);
        pressure_side_curve.insertKnot(blade_param2, pressure_side_curve.degree() + 1);
        gsInfo << suction_side_curve.knots() << "\n";

        gsKnotVector<T> kvnew = suction_side_curve.knots();
        std::vector<T> knots_for_B, knots_for_C, knots_for_D;
        for (size_t i = 0; i < kvnew.size(); i++) {
            if (kvnew[i] <= blade_param1) knots_for_B.push_back(kvnew[i]);
            if ((kvnew[i] >= blade_param1) && (kvnew[i] <= blade_param2)) knots_for_C.push_back(kvnew[i]);
            if (kvnew[i] >= blade_param2) knots_for_D.push_back(kvnew[i]);
        }
        gsKnotVector<T> kv_B = gsKnotVector<T> (knots_for_B); kv_B.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_B(kv_B.size()-kv_B.degree()-1, 2), pressure_coefs_for_B(kv_B.size()-kv_B.degree()-1, 2);
        for (size_t i = 0; i < kv_B.size()-kv_B.degree()-1; i++){
            suction_coefs_for_B.row(i) = suction_side_curve.coefs().row(i);
            pressure_coefs_for_B.row(i) = pressure_side_curve.coefs().row(i);
        }
        gsBSpline<T> suction_side_curve_B = gsBSpline<T> (kv_B, suction_coefs_for_B);
        gsBSpline<T> pressure_side_curve_B = gsBSpline<T> (kv_B, pressure_coefs_for_B);
        gsInfo << "kv_B: " << kv_B << "\n";
        gsInfo << suction_coefs_for_B << "\n";
        gsKnotVector<T> kv_C = gsKnotVector<T> (knots_for_C); kv_C.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_C(kv_C.size()-kv_C.degree()-1, 2), pressure_coefs_for_C(kv_C.size()-kv_C.degree()-1, 2);
        for (size_t i = 0; i < kv_C.size()-kv_C.degree()-1; i++){
            suction_coefs_for_C.row(i) = suction_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+i);
            pressure_coefs_for_C.row(i) = pressure_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+i);
        }
        gsBSpline<T> suction_side_curve_C = gsBSpline<T> (kv_C, suction_coefs_for_C);
        gsBSpline<T> pressure_side_curve_C = gsBSpline<T> (kv_C, pressure_coefs_for_C);
        gsInfo << "kv_C: " << kv_C << "\n";
        gsInfo << suction_coefs_for_C << "\n";
        gsKnotVector<T> kv_D = gsKnotVector<T> (knots_for_D); kv_D.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_D(kv_D.size()-kv_D.degree()-1, 2), pressure_coefs_for_D(kv_D.size()-kv_D.degree()-1, 2);
        for (size_t i = 0; i < kv_D.size()-kv_D.degree()-1; i++){
            suction_coefs_for_D.row(i) = suction_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+kv_C.size()-kv_C.degree()-1+i);
            pressure_coefs_for_D.row(i) = pressure_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+kv_C.size()-kv_C.degree()-1+i);
        }
        gsBSpline<T> suction_side_curve_D = gsBSpline<T> (kv_D, suction_coefs_for_D);
        gsBSpline<T> pressure_side_curve_D = gsBSpline<T> (kv_D, pressure_coefs_for_D);
        gsInfo << "kv_D: " << kv_D << "\n";
        gsInfo << suction_coefs_for_D << "\n";

        gsKnotVector<T> kv_BC = kv_B.knotUnion(kv_C);
        gsKnotVector<T> kv_BCD = kv_BC.knotUnion(kv_D);
        gsInfo << "kv_BCD: " << kv_BCD << "\n";

        std::vector<T> res;
        kv_BCD.difference(suction_side_curve_B.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_B.insertKnot(res[i]);
            pressure_side_curve_B.insertKnot(res[i]);
        }
        res.clear();
        kv_BCD.difference(suction_side_curve_C.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_C.insertKnot(res[i]);
            pressure_side_curve_C.insertKnot(res[i]);
        }
        res.clear();
        kv_BCD.difference(suction_side_curve_D.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_D.insertKnot(res[i]);
            pressure_side_curve_D.insertKnot(res[i]);
        }
        res.clear();
        suction_coefs_for_B.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_B.setZero();
        suction_coefs_for_C.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_C.setZero();
        suction_coefs_for_D.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_D.setZero();
        pressure_coefs_for_B.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_B.setZero();
        pressure_coefs_for_C.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_C.setZero();
        pressure_coefs_for_D.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_D.setZero();
        suction_coefs_for_B = suction_side_curve_B.coefs();
        suction_coefs_for_C = suction_side_curve_C.coefs();
        suction_coefs_for_D = suction_side_curve_D.coefs();
        pressure_coefs_for_B = pressure_side_curve_B.coefs();
        pressure_coefs_for_C = pressure_side_curve_C.coefs();
        pressure_coefs_for_D = pressure_side_curve_D.coefs();
        gsInfo << suction_side_curve_B << "\n";
        gsInfo << suction_side_curve_C << "\n";
        gsInfo << suction_side_curve_D << "\n";

        // --------------------inner joining points definition------------------------------------------------------------
        gsMatrix<T> X(1, 2), Y(1, 2), Z(1, 2), W(1, 2);
        X.row(0) = pressure_coefs_for_B.row(0) / 2 + suction_coefs_for_B.row(0) / 2;
        Y.row(0) = pressure_coefs_for_D.row(0) / 2 + suction_coefs_for_C.row(0) / 2;
        Z.row(0) = pressure_coefs_for_D.row(pressure_coefs_for_D.rows() - 1) / 2 + suction_coefs_for_D.row(0) / 2;
        W.row(0) = pressure_coefs_for_D.row(pressure_coefs_for_D.rows() - 1) / 3 + suction_coefs_for_D.row(suction_coefs_for_D.rows() - 1) / 3 + bottom_boundary_curve_E_cp.row(bottom_boundary_curve_E_cp.rows() - 1) / 6 +
                top_boundary_curve_E_cp.row(top_boundary_curve_E_cp.rows()-1) / 6;

        // --------------------cross-section curves-----------------------------------------------------------------------
        gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
        p0 << X.row(0);
        t0 << suction_coefs_for_B.row(0) - pressure_coefs_for_C.row(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        dir2 << suction_coefs_for_B.row(1) - suction_coefs_for_B.row(0);
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << suction_coefs_for_B.row(0);
        gsMatrix<T> cs_curve1_cp(4,2);
        cs_curve1_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve1_cp << "\n";
        gsBSpline<T> cs_curve1 = gsBSpline<T> ( kvcub, cs_curve1_cp );

        p3 << pressure_coefs_for_C.row(0);
        //t0 << suction_coefs_for_B.row(0) - pressure_coefs_for_C.row(0); t0 = t0/t0.norm();
        dir1 << X.row(0) - pressure_coefs_for_C.row(0);
        dir2 << top_boundary_curve_A_cp.row(0) - pressure_coefs_for_C.row(0);
        t1 = AxisDirectionNormed(dir1, dir2);
        t0 << pressure_coefs_for_C.row(0) - suction_coefs_for_B.row(0); t0 = t0/t0.norm();
        //dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        //dir2 << suction_coefs_for_B.row(1) - suction_coefs_for_B.row(0);
        //t1 = AxisDirectionNormed(dir1, dir2);
        p0 << X.row(0);
        gsMatrix<T> cs_curve2_cp(4,2);
        cs_curve2_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve2_cp << "\n";
        gsBSpline<T> cs_curve2 = gsBSpline<T> ( kvcub, cs_curve2_cp );

        p0 << Y.row(0);
        t0 << suction_coefs_for_C.row(0) - pressure_coefs_for_D.row(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << suction_coefs_for_C.row(1) - suction_coefs_for_C.row(0);
        t1(0) = dir1(1); t1(1) = - dir1(0); t1 = t1/t1.norm();
        p3 << suction_coefs_for_C.row(0);
        gsMatrix<T> cs_curve3_cp(4,2);
        cs_curve3_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve3_cp << "\n";
        gsBSpline<T> cs_curve3 = gsBSpline<T> ( kvcub, cs_curve3_cp );

        p3 << pressure_coefs_for_D.row(0);
        dir1 << pressure_coefs_for_D.row(1) - pressure_coefs_for_D.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        t0 << pressure_coefs_for_D.row(0) - suction_coefs_for_C.row(0); t0 = t0/t0.norm();
        p0 << Y.row(0);
        gsMatrix<T> cs_curve4_cp(4,2);
        cs_curve4_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve4_cp << "\n";
        gsBSpline<T> cs_curve4 = gsBSpline<T> ( kvcub, cs_curve4_cp );

        p0 << Z.row(0);
        t0 << suction_coefs_for_D.row(0) - pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << suction_coefs_for_D.row(1) - suction_coefs_for_D.row(0);
        t1(0) = dir1(1); t1(1) = - dir1(0); t1 = t1/t1.norm();
        p3 << suction_coefs_for_D.row(0);
        gsMatrix<T> cs_curve5_cp(4,2);
        cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve5_cp << "\n";
        gsBSpline<T> cs_curve5 = gsBSpline<T> ( kvcub, cs_curve5_cp );

        p3 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1);
        dir1 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1) - pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-2);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        t0 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1) - suction_coefs_for_D.row(0); t0 = t0/t0.norm();
        p0 << Z.row(0);
        gsMatrix<T> cs_curve6_cp(4,2);
        cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve6_cp << "\n";
        gsBSpline<T> cs_curve6 = gsBSpline<T> ( kvcub, cs_curve6_cp );

        p0 << X.row(0);
        dir1 << cs_curve1_cp.row(1) - cs_curve1_cp.row(0);
        t0(0) = dir1(1); t0(1) = - dir1(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << cs_curve3_cp.row(1) - cs_curve3_cp.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        p3 << Y.row(0);
        gsMatrix<T> cs_curve7_cp(4,2);
        cs_curve7_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve7_cp << "\n";
        gsBSpline<T> cs_curve7 = gsBSpline<T> ( kvcub, cs_curve7_cp );

        p0 << Y.row(0);
        dir1 << cs_curve3_cp.row(1) - cs_curve3_cp.row(0);
        t0(0) = dir1(1); t0(1) = - dir1(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << cs_curve5_cp.row(1) - cs_curve5_cp.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        p3 << Z.row(0);
        gsMatrix<T> cs_curve8_cp(4,2);
        cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve8_cp << "\n";
        gsBSpline<T> cs_curve8 = gsBSpline<T> ( kvcub, cs_curve8_cp );

        p0 << bottom_boundary_curve_A_cp.row(0);
        dir1 << top_boundary_curve_A_cp.row(0) - bottom_boundary_curve_A_cp.row(0);
        dir2 << pressure_coefs_for_B.row(0) - bottom_boundary_curve_A_cp.row(0);
        t0 = AxisDirectionNormed(dir1, dir2);
        t0 = 1.5 * t0;
        dir1 << cs_curve1_cp.row(1) - cs_curve1_cp.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        p3 << X.row(0);
        gsMatrix<T> cs_curve9_cp(4,2);
        cs_curve9_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve9_cp << "\n";
        gsBSpline<T> cs_curve9 = gsBSpline<T> ( kvcub, cs_curve9_cp );

        p0 << Z.row(0);
        dir1 << cs_curve5_cp.row(1) - cs_curve5_cp.row(0);
        t0(0) = dir1(1); t0(1) = - dir1(0); t0 = t0/t0.norm();
        t1 << length_x1 - length_x2, ystart_coor - yend_coor; t1 = t1/t1.norm();
        p3 << W.row(0);
        gsMatrix<T> cs_curve10_cp(4,2);
        cs_curve10_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve10_cp << "\n";
        gsBSpline<T> cs_curve10 = gsBSpline<T> ( kvcub, cs_curve10_cp );

        aux_cp << W.row(0),
                  suction_coefs_for_D.row(suction_coefs_for_D.rows()-1);
        gsBSpline<T> cs_curve11 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve11.degreeElevate(2);
        gsMatrix<T> cs_curve11_cp = cs_curve11.coefs();

        aux_cp << W.row(0),
                  bottom_boundary_curve_E_cp.row(bottom_boundary_curve_E_cp.rows()-1);
        gsBSpline<T> cs_curve12 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve12.degreeElevate(2);
        gsMatrix<T> cs_curve12_cp = cs_curve12.coefs();

        UnifyKnotVectors(cs_curve7, suction_side_curve_B);
        cs_curve7_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve7_cp.setZero();
        cs_curve7_cp = cs_curve7.coefs();
        UnifyKnotVectors(cs_curve8, suction_side_curve_C);
        cs_curve8_cp.resize(suction_coefs_for_C.rows(), 2); cs_curve8_cp.setZero();
        cs_curve8_cp = cs_curve8.coefs();
        UnifyKnotVectors(cs_curve9, pressure_side_curve_B);
        cs_curve9_cp.resize(pressure_coefs_for_B.rows(), 2); cs_curve9_cp.setZero();
        cs_curve9_cp = cs_curve9.coefs();
        UnifyKnotVectors(cs_curve10, suction_side_curve_D);
        cs_curve10_cp.resize(suction_coefs_for_D.rows(), 2); cs_curve10_cp.setZero();
        cs_curve10_cp = cs_curve10.coefs();
        UnifyKnotVectors(top_boundary_curve_A, pressure_side_curve_B);
        top_boundary_curve_A_cp.resize(pressure_coefs_for_B.rows(), 2); top_boundary_curve_A_cp.setZero();
        top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
        //UnifyKnotVectors(bottom_boundary_curve_A, pressure_side_curve_B);
        //bottom_boundary_curve_A_cp.resize(pressure_coefs_for_B.rows(), 2); bottom_boundary_curve_A_cp.setZero();
        //bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
        UnifyKnotVectors(bottom_boundary_curve_E, pressure_side_curve_B);
        bottom_boundary_curve_E_cp.resize(pressure_coefs_for_B.rows(), 2); bottom_boundary_curve_E_cp.setZero();
        bottom_boundary_curve_E_cp = bottom_boundary_curve_E.coefs();
        /*UnifyKnotVectors(top_boundary_curve_E, pressure_side_curve_B);
        top_boundary_curve_E_cp.resize(pressure_coefs_for_B.rows(), 2); top_boundary_curve_E_cp.setZero();
        top_boundary_curve_E_cp = top_boundary_curve_E.coefs();
        UnifyKnotVectors(cs_curve2, suction_side_curve_B);
        cs_curve2_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve2_cp.setZero();
        cs_curve2_cp = cs_curve2.coefs();
        UnifyKnotVectors(cs_curve4, suction_side_curve_B);
        cs_curve4_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve4_cp.setZero();
        cs_curve4_cp = cs_curve4.coefs();
        UnifyKnotVectors(cs_curve6, suction_side_curve_B);
        cs_curve6_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve6_cp.setZero();
        cs_curve6_cp = cs_curve6.coefs();
        UnifyKnotVectors(cs_curve12, suction_side_curve_B);
        cs_curve12_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve12_cp.setZero();
        cs_curve12_cp = cs_curve12.coefs();*/


        // -----------------------creating multipatch structure ----------------------------------------------------------------
        gsMatrix<T> patch1_cp(pressure_coefs_for_B.rows() * bottom_boundary_curve_A_cp.rows(), 2);
        discreteCoonsPatch(bottom_boundary_curve_A_cp, cs_curve2_cp, cs_curve9_cp, pressure_coefs_for_B, patch1_cp, true);
        gsTensorBSpline<2, T> patch1 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch1_cp);

        gsMatrix<T> patch2_cp(cs_curve9_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(cs_curve9_cp, top_boundary_curve_A_cp, left_boundary_curve_cp, cs_curve1_cp, patch2_cp, true);
        gsTensorBSpline<2, T> patch2 = gsTensorBSpline<2, T> (kv_BCD, kvcub, patch2_cp);

        gsMatrix<T> patch3_cp(pressure_coefs_for_C.rows() * cs_curve2_cp.rows(), 2);
        discreteCoonsPatch(cs_curve2_cp, cs_curve4_cp, cs_curve7_cp, pressure_coefs_for_C, patch3_cp, true);
        gsTensorBSpline<2, T> patch3 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch3_cp);

        gsMatrix<T> patch4_cp(cs_curve7_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(cs_curve7_cp, suction_coefs_for_B, cs_curve1_cp, cs_curve3_cp, patch4_cp, true);
        gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kv_BCD, kvcub, patch4_cp);

        gsMatrix<T> patch5_cp(pressure_coefs_for_D.rows() * cs_curve4_cp.rows(), 2);
        discreteCoonsPatch(cs_curve4_cp, cs_curve6_cp, cs_curve8_cp, pressure_coefs_for_D, patch5_cp, true);
        gsTensorBSpline<2, T> patch5 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch5_cp);

        gsMatrix<T> patch6_cp(cs_curve8_cp.rows() * cs_curve3_cp.rows(), 2);
        discreteCoonsPatch(cs_curve8_cp, suction_coefs_for_C, cs_curve3_cp, cs_curve5_cp, patch6_cp, true);
        gsTensorBSpline<2, T> patch6 = gsTensorBSpline<2, T> (kv_BCD, kvcub, patch6_cp);

        gsMatrix<T> patch7_cp(bottom_boundary_curve_E_cp.rows() * cs_curve6_cp.rows(), 2);
        discreteCoonsPatch(cs_curve6_cp, cs_curve12_cp, cs_curve10_cp, bottom_boundary_curve_E_cp, patch7_cp, true);
        gsTensorBSpline<2, T> patch7 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch7_cp);

        gsMatrix<T> patch8_cp(cs_curve10_cp.rows() * cs_curve5_cp.rows(), 2);
        discreteCoonsPatch(cs_curve10_cp, suction_coefs_for_D, cs_curve5_cp, cs_curve11_cp, patch8_cp, true);
        gsTensorBSpline<2, T> patch8 = gsTensorBSpline<2, T> (kv_BCD, kvcub, patch8_cp);

        gsMatrix<T> patch9_cp(cs_curve12_cp.rows() * cs_curve11_cp.rows(), 2);
        discreteCoonsPatch(cs_curve12_cp, top_boundary_curve_E_cp, cs_curve11_cp, right_boundary_curve_cp, patch9_cp, true);
        gsTensorBSpline<2, T> patch9 = gsTensorBSpline<2, T> (kvcub, kvcub, patch9_cp);

        // adding knots to fulfil periodic boundary conditions
        kv_BCD.difference(kvcub, res);
        for (unsigned i = 0; i < res.size(); i++) {
            patch1.insertKnot(res[i], 0);
            patch3.insertKnot(res[i], 0);
            patch5.insertKnot(res[i], 0);
            patch7.insertKnot(res[i], 0);
            patch9.insertKnot(res[i], 0);
        }
        // adding knots to other patches because of symmetry
        for (unsigned i = 0; i < res.size(); i++) {
            patch2.insertKnot(res[i], 1);
            patch4.insertKnot(res[i], 1);
            patch6.insertKnot(res[i], 1);
            patch8.insertKnot(res[i], 1);
            patch9.insertKnot(res[i], 1);
        }

        if(!coarse)
       {
           gsInfo << "Special refinement. \n";

           gsInfo << kv_BCD << "\n";
           gsInfo << patch1.knots(1) << "\n";
           std::vector<real_t> kv_BCD_unique = kv_BCD.unique();
           for (size_t i = 0; i < kv_BCD_unique.size(); i++)
               gsInfo << kv_BCD_unique[i] << ",";
           gsInfo << "\n";

           std::vector<real_t> kv_inserted_knots = smartKnotIdentification2(kv_BCD_unique);

           for (size_t i = 0; i < kv_inserted_knots.size(); i++) {
               patch1.insertKnot(kv_inserted_knots[i], 0);
               patch1.insertKnot(kv_inserted_knots[i], 1);
               patch2.insertKnot(kv_inserted_knots[i], 0);
               patch2.insertKnot(kv_inserted_knots[i], 1);
               patch3.insertKnot(kv_inserted_knots[i], 0);
               patch3.insertKnot(kv_inserted_knots[i], 1);
               patch4.insertKnot(kv_inserted_knots[i], 0);
               patch4.insertKnot(kv_inserted_knots[i], 1);
               patch5.insertKnot(kv_inserted_knots[i], 0);
               patch5.insertKnot(kv_inserted_knots[i], 1);
               patch6.insertKnot(kv_inserted_knots[i], 0);
               patch6.insertKnot(kv_inserted_knots[i], 1);
               patch7.insertKnot(kv_inserted_knots[i], 0);
               patch7.insertKnot(kv_inserted_knots[i], 1);
               patch8.insertKnot(kv_inserted_knots[i], 0);
               patch8.insertKnot(kv_inserted_knots[i], 1);
               patch9.insertKnot(kv_inserted_knots[i], 0);
               patch9.insertKnot(kv_inserted_knots[i], 1);
           }
       }



        gsMultiPatch<T> mpFinal;
        mpFinal.addPatch(patch1);
        mpFinal.addPatch(patch2);
        mpFinal.addPatch(patch3);
        mpFinal.addPatch(patch4);
        mpFinal.addPatch(patch5);
        mpFinal.addPatch(patch6);
        mpFinal.addPatch(patch7);
        mpFinal.addPatch(patch8);
        mpFinal.addPatch(patch9);

        mpFinal.addInterface(0, boundary::west, 1, boundary::south);
        mpFinal.addInterface(0, boundary::north, 2, boundary::south);
        mpFinal.addInterface(1, boundary::east, 3, boundary::west);
        mpFinal.addInterface(2, boundary::north, 4, boundary::south);
        mpFinal.addInterface(3, boundary::east, 5, boundary::west);
        mpFinal.addInterface(4, boundary::north, 6, boundary::south);
        mpFinal.addInterface(5, boundary::east, 7, boundary::west);
        mpFinal.addInterface(6, boundary::north, 8, boundary::south);
        mpFinal.addInterface(7, boundary::east, 8, boundary::west);
        mpFinal.addInterface(2, boundary::west, 3, boundary::south);
        mpFinal.addInterface(4, boundary::west, 5, boundary::south);
        mpFinal.addInterface(6, boundary::west, 7, boundary::south);

        // periodic interfaces
        mpFinal.addInterface(0, boundary::south, 1, boundary::north);
        mpFinal.addInterface(6, boundary::east, 8, boundary::north);

        mpFinal.addAutoBoundaries();


        // ---------------- plotting ---------------------------------------------------------------------------------------------------------------
        if (plot) {
            std::vector<gsGeometry<>*> curves;
            curves.clear();
            curves.push_back(&suction_side_curve_B);
            curves.push_back(&suction_side_curve_C);
            curves.push_back(&suction_side_curve_D);
            curves.push_back(&pressure_side_curve_B);
            curves.push_back(&pressure_side_curve_C);
            curves.push_back(&pressure_side_curve_D);
            curves.push_back(&cs_curve1);
            curves.push_back(&cs_curve2);
            curves.push_back(&cs_curve3);
            curves.push_back(&cs_curve4);
            curves.push_back(&cs_curve5);
            curves.push_back(&cs_curve6);
            curves.push_back(&cs_curve7);
            curves.push_back(&cs_curve8);
            curves.push_back(&cs_curve9);
            curves.push_back(&cs_curve10);
            curves.push_back(&cs_curve11);
            curves.push_back(&cs_curve12);
            curves.push_back(&left_boundary_curve);
            curves.push_back(&right_boundary_curve);
            curves.push_back(&top_boundary_curve_A);
            curves.push_back(&top_boundary_curve_E);
            curves.push_back(&bottom_boundary_curve_A);
            curves.push_back(&bottom_boundary_curve_E);


            gsWriteParaview( curves, "section_curves", 1000);

            //mpFinal.uniformRefine(); mpFinal.uniformRefine();
            gsWriteParaview( mpFinal, "patches", 50000, true);
        }



        return mpFinal;
    }

    gsMultiPatch<T> DomainBetweenBladeProfiles3c(T const & index,
        T const & length_x1,
        T const & length_x2,
        T const & pitch,
        T const & camberX,
        T const & camberY,
        T const & leadingAngle,
        T const & trailingAngle,
        T const & thicknessX,
        T const & thicknessY,
        T const & endingOffset,
        T const & outputAngle,
        T const & radius,
        T const & chordLength,
        T const & Angle,
        T const & rotationCenterX,
        T const & rotationCenterY,
        T const & uniformity_param,
        std::vector<T> const & kvfit_knots,
        bool const & coarse,
        gsVector<T> const & geom_Params)
    {
        T blade_param1 = geom_Params(0);
        T blade_param2 = geom_Params(1);

        gsKnotVector<T> kvfit = gsKnotVector<T> (kvfit_knots);
        gsKnotVector<T> kvcub(0, 1, 0, 4);
        gsKnotVector<T> kvlin(0, 1, 0, 2);

        bool plot = true;
        //bool plotMeshes = true;
        int num_samples = 100;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
               chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        BladeProfile<T> * pBladeProfile = 0;

        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
        pBladeProfile->compute2D(suction_side_curve, pressure_side_curve, kvfit, num_samples, 3, 2);

        //---------------transform given profile-----------------------------------------
        suction_side_curve.translate(-vec);
        pressure_side_curve.translate(-vec);
        pressure_side_curve.linearTransform(mat);
        suction_side_curve.linearTransform(mat);
        pBladeProfile->setPressureSide(pressure_side_curve);
        pBladeProfile->setSuctionSide(suction_side_curve);

        vec(0) = 0.0;
        vec(1) = pitch;
        suction_side_curve.translate(vec);
        pBladeProfile->setSuctionSide(suction_side_curve);

        gsMatrix<T> suction_side_cp = suction_side_curve.coefs();
        gsMatrix<T> pressure_side_cp = pressure_side_curve.coefs();

        // --------------outer boundary - part 1 -----------------------------------------
        gsMatrix<T> aux_cp(2, 2);
        gsMatrix<T> cp_bs = suction_side_curve.coefs();
        T ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
        T yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

        aux_cp << length_x1, ystart_coor,
                  length_x1, ystart_coor - pitch;
        gsBSpline<T> left_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve.degreeElevate(2);
        gsMatrix<T> left_boundary_curve_cp = left_boundary_curve.coefs();

        aux_cp << length_x2, yend_coor,
                  length_x2, yend_coor - pitch;
        gsBSpline<T> right_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_cp = right_boundary_curve.coefs();

        aux_cp << length_x1, ystart_coor,
                  suction_side_cp(0,0), suction_side_cp(0,1);
        gsBSpline<T> top_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();

        aux_cp << suction_side_cp(suction_side_curve.coefsSize()-1,0), suction_side_cp(suction_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor;
        gsBSpline<T> top_boundary_curve_E = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_E.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_E_cp = top_boundary_curve_E.coefs();

        aux_cp << length_x1, ystart_coor - pitch,
                  pressure_side_cp(0,0), pressure_side_cp(0,1);
        gsBSpline<T> bottom_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();

        aux_cp << pressure_side_cp(pressure_side_curve.coefsSize()-1,0), pressure_side_cp(pressure_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor - pitch;
        gsBSpline<T> bottom_boundary_curve_E = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_E.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_E_cp = bottom_boundary_curve_E.coefs();

        // --------------outer boundary - part 2 -----------------------------------------
        // separation of suction side and pressure side curves
        suction_side_curve.insertKnot(blade_param1, suction_side_curve.degree() + 1);
        pressure_side_curve.insertKnot(blade_param1, pressure_side_curve.degree() + 1);
        suction_side_curve.insertKnot(blade_param2, suction_side_curve.degree() + 1);
        pressure_side_curve.insertKnot(blade_param2, pressure_side_curve.degree() + 1);
        gsInfo << suction_side_curve.knots() << "\n";

        gsKnotVector<T> kvnew = suction_side_curve.knots();
        std::vector<T> knots_for_B, knots_for_C, knots_for_D;
        for (size_t i = 0; i < kvnew.size(); i++) {
            if (kvnew[i] <= blade_param1) knots_for_B.push_back(kvnew[i]);
            if ((kvnew[i] >= blade_param1) && (kvnew[i] <= blade_param2)) knots_for_C.push_back(kvnew[i]);
            if (kvnew[i] >= blade_param2) knots_for_D.push_back(kvnew[i]);
        }
        gsKnotVector<T> kv_B = gsKnotVector<T> (knots_for_B); kv_B.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_B(kv_B.size()-kv_B.degree()-1, 2), pressure_coefs_for_B(kv_B.size()-kv_B.degree()-1, 2);
        for (size_t i = 0; i < kv_B.size()-kv_B.degree()-1; i++){
            suction_coefs_for_B.row(i) = suction_side_curve.coefs().row(i);
            pressure_coefs_for_B.row(i) = pressure_side_curve.coefs().row(i);
        }
        gsBSpline<T> suction_side_curve_B = gsBSpline<T> (kv_B, suction_coefs_for_B);
        gsBSpline<T> pressure_side_curve_B = gsBSpline<T> (kv_B, pressure_coefs_for_B);
        gsInfo << "kv_B: " << kv_B << "\n";
        gsInfo << suction_coefs_for_B << "\n";
        gsKnotVector<T> kv_C = gsKnotVector<T> (knots_for_C); kv_C.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_C(kv_C.size()-kv_C.degree()-1, 2), pressure_coefs_for_C(kv_C.size()-kv_C.degree()-1, 2);
        for (size_t i = 0; i < kv_C.size()-kv_C.degree()-1; i++){
            suction_coefs_for_C.row(i) = suction_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+i);
            pressure_coefs_for_C.row(i) = pressure_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+i);
        }
        gsBSpline<T> suction_side_curve_C = gsBSpline<T> (kv_C, suction_coefs_for_C);
        gsBSpline<T> pressure_side_curve_C = gsBSpline<T> (kv_C, pressure_coefs_for_C);
        gsInfo << "kv_C: " << kv_C << "\n";
        gsInfo << suction_coefs_for_C << "\n";
        gsKnotVector<T> kv_D = gsKnotVector<T> (knots_for_D); kv_D.affineTransformTo(0.0, 1.0);
        gsMatrix<T> suction_coefs_for_D(kv_D.size()-kv_D.degree()-1, 2), pressure_coefs_for_D(kv_D.size()-kv_D.degree()-1, 2);
        for (size_t i = 0; i < kv_D.size()-kv_D.degree()-1; i++){
            suction_coefs_for_D.row(i) = suction_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+kv_C.size()-kv_C.degree()-1+i);
            pressure_coefs_for_D.row(i) = pressure_side_curve.coefs().row(kv_B.size()-kv_B.degree()-1+kv_C.size()-kv_C.degree()-1+i);
        }
        gsBSpline<T> suction_side_curve_D = gsBSpline<T> (kv_D, suction_coefs_for_D);
        gsBSpline<T> pressure_side_curve_D = gsBSpline<T> (kv_D, pressure_coefs_for_D);
        gsInfo << "kv_D: " << kv_D << "\n";
        gsInfo << suction_coefs_for_D << "\n";

        gsKnotVector<T> kv_BC = kv_B.knotUnion(kv_C);
        gsKnotVector<T> kv_BCD = kv_BC.knotUnion(kv_D);
        gsInfo << "kv_BCD: " << kv_BCD << "\n";

        std::vector<T> res;
        kv_BCD.difference(suction_side_curve_B.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_B.insertKnot(res[i]);
            pressure_side_curve_B.insertKnot(res[i]);
        }
        res.clear();
        kv_BCD.difference(suction_side_curve_C.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_C.insertKnot(res[i]);
            pressure_side_curve_C.insertKnot(res[i]);
        }
        res.clear();
        kv_BCD.difference(suction_side_curve_D.knots(), res);
        for (unsigned i = 0; i < res.size(); i++) {
            suction_side_curve_D.insertKnot(res[i]);
            pressure_side_curve_D.insertKnot(res[i]);
        }
        res.clear();
        suction_coefs_for_B.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_B.setZero();
        suction_coefs_for_C.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_C.setZero();
        suction_coefs_for_D.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); suction_coefs_for_D.setZero();
        pressure_coefs_for_B.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_B.setZero();
        pressure_coefs_for_C.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_C.setZero();
        pressure_coefs_for_D.resize(kv_BCD.size() - kv_BCD.degree() - 1, 2); pressure_coefs_for_D.setZero();
        suction_coefs_for_B = suction_side_curve_B.coefs();
        suction_coefs_for_C = suction_side_curve_C.coefs();
        suction_coefs_for_D = suction_side_curve_D.coefs();
        pressure_coefs_for_B = pressure_side_curve_B.coefs();
        pressure_coefs_for_C = pressure_side_curve_C.coefs();
        pressure_coefs_for_D = pressure_side_curve_D.coefs();
        gsInfo << suction_side_curve_B << "\n";
        gsInfo << suction_side_curve_C << "\n";
        gsInfo << suction_side_curve_D << "\n";

        // --------------------inner joining points definition------------------------------------------------------------
        gsMatrix<T> X(1, 2), Y(1, 2), Z(1, 2), W(1, 2), U(1, 2);
        X.row(0) = pressure_coefs_for_B.row(0) / 2 + suction_coefs_for_B.row(0) / 2;
        Y.row(0) = pressure_coefs_for_D.row(0) / 2 + suction_coefs_for_C.row(0) / 2;
        Z.row(0) = pressure_coefs_for_D.row(pressure_coefs_for_D.rows() - 1) / 2 + suction_coefs_for_D.row(0) / 2;
        W.row(0) = pressure_coefs_for_D.row(pressure_coefs_for_D.rows() - 1) / 3 + suction_coefs_for_D.row(suction_coefs_for_D.rows() - 1) / 3 + bottom_boundary_curve_E_cp.row(bottom_boundary_curve_E_cp.rows() - 1) / 6 +
                top_boundary_curve_E_cp.row(top_boundary_curve_E_cp.rows()-1) / 6;
        U.row(0) = pressure_coefs_for_B.row(0) / 3 + suction_coefs_for_B.row(0) / 3 + bottom_boundary_curve_A_cp.row(0) / 6 + top_boundary_curve_A_cp.row(0) / 6;

        // --------------------cross-section curves-----------------------------------------------------------------------
        gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
        p3 << X.row(0);
        t1 << suction_coefs_for_B.row(0) - pressure_coefs_for_C.row(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        dir2 << suction_coefs_for_B.row(1) - suction_coefs_for_B.row(0);
        t0 = AxisDirectionNormed(dir1, dir2);
        p0 << suction_coefs_for_B.row(0);
        gsMatrix<T> cs_curve1_cp(4,2);
        cs_curve1_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve1_cp << "\n";
        gsBSpline<T> cs_curve1 = gsBSpline<T> ( kvcub, cs_curve1_cp );

        p3 << pressure_coefs_for_C.row(0);
        //t0 << suction_coefs_for_B.row(0) - pressure_coefs_for_C.row(0); t0 = t0/t0.norm();
        dir1 << X.row(0) - pressure_coefs_for_C.row(0);
        dir2 << top_boundary_curve_A_cp.row(0) - pressure_coefs_for_C.row(0);
        t1 = AxisDirectionNormed(dir1, dir2);
        t0 << pressure_coefs_for_C.row(0) - suction_coefs_for_B.row(0); t0 = t0/t0.norm();
        //dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        //dir2 << suction_coefs_for_B.row(1) - suction_coefs_for_B.row(0);
        //t1 = AxisDirectionNormed(dir1, dir2);
        p0 << X.row(0);
        gsMatrix<T> cs_curve2_cp(4,2);
        cs_curve2_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve2_cp << "\n";
        gsBSpline<T> cs_curve2 = gsBSpline<T> ( kvcub, cs_curve2_cp );

        p3 << Y.row(0);
        t1 << suction_coefs_for_C.row(0) - pressure_coefs_for_D.row(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << suction_coefs_for_C.row(1) - suction_coefs_for_C.row(0);
        t0(0) = dir1(1); t0(1) = - dir1(0); t0 = t0/t0.norm();
        p0 << suction_coefs_for_C.row(0);
        gsMatrix<T> cs_curve3_cp(4,2);
        cs_curve3_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve3_cp << "\n";
        gsBSpline<T> cs_curve3 = gsBSpline<T> ( kvcub, cs_curve3_cp );

        p3 << pressure_coefs_for_D.row(0);
        dir1 << pressure_coefs_for_D.row(1) - pressure_coefs_for_D.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        t0 << pressure_coefs_for_D.row(0) - suction_coefs_for_C.row(0); t0 = t0/t0.norm();
        p0 << Y.row(0);
        gsMatrix<T> cs_curve4_cp(4,2);
        cs_curve4_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve4_cp << "\n";
        gsBSpline<T> cs_curve4 = gsBSpline<T> ( kvcub, cs_curve4_cp );

        p3 << Z.row(0);
        t1 << suction_coefs_for_D.row(0) - pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << suction_coefs_for_D.row(1) - suction_coefs_for_D.row(0);
        t0(0) = dir1(1); t0(1) = - dir1(0); t0 = t0/t0.norm();
        p0 << suction_coefs_for_D.row(0);
        gsMatrix<T> cs_curve5_cp(4,2);
        cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve5_cp << "\n";
        gsBSpline<T> cs_curve5 = gsBSpline<T> ( kvcub, cs_curve5_cp );

        p3 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1);
        dir1 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1) - pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-2);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        t0 << pressure_coefs_for_D.row(pressure_coefs_for_D.rows()-1) - suction_coefs_for_D.row(0); t0 = t0/t0.norm();
        p0 << Z.row(0);
        gsMatrix<T> cs_curve6_cp(4,2);
        cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve6_cp << "\n";
        gsBSpline<T> cs_curve6 = gsBSpline<T> ( kvcub, cs_curve6_cp );

        p0 << X.row(0);
        dir1 << cs_curve2_cp.row(1) - cs_curve2_cp.row(0);
        t0(0) = - dir1(1); t0(1) = dir1(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << cs_curve4_cp.row(1) - cs_curve4_cp.row(0);
        t1(0) = dir1(1); t1(1) = - dir1(0); t1 = t1/t1.norm();
        p3 << Y.row(0);
        gsMatrix<T> cs_curve7_cp(4,2);
        cs_curve7_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve7_cp << "\n";
        gsBSpline<T> cs_curve7 = gsBSpline<T> ( kvcub, cs_curve7_cp );

        p0 << Y.row(0);
        dir1 << cs_curve4_cp.row(1) - cs_curve4_cp.row(0);
        t0(0) = - dir1(1); t0(1) = dir1(0); t0 = t0/t0.norm();
        //dir1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0);
        //dir2 << bottomSplitPoint1 - suction_side_offset_trimmed_cp.row(0);
        //t0 = AxisDirectionNormed(dir1, dir2);
        //t1 << suction_side_offset_trimmed_cp.row(0) - suction_side_trimmed_cp.row(0); t1 = t1/t1.norm();
        dir1 << cs_curve6_cp.row(1) - cs_curve6_cp.row(0);
        t1(0) = dir1(1); t1(1) = - dir1(0); t1 = t1/t1.norm();
        p3 << Z.row(0);
        gsMatrix<T> cs_curve8_cp(4,2);
        cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve8_cp << "\n";
        gsBSpline<T> cs_curve8 = gsBSpline<T> ( kvcub, cs_curve8_cp );

        p0 << U.row(0);
        t0 << X.row(0) - U.row(0); t0 = t0/t0.norm();
        t1 << cs_curve7_cp.row(0) - cs_curve7_cp.row(1); t1 = t1/t1.norm();
        p3 << X.row(0);
        gsMatrix<T> cs_curve9_cp(4,2);
        cs_curve9_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve9_cp << "\n";
        gsBSpline<T> cs_curve9 = gsBSpline<T> ( kvcub, cs_curve9_cp );

        p0 << Z.row(0);
        dir1 << cs_curve6_cp.row(1) - cs_curve6_cp.row(0);
        t0(0) = - dir1(1); t0(1) = dir1(0); t0 = t0/t0.norm();
        t1 << length_x1 - length_x2, ystart_coor - yend_coor; t1 = t1/t1.norm();
        p3 << W.row(0);
        gsMatrix<T> cs_curve10_cp(4,2);
        cs_curve10_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve10_cp << "\n";
        gsBSpline<T> cs_curve10 = gsBSpline<T> ( kvcub, cs_curve10_cp );

        aux_cp << suction_coefs_for_D.row(suction_coefs_for_D.rows()-1),
                  W.row(0);
        gsBSpline<T> cs_curve11 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve11.degreeElevate(2);
        gsMatrix<T> cs_curve11_cp = cs_curve11.coefs();

        aux_cp << W.row(0),
                  bottom_boundary_curve_E_cp.row(bottom_boundary_curve_E_cp.rows()-1);
        gsBSpline<T> cs_curve12 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve12.degreeElevate(2);
        gsMatrix<T> cs_curve12_cp = cs_curve12.coefs();

        p0 << U.row(0);
        t0 << bottom_boundary_curve_A_cp.row(0) - top_boundary_curve_A_cp.row(0); t0 = t0/t0.norm();
        dir1 << cs_curve5_cp.row(1) - cs_curve5_cp.row(0);
        t1(0) = - dir1(1); t1(1) = dir1(0); t1 = t1/t1.norm();
        dir1 << pressure_coefs_for_B.row(1) - pressure_coefs_for_B.row(0);
        dir2 << bottom_boundary_curve_A_cp.row(0) - bottom_boundary_curve_A_cp.row(bottom_boundary_curve_A_cp.rows()-1);
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << pressure_coefs_for_B.row(0);
        gsMatrix<T> cs_curve13_cp(4,2);
        cs_curve13_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve13_cp << "\n";
        gsBSpline<T> cs_curve13 = gsBSpline<T> ( kvcub, cs_curve13_cp );

        aux_cp << top_boundary_curve_A_cp.row(0),
                  U.row(0);
        gsBSpline<T> cs_curve14 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve14.degreeElevate(2);
        gsMatrix<T> cs_curve14_cp = cs_curve14.coefs();

        UnifyKnotVectors(cs_curve7, suction_side_curve_B);
        cs_curve7_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve7_cp.setZero();
        cs_curve7_cp = cs_curve7.coefs();
        UnifyKnotVectors(cs_curve8, suction_side_curve_C);
        cs_curve8_cp.resize(suction_coefs_for_C.rows(), 2); cs_curve8_cp.setZero();
        cs_curve8_cp = cs_curve8.coefs();
        UnifyKnotVectors(cs_curve9, pressure_side_curve_B);
        cs_curve9_cp.resize(pressure_coefs_for_B.rows(), 2); cs_curve9_cp.setZero();
        cs_curve9_cp = cs_curve9.coefs();
        UnifyKnotVectors(cs_curve10, suction_side_curve_D);
        cs_curve10_cp.resize(suction_coefs_for_D.rows(), 2); cs_curve10_cp.setZero();
        cs_curve10_cp = cs_curve10.coefs();
        UnifyKnotVectors(top_boundary_curve_A, pressure_side_curve_B);
        top_boundary_curve_A_cp.resize(pressure_coefs_for_B.rows(), 2); top_boundary_curve_A_cp.setZero();
        top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
        //UnifyKnotVectors(bottom_boundary_curve_A, pressure_side_curve_B);
        //bottom_boundary_curve_A_cp.resize(pressure_coefs_for_B.rows(), 2); bottom_boundary_curve_A_cp.setZero();
        //bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
        UnifyKnotVectors(bottom_boundary_curve_E, pressure_side_curve_B);
        bottom_boundary_curve_E_cp.resize(pressure_coefs_for_B.rows(), 2); bottom_boundary_curve_E_cp.setZero();
        bottom_boundary_curve_E_cp = bottom_boundary_curve_E.coefs();
        /*UnifyKnotVectors(top_boundary_curve_E, pressure_side_curve_B);
        top_boundary_curve_E_cp.resize(pressure_coefs_for_B.rows(), 2); top_boundary_curve_E_cp.setZero();
        top_boundary_curve_E_cp = top_boundary_curve_E.coefs();
        UnifyKnotVectors(cs_curve2, suction_side_curve_B);
        cs_curve2_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve2_cp.setZero();
        cs_curve2_cp = cs_curve2.coefs();
        UnifyKnotVectors(cs_curve4, suction_side_curve_B);
        cs_curve4_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve4_cp.setZero();
        cs_curve4_cp = cs_curve4.coefs();
        UnifyKnotVectors(cs_curve6, suction_side_curve_B);
        cs_curve6_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve6_cp.setZero();
        cs_curve6_cp = cs_curve6.coefs();
        UnifyKnotVectors(cs_curve12, suction_side_curve_B);
        cs_curve12_cp.resize(suction_coefs_for_B.rows(), 2); cs_curve12_cp.setZero();
        cs_curve12_cp = cs_curve12.coefs();*/


        // -----------------------creating multipatch structure ----------------------------------------------------------------
        gsMatrix<T> patch1_cp(pressure_coefs_for_B.rows() * bottom_boundary_curve_A_cp.rows(), 2);
        discreteCoonsPatch(cs_curve13_cp, cs_curve2_cp, cs_curve9_cp, pressure_coefs_for_B, patch1_cp, true);
        gsTensorBSpline<2, T> patch1 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch1_cp);

        gsMatrix<T> patch2_cp(cs_curve9_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(cs_curve14_cp, cs_curve1_cp, top_boundary_curve_A_cp, cs_curve9_cp, patch2_cp, true);
        gsTensorBSpline<2, T> patch2 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch2_cp);

        gsMatrix<T> patch3_cp(pressure_coefs_for_C.rows() * cs_curve2_cp.rows(), 2);
        discreteCoonsPatch(cs_curve2_cp, cs_curve4_cp, cs_curve7_cp, pressure_coefs_for_C, patch3_cp, true);
        gsTensorBSpline<2, T> patch3 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch3_cp);

        gsMatrix<T> patch4_cp(cs_curve7_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(cs_curve1_cp, cs_curve3_cp, suction_coefs_for_B, cs_curve7_cp, patch4_cp, true);
        gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch4_cp);

        gsMatrix<T> patch5_cp(pressure_coefs_for_D.rows() * cs_curve4_cp.rows(), 2);
        discreteCoonsPatch(cs_curve4_cp, cs_curve6_cp, cs_curve8_cp, pressure_coefs_for_D, patch5_cp, true);
        gsTensorBSpline<2, T> patch5 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch5_cp);

        gsMatrix<T> patch6_cp(cs_curve8_cp.rows() * cs_curve3_cp.rows(), 2);
        discreteCoonsPatch(cs_curve3_cp, cs_curve5_cp, suction_coefs_for_C, cs_curve8_cp, patch6_cp, true);
        gsTensorBSpline<2, T> patch6 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch6_cp);

        gsMatrix<T> patch7_cp(bottom_boundary_curve_E_cp.rows() * cs_curve6_cp.rows(), 2);
        discreteCoonsPatch(cs_curve6_cp, cs_curve12_cp, cs_curve10_cp, bottom_boundary_curve_E_cp, patch7_cp, true);
        gsTensorBSpline<2, T> patch7 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch7_cp);

        gsMatrix<T> patch8_cp(cs_curve10_cp.rows() * cs_curve5_cp.rows(), 2);
        discreteCoonsPatch(cs_curve5_cp, cs_curve11_cp, suction_coefs_for_D, cs_curve10_cp, patch8_cp, true);
        gsTensorBSpline<2, T> patch8 = gsTensorBSpline<2, T> (kvcub, kv_BCD, patch8_cp);

        gsMatrix<T> patch9_cp(cs_curve12_cp.rows() * cs_curve11_cp.rows(), 2);
        discreteCoonsPatch(cs_curve11_cp, right_boundary_curve_cp, top_boundary_curve_E_cp, cs_curve12_cp, patch9_cp, true);
        gsTensorBSpline<2, T> patch9 = gsTensorBSpline<2, T> (kvcub, kvcub, patch9_cp);

        gsMatrix<T> patch10_cp(cs_curve14_cp.rows() * cs_curve13_cp.rows(), 2);
        discreteCoonsPatch(left_boundary_curve_cp, cs_curve13_cp, cs_curve14_cp, bottom_boundary_curve_A_cp, patch10_cp, true);
        gsTensorBSpline<2, T> patch10 = gsTensorBSpline<2, T> (kvcub, kvcub, patch10_cp);

        // adding knots to fulfil periodic boundary conditions
        kv_BCD.difference(kvcub, res);
        for (unsigned i = 0; i < res.size(); i++) {
            patch1.insertKnot(res[i], 0);
            patch3.insertKnot(res[i], 0);
            patch5.insertKnot(res[i], 0);
            patch7.insertKnot(res[i], 0);
            patch9.insertKnot(res[i], 0);
            patch10.insertKnot(res[i], 0);
        }
        // adding knots to other patches because of symmetry
        for (unsigned i = 0; i < res.size(); i++) {
            patch2.insertKnot(res[i], 0);
            patch4.insertKnot(res[i], 0);
            patch6.insertKnot(res[i], 0);
            patch8.insertKnot(res[i], 0);
            patch9.insertKnot(res[i], 1);
            patch10.insertKnot(res[i], 1);
        }

        if(!coarse)
       {
           gsInfo << "Special refinement. \n";

           gsInfo << kv_BCD << "\n";
           gsInfo << patch1.knots(1) << "\n";
           std::vector<real_t> kv_BCD_unique = kv_BCD.unique();
           for (size_t i = 0; i < kv_BCD_unique.size(); i++)
               gsInfo << kv_BCD_unique[i] << ",";
           gsInfo << "\n";

           std::vector<real_t> kv_inserted_knots = smartKnotIdentification2(kv_BCD_unique);

           for (size_t i = 0; i < kv_inserted_knots.size(); i++) {
               patch1.insertKnot(kv_inserted_knots[i], 0);
               patch1.insertKnot(kv_inserted_knots[i], 1);
               patch2.insertKnot(kv_inserted_knots[i], 0);
               patch2.insertKnot(kv_inserted_knots[i], 1);
               patch3.insertKnot(kv_inserted_knots[i], 0);
               patch3.insertKnot(kv_inserted_knots[i], 1);
               patch4.insertKnot(kv_inserted_knots[i], 0);
               patch4.insertKnot(kv_inserted_knots[i], 1);
               patch5.insertKnot(kv_inserted_knots[i], 0);
               patch5.insertKnot(kv_inserted_knots[i], 1);
               patch6.insertKnot(kv_inserted_knots[i], 0);
               patch6.insertKnot(kv_inserted_knots[i], 1);
               patch7.insertKnot(kv_inserted_knots[i], 0);
               patch7.insertKnot(kv_inserted_knots[i], 1);
               patch8.insertKnot(kv_inserted_knots[i], 0);
               patch8.insertKnot(kv_inserted_knots[i], 1);
               patch9.insertKnot(kv_inserted_knots[i], 0);
               patch9.insertKnot(kv_inserted_knots[i], 1);
               patch10.insertKnot(kv_inserted_knots[i], 0);
               patch10.insertKnot(kv_inserted_knots[i], 1);
           }
       }



        gsMultiPatch<T> mpFinal;
        mpFinal.addPatch(patch1);
        mpFinal.addPatch(patch2);
        mpFinal.addPatch(patch3);
        mpFinal.addPatch(patch4);
        mpFinal.addPatch(patch5);
        mpFinal.addPatch(patch6);
        mpFinal.addPatch(patch7);
        mpFinal.addPatch(patch8);
        mpFinal.addPatch(patch9);
        mpFinal.addPatch(patch10);

        mpFinal.addInterface(0, boundary::west, 1, boundary::east);
        mpFinal.addInterface(0, boundary::north, 2, boundary::south);
        mpFinal.addInterface(1, boundary::north, 3, boundary::south);
        mpFinal.addInterface(2, boundary::north, 4, boundary::south);
        mpFinal.addInterface(3, boundary::north, 5, boundary::south);
        mpFinal.addInterface(4, boundary::north, 6, boundary::south);
        mpFinal.addInterface(5, boundary::north, 7, boundary::south);
        mpFinal.addInterface(6, boundary::north, 8, boundary::east);
        mpFinal.addInterface(7, boundary::north, 8, boundary::south);
        mpFinal.addInterface(2, boundary::west, 3, boundary::east);
        mpFinal.addInterface(4, boundary::west, 5, boundary::east);
        mpFinal.addInterface(6, boundary::west, 7, boundary::east);
        mpFinal.addInterface(0, boundary::south, 9, boundary::north);
        mpFinal.addInterface(1, boundary::south, 9, boundary::west);

        // periodic interfaces
        mpFinal.addInterface(9, boundary::east, 1, boundary::west);
        mpFinal.addInterface(6, boundary::east, 8, boundary::west);

        mpFinal.addAutoBoundaries();


        // ---------------- plotting ---------------------------------------------------------------------------------------------------------------
        if (plot) {
            std::vector<gsGeometry<>*> curves;
            curves.clear();
            curves.push_back(&suction_side_curve_B);
            curves.push_back(&suction_side_curve_C);
            curves.push_back(&suction_side_curve_D);
            curves.push_back(&pressure_side_curve_B);
            curves.push_back(&pressure_side_curve_C);
            curves.push_back(&pressure_side_curve_D);
            curves.push_back(&cs_curve1);
            curves.push_back(&cs_curve2);
            curves.push_back(&cs_curve3);
            curves.push_back(&cs_curve4);
            curves.push_back(&cs_curve5);
            curves.push_back(&cs_curve6);
            curves.push_back(&cs_curve7);
            curves.push_back(&cs_curve8);
            curves.push_back(&cs_curve9);
            curves.push_back(&cs_curve10);
            curves.push_back(&cs_curve11);
            curves.push_back(&cs_curve12);
            curves.push_back(&cs_curve13);
            curves.push_back(&cs_curve14);
            curves.push_back(&left_boundary_curve);
            curves.push_back(&right_boundary_curve);
            curves.push_back(&top_boundary_curve_A);
            curves.push_back(&top_boundary_curve_E);
            curves.push_back(&bottom_boundary_curve_A);
            curves.push_back(&bottom_boundary_curve_E);


            gsWriteParaview( curves, "section_curves", 1000);

            //mpFinal.uniformRefine(); mpFinal.uniformRefine();
            gsWriteParaview( mpFinal, "patches", 50000, true);
        }



        return mpFinal;
    }


    gsMultiPatch<T> DomainAroundBladeProfile4(T const & index,
               T const & length_x1,
               T const & length_x2,
               T const & pitch,
               T const & camberX,
               T const & camberY,
               T const & leadingAngle,
               T const & trailingAngle,
               T const & thicknessX,
               T const & thicknessY,
               T const & endingOffset,
               T const & outputAngle,
               T const & radius,
               T const & chordLength,
               T const & Angle,
               T const & rotationCenterX,
               T const & rotationCenterY,
               T const & uniformity_param,
               std::vector<T> const & kvfit_knots,
               bool const & coarse,
               gsVector<T> const & geom_Params)
           {

               T blade_param = geom_Params(0);
               T boundary_param = geom_Params(1);
               //T fb = 1-(0.5 + camberY + thicknessY);
               T fb = 0.5 - camberY;

               gsKnotVector<T> kvfit(kvfit_knots);


               gsInfo << "degree of knot: " << kvfit.degree() << "\n";

                //----------------set parameters for blade profile----------------
               //bool plot = false;
               int num_samples = 30;
               gsVector<T> vec(2);
               //gsInfo << pitch << "\n ";
               vec(0) = rotationCenterX;
               vec(1) = rotationCenterY;
               gsMatrix<T> mat(2, 2);
               mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
                   chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
               //gsInfo << vec << "\n \n";
               //gsInfo << mat << "\n \n";
               gsBSpline<T> suction_side_curve;
               gsBSpline<T> pressure_side_curve;
               gsBSpline<T> suction_side_curve_transf;
               gsBSpline<T> pressure_side_curve_transf;


               //gsKnotVector<T> kvfit(0,1,4,4);

               BladeProfile<T> * pBladeProfile = 0;
               //unsigned degree = 3;
               //---------------compute blade profile for given parameters----------------------
               pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
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
               pBladeProfile->setSuctionSide(suction_side_curve_transf);
               pBladeProfile->setPressureSide(pressure_side_curve_transf);

               gsBSpline < T > bs = pBladeProfile->getPressureSide();
               gsBSpline < T > bp = pBladeProfile->getSuctionSide();

               T blade_param_new;
               size_t numinsert_blade_param_new;

               blade_param_new = smartKnotParameterInsert(blade_param,kvfit);
               if (kvfit.has(blade_param_new))
                   numinsert_blade_param_new = kvfit.degree();
               else
                   numinsert_blade_param_new = kvfit.degree()+1;


               bs.insertKnot(blade_param_new,numinsert_blade_param_new);
               bp.insertKnot(blade_param_new,numinsert_blade_param_new);
               kvfit.insert(blade_param_new,numinsert_blade_param_new);


               //---------------knot vectors of patches-----------------------

               gsKnotVector<T> kvuniform(0,1,0,kvfit.degree()+1);
               gsKnotVector<T> kvcubic(0,1,0,4);
               gsKnotVector<T> kvlinear(0,1,0,2);
               gsKnotVector<T> kvinlet = kvfit;
               gsKnotVector<T> kvperiodic = kvfit;

               unsigned num_knots_input = kvfit.degree()+1;
               unsigned num_knots_periodic = 0;

               unsigned k = 0;

               while (blade_param_new != kvfit.at(k))
               {
                   num_knots_input++;
                   k++;
               }
               num_knots_periodic = kvfit.size() - num_knots_input + kvfit.degree()+1;



               kvinlet.trimRight(kvfit.size() - num_knots_input);
               kvinlet.transform(0,1);
               kvperiodic.trimLeft(kvfit.size() - num_knots_periodic);
               kvperiodic.transform(0,1);

               //---------------some boundary coefs of patches-----------------------
               gsMultiPatch<T> mpFinal;

               unsigned num_cpblade = bs.coefsSize();
               unsigned num_cpinput = num_knots_input - kvfit.degree();
               unsigned num_cpperiodic = num_knots_periodic - kvfit.degree();



               gsMatrix < T > cp_bp(num_cpblade, 2);
               gsMatrix < T > cp_bs(num_cpblade, 2);
               //control points of suction,pressure side
               for (unsigned i = 0; i < num_cpblade; i++) {
                   cp_bp(i, 0) = bp.coef(i, 0);
                   cp_bp(i, 1) = bp.coef(i, 1);
                   cp_bs(i, 0) = bs.coef(i, 0);
                   cp_bs(i, 1) = bs.coef(i, 1);
               }

               int cpk = num_cpblade-1; //index of last coef for cp profile

               T ystart_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x1 + cp_bs(cpk, 1)*length_x1) / (
                   cp_bs(0, 0) - cp_bs(cpk, 0)));
               T yend_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x2 + cp_bs(cpk, 1)*length_x2) / (
                   cp_bs(0, 0) - cp_bs(cpk, 0)));
               T xperiodic_coor =(1-boundary_param)*length_x1 + boundary_param*length_x2;
               T yperiodicsuc_coor = (1-boundary_param)*(ystart_coor+ fb*pitch) + boundary_param*(yend_coor+fb*pitch);
               T yperiodicpress_coor = (1-boundary_param)*(ystart_coor- (1-fb)*pitch) + boundary_param*(yend_coor-(1-fb)*pitch);

               gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
               gsMatrix<T> par(1,1);

               //inlet

               gsKnotVector<T> kvfit3 = kvinlet;
               std::vector<T> knots_for_kvfit3(kvfit3.size() + kvfit3.degree() + (kvfit3.size() - 2 * kvfit3.degree() - 2));
               int j = 0;
               for (index_t i = kvfit3.size()-1; i > kvfit3.degree(); i--) { knots_for_kvfit3[j] = math::abs(1 - kvfit3[i])/2; j++; gsInfo << math::abs(1 - kvfit3[i])/2 << "\n"; }
               for (unsigned i = 1; i < kvfit3.size(); i++) { knots_for_kvfit3[j] = 0.5 + kvfit3[i]/2; j++; gsInfo << 0.5 + kvfit3[i]/2 << "\n"; }
               kvfit3 = gsKnotVector<T> (knots_for_kvfit3);
               gsInfo << kvfit3 << "\n";


               gsMatrix<T> inlet_suction_south_coef(kvfit3.size()-kvfit3.degree()-1, 2);

               T nPointsBefore = 0;
               while ((bs.knots())[nPointsBefore] < blade_param_new) { nPointsBefore++; }
               gsInfo << nPointsBefore << "\n";

               j = 0;
               for (index_t i = nPointsBefore-1; i > 0; i--) {
                   inlet_suction_south_coef.row(j) = bp.coef(i);
                   j++;
               }
               for (index_t i = 0; i < nPointsBefore; i++) {
                   inlet_suction_south_coef.row(j) = bs.coef(i);
                   j++;
               }

               gsMatrix<T> inlet_suction_north_coef(2, 2);
               gsMatrix<T> inlet_suction_west_coef(4, 2);
               gsMatrix<T> inlet_suction_east_coef(4, 2);
               //gsMatrix<T> inlet_suction_south_coef(2*(num_cpinput - 1)-1,2);
               inlet_suction_north_coef.setZero(2, 2);
               inlet_suction_west_coef.setZero(4, 2);
               //inlet_suction_south_coef.setZero(2*(num_cpinput - 1)-1,2);
               inlet_suction_east_coef.setZero(4, 2);

               gsBSpline<T> inlet_suction_south(kvfit3,inlet_suction_south_coef);

               inlet_suction_north_coef << length_x1, ystart_coor-(1-fb)*pitch,
               length_x1, ystart_coor+fb*pitch;
               gsBSpline<T> inlet_suction_north(kvlinear,inlet_suction_north_coef);
               inlet_suction_north.degreeElevate(kvfit3.degree()-1);
               for(unsigned i = kvfit3.degree()+1; i < kvfit3.size()- (kvfit3.degree()+1); i++)
               {
               inlet_suction_north.insertKnot(kvfit3.at(i));
               }

               p0 << inlet_suction_south_coef(kvfit3.size()-kvfit3.degree()-2,0), inlet_suction_south_coef(kvfit3.size()-kvfit3.degree()-2,1);
               par << 1.0;
               dir1 << inlet_suction_south.deriv(par).transpose();
               dir2 << -dir1(1), dir1(0);
               dir2 = dir2/dir2.norm();
               dir1=-dir1;
               t0 = 0.6*AxisDirectionNormed(dir1, dir2);
               p3 << length_x1, ystart_coor+fb*pitch;
               dir1 << length_x2-length_x1, yend_coor - ystart_coor;
               //dir2 << cp_bs (num_cpinput-2,0) - length_x1, cp_bs (num_cpinput-2,1) - (ystart_coor+fb*pitch);
               dir2 << 0, -1;
               dir2 = AxisDirectionNormed(dir1, dir2);
               t1 = 0.75*AxisDirectionNormed(dir1, dir2);

               inlet_suction_east_coef =  LSQFergusonShort(p0, t0, t1, p3);

               gsBSpline<T> inlet_suction_east(kvcubic,inlet_suction_east_coef);

               p0 << inlet_suction_south_coef(0,0), inlet_suction_south_coef(0,1);
               par << 1.0;
               dir1 << inlet_suction_south.deriv(par).transpose();
               dir2 << -dir1(1), -dir1(0);
               t0 = 0.5*dir2/dir2.norm();
               p3 << length_x1, ystart_coor-(1-fb)*pitch;
               dir1 << -length_x1+length_x2, -ystart_coor + yend_coor;
               dir1 = dir1/dir1.norm();
               dir2 << 0, 1;
               t1 = 0.75*AxisDirectionNormed(dir1, dir2);
               inlet_suction_west_coef = LSQFergusonShort(p0, t0, t1, p3);
               gsBSpline<T> inlet_suction_west(kvcubic,inlet_suction_west_coef);
               gsInfo << "inlet suction north " << inlet_suction_north.coefs() << " \n";
               gsInfo << "inlet suction south " << inlet_suction_south.coefs() << " \n";

               gsInfo << "inlet suction east " << inlet_suction_east.coefs() << " \n";

               gsInfo << "inlet suction west " << inlet_suction_west.coefs() << " \n";

               gsMatrix<T> inlet_suction_mat((kvfit3.size()-kvfit3.degree()-1)*(kvfit.degree()+1),2);
               inlet_suction_mat.setZero((kvfit3.size()-kvfit3.degree()-1)*(kvfit.degree()+1),2);
               discreteCoonsPatch(inlet_suction_south.coefs(),inlet_suction_north.coefs(),inlet_suction_west.coefs(),inlet_suction_east.coefs(),inlet_suction_mat,true);
               //springModelPatch(inlet_suction_south.coefs(),inlet_suction_north.coefs(),inlet_suction_west.coefs(),inlet_suction_east.coefs(),inlet_suction_mat,true);
               gsTensorBSplineBasis<2, T> inlet_suction_basis(kvfit3, kvuniform);
               gsTensorBSpline<2,T> inlet_suction(inlet_suction_basis,inlet_suction_mat);


              //---------------periodic_suction-------------------------------

               gsMatrix<T> periodic_suction_north_coef(2, 2);
               gsMatrix<T> periodic_suction_west_coef(4, 2);
               gsMatrix<T> periodic_suction_east_coef(4, 2);
               gsMatrix<T> periodic_suction_south_coef(num_cpperiodic - 1, 2);
               periodic_suction_north_coef.setZero(2, 2);
               periodic_suction_west_coef.setZero(4, 2);
               periodic_suction_south_coef.setZero(num_cpperiodic - 1, 2);
               periodic_suction_east_coef.setZero(4, 2);

               for(unsigned int i = num_cpinput - 1; i < num_cpinput-1+num_cpperiodic -1; i++)
               {
                   periodic_suction_south_coef(i - (num_cpinput-1), 0) =  cp_bs(i, 0);
                   periodic_suction_south_coef(i- (num_cpinput-1), 1) =  cp_bs(i, 1);
               }
               periodic_suction_west_coef <<  inlet_suction_east.coefs();

               periodic_suction_north_coef << length_x1, ystart_coor+fb*pitch,
                                           xperiodic_coor, yperiodicsuc_coor;

               gsBSpline<T> periodic_suction_south(kvperiodic,periodic_suction_south_coef);
               gsBSpline<T> periodic_suction_west(kvcubic,periodic_suction_west_coef);


               gsBSpline<T> periodic_suction_north(kvlinear,periodic_suction_north_coef);
               periodic_suction_north.degreeElevate(kvfit.degree()-1);
               for(unsigned int i = kvperiodic.degree()+1; i < kvperiodic.size()- (kvperiodic.degree()+1); i++)
               {
                   periodic_suction_north.insertKnot(kvperiodic.at(i));
               }

               p0 << cp_bs (cpk,0), cp_bs (cpk,1);
               par << 1.0;
               dir1 << bp.deriv(par).transpose();
               dir2 << -dir1(1),dir1(0);
               t0 = dir2/dir2.norm();
               p3 << xperiodic_coor, yperiodicsuc_coor;
               t1 << 0, -1;


               periodic_suction_east_coef =  LSQFergusonShort(p0, t0, t1, p3);

               gsBSpline<T> periodic_suction_east(kvcubic,periodic_suction_east_coef);

               gsInfo << "periodic suction north " << periodic_suction_north.coefs() << " \n";
               gsInfo << "periodic suction south " << periodic_suction_south.coefs() << " \n";

               gsInfo << "periodic suction east " << periodic_suction_east.coefs() << " \n";

               gsInfo << "periodic suction west " << periodic_suction_west.coefs() << " \n";


               gsMatrix<T> periodic_suction_mat((num_cpperiodic - 1)*(kvfit.degree()+1),2);
               discreteCoonsPatch(periodic_suction_south.coefs(),periodic_suction_north.coefs(),periodic_suction_west.coefs(),periodic_suction_east.coefs(),periodic_suction_mat,false);
               gsTensorBSplineBasis<2, T> periodic_suction_basis(kvperiodic, kvuniform);
               gsTensorBSpline<2,T> periodic_suction(periodic_suction_basis,periodic_suction_mat);




               //---------------outlet_suction-------------------------------

               gsMatrix<T> outlet_suction_north_coef(2, 2);
               gsMatrix<T> outlet_suction_west_coef(4, 2);
               gsMatrix<T> outlet_suction_east_coef(2, 2);
               gsMatrix<T> outlet_suction_south_coef(4, 2);
               outlet_suction_north_coef.setZero(2, 2);
               outlet_suction_west_coef.setZero(4, 2);
               outlet_suction_south_coef.setZero(4, 2);
               outlet_suction_east_coef.setZero(2, 2);

               p0 << cp_bs(cpk,0), cp_bs(cpk,1);
               par << 1.0;
               dir1 << bp.deriv(par).transpose();
               dir2 << bs.deriv(par).transpose();
                t0 =0.5*AxisDirectionNormed(dir1, dir2);
                p3 << length_x2, yend_coor;
               t1 << length_x1-length_x2, ystart_coor - yend_coor;

               outlet_suction_south_coef = LSQFergusonShort(p0, t0, t1, p3);



               outlet_suction_west_coef <<  periodic_suction_east.coefs();
               outlet_suction_east_coef <<  length_x2, yend_coor,
                                            length_x2, yend_coor + fb*pitch;
               outlet_suction_north_coef <<  xperiodic_coor, yperiodicsuc_coor,
                       length_x2, yend_coor + fb*pitch;

               gsBSpline<T> outlet_suction_south(kvcubic,outlet_suction_south_coef);

               gsBSpline<T> outlet_suction_west(kvcubic,outlet_suction_west_coef);

               gsBSpline<T> outlet_suction_east(kvlinear,outlet_suction_east_coef);
               outlet_suction_east.degreeElevate(kvfit.degree()-1);
               gsBSpline<T> outlet_suction_north(kvlinear,outlet_suction_north_coef);
               outlet_suction_north.degreeElevate(kvfit.degree()-1);

               gsMatrix<T> outlet_suction_mat((kvfit.degree()+1)*(kvfit.degree()+1),2);
               springModelPatch(outlet_suction_south.coefs(),outlet_suction_north.coefs(),outlet_suction_west.coefs(),outlet_suction_east.coefs(),outlet_suction_mat,false);
               gsTensorBSplineBasis<2, T> outlet_suction_basis(kvuniform, kvuniform);
               gsTensorBSpline<2,T> outlet_suction(outlet_suction_basis, outlet_suction_mat);


               //---------------periodic_pressure-------------------------------

               /*gsMatrix<T> periodic_pressure_south_coef(num_cpperiodic - 1, 2);
               gsMatrix<T> periodic_pressure_west_coef(4, 2);
               gsMatrix<T> periodic_pressure_east_coef(4, 2);
               gsMatrix<T> periodic_pressure_north_coef(2, 2);
               periodic_pressure_south_coef.setZero(num_cpperiodic - 1, 2);
               periodic_pressure_west_coef.setZero(4, 2);
               periodic_pressure_north_coef.setZero(2, 2);
               periodic_pressure_east_coef.setZero(4, 2);

               for(unsigned int i = num_cpinput - 1; i < num_cpinput-1+num_cpperiodic -1; i++)
               {
                   periodic_pressure_south_coef(num_cpperiodic - 2 - i + num_cpinput - 1, 0) =  cp_bp(i, 0);
                   periodic_pressure_south_coef(num_cpperiodic - 2 - i + num_cpinput - 1, 1) =  cp_bp(i, 1);
               }
               gsKnotVector<T> kvperiodicrev = kvperiodic;
               kvperiodicrev.reverse();
               gsBSpline<T> periodic_pressure_south(kvperiodicrev,periodic_pressure_south_coef);

                periodic_pressure_east_coef = inlet_suction_west.coefs();
                gsBSpline<T> periodic_pressure_east(kvcubic,periodic_pressure_east_coef);

                gsInfo << "periodic pressure south " << periodic_pressure_south.coefs() << " \n";

                gsInfo << "periodic pressure east " << periodic_pressure_east.coefs() << " \n";

               periodic_pressure_north_coef << xperiodic_coor, yperiodicpress_coor,
               length_x1, ystart_coor-(1-fb)*pitch;
               gsBSpline<T> periodic_pressure_north(kvlinear,periodic_pressure_north_coef);
               periodic_pressure_north.degreeElevate(kvfit.degree()-1);
               for(unsigned i = kvperiodic.degree()+1; i < kvperiodic.size()- (kvperiodic.degree()+1); i++)
               {
                   gsInfo << kvperiodic.at(i) << "\n";
                   periodic_pressure_north.insertKnot(kvperiodicrev.at(i));
               }

               gsInfo << "periodic pressure north " << periodic_pressure_north.coefs() << " \n";


               p0 << cp_bp(cpk,0), cp_bp(cpk,1);
               par << 1.0;
               dir1 << bs.deriv(par).transpose();
               dir2 << dir1(0),-dir1(1);
               t0 << dir2/dir2.norm();
               p3 << xperiodic_coor, yperiodicpress_coor;
               t1 << 0, 0.8;

               periodic_pressure_west_coef =  LSQFergusonShort(p0, t0, t1, p3);
               gsBSpline<T> periodic_pressure_west(kvcubic,periodic_pressure_west_coef);




               gsInfo << "periodic pressure west " << periodic_pressure_west.coefs() << " \n";

               gsMatrix<T> periodic_pressure_mat((num_cpperiodic - 1)*(kvfit.degree()+1),2);
               discreteCoonsPatch(periodic_pressure_south.coefs(),periodic_pressure_north.coefs(),periodic_pressure_west.coefs(),periodic_pressure_east.coefs(),periodic_pressure_mat,false);
               gsTensorBSplineBasis<2, T> periodic_pressure_basis(kvperiodicrev, kvuniform);
               gsTensorBSpline<2,T> periodic_pressure(periodic_pressure_basis,periodic_pressure_mat);
               */
               gsMatrix<T> periodic_pressure_south_coef(num_cpperiodic - 1, 2);
               gsMatrix<T> periodic_pressure_west_coef(4, 2);
               gsMatrix<T> periodic_pressure_east_coef(4, 2);
               gsMatrix<T> periodic_pressure_north_coef(2, 2);
               periodic_pressure_south_coef.setZero(num_cpperiodic - 1, 2);
               periodic_pressure_west_coef.setZero(4, 2);
               periodic_pressure_north_coef.setZero(2, 2);
               periodic_pressure_east_coef.setZero(4, 2);

               for(unsigned int i = num_cpinput - 1; i < num_cpinput-1+num_cpperiodic -1; i++)
               {
                   periodic_pressure_south_coef(i - num_cpinput + 1, 0) =  cp_bp(i, 0);
                   periodic_pressure_south_coef(i - num_cpinput + 1, 1) =  cp_bp(i, 1);
               }
               gsKnotVector<T> kvperiodicrev = kvperiodic;
               //kvperiodicrev.reverse();
               gsBSpline<T> periodic_pressure_south(kvperiodicrev,periodic_pressure_south_coef);

                periodic_pressure_east_coef = inlet_suction_west.coefs();
                gsBSpline<T> periodic_pressure_east(kvcubic,periodic_pressure_east_coef);

                gsInfo << "periodic pressure south " << periodic_pressure_south.coefs() << " \n";

                gsInfo << "periodic pressure east " << periodic_pressure_east.coefs() << " \n";

               periodic_pressure_north_coef << length_x1, ystart_coor-(1-fb)*pitch,
               xperiodic_coor, yperiodicpress_coor;
               gsBSpline<T> periodic_pressure_north(kvlinear,periodic_pressure_north_coef);
               periodic_pressure_north.degreeElevate(kvfit.degree()-1);
               for(unsigned i = kvperiodic.degree()+1; i < kvperiodic.size()- (kvperiodic.degree()+1); i++)
               {
                   gsInfo << kvperiodic.at(i) << "\n";
                   periodic_pressure_north.insertKnot(kvperiodicrev.at(i));
               }

               gsInfo << "periodic pressure north " << periodic_pressure_north.coefs() << " \n";


               p0 << cp_bp(cpk,0), cp_bp(cpk,1);
               par << 1.0;
               dir1 << bs.deriv(par).transpose();
               dir2 << dir1(0),-dir1(1);
               t0 << dir2/dir2.norm();
               p3 << xperiodic_coor, yperiodicpress_coor;
               t1 << 0, 0.8;

               periodic_pressure_west_coef =  LSQFergusonShort(p0, t0, t1, p3);
               gsBSpline<T> periodic_pressure_west(kvcubic,periodic_pressure_west_coef);




               gsInfo << "periodic pressure west " << periodic_pressure_west.coefs() << " \n";

               gsMatrix<T> periodic_pressure_mat((num_cpperiodic - 1)*(kvfit.degree()+1),2);
               //discreteCoonsPatch(periodic_pressure_south.coefs(),periodic_pressure_north.coefs(),periodic_pressure_west.coefs(),periodic_pressure_east.coefs(),periodic_pressure_mat,false);
               discreteCoonsPatch(periodic_pressure_east.coefs(), periodic_pressure_west.coefs(), periodic_pressure_south.coefs(),periodic_pressure_north.coefs(),periodic_pressure_mat,false);
               gsTensorBSplineBasis<2, T> periodic_pressure_basis(kvuniform, kvperiodicrev);
               gsTensorBSpline<2,T> periodic_pressure(periodic_pressure_basis,periodic_pressure_mat);


               //---------------outlet_pressure-------------------------------

               gsMatrix<T> outlet_pressure_north_coef(2, 2);
               gsMatrix<T> outlet_pressure_west_coef(4, 2);
               gsMatrix<T> outlet_pressure_east_coef(2, 2);
               gsMatrix<T> outlet_pressure_south_coef(4, 2);
               outlet_pressure_north_coef.setZero(2, 2);
               outlet_pressure_west_coef.setZero(4, 2);
               outlet_pressure_south_coef.setZero(4, 2);
               outlet_pressure_east_coef.setZero(2, 2);

               outlet_pressure_south_coef << outlet_suction_south.coefs();
               outlet_pressure_west_coef <<  periodic_pressure_west.coefs();
               outlet_pressure_east_coef <<  length_x2, yend_coor, length_x2, yend_coor - (1-fb)*pitch;

               outlet_pressure_north_coef <<  xperiodic_coor, yperiodicpress_coor,
                       length_x2, yend_coor - (1-fb)*pitch;

               gsBSpline<T> outlet_pressure_south(kvcubic,outlet_pressure_south_coef);

               gsBSpline<T> outlet_pressure_west(kvcubic,outlet_pressure_west_coef);

               gsBSpline<T> outlet_pressure_east(kvlinear,outlet_pressure_east_coef);
               outlet_pressure_east.degreeElevate(kvfit.degree()-1);
               gsBSpline<T> outlet_pressure_north(kvlinear,outlet_pressure_north_coef);
               outlet_pressure_north.degreeElevate(kvfit.degree()-1);

               gsInfo << "outlet pressure north " << outlet_pressure_north.coefs() << " \n";
               gsInfo << "outlet pressure south " << outlet_pressure_south.coefs() << " \n";

               gsInfo << "outlet pressure east " << outlet_pressure_east.coefs() << " \n";

               gsInfo << "outlet pressure west " << outlet_pressure_west.coefs() << " \n";

               gsMatrix<T> outlet_pressure_mat((kvfit.degree()+1)*(kvfit.degree()+1),2);

               discreteCoonsPatch(outlet_pressure_west.coefs(),outlet_pressure_east.coefs(),outlet_pressure_south.coefs(),outlet_pressure_north.coefs(),outlet_pressure_mat,false);
               gsTensorBSplineBasis<2, T> outlet_pressure_basis(kvuniform, kvuniform);
               gsTensorBSpline<2,T> outlet_pressure(outlet_pressure_basis, outlet_pressure_mat);


               if(!coarse)
              {

                  std::vector<real_t> inlet_suction_kv0_unique = inlet_suction.knots(0).unique();
                  std::vector<real_t> inlet_suction_inserted_knots = smartKnotIdentification2(inlet_suction_kv0_unique);

                  for (size_t i = 0; i < inlet_suction_inserted_knots.size(); i++) {
                      inlet_suction.insertKnot(inlet_suction_inserted_knots[i], 0);
                  }


                  std::vector<real_t> periodic_pressure_kv0_unique = periodic_pressure.knots(1).unique();
                  std::vector<real_t> periodic_pressure_inserted_knots = smartKnotIdentification2(periodic_pressure_kv0_unique);

                  for (size_t i = 0; i < periodic_pressure_inserted_knots.size(); i++) {
                      periodic_pressure.insertKnot(periodic_pressure_inserted_knots[i], 1);

                  }


                  std::vector<real_t> periodic_suction_kv0_unique = periodic_suction.knots(0).unique();
                  std::vector<real_t> periodic_suction_inserted_knots = smartKnotIdentification2(periodic_suction_kv0_unique);

                  for (size_t i = 0; i < periodic_suction_inserted_knots.size(); i++) {
                      periodic_suction.insertKnot(periodic_suction_inserted_knots[i], 0);

                  }

                   gsInfo << periodic_suction.knots(0) << "\n";
                   gsInfo << periodic_suction.knots(1) << "\n";

                   outlet_suction.insertKnot(0.25, 0);
                   outlet_suction.insertKnot(0.5, 0);
                   outlet_suction.insertKnot(0.75, 0);
                   outlet_pressure.insertKnot(0.25, 0);
                   outlet_pressure.insertKnot(0.5, 0);
                    outlet_pressure.insertKnot(0.75, 0);


                   inlet_suction.insertKnot(0.25, 1);
                   inlet_suction.insertKnot(0.5, 1);
                   inlet_suction.insertKnot(0.75, 1);
                   outlet_suction.insertKnot(0.25, 1);
                   outlet_suction.insertKnot(0.5, 1);
                   outlet_suction.insertKnot(0.75, 1);
                   periodic_pressure.insertKnot(0.25, 0);
                   periodic_pressure.insertKnot(0.5, 0);
                   periodic_pressure.insertKnot(0.75, 0);
                   periodic_suction.insertKnot(0.25, 1);
                   periodic_suction.insertKnot(0.5, 1);
                   periodic_suction.insertKnot(0.75, 1);
                   outlet_pressure.insertKnot(0.25, 1);
                   outlet_pressure.insertKnot(0.5, 1);
                   outlet_pressure.insertKnot(0.75, 1);

                   inlet_suction.insertKnot(0.15, 0);
                   inlet_suction.insertKnot(0.85, 0);
              }

                   mpFinal.addPatch(inlet_suction);
                   mpFinal.addPatch(periodic_suction);
                   mpFinal.addPatch(outlet_suction);
                   mpFinal.addPatch(periodic_pressure);
                   mpFinal.addPatch(outlet_pressure);

                   //=================================optimization===========================================

                   gsInfo << "Optimalizace 4. platu (mezi offsety) ...\n";
                   gsInfo << "----------------------------------------\n";

                   T orthogonality = 0.0;
                   T skewness = 0.005;
                   T eccentricity = 0.0;
                   T intersection = 0.0;
                   T uniformity = 0.01;
                   T area = 1-skewness-uniformity;
                   T length = 0.0;
                   T epsilon = 1e-7;

                   gsQualityMeasure<T> optimization(mpFinal.patch(0));
                   //    T opt_val = optimization.functional(orthogonality, skewness,
                   //                                             eccentricity, uniformity,
                   //                                             length, area,
                   //                                             intersection, epsilon);
                   optimization.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);
                   gsQualityMeasure<T> optimization2(mpFinal.patch(1));
                   optimization2.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);

                   //    gsInfo << "Value of functional: "
                   //           << opt_val
                   //           << "\n";
                   gsInfo << "Konec optimalizace.\n";

                   mpFinal.addInterface(0, boundary::east, 1, boundary::west);
                   mpFinal.addInterface(1, boundary::east, 2, boundary::west);
                   mpFinal.addInterface(0, boundary::west, 3, boundary::south);
                   mpFinal.addInterface(3, boundary::north, 4, boundary::south);
                   // mpFinal.addInterface(0, boundary::west, 3, boundary::east);
                   mpFinal.addInterface(2, boundary::south, 4, boundary::west);

                   //periodic
                   mpFinal.addInterface(1, boundary::north, 3, boundary::east);
                   mpFinal.addInterface(2, boundary::north, 4, boundary::east);

               mpFinal.addAutoBoundaries();

               gsWriteParaview(mpFinal,"patch",10000,true,true);
               return mpFinal;
           }

        gsMultiPatch<T> DomainAroundBladeProfile4_Linear(T const & index,
            T const & length_x1,
            T const & length_x2,
            T const & pitch,
            T const & camberX,
            T const & camberY,
            T const & leadingAngle,
            T const & trailingAngle,
            T const & thicknessX,
            T const & thicknessY,
            T const & endingOffset,
            T const & outputAngle,
            T const & radius,
            T const & chordLength,
            T const & Angle,
            T const & rotationCenterX,
            T const & rotationCenterY,
            T const & uniformity_param,
            std::vector<T> const & kvfit_knots,
            bool const & coarse,
            gsVector<T> const & geom_Params)
        {
            T blade_param = 0.01;
            T boundary_param = 0.6;
            T fb = 1-(0.5 + camberY + thicknessY);
            //T fb = 0.5;
             //----------------knot vector seTing----------------
            //std::vector<T> kvfit_knots = {0.,0.,0.,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.,1.,1.};
            //std::vector<T> kvfit_knots = {0.,0.,0.,0.,0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617,1.,1.,1.,1.};
            gsKnotVector<T> kvfit(kvfit_knots);

            //gsKnotVector<T> kvfit(0,1,4,4);
            gsInfo << "degree of knot: " << kvfit.degree() << "\n";

             //----------------set parameters for blade profile----------------
            //bool plot = false;
            int num_samples = 30;
            gsVector<T> vec(2);
            //gsInfo << pitch << "\n ";
            vec(0) = rotationCenterX;
            vec(1) = rotationCenterY;
            gsMatrix<T> mat(2, 2);
            mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
                chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
            //gsInfo << vec << "\n \n";
            //gsInfo << mat << "\n \n";
            gsBSpline<T> suction_side_curve;
            gsBSpline<T> pressure_side_curve;
            gsBSpline<T> suction_side_curve_transf;
            gsBSpline<T> pressure_side_curve_transf;


            //gsKnotVector<T> kvfit(0,1,4,4);

            BladeProfile<T> * pBladeProfile = 0;
            //unsigned degree = 3;
            //---------------compute blade profile for given parameters----------------------
            pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX,
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
            pBladeProfile->setSuctionSide(suction_side_curve_transf);
            pBladeProfile->setPressureSide(pressure_side_curve_transf);

            gsBSpline < T > bs = pBladeProfile->getPressureSide();
            gsBSpline < T > bp = pBladeProfile->getSuctionSide();



    //        gsWriteParaview(bs,"blade_suction_spline");
    //        gsWriteParaview(bp,"blade_pressure_spline");


            GISMO_ASSERT(kvfit.inDomain(blade_param), "blade_param is not in the parametric domain");

            std::vector<real_t>::const_iterator span;
            span = kvfit.iFind(blade_param);

            //insert blade_param as the boundary of span or in the middle
            T middle_span = (span[0]+span[1])/2.;
            std::vector<T> diff_span(3);
            diff_span[0] = math::abs(span[0]-blade_param);
            diff_span[1] = math::abs(span[1]-blade_param);
            diff_span[2] = math::abs(middle_span-blade_param);

            T blade_param_new;

            auto min_diff_span = std::min_element( diff_span.begin(), diff_span.end());
            if (*min_diff_span == diff_span[2] || (span[0]==0. || span[1]==1.))
            {
                blade_param_new = middle_span;
                bs.insertKnot(middle_span,kvfit.degree()+1);
                bp.insertKnot(middle_span,kvfit.degree()+1);
                kvfit.insert(middle_span,kvfit.degree()+1);
            }
            else
            {
                blade_param_new = span[std::distance(diff_span.begin(), min_diff_span)];
                bs.insertKnot(blade_param_new,kvfit.degree());
                bp.insertKnot(blade_param_new,kvfit.degree());
                kvfit.insert(blade_param_new,kvfit.degree());
            }




            //---------------knot vectors of patches-----------------------

            gsKnotVector<T> kvuniform(0,1,0,kvfit.degree()+1);
            gsKnotVector<T> kvlinear(0,1,0,2);
            gsKnotVector<T> kvinlet = kvfit;
            gsKnotVector<T> kvperiodic = kvfit;

    //        gsInfo << "kvinlet" << kvinlet << "\n";
    //        gsInfo << "kvperiodic" << kvperiodic << "\n";

            unsigned num_knots_input = kvfit.degree()+1;
            unsigned num_knots_periodic = 0;

            unsigned k = 0;

            while (blade_param_new != kvfit.at(k))
            {
                num_knots_input++;
                k++;
            }
            num_knots_periodic = kvfit.size() - num_knots_input + kvfit.degree()+1;

    //        gsInfo << "num_knots_input" << num_knots_input << "\n";
    //        gsInfo << "num_knots_periodic" << num_knots_periodic << "\n";

            kvinlet.trimRight(kvfit.size() - num_knots_input);
            kvinlet.transform(0,1);
            kvperiodic.trimLeft(kvfit.size() - num_knots_periodic);
            kvperiodic.transform(0,1);

    //        gsInfo << "kvinlet" << kvinlet << "\n";
    //        gsInfo << "kvperiodic" << kvperiodic << "\n";

            //---------------some boundary coefs of patches-----------------------
            gsMultiPatch<T> mpFinal;

            unsigned num_cpblade = bs.coefsSize();
            unsigned num_cpinput = num_knots_input - kvfit.degree();
            unsigned num_cpperiodic = num_knots_periodic - kvfit.degree();



            gsMatrix < T > cp_bp(num_cpblade, 2);
            gsMatrix < T > cp_bs(num_cpblade, 2);
            //control points of suction,pressure side
            for (unsigned i = 0; i < num_cpblade; i++) {
                cp_bp(i, 0) = bp.coef(i, 0);
                cp_bp(i, 1) = bp.coef(i, 1);
                cp_bs(i, 0) = bs.coef(i, 0);
                cp_bs(i, 1) = bs.coef(i, 1);
            }

    //        gsInfo << "cp_bp" << cp_bp << "\n";
    //        gsInfo << "cp_bs" << cp_bs << "\n";

            int cpk = num_cpblade-1; //index of last coef for cp profile

            T ystart_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x1 + cp_bs(cpk, 1)*length_x1) / (
                cp_bs(0, 0) - cp_bs(cpk, 0)));
            T yend_coor = -((cp_bs(0, 1)*cp_bs(cpk, 0) - cp_bs(0, 0)*cp_bs(cpk, 1) - cp_bs(0, 1)*length_x2 + cp_bs(cpk, 1)*length_x2) / (
                cp_bs(0, 0) - cp_bs(cpk, 0)));
            T xperiodic_coor =(1-boundary_param)*length_x1 + boundary_param*length_x2;
            T yperiodicsuc_coor = (1-boundary_param)*(ystart_coor+ fb*pitch) + boundary_param*(yend_coor+fb*pitch);
            T yperiodicpress_coor = (1-boundary_param)*(ystart_coor- (1-fb)*pitch) + boundary_param*(yend_coor-(1-fb)*pitch);

            //---------------inlet_suction-------------------------------

            gsMatrix<T> inlet_suction_north_coef(2, 2);
            gsMatrix<T> inlet_suction_west_coef(2, 2);
            gsMatrix<T> inlet_suction_east_coef(2, 2);
            gsMatrix<T> inlet_suction_south_coef(num_cpinput - 1, 2);
            inlet_suction_north_coef.setZero(2, 2);
            inlet_suction_west_coef.setZero(2, 2);
            inlet_suction_south_coef.setZero(num_cpinput - 1, 2);
            inlet_suction_east_coef.setZero(2, 2);

            for(unsigned i = 0; i < num_cpinput - 1; i++)
            {
                inlet_suction_south_coef(i, 0) =  cp_bs(i, 0);
                inlet_suction_south_coef(i, 1) =  cp_bs(i, 1);
            }
            inlet_suction_west_coef <<  cp_bs (0,0), cp_bs (0,1),
                                        length_x1, ystart_coor;
            inlet_suction_east_coef <<  cp_bs (num_cpinput-1,0), cp_bs (num_cpinput-1,1),
                                        length_x1, ystart_coor+fb*pitch;
            inlet_suction_north_coef << length_x1, ystart_coor,
                                        length_x1, ystart_coor+fb*pitch;

            gsBSpline<T> inlet_suction_south(kvinlet,inlet_suction_south_coef);
            gsBSpline<T> inlet_suction_west(kvlinear,inlet_suction_west_coef);
            inlet_suction_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> inlet_suction_east(kvlinear,inlet_suction_east_coef);
            inlet_suction_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> inlet_suction_north(kvlinear,inlet_suction_north_coef);
            inlet_suction_north.degreeElevate(kvfit.degree()-1);
            for(unsigned i = kvinlet.degree()+1; i < kvinlet.size()- (kvinlet.degree()+1); i++)
            {
                inlet_suction_north.insertKnot(kvinlet.at(i));
            }


    //        gsInfo << "inlet suction north " << inlet_suction_north.coefs() << " \n";
    //        gsInfo << "inlet suction south " << inlet_suction_south.coefs() << " \n";

    //        gsInfo << "inlet suction east " << inlet_suction_east.coefs() << " \n";

    //        gsInfo << "inlet suction west " << inlet_suction_west.coefs() << " \n";

            gsMatrix<T> inlet_suction_mat((num_cpinput - 1)*(kvfit.degree()+1),2);
            inlet_suction_mat.setZero((num_cpinput - 1)*(kvfit.degree()+1),2);
            discreteCoonsPatch(inlet_suction_south.coefs(),inlet_suction_north.coefs(),inlet_suction_west.coefs(),inlet_suction_east.coefs(),inlet_suction_mat,true);
            gsTensorBSplineBasis<2, T> inlet_suction_basis(kvinlet, kvuniform);
            gsTensorBSpline<2,T> inlet_suction(inlet_suction_basis,inlet_suction_mat);


        //    gsMultiPatch<T> boundaries_inlet_suction;
        //    boundaries_inlet_suction.addPatch(inlet_suction_west);
        //    boundaries_inlet_suction.addPatch(inlet_suction_south);
        //    boundaries_inlet_suction.addPatch(inlet_suction_east);
        //    boundaries_inlet_suction.addPatch(inlet_suction_north);
        //    gsCoonsPatch<T> patch_inlet_suction = coonsPatch(boundaries_inlet_suction);
        //    patch_inlet_suction.compute();

            mpFinal.addPatch(inlet_suction);

          //gsInfo << "inlet computed \n";

            //---------------periodic_suction-------------------------------

            gsMatrix<T> periodic_suction_north_coef(2, 2);
            gsMatrix<T> periodic_suction_west_coef(2, 2);
            gsMatrix<T> periodic_suction_east_coef(2, 2);
            gsMatrix<T> periodic_suction_south_coef(num_cpperiodic - 1, 2);
            periodic_suction_north_coef.setZero(2, 2);
            periodic_suction_west_coef.setZero(2, 2);
            periodic_suction_south_coef.setZero(num_cpperiodic - 1, 2);
            periodic_suction_east_coef.setZero(2, 2);

            for(unsigned int i = num_cpinput - 1; i < num_cpinput-1+num_cpperiodic -1; i++)
            {
                periodic_suction_south_coef(i - (num_cpinput-1), 0) =  cp_bs(i, 0);
                periodic_suction_south_coef(i- (num_cpinput-1), 1) =  cp_bs(i, 1);
            }
            periodic_suction_west_coef <<  inlet_suction_east_coef;
            periodic_suction_east_coef <<  cp_bs (cpk,0), cp_bs (cpk,1),
                                        xperiodic_coor, yperiodicsuc_coor;
            periodic_suction_north_coef << length_x1, ystart_coor+fb*pitch,
                                        xperiodic_coor, yperiodicsuc_coor;

            gsBSpline<T> periodic_suction_south(kvperiodic,periodic_suction_south_coef);
            gsBSpline<T> periodic_suction_west(kvlinear,periodic_suction_west_coef);
            periodic_suction_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> periodic_suction_east(kvlinear,periodic_suction_east_coef);
            periodic_suction_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> periodic_suction_north(kvlinear,periodic_suction_north_coef);
            periodic_suction_north.degreeElevate(kvfit.degree()-1);
            for(unsigned int i = kvperiodic.degree()+1; i < kvperiodic.size()- (kvperiodic.degree()+1); i++)
            {
                periodic_suction_north.insertKnot(kvperiodic.at(i));
            }

    //        gsInfo << "periodic suction north " << periodic_suction_north.coefs() << " \n";
    //        gsInfo << "periodic suction south " << periodic_suction_south.coefs() << " \n";

    //        gsInfo << "periodic suction east " << periodic_suction_east.coefs() << " \n";

    //        gsInfo << "periodic suction west " << periodic_suction_west.coefs() << " \n";


        //    gsMultiPatch<T> boundaries_periodic_suction;
        //    boundaries_periodic_suction.addPatch(periodic_suction_west);
        //    boundaries_periodic_suction.addPatch(periodic_suction_south);
        //    boundaries_periodic_suction.addPatch(periodic_suction_east);
        //    boundaries_periodic_suction.addPatch(periodic_suction_north);
        //    gsCoonsPatch<T> patch_periodic_suction = coonsPatch(boundaries_periodic_suction);
        //    patch_periodic_suction.compute();

            gsMatrix<T> periodic_suction_mat((num_cpperiodic - 1)*(kvfit.degree()+1),2);
            discreteCoonsPatch(periodic_suction_south.coefs(),periodic_suction_north.coefs(),periodic_suction_west.coefs(),periodic_suction_east.coefs(),periodic_suction_mat,false);
            gsTensorBSplineBasis<2, T> periodic_suction_basis(kvperiodic, kvuniform);
            gsTensorBSpline<2,T> periodic_suction(periodic_suction_basis,periodic_suction_mat);

            mpFinal.addPatch(periodic_suction);


            //---------------outlet_suction-------------------------------

            gsMatrix<T> outlet_suction_north_coef(2, 2);
            gsMatrix<T> outlet_suction_west_coef(2, 2);
            gsMatrix<T> outlet_suction_east_coef(2, 2);
            gsMatrix<T> outlet_suction_south_coef(2, 2);
            outlet_suction_north_coef.setZero(2, 2);
            outlet_suction_west_coef.setZero(2, 2);
            outlet_suction_south_coef.setZero(2, 2);
            outlet_suction_east_coef.setZero(2, 2);

            outlet_suction_south_coef << cp_bs(cpk,0), cp_bs(cpk,1),
                                        length_x2, yend_coor;
            outlet_suction_west_coef <<  periodic_suction_east_coef;
            outlet_suction_east_coef <<  length_x2, yend_coor,
                                         length_x2, yend_coor + fb*pitch;
            outlet_suction_north_coef <<  xperiodic_coor, yperiodicsuc_coor,
                    length_x2, yend_coor + fb*pitch;

            gsBSpline<T> outlet_suction_south(kvlinear,outlet_suction_south_coef);
            outlet_suction_south.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_suction_west(kvlinear,outlet_suction_west_coef);
            outlet_suction_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_suction_east(kvlinear,outlet_suction_east_coef);
            outlet_suction_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_suction_north(kvlinear,outlet_suction_north_coef);
            outlet_suction_north.degreeElevate(kvfit.degree()-1);

        //    gsMultiPatch<T> boundaries_outlet_suction;
        //    boundaries_outlet_suction.addPatch(outlet_suction_west);
        //    boundaries_outlet_suction.addPatch(outlet_suction_south);
        //    boundaries_outlet_suction.addPatch(outlet_suction_east);
        //    boundaries_outlet_suction.addPatch(outlet_suction_north);
        //    gsCoonsPatch<T> patch_outlet_suction = coonsPatch(boundaries_outlet_suction);
        //    patch_outlet_suction.compute();

            gsMatrix<T> outlet_suction_mat((kvfit.degree()+1)*(kvfit.degree()+1),2);
            discreteCoonsPatch(outlet_suction_south.coefs(),outlet_suction_north.coefs(),outlet_suction_west.coefs(),outlet_suction_east.coefs(),outlet_suction_mat,false);
            gsTensorBSplineBasis<2, T> outlet_suction_basis(kvuniform, kvuniform);
            gsTensorBSpline<2,T> outlet_suction(outlet_suction_basis, outlet_suction_mat);


            mpFinal.addPatch(outlet_suction);


            //---------------inlet_pressure-------------------------------

            //gsInfo << "kvinlet" << kvinlet << "\n";
            gsMatrix<T> inlet_pressure_north_coef(num_cpinput - 1, 2);
            gsMatrix<T> inlet_pressure_west_coef(2, 2);
            gsMatrix<T> inlet_pressure_east_coef(2, 2);
            gsMatrix<T> inlet_pressure_south_coef(2, 2);
            inlet_pressure_north_coef.setZero(num_cpinput - 1, 2);
            inlet_pressure_west_coef.setZero(2, 2);
            inlet_pressure_south_coef.setZero(2, 2);
            inlet_pressure_east_coef.setZero(2, 2);

            for(unsigned i = 0; i < num_cpinput - 1; i++)
            {
                inlet_pressure_north_coef(i , 0) =  cp_bp(i, 0);
                inlet_pressure_north_coef(i, 1) =  cp_bp(i, 1);
            }
            inlet_pressure_west_coef <<  length_x1, ystart_coor,
                                        cp_bp (0,0), cp_bp (0,1),
            inlet_pressure_east_coef <<  length_x1, ystart_coor-(1-fb)*pitch,
                    cp_bp (num_cpinput-1,0), cp_bp (num_cpinput-1,1);
                                        ;
            inlet_pressure_south_coef << length_x1, ystart_coor,
                                    length_x1, ystart_coor-(1-fb)*pitch;


            gsBSpline<T> inlet_pressure_west(kvlinear,inlet_pressure_west_coef);
            inlet_pressure_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> inlet_pressure_east(kvlinear,inlet_pressure_east_coef);
            inlet_pressure_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> inlet_pressure_north(kvinlet,inlet_pressure_north_coef);

            gsBSpline<T> inlet_pressure_south(kvlinear,inlet_pressure_south_coef);
            inlet_pressure_south.degreeElevate(kvfit.degree()-1);

            for(unsigned i = kvinlet.degree()+1; i < kvinlet.size()- (kvinlet.degree()+1); i++)
            {
                inlet_pressure_south.insertKnot(kvinlet.at(i));
            }

    //        gsInfo << "inlet pressure north " << inlet_pressure_north.coefs() << " \n";
    //        gsInfo << "inlet pressure south " << inlet_pressure_south.coefs() << " \n";

    //        gsInfo << "inlet pressure east " << inlet_pressure_east.coefs() << " \n";

    //        gsInfo << "inlet pressure west " << inlet_pressure_west.coefs() << " \n";

    //        gsInfo << inlet_pressure_west.basis();
    //        gsInfo << inlet_pressure_south.basis();
    //        gsInfo << inlet_pressure_north.basis();
    //        gsInfo << inlet_pressure_east.basis();




            gsMatrix<T> inlet_pressure_mat((num_cpinput - 1)*(kvfit.degree()+1),2);
            discreteCoonsPatch(inlet_pressure_south.coefs(),inlet_pressure_north.coefs(),inlet_pressure_west.coefs(),inlet_pressure_east.coefs(),inlet_pressure_mat,false);
            gsTensorBSplineBasis<2, T> inlet_pressure_basis(kvinlet, kvuniform);
            gsTensorBSpline<2,T> inlet_pressure(inlet_pressure_basis,inlet_pressure_mat);

            mpFinal.addPatch(inlet_pressure);



     //     gsInfo << "inlet computed \n";

            //---------------periodic_pressure-------------------------------

            gsMatrix<T> periodic_pressure_north_coef(num_cpperiodic - 1, 2);
            gsMatrix<T> periodic_pressure_west_coef(2, 2);
            gsMatrix<T> periodic_pressure_east_coef(2, 2);
            gsMatrix<T> periodic_pressure_south_coef(2, 2);
            periodic_pressure_north_coef.setZero(num_cpperiodic - 1, 2);
            periodic_pressure_west_coef.setZero(2, 2);
            periodic_pressure_south_coef.setZero(2, 2);
            periodic_pressure_east_coef.setZero(2, 2);

            for(unsigned int i = num_cpinput - 1; i < num_cpinput-1+num_cpperiodic -1; i++)
            {
                periodic_pressure_north_coef(i - (num_cpinput-1), 0) =  cp_bp(i, 0);
                periodic_pressure_north_coef(i- (num_cpinput-1), 1) =  cp_bp(i, 1);
            }
            periodic_pressure_west_coef <<  inlet_pressure_east_coef;
            periodic_pressure_east_coef <<  xperiodic_coor, yperiodicpress_coor,
                                            cp_bp (cpk,0), cp_bp (cpk,1);
            periodic_pressure_south_coef << length_x1, ystart_coor-(1-fb)*pitch,
                                            xperiodic_coor, yperiodicpress_coor;

            gsBSpline<T> periodic_pressure_north(kvperiodic,periodic_pressure_north_coef);
            gsBSpline<T> periodic_pressure_west(kvlinear,periodic_pressure_west_coef);
            periodic_pressure_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> periodic_pressure_east(kvlinear,periodic_pressure_east_coef);
            periodic_pressure_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> periodic_pressure_south(kvlinear,periodic_pressure_south_coef);
            periodic_pressure_south.degreeElevate(kvfit.degree()-1);
            for(unsigned int i = kvperiodic.degree()+1; i < kvperiodic.size()- (kvperiodic.degree()+1); i++)
            {
                periodic_pressure_south.insertKnot(kvperiodic.at(i));
            }

    //        gsInfo << "periodic pressure north " << periodic_pressure_north.coefs() << " \n";
    //        gsInfo << "periodic pressure south " << periodic_pressure_south.coefs() << " \n";

    //        gsInfo << "periodic pressure east " << periodic_pressure_east.coefs() << " \n";

    //        gsInfo << "periodic pressure west " << periodic_pressure_west.coefs() << " \n";


        //    gsMultiPatch<T> boundaries_periodic_pressure;
        //    boundaries_periodic_pressure.addPatch(periodic_pressure_west);
        //    boundaries_periodic_pressure.addPatch(periodic_pressure_south);
        //    boundaries_periodic_pressure.addPatch(periodic_pressure_east);
        //    boundaries_periodic_pressure.addPatch(periodic_pressure_north);
        //    gsCoonsPatch<T> patch_periodic_pressure = coonsPatch(boundaries_periodic_pressure);
        //    patch_periodic_pressure.compute();

            gsMatrix<T> periodic_pressure_mat((num_cpperiodic - 1)*(kvfit.degree()+1),2);
            discreteCoonsPatch(periodic_pressure_south.coefs(),periodic_pressure_north.coefs(),periodic_pressure_west.coefs(),periodic_pressure_east.coefs(),periodic_pressure_mat,false);
            gsTensorBSplineBasis<2, T> periodic_pressure_basis(kvperiodic, kvuniform);
            gsTensorBSpline<2,T> periodic_pressure(periodic_pressure_basis,periodic_pressure_mat);


            mpFinal.addPatch(periodic_pressure);


            //---------------outlet_pressure-------------------------------

            gsMatrix<T> outlet_pressure_north_coef(2, 2);
            gsMatrix<T> outlet_pressure_west_coef(2, 2);
            gsMatrix<T> outlet_pressure_east_coef(2, 2);
            gsMatrix<T> outlet_pressure_south_coef(2, 2);
            outlet_pressure_north_coef.setZero(2, 2);
            outlet_pressure_west_coef.setZero(2, 2);
            outlet_pressure_south_coef.setZero(2, 2);
            outlet_pressure_east_coef.setZero(2, 2);

            outlet_pressure_north_coef << cp_bp(cpk,0), cp_bp(cpk,1),
                                        length_x2, yend_coor;
            outlet_pressure_west_coef <<  periodic_pressure_east_coef;
            outlet_pressure_east_coef <<  length_x2, yend_coor - (1-fb)*pitch,
                    length_x2, yend_coor;
            outlet_pressure_south_coef <<  xperiodic_coor, yperiodicpress_coor,
                    length_x2, yend_coor - (1-fb)*pitch;

            gsBSpline<T> outlet_pressure_south(kvlinear,outlet_pressure_south_coef);
            outlet_pressure_south.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_pressure_west(kvlinear,outlet_pressure_west_coef);
            outlet_pressure_west.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_pressure_east(kvlinear,outlet_pressure_east_coef);
            outlet_pressure_east.degreeElevate(kvfit.degree()-1);
            gsBSpline<T> outlet_pressure_north(kvlinear,outlet_pressure_north_coef);
            outlet_pressure_north.degreeElevate(kvfit.degree()-1);


        //    gsMultiPatch<T> boundaries_outlet_pressure;
        //    boundaries_outlet_pressure.addPatch(outlet_pressure_west);
        //    boundaries_outlet_pressure.addPatch(outlet_pressure_south);
        //    boundaries_outlet_pressure.addPatch(outlet_pressure_east);
        //    boundaries_outlet_pressure.addPatch(outlet_pressure_north);
        //    gsCoonsPatch<T> patch_outlet_pressure = coonsPatch(boundaries_outlet_pressure);
        //    patch_outlet_pressure.compute();

            gsMatrix<T> outlet_pressure_mat((kvfit.degree()+1)*(kvfit.degree()+1),2);

            discreteCoonsPatch(outlet_pressure_south.coefs(),outlet_pressure_north.coefs(),outlet_pressure_west.coefs(),outlet_pressure_east.coefs(),outlet_pressure_mat,false);
            gsTensorBSplineBasis<2, T> outlet_pressure_basis(kvuniform, kvuniform);
            gsTensorBSpline<2,T> outlet_pressure(outlet_pressure_basis, outlet_pressure_mat);

            mpFinal.addPatch(outlet_pressure);







            mpFinal.addInterface(0, boundary::east, 1, boundary::west);
            mpFinal.addInterface(1, boundary::east, 2, boundary::west);
             mpFinal.addInterface(3, boundary::east, 4, boundary::west);
              mpFinal.addInterface(4, boundary::east, 5, boundary::west);
               mpFinal.addInterface(0, boundary::west, 3, boundary::east);
               mpFinal.addInterface(2, boundary::south, 5, boundary::north);

               //periodic
               mpFinal.addInterface(1, boundary::north, 4, boundary::south);
               mpFinal.addInterface(2, boundary::north, 5, boundary::south);

            mpFinal.addAutoBoundaries();


            return mpFinal;
        }

    gsMultiPatch<T> DomainBetweenBladeProfiles5(T const & index,
                                                T const & length_x1,
                                                T const & length_x2,
                                                T const & pitch,
                                                T const & camberX,
                                                T const & camberY,
                                                T const & leadingAngle,
                                                T const & trailingAngle,
                                                T const & thicknessX,
                                                T const & thicknessY,
                                                T const & endingOffset,
                                                T const & outputAngle,
                                                T const & radius,
                                                T const & chordLength,
                                                T const & Angle,
                                                T const & rotationCenterX,
                                                T const & rotationCenterY,
                                                T const & uniformity_param,
                                                std::vector<T> const & kvfit_knots,
                                                bool const & coarse,
                                                gsVector<T> const & geom_Params) {

        //----------------set parameters for blade profile----------------
        T offset_distance = geom_Params(0);
        T wake_width = geom_Params(1);
        //gsKnotVector<T> kvfit(0, 1, 4, 4);
        //std::vector<T> kvfit_knots = {0.0, 0.0, 0.0, 0.0, 0.0125335, 0.02542396, 0.1074354, 0.1894468, 0.3920851, 0.5947234, 0.7973617, 1.0, 1.0, 1.0, 1.0};
        gsKnotVector<T> kvfit = gsKnotVector<T> (kvfit_knots);
        gsKnotVector<T> kvcub(0, 1, 0, 4);
        gsKnotVector<T> kvlin(0, 1, 0, 2);

        bool plot = true;
        //bool plotMeshes = true;
        int num_samples = 100;
        gsVector<T> vec(2);
        //gsInfo << pitch << "\n ";
        vec(0) = rotationCenterX;
        vec(1) = rotationCenterY;
        gsMatrix<T> mat(2, 2);
        mat << chordLength * math::cos(Angle), chordLength * math::sin(Angle),
               chordLength * math::sin(Angle), -chordLength * math::cos(Angle);
        gsBSpline<T> suction_side_curve;
        gsBSpline<T> pressure_side_curve;
        gsBSpline<T> suction_side_offset_curve;
        gsBSpline<T> pressure_side_offset_curve;
        BladeProfile<T> * pBladeProfile = 0;

        //---------------compute blade profile for given parameters----------------------
        pBladeProfile = new BladeProfile<T>(camberX, camberY, leadingAngle, trailingAngle, thicknessX, thicknessY, endingOffset, outputAngle, radius, chordLength, Angle, rotationCenterX, rotationCenterY, 0.0);
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

        gsMatrix<T> suction_side_cp = suction_side_curve.coefs();
        gsMatrix<T> pressure_side_cp = pressure_side_curve.coefs();
        gsMatrix<T> suction_side_offset_cp = suction_side_offset_curve.coefs();
        gsMatrix<T> pressure_side_offset_cp = pressure_side_offset_curve.coefs();

        // --------------- cross-section curves ---------------------------------------------------------------------------------------------------
        gsMatrix<T> aux_cp(2, 2);
        aux_cp << suction_side_offset_cp.row(0),
                  suction_side_cp.row(0);
        gsBSpline<T> cs_curve1 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve1.degreeElevate(2);
        gsMatrix<T> cs_curve1_cp = cs_curve1.coefs();
        aux_cp << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1),
                  suction_side_cp.row(suction_side_curve.coefsSize()-1);
        gsBSpline<T> cs_curve2 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve2.degreeElevate(2);
        gsMatrix<T> cs_curve2_cp = cs_curve2.coefs();
        aux_cp << pressure_side_offset_cp.row(0),
                  pressure_side_cp.row(0);
        gsBSpline<T> cs_curve3 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve3.degreeElevate(2);
        gsMatrix<T> cs_curve3_cp = cs_curve3.coefs();
        aux_cp << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1),
                  pressure_side_cp.row(pressure_side_curve.coefsSize()-1);
        gsBSpline<T> cs_curve4 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve4.degreeElevate(2);
        gsMatrix<T> cs_curve4_cp = cs_curve4.coefs();

        // ---------------- outer boundary -------------------------------------------------------------------------------------------------------
        gsMatrix<T> cp_bs = suction_side_curve.coefs();
        T ystart_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x1 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x1) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));
        T yend_coor = -((cp_bs(0, 1)*cp_bs(suction_side_curve.coefsSize()-1, 0) - cp_bs(0, 0)*cp_bs(suction_side_curve.coefsSize()-1, 1) - cp_bs(0, 1)*length_x2 + cp_bs(suction_side_curve.coefsSize()-1, 1)*length_x2) / (cp_bs(0, 0) - cp_bs(suction_side_curve.coefsSize()-1, 0)));

        aux_cp << length_x1, ystart_coor - pitch,
                  length_x1, ystart_coor;
        gsBSpline<T> left_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        left_boundary_curve.degreeElevate(2);
        gsMatrix<T> left_boundary_curve_cp = left_boundary_curve.coefs();

        /*aux_cp << length_x2, yend_coor + ft*pitch,
                  length_x2, yend_coor - fb*pitch;
        gsBSpline<T> right_boundary_curve = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve.degreeElevate(2);*/
        gsMatrix<T> rightSplitPoint1(1,2), rightSplitPoint2(1,2);
        T inserted_knot_right = wake_width*offset_distance/pitch;
        rightSplitPoint1(0,0) = length_x2;
        rightSplitPoint1(0,1) = (1 - inserted_knot_right) * (yend_coor - pitch) + inserted_knot_right * yend_coor;
        rightSplitPoint2(0,0) = length_x2;
        rightSplitPoint2(0,1) = inserted_knot_right * (yend_coor - pitch) + (1 - inserted_knot_right) * yend_coor;
        aux_cp << rightSplitPoint1(0,0), rightSplitPoint1(0,1),
                  length_x2, yend_coor - pitch;
        gsBSpline<T> right_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_A_cp = right_boundary_curve_A.coefs();
        aux_cp << rightSplitPoint1(0,0), rightSplitPoint1(0,1),
                  rightSplitPoint2(0,0), rightSplitPoint2(0,1);
        gsBSpline<T> right_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve_B.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_B_cp = right_boundary_curve_B.coefs();
        aux_cp << rightSplitPoint2(0,0), rightSplitPoint2(0,1),
                  length_x2, yend_coor;
        gsBSpline<T> right_boundary_curve_C = gsBSpline<T> ( kvlin, aux_cp );
        right_boundary_curve_C.degreeElevate(2);
        gsMatrix<T> right_boundary_curve_C_cp = right_boundary_curve_C.coefs();

        aux_cp << length_x1, ystart_coor,
                  suction_side_offset_cp(0,0), suction_side_offset_cp(0,1);
        gsBSpline<T> top_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_A_cp = top_boundary_curve_A.coefs();
        aux_cp << suction_side_cp(suction_side_curve.coefsSize()-1,0), suction_side_cp(suction_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor;
        gsBSpline<T> top_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        top_boundary_curve_B.degreeElevate(2);
        gsMatrix<T> top_boundary_curve_B_cp = top_boundary_curve_B.coefs();

        aux_cp << length_x1, ystart_coor - pitch,
                  pressure_side_offset_cp(0,0), pressure_side_offset_cp(0,1);
        gsBSpline<T> bottom_boundary_curve_A = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_A.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_A_cp = bottom_boundary_curve_A.coefs();
        aux_cp << pressure_side_cp(pressure_side_curve.coefsSize()-1,0), pressure_side_cp(pressure_side_curve.coefsSize()-1,1),
                  length_x2, yend_coor - pitch;
        gsBSpline<T> bottom_boundary_curve_B = gsBSpline<T> ( kvlin, aux_cp );
        bottom_boundary_curve_B.degreeElevate(2);
        gsMatrix<T> bottom_boundary_curve_B_cp = bottom_boundary_curve_B.coefs();

        // ---------------- cross-section curves 2 -------------------------------------------------------------------------------------------------
        gsMatrix<T> p0(1,2), t0(1,2), t1 (1,2), p3(1,2), dir1(1,2), dir2(1,2);
        /*p0 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
        t0 = suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-2); t0 = t0/t0.norm();
        t1 << - length_x2 + length_x1, - yend_coor + ystart_coor; t1 = t1/t1.norm();
        p3 << rightSplitPoint2(0,0), rightSplitPoint2(0,1);
        gsMatrix<T> cs_curve5_cp(4,2);
        cs_curve5_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve5_cp << "\n";
        gsBSpline<T> cs_curve5 = gsBSpline<T> ( kvcub, cs_curve5_cp );*/
        aux_cp << suction_side_offset_cp(suction_side_offset_curve.coefsSize()-1,0), suction_side_offset_cp(suction_side_offset_curve.coefsSize()-1,1),
                  rightSplitPoint2(0,0), rightSplitPoint2(0,1);
        gsBSpline<T> cs_curve5 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve5.degreeElevate(2);
        gsMatrix<T> cs_curve5_cp = cs_curve5.coefs();

        /*p0 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
        t0 = pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-2); t0 = t0/t0.norm();
        t1 << - length_x2 + length_x1, - yend_coor + ystart_coor; t1 = t1/t1.norm();
        p3 << rightSplitPoint1(0,0), rightSplitPoint1(0,1);
        gsMatrix<T> cs_curve6_cp(4,2);
        cs_curve6_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve6_cp << "\n";
        gsBSpline<T> cs_curve6 = gsBSpline<T> ( kvcub, cs_curve6_cp );*/
        aux_cp << pressure_side_offset_cp(pressure_side_offset_curve.coefsSize()-1,0), pressure_side_offset_cp(pressure_side_offset_curve.coefsSize()-1,1),
                  rightSplitPoint1(0,0), rightSplitPoint1(0,1);
        gsBSpline<T> cs_curve6 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve6.degreeElevate(2);
        gsMatrix<T> cs_curve6_cp = cs_curve6.coefs();

        p0 << pressure_side_offset_cp.row(0);
        dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        //dir2 << ystart_coor - yend_coor, length_x2 - length_x1;
        dir2 << pressure_side_offset_cp.row(1) - pressure_side_offset_cp.row(0);
        t0 = AxisDirectionNormed(dir1, dir2);
        dir1 << length_x1 - length_x2, ystart_coor - yend_coor;
        //dir2 << yend_coor - ystart_coor, length_x1 - length_x2;
        dir2 << suction_side_offset_cp.row(1) - suction_side_offset_cp.row(0);
        t1 = AxisDirectionNormed(dir1, dir2);
        //t1 << length_x1 - suction_side_offset_cp(0,0), ystart_coor - pitch - suction_side_offset_cp(0,1); t1 = t1/t1.norm();
        p3 << suction_side_offset_cp.row(0);
        gsMatrix<T> cs_curve7_cp(4,2);
        cs_curve7_cp = LSQFergusonSmooth(p0, t0, t1, p3);
        gsInfo << cs_curve7_cp << "\n";
        gsBSpline<T> cs_curve7 = gsBSpline<T> ( kvcub, cs_curve7_cp );

        p0 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
        //t0 << (suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1))/(suction_side_offset_trimmed_cp.row(suction_side_trimmed.coefsSize()-1) - suction_side_trimmed_cp.row(suction_side_trimmed.coefsSize()-1)).norm();
        //dir1 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - pressure_side_cp.row(pressure_side_curve.coefsSize()-1);
        //dir2 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1);
        //t0 = AxisDirectionNormed(dir1, dir2);
        t0 << (suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1))/((suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1)).norm());
        t0 = 3 * t0 / 4;
        //t1 << - yend_coor + ystart_coor, length_x2 - length_x1; t1 = t1/t1.norm();
        dir1 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1) - suction_side_cp.row(suction_side_curve.coefsSize()-1);
        dir2 << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1) - suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
        t1 = AxisDirectionNormed(dir1, dir2);
        p3 << suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
        gsMatrix<T> cs_curve8_cp(4,2);
        cs_curve8_cp = LSQFergusonShort(p0, t0, t1, p3);
        gsInfo << cs_curve8_cp << "\n";
        gsBSpline<T> cs_curve8 = gsBSpline<T> ( kvcub, cs_curve8_cp );

        /*
        aux_cp << pressure_side_offset_cp.row(pressure_side_offset_curve.coefsSize()-1),
                  suction_side_offset_cp.row(suction_side_offset_curve.coefsSize()-1);
        gsBSpline<T> cs_curve8 = gsBSpline<T> ( kvlin, aux_cp );
        cs_curve8.degreeElevate(2);
        gsMatrix<T> cs_curve8_cp = cs_curve8.coefs();
        */

        // ---------------- construction of patches ------------------------------------------------------------------------------------------------
        gsMatrix<T> patch1_cp(suction_side_cp.rows() * cs_curve1_cp.rows(), 2);
        discreteCoonsPatch(suction_side_offset_cp, suction_side_cp, cs_curve1_cp, cs_curve2_cp, patch1_cp, true);
        gsTensorBSpline<2, T> patch1 = gsTensorBSpline<2, T> (kvfit, kvcub, patch1_cp);

        gsMatrix<T> patch2_cp(pressure_side_cp.rows() * cs_curve3_cp.rows(), 2);
        discreteCoonsPatch(cs_curve3_cp, cs_curve4_cp, pressure_side_offset_cp, pressure_side_cp, patch2_cp, true);
        gsTensorBSpline<2, T> patch2 = gsTensorBSpline<2, T> (kvcub, kvfit, patch2_cp);

        gsMatrix<T> patch3_cp(top_boundary_curve_A_cp.rows() * left_boundary_curve_cp.rows(), 2);
        discreteCoonsPatch(bottom_boundary_curve_A_cp, top_boundary_curve_A_cp, left_boundary_curve_cp, cs_curve7_cp, patch3_cp, true);
        gsTensorBSpline<2, T> patch3 = gsTensorBSpline<2, T> (kvcub, kvcub, patch3_cp);

        //gsMatrix<T> patch4_cp(suction_side_offset_cp.rows() * cs_curve7_cp.rows(), 2);
        //discreteCoonsPatch(pressure_side_offset_cp, suction_side_offset_cp, cs_curve7_cp, cs_curve8_cp, patch4_cp, true);
        //gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvfit, kvcub, patch4_cp);
        gsMatrix<T> patch4_cp(suction_side_offset_cp.rows() * cs_curve7_cp.rows(), 2);
        gsMatrix<T> patch4_cp1(suction_side_offset_cp.rows() * cs_curve7_cp.rows(), 2);
        gsMatrix<T> patch4_cp2(suction_side_offset_cp.rows() * cs_curve7_cp.rows(), 2);
        springModelPatch(pressure_side_offset_cp, suction_side_offset_cp, cs_curve7_cp, cs_curve8_cp, patch4_cp1, true);
        discreteCoonsPatch(pressure_side_offset_cp, suction_side_offset_cp, cs_curve7_cp, cs_curve8_cp, patch4_cp2, true);
        patch4_cp = (patch4_cp1 + patch4_cp2) / 2;
        //gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvfit2, kvcub, patch4_cp);
        gsTensorBSpline<2, T> patch4 = gsTensorBSpline<2, T> (kvfit, kvcub, patch4_cp);


        gsMatrix<T> patch5_cp(cs_curve5_cp.rows() * cs_curve8_cp.rows(), 2);
        discreteCoonsPatch(cs_curve6_cp, cs_curve5_cp, cs_curve8_cp, right_boundary_curve_B_cp, patch5_cp, true);
        gsTensorBSpline<2, T> patch5 = gsTensorBSpline<2, T> (kvcub, kvcub, patch5_cp);

        gsMatrix<T> patch6_cp(cs_curve5_cp.rows() * cs_curve2_cp.rows(), 2);
        discreteCoonsPatch(cs_curve5_cp, top_boundary_curve_B_cp, cs_curve2_cp, right_boundary_curve_C_cp, patch6_cp, true);
        gsTensorBSpline<2, T> patch6 = gsTensorBSpline<2, T> (kvcub, kvcub, patch6_cp);

        gsMatrix<T> patch7_cp(cs_curve6_cp.rows() * cs_curve4_cp.rows(), 2);
        discreteCoonsPatch(cs_curve4_cp, right_boundary_curve_A_cp, cs_curve6_cp, bottom_boundary_curve_B_cp, patch7_cp, true);
        gsTensorBSpline<2, T> patch7 = gsTensorBSpline<2, T> (kvcub, kvcub, patch7_cp);

        // ---------------- construction of multipatch ---------------------------------------------------------------------------------------------
        if(!coarse)
       {
           gsInfo << "Special refinement. \n";

           gsInfo << patch1.knots(0) << "\n";
           gsInfo << patch1.knots(1) << "\n";
           std::vector<real_t> patch1_kv0_unique = patch1.knots(0).unique();
           for (size_t i = 0; i < patch1_kv0_unique.size(); i++)
               gsInfo << patch1_kv0_unique[i] << ",";
           gsInfo << "\n";

           std::vector<real_t> patch1_inserted_knots = smartKnotIdentification(patch1_kv0_unique);

           for (size_t i = 0; i < patch1_inserted_knots.size(); i++) {
               patch1.insertKnot(patch1_inserted_knots[i], 0);
               patch2.insertKnot(patch1_inserted_knots[i], 1);
               patch4.insertKnot(patch1_inserted_knots[i], 0);
           }

           // patch 4 + patch 3 + patch 5
           patch4.insertKnot(0.33, 1);
           patch4.insertKnot(0.66, 1);
           patch3.insertKnot(0.33, 1);
           patch3.insertKnot(0.66, 1);
           patch5.insertKnot(0.33, 1);
           patch5.insertKnot(0.66, 1);

           // patch 3
           //patch3.insertKnot(0.5, 0);
           std::vector<real_t> patch3_kv0_unique_aux = patch3.knots(0).unique();
           std::vector<real_t> patch3_kv0_unique(patch3.knots(0).unique().size() + 1);
           patch3_kv0_unique[0] = patch3_kv0_unique_aux[0];
           patch3_kv0_unique[1] = patch1_kv0_unique[1];
           for (size_t i = 1; i < patch3_kv0_unique_aux.size(); i++)
               patch3_kv0_unique[i+1] = patch3_kv0_unique_aux[i];
           std::vector<real_t> patch3_inserted_knots = smartKnotIdentification(patch3_kv0_unique);
           for (size_t i = 0; i < patch3_inserted_knots.size(); i++) {
               patch3.insertKnot(1-patch3_inserted_knots[i], 0);
           }

           // patch 5 + patch 6 + patch 7
           patch5.insertKnot(0.25, 0);
           patch5.insertKnot(0.5, 0);
           patch5.insertKnot(0.75, 0);
           patch6.insertKnot(0.25, 0);
           patch6.insertKnot(0.5, 0);
           patch6.insertKnot(0.75, 0);
           patch7.insertKnot(0.25, 1);
           patch7.insertKnot(0.5, 1);
           patch7.insertKnot(0.75, 1);
       }


        gsMultiPatch<T> mpFinal;
        mpFinal.addPatch(patch1);
        mpFinal.addPatch(patch2);
        mpFinal.addPatch(patch3);
        mpFinal.addPatch(patch4);
        mpFinal.addPatch(patch5);
        mpFinal.addPatch(patch6);
        mpFinal.addPatch(patch7);

        //=================================optimization===========================================

        gsInfo << "Optimalizace 4. platu (mezi offsety) ...\n";
        gsInfo << "----------------------------------------\n";

        T orthogonality = 0.05;
        T skewness = 0.0;
        T eccentricity = 0.0;
        T intersection = 0.0;
        T uniformity = 0.01;
        T area = 1-skewness-uniformity;
        T length = 0.0;
        T epsilon = 1e-7;

        gsQualityMeasure<T> optimization(mpFinal.patch(3));
        //    T opt_val = optimization.functional(orthogonality, skewness,
        //                                             eccentricity, uniformity,
        //                                             length, area,
        //                                             intersection, epsilon);
        optimization.optimize(orthogonality, skewness, eccentricity, uniformity, length, area, intersection, epsilon);


        //    gsInfo << "Value of functional: "
        //           << opt_val
        //           << "\n";
        gsInfo << "Konec optimalizace.\n";

        mpFinal.addInterface(2, boundary::east, 3, boundary::west);
        mpFinal.addInterface(3, boundary::south, 1, boundary::west);
        mpFinal.addInterface(3, boundary::north, 0, boundary::south);
        mpFinal.addInterface(3, boundary::east, 4, boundary::west);
        mpFinal.addInterface(4, boundary::south, 6, boundary::west);
        mpFinal.addInterface(4, boundary::north, 5, boundary::south);
        mpFinal.addInterface(0, boundary::east, 5, boundary::west);
        mpFinal.addInterface(1, boundary::north, 6, boundary::south);
        // periodic interfaces
        mpFinal.addInterface(2, boundary::south, 2, boundary::north);
        mpFinal.addInterface(0, boundary::west, 1, boundary::south);
        mpFinal.addInterface(5, boundary::north, 6, boundary::east);
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

            //mpFinal.uniformRefine(); mpFinal.uniformRefine();
            gsWriteParaview( mpFinal, "patches", 50000, true);
        }

        return mpFinal;
    }

    // cubic B-spline approximation of a sector of an annular cylinder (single patch)
    // the annulus sector (in yz plane) is described such that its center lies in the origin and it is symmetrical about y axis
    // r1, r2 - radii of the annulus
    // L - length of the patch in the x-direction
    // phi - the angle of the sector
    // xPos - position in the x-direction
    
    gsTensorBSpline<3, T> cubicBSplineAnnulusSector3D_patch(const T r1, const T r2, const T L, const T phi, const T xPos = 0.0)
    {
        gsKnotVector<T> kv(0, 1, 0, 4); // first, last, inter, mult_end

        T R = r2 - r1; // width of the annulus
        T phi1 = ((EIGEN_PI - phi) / 2);
        T sin1 = math::sin(phi1);
        T cos1 = math::cos(phi1);

        // z-coordinates of the "outer" control points (right boundary)
        T z0 = cos1 * r1;
        T z1 = cos1 * (r1 + (1. / 3)*R);
        T z2 = cos1 * (r2 - (1. / 3)*R);
        T z3 = cos1 * r2;

        // y-coordinates of the "outer" control points right boundary
        T y0 = sin1 * r1;
        T y1 = sin1 * (r1 + (1. / 3)*R);
        T y2 = sin1 * (r2 - (1. / 3)*R);
        T y3 = sin1 * r2;

        // z-coordinates of the "inner" control points
        T z00 = z0 - (4. / 3)*math::cos(phi / 2)*math::tan(phi / 4)*r1; // cubic B-spline approximation such that the midpoint lies on the circle
        T z11 = z1 - (4. / 3)*math::cos(phi / 2)*math::tan(phi / 4)*(r1 + (1. / 3)*R);
        T z22 = z2 - (4. / 3)*math::cos(phi / 2)*math::tan(phi / 4)*(r2 - (1. / 3)*R);
        T z33 = z3 - (4. / 3)*math::cos(phi / 2)*math::tan(phi / 4)*r2;

        // y-coordinates of the "inner" control points
        T y00 = y0 + (4. / 3)*math::sin(phi / 2)*math::tan(phi / 4)*r1;
        T y11 = y1 + (4. / 3)*math::sin(phi / 2)*math::tan(phi / 4)*(r1 + (1. / 3)*R);
        T y22 = y2 + (4. / 3)*math::sin(phi / 2)*math::tan(phi / 4)*(r2 - (1. / 3)*R);
        T y33 = y3 + (4. / 3)*math::sin(phi / 2)*math::tan(phi / 4)*r2;

        gsMatrix<T> coef(64, 3);
        coef << xPos, y0, -z0,
            xPos, y1, -z1,
            xPos, y2, -z2,
            xPos, y3, -z3,
            xPos, y00, -z00,
            xPos, y11, -z11,
            xPos, y22, -z22,
            xPos, y33, -z33,
            xPos, y00, z00,
            xPos, y11, z11,
            xPos, y22, z22,
            xPos, y33, z33,
            xPos, y0, z0,
            xPos, y1, z1,
            xPos, y2, z2,
            xPos, y3, z3,

            xPos + L / 3, y0, -z0,
            xPos + L / 3, y1, -z1,
            xPos + L / 3, y2, -z2,
            xPos + L / 3, y3, -z3,
            xPos + L / 3, y00, -z00,
            xPos + L / 3, y11, -z11,
            xPos + L / 3, y22, -z22,
            xPos + L / 3, y33, -z33,
            xPos + L / 3, y00, z00,
            xPos + L / 3, y11, z11,
            xPos + L / 3, y22, z22,
            xPos + L / 3, y33, z33,
            xPos + L / 3, y0, z0,
            xPos + L / 3, y1, z1,
            xPos + L / 3, y2, z2,
            xPos + L / 3, y3, z3,

            xPos + (2. / 3)*L, y0, -z0,
            xPos + (2. / 3)*L, y1, -z1,
            xPos + (2. / 3)*L, y2, -z2,
            xPos + (2. / 3)*L, y3, -z3,
            xPos + (2. / 3)*L, y00, -z00,
            xPos + (2. / 3)*L, y11, -z11,
            xPos + (2. / 3)*L, y22, -z22,
            xPos + (2. / 3)*L, y33, -z33,
            xPos + (2. / 3)*L, y00, z00,
            xPos + (2. / 3)*L, y11, z11,
            xPos + (2. / 3)*L, y22, z22,
            xPos + (2. / 3)*L, y33, z33,
            xPos + (2. / 3)*L, y0, z0,
            xPos + (2. / 3)*L, y1, z1,
            xPos + (2. / 3)*L, y2, z2,
            xPos + (2. / 3)*L, y3, z3,

            xPos + L, y0, -z0,
            xPos + L, y1, -z1,
            xPos + L, y2, -z2,
            xPos + L, y3, -z3,
            xPos + L, y00, -z00,
            xPos + L, y11, -z11,
            xPos + L, y22, -z22,
            xPos + L, y33, -z33,
            xPos + L, y00, z00,
            xPos + L, y11, z11,
            xPos + L, y22, z22,
            xPos + L, y33, z33,
            xPos + L, y0, z0,
            xPos + L, y1, z1,
            xPos + L, y2, z2,
            xPos + L, y3, z3;

        return gsTensorBSpline<3, T>(kv, kv, kv, coef);
    }

    // cubic B-spline approximation of a sector of an annular cylinder
    // the annulus sector (in yz plane) is described such that its center lies in the origin and it is symmetrical about y axis
    // r1, r2 - radii of the annulus
    // L - length of the domain in the x-direction
    // phi - the angle of the sector
    // np - number of patches in the x-direction
    gsMultiPatch<T> cubicBSplineAnnulusSector3D(const T r1, const T r2, const T L, const T phi, const int np = 1)
    {
        gsMultiPatch<T> mp;

        T Lp = L / np;

        for (int i = 0; i < np; i++)
            mp.addPatch(cubicBSplineAnnulusSector3D_patch(r1, r2, Lp, phi, i*Lp));

        return mp;
    }

}; //uwbGeometryCreators

} //namespace gismo
