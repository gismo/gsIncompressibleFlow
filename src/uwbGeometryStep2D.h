#ifdef _MSC_VER // to be removed
#define _USE_MATH_DEFINES
#endif

#ifndef UWBGEOMETRYSTEP2D_H
#define UWBGEOMETRYSTEP2D_H

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <gismo.h>

#include <math.h>
#include <string.h>
#include <sstream>

using namespace gismo;

template<class T> gsTensorBSpline<2, T> BSplineRect(int deg, const T llx, const T lly, const T a, const T b) // llx - lower left x, lly - lower left y
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
    default:
        GISMO_ERROR("Degree not implemented.");
        break;
    }

    return gsTensorBSpline<2, T>(kv, kv, coef);
}


template<class T> gsMultiPatch<T> BSplineStep2D(int deg, T const & a, T const & b, T const & a_in)
{
    gsMultiPatch<T> mp;

    mp.addPatch(BSplineRect(deg, 0.0, 0.0, a, b / 2));
    mp.addPatch(BSplineRect(deg, 0.0, b / 2, a, b / 2));
    mp.addPatch(BSplineRect(deg, -a_in, b / 2, a_in, b / 2));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);

    //mp.addInterface(0, boundary::south, 1, boundary::north);

    mp.addAutoBoundaries();

    return mp;
}

template<class T> void refineBasis(gsMultiBasis<T>& basis, int numRefine, int numRefineU, int numRefineLocal, real_t addRefPart, int dim, real_t a, real_t b)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    int bRefine = math::floor((a / b) - 1) + numRefineU;

    const gsTensorBSplineBasis<2, T>*  basis0 = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&(basis.basis(0))); //basis of the patch 0

    gsMatrix<> box(2, 2);
    box << 0, 1, 0, 0;
    for (int i = 0; i < bRefine; i++){
        basis.refine(0, box);
        basis.refine(1, box);
    }

    box << 0, addRefPart, 0, 0;
    basis.refine(0, box);
    basis.refine(1, box);

    for (int i = 0; i < numRefineLocal; i++){
        int sizeKnots_v = basis0->knots(1).size() - 1;
        real_t lastKnot_v = basis0->knot(1, sizeKnots_v);
        real_t vKnot_bottom = basis0->knot(1, basis0->degree(1) + 1);
        real_t vKnot_upper = basis0->knot(1, sizeKnots_v - (basis0->degree(1) + 1));

        box << 0, 0, 0, vKnot_bottom;
        basis.refine(0, box);
        basis.refine(1, box);
        basis.refine(2, box);

        box << 0, 0, vKnot_upper, lastKnot_v;
        basis.refine(0, box);
        basis.refine(1, box);
        basis.refine(2, box);
    }
}

template<class T> void geometryStep2D(std::map<std::string, gsVector<std::string>> parList, std::string outputDIR)
{
    int numRefine, numRefineU, numRefineLocal, addRefPart, deg;
    int dim = 2;
    T a, b, a_in;

    std::string GSF;
    get_parameter(parList,GSF,"geometry_settings_file");
    std::map<std::string, gsVector<std::string>> paramList = readInputFile(GSF);

    get_parameter(paramList,numRefine,"numRefine");
    get_parameter(paramList,numRefineU,"numRefineU");
    get_parameter(paramList,numRefineLocal,"numRefineLocal");
    get_parameter(paramList,deg,"deg");
    get_parameter(paramList,a,"a");
    get_parameter(paramList,a_in,"a_in");
    get_parameter(paramList,b,"b");
    get_parameter(paramList,addRefPart,"addRefPart");

    gsMultiPatch<T> patches;
    patches = BSplineStep2D<T>(deg, a, b, a_in);
    gsMultiBasis<T> tbasis(patches);

    refineBasis(tbasis, numRefine, numRefineU, numRefineLocal, addRefPart, dim, a, b);

    std::string stdeg;
    switch(deg){
    case 1:{
        stdeg = "_lin";
        break;}
    case 2:{
        stdeg = "_kva";
        break;}
    case 3:{
        stdeg = "kub";
        break;}
    default:{
        GISMO_ERROR("Wrong degree!");
        break;}
    }

    std::string fileName = "step2D_ref_" + util::to_string(numRefine) + "_" +
                            util::to_string(numRefineU) + "_" + util::to_string(numRefineLocal) + stdeg;
    gsFileData<> fdpatches;
    fdpatches << patches;
    fdpatches.save(outputDIR + fileName + "_patches.xml");

    gsFileData<> fdtbasis;
    fdtbasis << tbasis;
    fdtbasis.save(outputDIR + fileName + "_tbasis.xml");
 }

#endif
