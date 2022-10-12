/** @file uwbGeometryUtils.h

Author(s): H. Hornikova
*/

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <gismo.h>


using namespace gismo;

// deg - basis degree
// llx - x-coordinate of the lower left corner
// lly - y-coordinate of the lower left corner
// a - width
// b - height
// numSep - number of C^0 separators (for rIgA)

template<class T>
gsTensorBSpline<2, T> BSplineRectangle(int deg, const T llx, const T lly, const T a, const T b, int numSep = 0)
{
    gsKnotVector<T> kv(0, 1, numSep, deg + 1, deg); // first, last, num_inter, mult_end, mult_inter

    int n = kv.size() - deg - 1;
    gsMatrix<T> coef(n*n, 2);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            coef.row(i + j*n) << llx + (i*a) / (deg*(numSep + 1)), lly + (j*b) / (deg*(numSep + 1));

    return gsTensorBSpline<2, T>(kv, kv, coef);
}


template<class T>
gsTensorBSpline<3, T> BSplinePrism(int deg, const T llx, const T lly, const T llz, const T a, const T b, const T c, int numSep = 0)
{
    gsKnotVector<T> kv(0, 1, numSep, deg + 1, deg);

    int n = kv.size() - deg - 1;
    gsMatrix<T> coef(n*n*n, 3);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                coef.row(i + j*n + k*n*n) << llx + (i*a) / (deg*(numSep + 1)), lly + (j*b) / (deg*(numSep + 1)), llz + (k*c) / (deg*(numSep + 1));

    return gsTensorBSpline<3, T>(kv, kv, kv, coef);
}


template<class T> gsMultiPatch<T> BSplineCavity2D(int deg, const T a, const T b, const int np = 1, int numSep = 0)
{
    gsMultiPatch<T> mp;

    T aPatch = a / np;
    T bPatch = b / np;

    for (int j = 0; j < np; j++)
        for (int i = 0; i < np; i++)
            mp.addPatch(BSplineRectangle(deg, i*aPatch, j*bPatch, aPatch, bPatch, numSep));

    mp.computeTopology();

    return mp;
}


template<class T> gsMultiPatch<T> BSplineStep2D(int deg, T const & a, T const & b, T const & a_in, bool periodic = false)
{
    gsMultiPatch<T> mp;

    mp.addPatch(BSplineRectangle(deg, 0.0, 0.0, a, b / 2));
    mp.addPatch(BSplineRectangle(deg, 0.0, b / 2, a, b / 2));
    mp.addPatch(BSplineRectangle(deg, -a_in, b / 2, a_in, b / 2));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);

    if(periodic)
        mp.addInterface(0, boundary::south, 1, boundary::north);

    mp.addAutoBoundaries();

    return mp;
}


template<class T> gsMultiPatch<T> BSplineStep3D(int deg, T const & a, T const & b, T const & c, T const & a_in, bool periodic = false)
{
    gsMultiPatch<T> mp;

    mp.addPatch(BSplinePrism<T>(deg, 0.0, 0.0, 0.0, a, b / 2, c));
    mp.addPatch(BSplinePrism<T>(deg, 0.0, b / 2, 0.0, a, b / 2, c));
    mp.addPatch(BSplinePrism<T>(deg, -a_in, b / 2, 0.0, a_in, b / 2, c));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(2, boundary::east, 1, boundary::west);

    if (periodic)
    {
        mp.addInterface(0, boundary::front, (size_t)0, boundary::back);
        mp.addInterface(1, boundary::front, 1, boundary::back);
        mp.addInterface(2, boundary::front, 2, boundary::back);
    }
    
    mp.addAutoBoundaries();

    return mp;
}


template<class T>
void defineBCs_cavity2D(gsBoundaryConditions<T>& bcInfo, const int np, std::vector<std::pair<int, boxSide> >& bndWall, std::string lidVel = "1")
{
    gsFunctionExpr<> Uwall("0", "0", 2);
    gsFunctionExpr<> Ulid(lidVel, "0", 2);

    for (int i = 1; i <= np; i++)
    {
        bcInfo.addCondition(np*np - i, boundary::north, condition_type::dirichlet, Ulid, 0);
        bndWall.push_back(std::make_pair(np*np - i, boundary::north));
    }

    for (int i = 0; i < np; i++)
    {
        bcInfo.addCondition(i, boundary::south, condition_type::dirichlet, Uwall, 0);
        bndWall.push_back(std::make_pair(i, boundary::south));
    }

    for (int i = 0; i < np; i++)
    {
        bcInfo.addCondition(i * np, boundary::west, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition((i + 1)*np - 1, boundary::east, condition_type::dirichlet, Uwall, 0);
        bndWall.push_back(std::make_pair(i * np, boundary::west));
        bndWall.push_back(std::make_pair((i + 1)*np - 1, boundary::east));
    }
}

template<class T> void addBCs(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndWall, gsFunctionExpr<> Uin, gsFunctionExpr<> Uwall)
{
    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin);

    for (size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall);
}

template<class T> void defineBCs_step2D(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false, std::string inVel = "default")
{
    if (inVel == "default")
        inVel = "(-4*(y-1.5)^2 + 1)";

    gsFunctionExpr<> Uin, Uwall;
    Uin = gsFunctionExpr<>(inVel, "0", 2);
    Uwall = gsFunctionExpr<>("0", "0", 2);

    if (!periodic)
    {
        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::south));
        bndWall.push_back(std::make_pair(1, boundary::north));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));
    }
    else
    {
        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));
    }

    addBCs(bcInfo, bndIn, bndWall, Uin, Uwall);
}


template<class T> void defineBCs_step3D(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false, std::string inVel = "default")
{
    gsFunctionExpr<> Uin, Uwall;

    if (!periodic)
    {
        if (inVel == "default")
            inVel = "(-4*(y-1.5)^2 + 1)*(-(z-1)^2 + 1)";

        Uin = gsFunctionExpr<>(inVel, "0", "0", 3);
        Uwall = gsFunctionExpr<>("0", "0", "0", 3);

        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::south));
        bndWall.push_back(std::make_pair(0, boundary::front));
        bndWall.push_back(std::make_pair(0, boundary::back));
        bndWall.push_back(std::make_pair(1, boundary::north));
        bndWall.push_back(std::make_pair(1, boundary::front));
        bndWall.push_back(std::make_pair(1, boundary::back));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndWall.push_back(std::make_pair(2, boundary::front));
        bndWall.push_back(std::make_pair(2, boundary::back));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));
    }
    else
    {
        if (inVel == "default")
            inVel = "(-4*(y-1.5)^2 + 1)";

        Uin = gsFunctionExpr<>(inVel, "0", "0", 3);
        Uwall = gsFunctionExpr<>("0", "0", "0", 3);

        bndIn.push_back(std::make_pair(2, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::west));
        bndWall.push_back(std::make_pair(0, boundary::south));
        bndWall.push_back(std::make_pair(1, boundary::north));
        bndWall.push_back(std::make_pair(2, boundary::south));
        bndWall.push_back(std::make_pair(2, boundary::north));
        bndOut.push_back(std::make_pair(0, boundary::east));
        bndOut.push_back(std::make_pair(1, boundary::east));
    }

    addBCs(bcInfo, bndIn, bndWall, Uin, Uwall);
}


template<class T> void defineBCs_step(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, int dim, bool periodic = false, std::string inVel = "default")
{
    switch (dim)
    {
    case 2:
        defineBCs_step2D(bcInfo, bndIn, bndOut, bndWall, periodic, inVel);
        break;
    case 3:
        defineBCs_step3D(bcInfo, bndIn, bndOut, bndWall, periodic, inVel);
        break;
    default:
        GISMO_ERROR("Wrong dimension!");
        break;
    }
}


template< int d, class T> void refineFirstKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
{
    const gsTensorBSplineBasis<d, T>*  patchBasis = dynamic_cast<const gsTensorBSplineBasis<d, T>*>(&(basis.basis(patch)));
    gsMatrix<> box(d, 2);

    for (int i = 0; i < numRefine; i++)
    {
        box.setZero();
        T knot = patchBasis->knot(dir, patchBasis->degree(dir) + 1);
        box(dir, 1) = knot;
        basis.refine(patch, box);
    }
}


template< int d, class T> void refineLastKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
{
    const gsTensorBSplineBasis<d, T>*  patchBasis = dynamic_cast<const gsTensorBSplineBasis<d, T>*>(&(basis.basis(patch)));
    gsMatrix<> box(d, 2);

    for (int i = 0; i < numRefine; i++)
    {
        box.setZero();
        int sizeKnots = patchBasis->knots(dir).size() - 1;
        T lastKnot = patchBasis->knot(dir, sizeKnots);
        T knot = patchBasis->knot(dir, sizeKnots - (patchBasis->degree(dir) + 1));
        box(dir, 0) = knot;
        box(dir, 1) = lastKnot;
        basis.refine(patch, box);
    }
}


template<class T> void refineBasis_cavity2D(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, int numSep = 0)
{
    int npDir = math::sqrt(basis.nPieces());
    int npRef = std::log2(npDir);
    int sepRef = std::log2(numSep + 1);

    for (int i = 0; i < numRefine - npRef - sepRef; ++i)
        basis.uniformRefine();

    if (numRefineLocal && npDir == 1)
    {
        for (int dir = 0; dir <= 1; dir++)
        {
            refineFirstKnotSpan<2, T>(basis, numRefineLocal, 0, dir);
            refineLastKnotSpan<2, T>(basis, numRefineLocal, 0, dir);
        }
    }
}


template<int d, class T> void refineLocal_step(gsMultiBasis<T>& basis, int numRefineWalls, int numRefineCorner)
{
    // refinement near upper and bottom walls
    refineFirstKnotSpan<d, T>(basis, numRefineWalls, 0, 1);
    refineLastKnotSpan<d, T>(basis, numRefineWalls, 1, 1);
    refineLastKnotSpan<d, T>(basis, numRefineWalls, 2, 1);

    // refinement near the corner 
    refineFirstKnotSpan<d, T>(basis, numRefineCorner, 0, 0);
    refineFirstKnotSpan<d, T>(basis, numRefineCorner, 1, 0);
    refineLastKnotSpan<d, T>(basis, numRefineCorner, 2, 0);
    refineLastKnotSpan<d, T>(basis, numRefineCorner, 0, 1);
    refineFirstKnotSpan<d, T>(basis, numRefineCorner, 1, 1);
    refineFirstKnotSpan<d, T>(basis, numRefineCorner, 2, 1);
}


template<class T> void refineBasis_step(gsMultiBasis<T>& basis, int numRefine, int numRefineWalls, int numRefineCorner, int numRefineU, real_t addRefPart, int dim, real_t a, real_t b, real_t c = 0.0)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    gsMatrix<> box(dim, 2);

    int uRefine = math::floor(std::log2(a / b)) + 1 + numRefineU;
    box.setZero();
    box(0, 1) = 1; 
    for (int i = 0; i < uRefine; i++)
        for (int p = 0; p < 2; p++)
            basis.refine(p, box);

    if (dim == 3)
    {
        int wRefine = math::floor(std::log2(c / b)) + 1;
        box.setZero();
        box(2, 1) = 1;
        for (int i = 0; i < wRefine; i++)
            for (int p = 0; p < 3; p++)
                basis.refine(p, box);
    }

    box.setZero();
    box(0,1) = addRefPart;
    for (int p = 0; p < 2; p++)
        basis.refine(p, box);

    switch (dim)
    {
    case 2:
        refineLocal_step<2, T>(basis, numRefineWalls, numRefineCorner);
        break;
    case 3:
        refineLocal_step<3, T>(basis, numRefineWalls, numRefineCorner);
        break;
    default:
        GISMO_ERROR("Wrong dimension!");
        break;
    }
}