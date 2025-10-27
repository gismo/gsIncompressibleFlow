/** @file gsFlowUtils.h

    Miscellaneous useful functions for the incompressible flow solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include <gismo.h>

namespace gismo
{

template<class T>
class gsFlowSolverParams;

#ifdef GISMO_WITH_PARDISO
/// @brief Setup the pardiso solver.
/// @param[in,out] solver a reference to the pardiso solver
template <class T>
inline void pardisoSetup(typename gsSparseSolver<T>::PardisoLU& solver)
{
    solver.setParam(7, 15);
    solver.setParam(9, 10);
    solver.setParam(12, 0);
}
#endif


inline void startAnimationFile(std::ofstream& file)
{
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"0.1\">";
    file << "<Collection>\n";
}


inline void endAnimationFile(std::ofstream& file)
{
    file << "</Collection>\n";
    file << "</VTKFile>\n";
    file.close();
}


/// @brief Creates a one-column matrix (vector) of unique values from the input matrix (useful for creating a unique list of active basis functions in several quadrature points).
/// @param[in]  mat  a const reference to the input matrix
template <class T>
gsMatrix<T> createVectorOfUniqueIndices(const gsMatrix<T>& mat)
{
    std::map<T, bool> map;
    std::vector<T> vecTmp;

    for (index_t j = 0; j < mat.cols(); j++)
    {
        for (index_t i = 0; i < mat.rows(); i++)
        {
            if(map.count(mat(i,j)) == 0)
            {
                map[mat(i,j)] = true;
                vecTmp.push_back(mat(i,j));
            }
        }
    }

    gsMatrix<T> vec(vecTmp.size(), 1);
    for (size_t i = 0; i < vecTmp.size(); i++)
        vec(i) = vecTmp[i];


    return vec;
}


/// @brief Get a vector of nonzero entries per outer index (row or column depending on the matrix storage order).
/// @tparam T           real number type
/// @tparam MatOrder    matrix storage order (RowMajor/ColMajor)
/// @param mat[in]      a const reference to the matrix
template <class T, int MatOrder>
gsVector<index_t> getNnzVectorPerOuter(const gsSparseMatrix<T, MatOrder>& mat)
{
    gsVector<index_t> nnzPerOuter(mat.outerSize());
    nnzPerOuter.setZero();

    for (index_t outer = 0; outer < mat.outerSize(); outer++)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(mat, outer); it; ++it)
            ++nnzPerOuter[outer];

    return nnzPerOuter;
}

template <class T, int MatOrder>
int getMaxNnzPerOuter(const gsSparseMatrix<T, MatOrder>& mat)
{
    int maxNnzInOuter = 0;

    for (index_t outer = 0; outer < mat.outerSize(); outer++)
    {
        int nnzInOuter = 0;

        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(mat, outer); it; ++it)
            nnzInOuter++;

        maxNnzInOuter = math::max(maxNnzInOuter, nnzInOuter);
    }

    return maxNnzInOuter;
}



/// @brief Fill a diagonal approximation of an inverse matrix.
/// @tparam T           real number type
/// @tparam MatOrder    matrix storage order (RowMajor/ColMajor)
/// @param[in]  mat     a const reference to the matrix of which the inverse is approximated
/// @param[out] diagInv a reference to the resulting inverse approximation
/// @param[in]  repeat  number of the diagonal block repetition (e.g. for velocity components)
/// @param[in]  lumping use lumping to define the diagonal approximation
template <class T, int MatOrder>
void diagInvMatrix_into(const gsSparseMatrix<T, MatOrder>& mat, gsSparseMatrix<T, MatOrder>& diagInv, int repeat, bool lumping = false)
{
    GISMO_ENSURE(mat.nonZeros() != 0, "diagInvMatrix_into(): The matrix is empty!");

    int varDofs = mat.rows();

    diagInv.resize(repeat*varDofs, repeat*varDofs);
    diagInv.reserve(gsVector<int>::Constant(diagInv.cols(), 1));

    const gsSparseMatrix<T, MatOrder>* matPtr = &mat;
    gsSparseMatrix<T, MatOrder> lumped(varDofs, varDofs);

    if (lumping)
    {
        lumped.reserve(gsVector<int>::Constant(varDofs, 1));

        for (int j = 0; j < varDofs; j++)
            lumped.insert(j, j) = math::abs(mat.at(j, j)); // abs value because of "lumping" diag(A) in SIMPLE-type prec., does not change lumped mass matrix in IgA

        for (int j = 0; j < varDofs; j++)
        {
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(mat, j); it; ++it)
            {
                int i = it.row();

                if (i != j)
                    lumped.coeffRef(i, i) += math::abs(it.value());
            }
        }

        matPtr = &lumped;
    }

    for (int i = 0; i < varDofs; i++)
    {
        T tmp = 1 / matPtr->coeff(i, i);

        for (int s = 0; s < repeat; s++)
            diagInv.coeffRef(i + s * varDofs, i + s * varDofs) = tmp;
    }
}


/// @brief Returns a B-spline parametrization of a rectangle of a given degree in both directions.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param llx      \a x coordinate of the lower left corner
/// @param lly      \a y coordinate of the lower left corner
/// @param a        width of the rectangle
/// @param b        height of the rectangle
/// @param numSep   number of \f C^0 \f separators in the parametrization (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template <class T>
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

// TODO: Leave out?
template<class T>
gsTensorBSpline<2, T> BSplineRectangle2(int deg, const T llx, const T lly, const T a, const T b, int numSep = 0)
{
    gsKnotVector<T> kv1(0, 1, a-1, deg + 1, 1); // first, last, num_inter, mult_end, mult_inter
    gsKnotVector<T> kv2(0, 1, b-1, deg + 1, 1); // first, last, num_inter, mult_end, mult_inter

    int m = kv1.size() - deg - 1;
    int n = kv2.size() - deg - 1;
    gsMatrix<T> coef(m*n, 2);

    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            coef.row(i + j*n) << llx + (i*a) / (deg+a-1), lly + (j*b) / (deg+b-1);

    return gsTensorBSpline<2, T>(kv1, kv2, coef);
}


/// @brief Returns a B-spline parametrization of a 3D block of a given degree in all directions.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param llx      \a x coordinate of a corner
/// @param lly      \a y coordinate of a corner
/// @param llz      \a z coordinate of a corner
/// @param a        width of the block
/// @param b        height of the block
/// @param c        depth of the block
/// @param numSep   number of \f C^0 \f separators in the parametrization (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template <class T>
gsTensorBSpline<3, T> BSplineBlock(int deg, const T llx, const T lly, const T llz, const T a, const T b, const T c, int numSep = 0)
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


/// @brief Returns a B-spline multipatch domain for 2D problems of flow in a cavity.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        width of the cavity
/// @param b        height of the cavity
/// @param np       number of patches in each direction
/// @param numSep   number of \f C^0 \f separators in each patch (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template <class T>
gsMultiPatch<T> BSplineCavity2D(int deg, const T a, const T b, const int np = 1, int numSep = 0)
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


/// @brief Returns a B-spline multipatch domain for 3D problems of flow in a cavity.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        width of the cavity
/// @param b        height of the cavity
/// @param c        depth of the cavity
/// @param numSep   number of \f C^0 \f separators in each patch (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template <class T>
gsMultiPatch<T> BSplineCavity3D(int deg, const T a, const T b, const T c, int numSep = 0)
{
    gsMultiPatch<T> mp;
    mp.addPatch(BSplineBlock(deg, 0.0, 0.0, 0.0, a, b, c, numSep));
    mp.computeTopology();
    return mp;
}


/// @brief Returns a B-spline multipatch domain for 2D problems of flow over a backward facing step.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        length of the domain behind the step
/// @param b        total height of the domain
/// @param a_in     length of the inflow part (before the step)
/// @param h        height of the step
/// @param periodic periodic domain (if true, the bottom and top boundaries of the channel behind step are defined as interface)
template <class T>
gsMultiPatch<T> BSplineStep2D(int deg, const T a, const T b, const T a_in, T h = 0, bool periodic = false)
{
    gsMultiPatch<T> mp;

    if (h == 0)
        h = b / 2;

    mp.addPatch(BSplineRectangle(deg, 0.0, 0.0, a, h));
    mp.addPatch(BSplineRectangle(deg, 0.0, h, a, b - h));
    mp.addPatch(BSplineRectangle(deg, -a_in, h, a_in, b - h));

    // basic refinement
    index_t aNumElem = std::floor(2*a/b);
    index_t ainNumElem = std::floor(2*a_in/b);
    T aStep = 1.0 / aNumElem;
    T ainStep = 1.0 / ainNumElem;

    for (index_t p = 0; p < 2; p++)
        for (index_t i = 1; i < aNumElem; i++)
            mp.patch(p).insertKnot(i * aStep, 0);

    for (index_t i = 1; i < ainNumElem; i++)
        mp.patch(2).insertKnot(i * ainStep, 0);

    // topology
    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);

    if(periodic)
        mp.addInterface(0, boundary::south, 1, boundary::north);

    mp.addAutoBoundaries();

    return mp;
}


/// @brief Returns a B-spline multipatch domain for 3D problems of flow over a backward facing step.
/// @tparam T       real number type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        length of the domain behind the step
/// @param b        total height of the domain
/// @param c        width of the domain
/// @param a_in     length of the inflow part (before the step)
/// @param h        height of the step
/// @param periodic periodic domain (if true, the bottom and top boundaries of the channel behind step are defined as interface)
template <class T>
gsMultiPatch<T> BSplineStep3D(int deg, const T a, const T b, const T c, const T a_in, T h = 0, bool periodic = false)
{
    gsMultiPatch<T> mp;

    if (h == 0)
        h = b / 2;

    mp.addPatch(BSplineBlock<T>(deg, 0.0, 0.0, 0.0, a, h, c));
    mp.addPatch(BSplineBlock<T>(deg, 0.0, h, 0.0, a, b - h, c));
    mp.addPatch(BSplineBlock<T>(deg, -a_in, h, 0.0, a_in, b - h, c));

    // basic refinement
    index_t aNumElem = std::floor(2*a/b);
    index_t ainNumElem = std::floor(2*a_in/b);
    index_t cNumElem = std::floor(2*c/b);
    T aStep = 1.0 / aNumElem;
    T ainStep = 1.0 / ainNumElem;
    T cStep = 1.0 / cNumElem;

    for (index_t p = 0; p < 2; p++)
        for (index_t i = 1; i < aNumElem; i++)
            mp.patch(p).insertKnot(i * aStep, 0);

    for (index_t i = 1; i < ainNumElem; i++)
        mp.patch(2).insertKnot(i * ainStep, 0);

    for (index_t p = 0; p < 3; p++)
        for (index_t i = 1; i < cNumElem; i++)
            mp.patch(p).insertKnot(i * cStep, 2);

    // topology
    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(2, boundary::east, 1, boundary::west);

    if (periodic)
    {
        mp.addInterface(0, boundary::front, 0, boundary::back);
        mp.addInterface(1, boundary::front, 1, boundary::back);
        mp.addInterface(2, boundary::front, 2, boundary::back);
    }
    
    mp.addAutoBoundaries();

    return mp;
}


inline void defineBndParts_step(index_t dim, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false)
{
    switch(dim)
    {
        case 2:
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
            break;

        case 3:
            if (!periodic)
            {
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
                bndIn.push_back(std::make_pair(2, boundary::west));
                bndWall.push_back(std::make_pair(0, boundary::west));
                bndWall.push_back(std::make_pair(0, boundary::south));
                bndWall.push_back(std::make_pair(1, boundary::north));
                bndWall.push_back(std::make_pair(2, boundary::south));
                bndWall.push_back(std::make_pair(2, boundary::north));
                bndOut.push_back(std::make_pair(0, boundary::east));
                bndOut.push_back(std::make_pair(1, boundary::east));
            }
            break;

        default:
            GISMO_ERROR("Wrong dimension!");
    }
}

inline void defineBndParts_cavity(index_t dim, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, const int np = 1)
{
    switch(dim)
    {
        case 2:
            for (int i = 1; i <= np; i++)
                bndWall.push_back(std::make_pair(np*np - i, boundary::north));

            for (int i = 0; i < np; i++)
            {
                bndWall.push_back(std::make_pair(i, boundary::south));
                bndWall.push_back(std::make_pair(i * np, boundary::west));
                bndWall.push_back(std::make_pair((i + 1)*np - 1, boundary::east));
            }
            break;

        case 3:
            GISMO_ASSERT(np == 1, "Number of patches in each direction must be 1 for 3D cavity.");
            bndWall.push_back(std::make_pair(0, boundary::north));
            bndWall.push_back(std::make_pair(0, boundary::south));
            bndWall.push_back(std::make_pair(0, boundary::west));
            bndWall.push_back(std::make_pair(0, boundary::east));
            bndWall.push_back(std::make_pair(0, boundary::front));
            bndWall.push_back(std::make_pair(0, boundary::back));
            break;

        default:
            GISMO_ERROR("Wrong dimension!");
    }
}

inline void defineBndParts_profile2D(std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall)
{
    bndIn.push_back(std::make_pair(0, boundary::west));
    bndWall.push_back(std::make_pair(1, boundary::north));
    bndWall.push_back(std::make_pair(1, boundary::south));
    bndOut.push_back(std::make_pair(2, boundary::east));
}

/// @brief Define which patch sides belong to the inflow/outflow/wall boundary part for predefined domain geometries.
/// @param[in]  geo      geometry ID (1 - backward facing step, 2 - cavity, 3 - blade profile) 
/// @param[out] bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[out] bndOut   reference to a container of patch sides corresponding to outflow boundaries
/// @param[out] bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  periodic periodic domain (true/false, applicable for step)
/// @param[in]  np       number of patches in each direction (applicable for 2D cavity)
inline void defineBndParts(index_t geo, index_t dim, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false, int np = 1)
{
    bndIn.clear();
    bndOut.clear();
    bndWall.clear();

    switch(geo)
    {
        default:
            gsWarn << "Unknown geometry ID in function defineBndParts(), using geo = 1 (backward facing step).";

        case 1:
            defineBndParts_step(dim, bndIn, bndOut, bndWall, periodic);
            break;

        case 2:
            defineBndParts_cavity(dim, bndIn, bndOut, bndWall, np);
            break;

        case 3:
            if (dim == 3)
                gsWarn << "Geometry 3 is only 2D!\n";
            defineBndParts_profile2D(bndIn, bndOut, bndWall);
            break;
    }
}


/// @brief Define boundary conditions for the corresponding boundary parts.
/// @tparam T            real number type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[in]  bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[in]  bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  Uin      the inflow velocity as gsFunctionExpr
/// @param[in]  Uwall    the wall velocity as gsFunctionExpr
template <class T>
void addBCs(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndWall, gsFunctionExpr<T> Uin, gsFunctionExpr<T> Uwall)
{
    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin);

    for (size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall);
}


/// @brief Define boundary conditions for the corresponding boundary parts.
/// @tparam T            real number type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[in]  bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[in]  bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  Uin      the inflow velocity as gsFunctionExpr
/// @param[in]  Uwall    the wall velocity as gsFunctionExpr
/// @param[in]  unk      specifies which unknown variable the boundary condition refers to
template<class T>
void addBCs(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndWall, gsFunctionExpr<T> Uin, gsFunctionExpr<T> Uwall, short_t unk)
{
    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin, unk);

    for (size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall, unk);
}


/// @brief Define boundary conditions for the backward-facing step problem.
/// @tparam T            real number type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param dim           space dimension
/// @param[in]  periodic periodic domain (true/false)
/// @param[in]  inVel    the \a x component of the inflow velocity as string
template <class T>
void defineBCs_step(gsBoundaryConditions<T>& bcInfo, int dim, bool periodic = false, std::string inVel = "default")
{
    gsFunctionExpr<T> Uin, Uwall;

    switch (dim)
    {
    case 2:
    {
        if (inVel == "default")
            inVel = "(-4*(y-1.5)^2 + 1)";
        Uin = gsFunctionExpr<T>(inVel, "0", 2);
        Uwall = gsFunctionExpr<T>("0", "0", 2);
        break;
    }

    case 3:
    {
        if (!periodic)
        {
            if (inVel == "default")
                inVel = "(-4*(y-1.5)^2 + 1)*(-(z-1)^2 + 1)";

            Uin = gsFunctionExpr<T>(inVel, "0", "0", 3);
            Uwall = gsFunctionExpr<T>("0", "0", "0", 3);
        }
        else
        {
            if (inVel == "default")
                inVel = "(-4*(y-1.5)^2 + 1)";

            Uin = gsFunctionExpr<T>(inVel, "0", "0", 3);
            Uwall = gsFunctionExpr<T>("0", "0", "0", 3);
        }
        break;
    }
        
    default:
        GISMO_ERROR("Wrong dimension!");
        break;
    }

    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;
    defineBndParts_step(dim, bndIn, bndOut, bndWall, periodic);
    addBCs(bcInfo, bndIn, bndWall, Uin, Uwall);
}


/// @brief Define boundary conditions for the lid-driven cavity problem.
/// @tparam T            real number type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions
/// @param dim           space dimension
/// @param[in]  np       number of patches in each direction
/// @param[in]  lidVel   the \a x component of the lid velocity
/// @param[in]  fixPressure   set pressure to zero in one corner
template <class T>
void defineBCs_cavity(gsBoundaryConditions<T>& bcInfo, int dim, const int np = 1, std::string lidVelX = "1", std::string lidVelZ = "0", bool fixPressure = true)
{
    gsFunctionExpr<T> Uwall, Ulid;

    switch (dim)
    {
    case 2:
    {
        Uwall = gsFunctionExpr<T>("0", "0", 2);
        Ulid = gsFunctionExpr<T>(lidVelX, "0", 2);

        for (int i = 1; i <= np; i++)
            bcInfo.addCondition(np*np - i, boundary::north, condition_type::dirichlet, Ulid, 0);

        for (int i = 0; i < np; i++)
        {
            bcInfo.addCondition(i, boundary::south, condition_type::dirichlet, Uwall, 0);
            bcInfo.addCondition(i * np, boundary::west, condition_type::dirichlet, Uwall, 0);
            bcInfo.addCondition((i + 1)*np - 1, boundary::east, condition_type::dirichlet, Uwall, 0);
        }

        if (fixPressure)
            bcInfo.addCornerValue(boundary::southwest, 0.0, 0, 1); // corner, value, patch, unknown

        break;
    }

    case 3:
    {
        GISMO_ASSERT(np == 1, "Number of patches in each direction must be 1 for 3D cavity.");
        Uwall = gsFunctionExpr<T>("0", "0", "0", 3);
        Ulid = gsFunctionExpr<T>(lidVelX, "0", lidVelZ, 3);

        bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, Ulid, 0);
        bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, Uwall, 0);
        bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, Uwall, 0);

        if (fixPressure)
            bcInfo.addCornerValue(boundary::southwestfront, 0.0, 0, 1); // corner, value, patch, unknown

        break;
    }

    default:
        GISMO_ERROR("Wrong dimension!");
    }
}


/// @brief Define boundary conditions for the 2D blade profile problem.
/// @tparam T            real number type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[in]  inVel    the \a x component of the inflow velocity as string
/// @param[in] inVelX    x-coordinate of inflow velocity
/// @param[in] inVelY    y-coordinate of inflow velocity
template <class T>
void defineBCs_profile2D(gsBoundaryConditions<T>& bcInfo, T inVelX, T inVelY)
{
    gsFunctionExpr<T> Uin = gsFunctionExpr<T>(util::to_string(inVelX), util::to_string(inVelY), 2);
    gsFunctionExpr<T> Uwall = gsFunctionExpr<T>("0", "0", 2);

    std::vector<std::pair<int, boxSide> > bndIn, bndOut, bndWall;
    defineBndParts_profile2D(bndIn, bndOut, bndWall);
    addBCs(bcInfo, bndIn, bndWall, Uin, Uwall);
}


/// @brief Refine basis near wall (the first knot span).
/// @tparam T           real number type
/// @tparam d           space dimension
/// @param basis        reference to the basis to be refined
/// @param numRefine    number of recursive refinements
/// @param patch        patch number
/// @param dir          direction in which the basis will be refined
template< int d, class T>
void refineFirstKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
{
    const gsTensorBSplineBasis<d, T>*  patchBasis = dynamic_cast<const gsTensorBSplineBasis<d, T>*>(&(basis.basis(patch)));
    gsMatrix<T> box(d, 2);

    for (int i = 0; i < numRefine; i++)
    {
        box.setZero();
        T knot = patchBasis->knot(dir, patchBasis->degree(dir) + 1);
        box(dir, 1) = knot;
        basis.refine(patch, box);
    }
}


/// @brief Refine basis near wall (the last knot span).
/// @tparam T           real number type
/// @tparam d           space dimension
/// @param basis        reference to the basis to be refined
/// @param numRefine    number of recursive refinements
/// @param patch        patch number
/// @param dir          direction in which the basis will be refined
template< int d, class T>
void refineLastKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
{
    const gsTensorBSplineBasis<d, T>*  patchBasis = dynamic_cast<const gsTensorBSplineBasis<d, T>*>(&(basis.basis(patch)));
    gsMatrix<T> box(d, 2);

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


/// @brief Refine basis for the 2D lid-driven cavity problem.
/// @tparam T               real number type
/// @param basis            reference to the basis to be refined
/// @param numRefine        number of uniform refinements
/// @param numRefineLocal   number of near-wall refinements
/// @param dim              space dimension
/// @param numSep           number of \f C^0 \f separators in each patch (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template <class T>
void refineBasis_cavity(gsMultiBasis<T>& basis, int numRefine, int numRefineLocal, int dim, int numSep = 0)
{
    switch (dim)
    {
        case 2:
        {
            int npDir = math::sqrt(basis.nPieces());
            int npRef = std::log2(npDir);
            int sepRef = std::log2(numSep + 1);

            for (int i = 0; i < numRefine - npRef - sepRef; ++i)
                basis.uniformRefine();

            if (numRefineLocal && npDir == 1)
            {
                for (int dir = 0; dir < dim; dir++)
                {
                    refineFirstKnotSpan<2, T>(basis, numRefineLocal, 0, dir);
                    refineLastKnotSpan<2, T>(basis, numRefineLocal, 0, dir);
                }
            }

            break;
        }
        case 3:
        {
            for (int i = 0; i < numRefine; ++i)
                basis.uniformRefine();

            for (int dir = 0; dir < dim; dir++)
            {
                refineFirstKnotSpan<3, T>(basis, numRefineLocal, 0, dir);
                refineLastKnotSpan<3, T>(basis, numRefineLocal, 0, dir);
            }

            break;
        }
        default:
            GISMO_ERROR("Wrong dimension!");
            break;
    }
}


/// @brief Refine basis for the backward-facing step problem near walls.
/// @tparam T               real number type
/// @tparam d               space dimension
/// @param basis            reference to the basis to be refined
/// @param numRefineWalls   number of refinements near top and bottom walls
/// @param numRefineCorner  number of refinements near the step corner
template<int d, class T>
void refineLocal_step(gsMultiBasis<T>& basis, int numRefineWalls, int numRefineCorner)
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


/// @brief Refine basis for the backward-facing step problem.
/// @tparam T               real number type
/// @param basis            reference to the basis to be refined
/// @param numRefine        number of uniform refinements
/// @param numRefineWalls   number of refinements near top and bottom walls
/// @param numRefineCorner  number of refinements near the step corner
/// @param numRefineU       number of uniform refinements of patches 0 and 1 in \a u direction
/// @param addRefPart       a value of the \a u parameter for additional refinement behind the step 
/// @param dim              space dimension
template <class T>
void refineBasis_step(gsMultiBasis<T>& basis, int numRefine, int numRefineWalls, int numRefineCorner, int numRefineU, real_t addRefPart, int dim)
{
    gsMatrix<T> box(dim, 2);

    box.setZero();
    box(0,1) = addRefPart;
    for (int p = 0; p < 2; p++)
        basis.refine(p, box);

    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();
        
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


/// @brief Refine basis for the 2D profile problem.
/// @tparam T                   real number type
/// @param[in, out] basis       reference to the basis to be refined
/// @param[in] numRefine        number of uniform refinements
/// @param[in] numRefineBlade   number of refinements near the blade
/// @param[in] numRefineLead    number of refinements near the beginning of the blade
template <class T>
void refineBasis_profile2D(gsMultiBasis<T>& basis, int numRefine, int numRefineBlade, int numRefineLead)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    if (numRefineBlade)
    {
        for (int p = 0; p < basis.nPieces(); p++)
        {
            refineFirstKnotSpan<2, T>(basis, numRefineBlade, p, 1);
            refineLastKnotSpan<2, T>(basis, numRefineBlade, p, 1);
        }
    }

    if (numRefineLead)
    {
        refineLastKnotSpan<2, T>(basis, numRefineLead, 0, 0);
        refineFirstKnotSpan<2, T>(basis, numRefineLead, 1, 0);
    }
}


template <class T>
void plotQuantityFromSolution(std::string quantity, const gsField<T>& solField, std::string const & fn, unsigned npts)
{
    gsParaviewCollection pwCollection(fn);
    std::string fileName;

    for (unsigned int p = 0; p < solField.patches().nPatches(); p++)
    {
        fileName = fn + util::to_string(p);

        gsMatrix<T> geoVals, quantityVals;
        gsVector<unsigned> np;

        gsGeometry<T>* patch = &solField.patches().patch(p);
        const int parDim = patch->domainDim();
        const int tarDim = patch->targetDim();

        gsMatrix<T> ab = patch->support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);
        np = uniformSampleCount(a, b, npts);
        gsMatrix<T> pts = gsPointGrid(a, b, np);

        // evaluate geometry at pts
        geoVals = patch->eval(pts);

        if (3 - parDim > 0)
        {
            np.conservativeResize(3);
            np.bottomRows(3 - parDim).setOnes();
        }
        if (3 - tarDim > 0)
        {
            geoVals.conservativeResize(3, geoVals.cols());
            geoVals.bottomRows(3 - tarDim).setZero();
        }

        // evaluate quantity at pts
        if (quantity == "divergence")
            solField.function(p).div_into(pts, quantityVals); // incorrect?? (derivatives in parametric space?)
        else
            GISMO_ERROR("Function plotQuantityFromSolution for " + quantity + " not implemented.");

        // paraview export
        gsWriteParaviewTPgrid(geoVals, quantityVals, np.template cast<index_t>(), fileName);

        pwCollection.addPart(fileName + ".vts");
    }
    pwCollection.save();
}


template<class T>    
gsField<T> computeDistanceField(typename gsFlowSolverParams<T>::Ptr paramsPtr)
{
    gsInfo << "Computing required distance field ... ";

    gsMultiPatch<T> patches = paramsPtr->getPde().patches();    // multipatch representing the computational domain
    gsMultiBasis<T> basis = paramsPtr->getBasis(1);             // pressure basis as base basis for distance computation

    std::vector<std::pair<int, boxSide> > bndIn = paramsPtr->getBndIn();
    std::vector<std::pair<int, boxSide> > bndWall = paramsPtr->getBndWall();
    index_t numRefs = paramsPtr->options().getInt("TM.addRefsDF");

    for (int i = 0; i < numRefs; ++i)       // additional refinements for distance computation
        basis.uniformRefine();

    gsFunctionExpr<real_t> fw("1", "0", "0", patches.targetDim());
    gsFunctionExpr<real_t> gw("0", "0", "0", patches.targetDim());
    gsFunctionExpr<real_t> wallw("0.0", "0.0", "0.0", patches.targetDim());

    // boundary conditions for related Poisson problem
    gsBoundaryConditions<T> bcDF;
    addBCs(bcDF, bndIn, bndWall, gw, wallw, 0);

    // solving Poisson problem
    gsPoissonAssembler<T> assembler(patches, basis, bcDF, fw, dirichlet::elimination, iFace::glue);
    assembler.assemble();

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU solver;
    pardisoSetup<T>(solver);
    #else
    typename gsSparseSolver<T>::LU solver;
    #endif 

    gsMatrix<T> solVector = solver.compute( assembler.matrix() ).solve ( assembler.rhs() );
    for (int i = 0; i < solVector.rows(); i++)
    {
        if (solVector(i, 0) < 0)
            solVector(i, 0) = math::abs(solVector(i, 0));
    }

    gsField<T> solPoisson = assembler.constructSolution(solVector, 0);
    //gsWriteParaview<T>(solPoisson, "PoissonSolution", 10000);

    // evaluating distance values at a grid of points
    size_t np = patches.nPatches();
    gsMultiPatch<T>* wallDistanceMP = new gsMultiPatch<T>;
    for (size_t i = 0; i < np; i++)
    {
        index_t patchId = i;
        const gsBasis<T> & basisp = basis.piece(patchId);

        std::vector< gsVector<T> > rr;
        rr.reserve(patches.parDim());

        for (short_t j = 0; j < patches.parDim(); ++j)            // computing grid of point
        {
            rr.push_back(basisp.component(j).anchors().transpose());
        }
        gsMatrix<T> gridPts = gsPointGrid<T>(rr);

        gsMapData<T> mapData;
        unsigned geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        mapData.flags = geoFlags;
        mapData.patchId = patchId;
        mapData.points = gridPts;
        paramsPtr->getPde().patches().patch(patchId).computeMap(mapData);

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        basisp.deriv_into(gridPts, parGrads);
        
        gsMatrix<T> solPoissonCoeffVec = solPoisson.coefficientVector(patchId);
        std::vector<gsMatrix<T> > solPoissonGrads(gridPts.cols());
        for (index_t k = 0; k < gridPts.cols(); k++)
        {
            // eval solPoissonGrads at all pts
            basisp.active_into(gridPts.col(k), actives);
            int numActP = actives.rows();
            gsMatrix<T> solActPoissonCoeffs(1, numActP);
            for (int j = 0; j < numActP; j++)
                solActPoissonCoeffs(0, j) = solPoissonCoeffVec(actives(j, 0), 0);

            transformGradients(mapData, k, parGrads, physGrad);
            solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGrad.transpose();
        }

        gsMatrix<T> solPoissonVals = solPoisson.value(gridPts, patchId);
        
        // Evaluate wall distance at pts
        gsMatrix<T> wallDistanceVals(1, gridPts.cols());
        for (index_t k = 0; k < gridPts.cols(); k++)
        {
            wallDistanceVals(0, k) = math::sqrt(math::pow(solPoissonGrads[k].norm(), 2) + 2 * solPoissonVals(0, k)) - solPoissonGrads[k].norm();
            if (math::isnan(wallDistanceVals(0, k)))
            {
                wallDistanceVals(0, k) = 0.;
            }
        }

        typename gsGeometry<T>::uPtr geo = basisp.interpolateAtAnchors(wallDistanceVals);    // interpolating distances at grid points 
        const gsMatrix<T> & distanceCoeffs = geo->coefs();
        wallDistanceMP->addPatch(basisp.makeGeometry(distanceCoeffs));
    }
    
    gsInfo << "Done" << std::endl;
    gsField<T> result = gsField<T>(paramsPtr->getPde().patches(), typename gsFunctionSet<T>::Ptr(wallDistanceMP), true);
    return result;
}

template<class T>
T computeDimensionlessWallDistance(typename gsFlowSolverParams<T>::Ptr paramsPtr, /*gsMatrix<T> solution, std::vector<gsMultiBasis<T>> multiBasis, gsVector<int> distancePatches, std::vector<boxSide> distanceSides,*/ T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus, unsigned npts, bool print, bool estimate)
{
    std::ofstream ofile;
    if (print)
    {
        time_t t = std::time(0);   // get time now
        struct tm * now = std::localtime(&t);
        std::stringstream filename;
        filename << "wallDistance_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_" << now->tm_hour << "-" << now->tm_min << ".txt";

        ofile.open(filename.str().c_str());

        ofile << "\n             Dimensionless wall distance information                       \n";
        ofile << "=========================================================================================================\n";
        ofile << "patch index  |     side      |   number of unfit points (y+ > " << maxYplus << ") |   min y+   |   max y+   |   average y+   \n";
        ofile << "---------------------------------------------------------------------------------------------------------\n";
    }

    gsMultiPatch<T> patches = paramsPtr->getPde().patches();    // multipatch representing the computational domain
    std::vector<gsMultiBasis<T> > basis = paramsPtr->getBases();           // pressure basis as base basis for distance computation
    std::vector<std::pair<int, boxSide> > bndWall = paramsPtr->getBndWall();
    index_t tarDim = patches.targetDim();

    gsVector<T> minYPlusOverSides;
    minYPlusOverSides.setZero(bndWall.size());

    T skinFrCoeff = math::pow(2.0*math::log10(reynoldsNumber) - 0.65, -2.3);
    T tau_w = 0.5 * skinFrCoeff * math::pow(uFreeStream, 2);
    T u_tau = math::sqrt(tau_w);

    int unfitPointsNum;
    for (size_t numSide = 0; numSide < bndWall.size(); numSide++)
    {
        int patchIndex = bndWall[numSide].first;
        boxSide side = bndWall[numSide].second;

        typename gsGeometry<T>::uPtr geomBoundary = patches.patch(patchIndex).boundary(side);

        gsMatrix<T> ab = geomBoundary->support();
        gsVector<T> a = ab.col(0);
        gsVector<T> b = ab.col(1);

        gsVector<unsigned> np = uniformSampleCount(a, b, npts);
        gsMatrix<T> paramBoundaryNodes = gsPointGrid(a, b, np);
        int numOfPoints = paramBoundaryNodes.cols();

        gsMatrix<T> paramBoundaryNodesFull;
        paramBoundaryNodesFull.setZero(patches.dim(), paramBoundaryNodes.cols());
        gsMatrix<T> paramFirstNodesFull = paramBoundaryNodesFull;

        int param = side.parameter();
        int index = 0;

        const gsBasis<T>& basisPatch = (basis.at(1)).piece(patchIndex);
        bool THB = false;
        //----------------------------------------------
        if (tarDim == 2)
        {
            if (dynamic_cast<const gsTHBSplineBasis<2, T>*>(&basisPatch))
                THB = true;
        }
        else
        {
            if (dynamic_cast<const gsTHBSplineBasis<3, T>*>(&basisPatch))
                THB = true;
        }
        //----------------------------------------------

        for (int s = 0; s < patches.dim(); s++)
        {
            if (s == side.direction())
            {
                if (param == 1) // east or north or back
                {
                    paramBoundaryNodesFull.middleRows(s, 1).setOnes();

                    real_t nodeValue;
                    if (THB)
                    {
                        if (tarDim == 2)
                        {
                            const gsTHBSplineBasis<2, T>* pBasisTHB = dynamic_cast<const gsTHBSplineBasis<2, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasisTHB->knot(pBasisTHB->maxLevel(),s, pBasisTHB->numKnots(pBasisTHB->maxLevel(),s) - pBasisTHB->degree(s) - 2);
                        }
                        else
                        {
                            const gsTHBSplineBasis<3, T>* pBasisTHB = dynamic_cast<const gsTHBSplineBasis<3, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasisTHB->knot(pBasisTHB->maxLevel(),s, pBasisTHB->numKnots(pBasisTHB->maxLevel(),s) - pBasisTHB->degree(s) - 2);
                        }
                    }
                    else
                    {
                        if (tarDim == 2)
                        {
                            const gsTensorBSplineBasis<2, T>* pBasis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasis->knot(s, pBasis->knots(s).size() - pBasis->degree(s) - 2);
                        }
                        else
                        {
                            const gsTensorBSplineBasis<3, T>* pBasis = dynamic_cast<const gsTensorBSplineBasis<3, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasis->knot(s, pBasis->knots(s).size() - pBasis->degree(s) - 2);
                        }
                    }
                    //basis->knot(0, basis->knots(0).size() - basis->degree(0) - 2)//east
                    //basis->knot(1, basis->knots(1).size() - basis->degree(1) - 2)//north
                    //basis->knot(2, basis->knots(2).size() - basis->degree(2) - 2)//back

                    paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                }
                else // param = 0 ... i.e. west or south or front
                {
                    real_t nodeValue;
                    if (THB)
                    {
                        if (tarDim == 2)
                        {
                            const gsTHBSplineBasis<2, T>* pBasisTHB = dynamic_cast<const gsTHBSplineBasis<2, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasisTHB->knot(pBasisTHB->maxLevel(),s, pBasisTHB->degree(s) + 1);
                        }
                        else
                        {
                            const gsTHBSplineBasis<3, T>* pBasisTHB = dynamic_cast<const gsTHBSplineBasis<3, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasisTHB->knot(pBasisTHB->maxLevel(),s, pBasisTHB->degree(s) + 1);
                        }
                    }
                    else
                    {
                        if (tarDim == 2)
                        {
                            const gsTensorBSplineBasis<2, T>* pBasis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasis->knot(s, pBasis->degree(s) + 1);
                        }
                        else
                        {
                            const gsTensorBSplineBasis<3, T>* pBasis = dynamic_cast<const gsTensorBSplineBasis<3, T>*>(&basis.at(1).piece(patchIndex));
                            nodeValue = pBasis->knot(s, pBasis->degree(s) + 1);
                        }
                    }
                    //basis->knot(0, basis->degree(0) + 1)//west
                    //basis->knot(1, basis->degree(1) + 1)//south
                    //basis->knot(2, basis->degree(2) + 1)//front

                    paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                }

            }
            else
            {
                paramBoundaryNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                paramFirstNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                index++;
            }
        }

        gsMatrix<T> physBoundaryNodesFull;
        gsMatrix<T> physFirstNodesFull;
        patches.patch(patchIndex).eval_into(paramBoundaryNodesFull, physBoundaryNodesFull); // physical coordinates of paramBoundaryNodes
        patches.patch(patchIndex).eval_into(paramFirstNodesFull, physFirstNodesFull); // physical coordinates of paramBoundaryNodes

        T averageValue = 0.;
        unfitPointsNum = 0;
        gsVector<T> yPlusAtSide;
        yPlusAtSide.setZero(numOfPoints);
        for (int nodeIndex = 0; nodeIndex < numOfPoints; nodeIndex++)
        {
            T distanceToWall = (physFirstNodesFull.col(nodeIndex) - physBoundaryNodesFull.col(nodeIndex)).norm();
            yPlusAtSide[nodeIndex] = distanceToWall * u_tau / viscosity;

            averageValue += yPlusAtSide[nodeIndex];
            if (yPlusAtSide[nodeIndex] > maxYplus)
                unfitPointsNum++;
        }
        if (print)
            ofile << "      " << patchIndex << "      |   " << side << "    |                 " << unfitPointsNum << "                 |   " << yPlusAtSide.minCoeff() <<
            "  |   " << yPlusAtSide.maxCoeff() << "  |   " << averageValue / numOfPoints << "\n";

        minYPlusOverSides[numSide] = yPlusAtSide.minCoeff();
    }
    //ofile << "=========================================================================================================\n\n";
    if (print)
        ofile.close();

    return minYPlusOverSides.minCoeff() * viscosity / u_tau;
} 


template<class T>
gsMultiPatch<T> linearizeGeometry(gsMultiPatch<T> mp, int uRefine = 0)
{
    gsTensorBSpline<2, T>* check1 = dynamic_cast<gsTensorBSpline<2, T>*>(&(mp.patch(0)));
    gsTensorBSpline<3, T>* check2 = dynamic_cast<gsTensorBSpline<3, T>*>(&(mp.patch(0)));

    if (check1 == NULL && check2 == NULL)
        GISMO_ERROR("linearizeGeometry() can only be applied to tensor B-spline multi-patches.");

    gsMultiPatch<T> mp_new;
    int dim = mp.dim();
    gsMatrix<int> degs(mp.nPatches(), dim);

    for (index_t i = 0; i < mp.nPatches(); i++)
        for (int j = 0; j < dim; j++) 
            degs(i, j) = mp[i].degree(j);

    for (index_t i = 0; i < mp.nPatches(); i++)
    {
        if (dim == 2)
        {
            gsTensorBSpline<2, T>* tbspline = dynamic_cast<gsTensorBSpline<2, T>*>(&(mp.patch(i)));

            for (int j = 0; j < uRefine; j++)
                tbspline->uniformRefine();

            std::vector<T> knots1 = (tbspline->knots(0)).breaks();
            std::vector<T> knots2 = (tbspline->knots(1)).breaks();
            gsMatrix<T> parpoints(2, (knots1.size()) * (knots2.size()));
            gsMatrix<T> points(2, (knots1.size()) * (knots2.size()));

            for (index_t b = 0; b < knots2.size(); b++)
                for (index_t a = 0; a < knots1.size(); a++)
                {
                    parpoints(0, b * (knots1.size()) + a) = knots1[a];
                    parpoints(1, b * (knots1.size()) + a) = knots2[b];
                }

            tbspline->eval_into(parpoints, points);

            gsKnotVector<T> kv1((tbspline->knots(0)).unique(), 1, 0);
            gsKnotVector<T> kv2((tbspline->knots(1)).unique(), 1, 0);

            gsTensorBSpline<2, T> tbspline_new(kv1, kv2, points.transpose());
            mp_new.addPatch(tbspline_new);
        }
        else 
        {
            gsTensorBSpline<3, T>* tbspline = dynamic_cast<gsTensorBSpline<3, T>*>(&(mp.patch(i)));

            for (int j = 0; j < uRefine; j++)
                tbspline->uniformRefine();

            std::vector<T> knots1 = (tbspline->knots(0)).breaks();
            std::vector<T> knots2 = (tbspline->knots(1)).breaks();
            std::vector<T> knots3 = (tbspline->knots(2)).breaks();
            gsMatrix<T> parpoints(3, (knots1.size()) * (knots2.size()) * (knots3.size()));
            gsMatrix<T> points(3, (knots1.size()) * (knots2.size()) * (knots3.size()));

            for (index_t c = 0; c < knots3.size(); c++)
                for (index_t b = 0; b < knots2.size(); b++)
                    for (index_t a = 0; a < knots1.size(); a++)
                    {
                        parpoints(0, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots1[a];
                        parpoints(1, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots2[b];
                        parpoints(2, c * (knots1.size()) * (knots2.size()) + b * (knots1.size()) + a) = knots3[c];
                    }

            tbspline->eval_into(parpoints, points);

            gsKnotVector<T> kv1((tbspline->knots(0)).unique(), 1, 0);
            gsKnotVector<T> kv2((tbspline->knots(1)).unique(), 1, 0);
            gsKnotVector<T> kv3((tbspline->knots(2)).unique(), 1, 0);

            gsTensorBSpline<3, T> tbspline_new(kv1, kv2, kv3, points.transpose());
            mp_new.addPatch(tbspline_new);
        }
    }

    for (index_t i = 0; i < mp.nInterfaces(); i++)
        mp_new.addInterface(mp.bInterface(i));

    mp_new.addAutoBoundaries();

    return mp_new;
}


// -----------------------------------------------


template <class T>
struct gsVectorHash
{
    std::size_t operator()(const gsVector<T>& vec) const
    {
        std::size_t seed = vec.size();
        for (int i = 0; i < vec.size(); i++)
            seed ^= std::hash<T>()(vec(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};


inline index_t mapTensorID(index_t i, index_t j, index_t k, size_t size1, size_t size2)
{
    index_t result = (k*size1*size2)+(j*size1)+i;
    return result;
}

} //namespace gismo