/** @file gsINSUtils.h

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

/// @brief Fill a diagonal approximation of an inverse matrix.
/// @param[in]  mat     a const reference to the matrix of which the inverse is approximated
/// @param[out] diagInv a reference to the resulting inverse approximation
/// @param[in]  repeat  number of the diagonal block repetition (e.g. for velocity components)
/// @param[in]  lumping use lumping to define the diagonal approximation
template<class T>
void diagInvMatrix_into(const gsSparseMatrix<T>& mat, gsSparseMatrix<T>& diagInv, int repeat, bool lumping = false)
{
    GISMO_ENSURE(mat.nonZeros() != 0, "diagInvMatrix_into(): The matrix is empty!");

    int varDofs = mat.rows();

    diagInv.resize(repeat*varDofs, repeat*varDofs);
    diagInv.reserve(gsVector<int>::Constant(diagInv.cols(), 1));

    const gsSparseMatrix<T>* matPtr = &mat;
    gsSparseMatrix<T> lumped(varDofs, varDofs);

    if (lumping)
    {
        lumped.reserve(gsVector<int>::Constant(varDofs, 1));

        for (int j = 0; j < varDofs; j++)
            lumped.insert(j, j) = math::abs(mat.at(j, j)); // abs value because of "lumping" diag(A) in SIMPLE-type prec., does not change lumped mass matrix in IgA

        for (int j = 0; j < varDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(mat, j); it; ++it)
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
/// @tparam T       coefficient type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param llx      \a x coordinate of the lower left corner
/// @param lly      \a y coordinate of the lower left corner
/// @param a        width of the rectangle
/// @param b        height of the rectangle
/// @param numSep   number of \f C^0 \f separators in the parametrization (knots of multiplicity \a deg uniformly distributed in each parametric direction)
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


/// @brief Returns a B-spline parametrization of a 3D block of a given degree in all directions.
/// @tparam T       coefficient type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param llx      \a x coordinate of a corner
/// @param lly      \a y coordinate of a corner
/// @param llz      \a z coordinate of a corner
/// @param a        width of the block
/// @param b        height of the block
/// @param c        depth of the block
/// @param numSep   number of \f C^0 \f separators in the parametrization (knots of multiplicity \a deg uniformly distributed in each parametric direction)
template<class T>
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
/// @tparam T       coefficient type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        width of the cavity
/// @param b        height of the cavity
/// @param np       number of patches in each direction
/// @param numSep   number of \f C^0 \f separators in each patch (knots of multiplicity \a deg uniformly distributed in each parametric direction)
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


/// @brief Returns a B-spline multipatch domain for 2D problems of flow over a backward facing step.
/// @tparam T       coefficient type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        length of the domain behind the step
/// @param b        total height of the domain
/// @param a_in     length of the inflow part (before the step)
/// @param h        height of the step
/// @param periodic periodic domain (if true, the bottom and top boundaries of the channel behind step are defined as interface)
template<class T> gsMultiPatch<T> BSplineStep2D(int deg, const T a, const T b, const T a_in, T h = 0, bool periodic = false)
{
    gsMultiPatch<T> mp;

    if (h == 0)
        h = b / 2;

    mp.addPatch(BSplineRectangle(deg, 0.0, 0.0, a, h));
    mp.addPatch(BSplineRectangle(deg, 0.0, h, a, b - h));
    mp.addPatch(BSplineRectangle(deg, -a_in, h, a_in, b - h));

    mp.addInterface(0, boundary::north, 1, boundary::south);
    mp.addInterface(1, boundary::west, 2, boundary::east);

    if(periodic)
        mp.addInterface(0, boundary::south, 1, boundary::north);

    mp.addAutoBoundaries();

    return mp;
}


/// @brief Returns a B-spline multipatch domain for 3D problems of flow over a backward facing step.
/// @tparam T       coefficient type
/// @param deg      polynomial degree of the B-spline parametrization
/// @param a        length of the domain behind the step
/// @param b        total height of the domain
/// @param c        width of the domain
/// @param a_in     length of the inflow part (before the step)
/// @param h        height of the step
/// @param periodic periodic domain (if true, the bottom and top boundaries of the channel behind step are defined as interface)
template<class T> gsMultiPatch<T> BSplineStep3D(int deg, const T a, const T b, const T c, const T a_in, T h = 0, bool periodic = false)
{
    gsMultiPatch<T> mp;

    if (h == 0)
        h = b / 2;

    mp.addPatch(BSplineBlock<T>(deg, 0.0, 0.0, 0.0, a, h, c));
    mp.addPatch(BSplineBlock<T>(deg, 0.0, h, 0.0, a, b - h, c));
    mp.addPatch(BSplineBlock<T>(deg, -a_in, h, 0.0, a_in, b - h, c));

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


/// @brief Define boundary conditions for the corresponding boundary parts.
/// @tparam T            coefficient type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[in]  bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[in]  bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  Uin      the inflow velocity as gsFunctionExpr
/// @param[in]  Uwall    the wall velocity as gsFunctionExpr
template<class T> void addBCs(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndWall, gsFunctionExpr<T> Uin, gsFunctionExpr<T> Uwall)
{
    for (size_t i = 0; i < bndIn.size(); i++)
        bcInfo.addCondition(bndIn[i].first, bndIn[i].second, condition_type::dirichlet, Uin);

    for (size_t i = 0; i < bndWall.size(); i++)
        bcInfo.addCondition(bndWall[i].first, bndWall[i].second, condition_type::dirichlet, Uwall);
}


/// @brief Define boundary conditions for the 2D lid-driven cavity problem.
/// @tparam T           coefficient type
/// @param[out] bcInfo  reference to the boundary conditions as gsBoundaryConditions 
/// @param[in]  np      number of patches in each direction
/// @param[out] bndWall reference to a container of patch sides corresponding to solid walls
/// @param[in]  lidVel  the lid velocity
template<class T>
void defineBCs_cavity2D(gsBoundaryConditions<T>& bcInfo, const int np, std::vector<std::pair<int, boxSide> >& bndWall, std::string lidVel = "1")
{
    gsFunctionExpr<T> Uwall("0", "0", 2);
    gsFunctionExpr<T> Ulid(lidVel, "0", 2);

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


/// @brief Define boundary conditions for the 2D backward-facing step problem.
/// @tparam T            coefficient type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[out] bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[out] bndOut   reference to a container of patch sides corresponding to outflow boundaries
/// @param[out] bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  periodic periodic domain (true/false)
/// @param[in]  inVel    the \a x component of the inflow velocity as string
template<class T> void defineBCs_step2D(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false, std::string inVel = "default")
{
    if (inVel == "default")
        inVel = "(-4*(y-1.5)^2 + 1)";

    gsFunctionExpr<T> Uin, Uwall;
    Uin = gsFunctionExpr<T>(inVel, "0", 2);
    Uwall = gsFunctionExpr<T>("0", "0", 2);

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


/// @brief Define boundary conditions for the 3D backward-facing step problem.
/// @tparam T            coefficient type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[out] bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[out] bndOut   reference to a container of patch sides corresponding to outflow boundaries
/// @param[out] bndWall  reference to a container of patch sides corresponding to solid walls
/// @param[in]  periodic periodic domain (true/false)
/// @param[in]  inVel    the \a x component of the inflow velocity as string
template<class T> void defineBCs_step3D(gsBoundaryConditions<T>& bcInfo, std::vector<std::pair<int, boxSide> >& bndIn, std::vector<std::pair<int, boxSide> >& bndOut, std::vector<std::pair<int, boxSide> >& bndWall, bool periodic = false, std::string inVel = "default")
{
    gsFunctionExpr<T> Uin, Uwall;

    if (!periodic)
    {
        if (inVel == "default")
            inVel = "(-4*(y-1.5)^2 + 1)*(-(z-1)^2 + 1)";

        Uin = gsFunctionExpr<T>(inVel, "0", "0", 3);
        Uwall = gsFunctionExpr<T>("0", "0", "0", 3);

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

        Uin = gsFunctionExpr<T>(inVel, "0", "0", 3);
        Uwall = gsFunctionExpr<T>("0", "0", "0", 3);

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


/// @brief Define boundary conditions for the backward-facing step problem.
/// @tparam T            coefficient type
/// @param[out] bcInfo   reference to the boundary conditions as gsBoundaryConditions 
/// @param[out] bndIn    reference to a container of patch sides corresponding to inflow boundaries
/// @param[out] bndOut   reference to a container of patch sides corresponding to outflow boundaries
/// @param[out] bndWall  reference to a container of patch sides corresponding to solid walls
/// @param dim           space dimension
/// @param[in]  periodic periodic domain (true/false)
/// @param[in]  inVel    the \a x component of the inflow velocity as string
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


/// @brief Refine basis near wall (the first knot span).
/// @tparam T           coefficient type
/// @tparam d           space dimension
/// @param basis        reference to the basis to be refined
/// @param numRefine    number of recursive refinements
/// @param patch        patch number
/// @param dir          direction in which the basis will be refined
template< int d, class T> void refineFirstKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
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
/// @tparam T           coefficient type
/// @tparam d           space dimension
/// @param basis        reference to the basis to be refined
/// @param numRefine    number of recursive refinements
/// @param patch        patch number
/// @param dir          direction in which the basis will be refined
template< int d, class T> void refineLastKnotSpan(gsMultiBasis<T>& basis, int numRefine, int patch, int dir)
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
/// @tparam T               coefficient type
/// @param basis            reference to the basis to be refined
/// @param numRefine        number of uniform refinements
/// @param numRefineLocal   number of near-wall refinements
/// @param numSep           number of \f C^0 \f separators in each patch (knots of multiplicity \a deg uniformly distributed in each parametric direction)
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


/// @brief Refine basis for the backward-facing step problem near walls.
/// @tparam T               coefficient type
/// @tparam d               space dimension
/// @param basis            reference to the basis to be refined
/// @param numRefineWalls   number of refinements near top and bottom walls
/// @param numRefineCorner  number of refinements near the step corner
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


/// @brief Refine basis for the backward-facing step problem.
/// @tparam T               coefficient type
/// @param basis            reference to the basis to be refined
/// @param numRefine        number of uniform refinements
/// @param numRefineWalls   number of refinements near top and bottom walls
/// @param numRefineCorner  number of refinements near the step corner
/// @param numRefineU       number of uniform refinements of patches 0 and 1 in \a u direction
/// @param addRefPart       a value of the \a u parameter for additional refinement behind the step 
/// @param dim              space dimension
/// @param a                length of the domain behind the step
/// @param b                total height of the domain
/// @param c                width of the domain (3D case)
template<class T> void refineBasis_step(gsMultiBasis<T>& basis, int numRefine, int numRefineWalls, int numRefineCorner, int numRefineU, real_t addRefPart, int dim, real_t a, real_t b, real_t c = 0.0)
{
    for (int i = 0; i < numRefine; ++i)
        basis.uniformRefine();

    gsMatrix<T> box(dim, 2);

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

} //namespace gismo