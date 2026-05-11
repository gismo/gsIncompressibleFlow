/** @file gsINSMeshDeformer.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li (2025)
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsIO/gsOptionList.h>
#include <gsMatrix/gsSparseSolver.h>

#include <algorithm>
#include <string>
#include <vector>

namespace gismo
{

namespace mesh_deform_method
{
    enum method
    {
        HE = 0,
        IHE = 1,
        LE = 2,
        ILE = 3,
        TINE = 4,
        TINE_StVK = 5,
        BHE = 6,
        IBHE = 7
    };
}

/// Mesh deformer used by incompressible-flow ALE examples.
///
/// The class owns the mesh-displacement state and updates the referenced
/// geometry. Boundary control points are constrained strongly. Interior control
/// points are obtained from Galerkin systems on the geometry basis. Harmonic
/// methods solve a Laplace extension, linear elasticity solves the isotropic
/// small-strain elasticity system, and biharmonic methods solve a scalar
/// biharmonic extension component-wise. Incremental methods solve for the
/// boundary residual and accumulate the result.
template <class T>
class gsINSMeshDeformer
{
public:
    gsINSMeshDeformer(gsMultiPatch<T>& fluidGeometry,
                      mesh_deform_method::method method = mesh_deform_method::IHE)
    : m_geometry(fluidGeometry),
      m_originalGeometry(fluidGeometry),
      m_method(method),
      m_options(defaultOptions()),
      m_solverVelBasis(nullptr),
      m_isInitialized(false),
      m_hasSavedState(false),
      m_lastCheckResult(-1)
    {}

    void setMovingBoundaries(const std::vector<patchSide>& movingBoundaries)
    {
        m_movingBoundaries = movingBoundaries;
    }

    void initialize()
    {
        if (m_isInitialized)
            return;

        makeZeroLike(m_geometry, m_boundaryDisplacement);
        makeZeroLike(m_geometry, m_meshDisplacement);
        makeZeroLike(m_geometry, m_savedMeshDisplacement);

        m_isInitialized = true;
        gsInfo << "gsINSMeshDeformer initialized with method: " << methodName() << "\n";
    }

    void updateMesh(const gsMultiPatch<T>& boundaryDisplacementValues)
    {
        if (!m_isInitialized)
            initialize();

        GISMO_ENSURE(isSupportedMethod(m_method),
                     "Unsupported mesh deformation method: " << static_cast<int>(m_method));

        copyBoundaryDisplacement(boundaryDisplacementValues);

        if (isIncrementalMethod())
        {
            gsMultiPatch<T> update;
            solveAllPatches(true, update);
            addToDisplacement(update);
        }
        else
        {
            solveAllPatches(false, m_meshDisplacement);
        }

        applyDisplacementToGeometry();
        m_lastCheckResult = m_options.getSwitch("Check") ? checkCurrentGeometry() : -1;
        if (m_lastCheckResult >= 0)
            gsWarn << "[gsINSMeshDeformer] Mesh check found a non-positive Jacobian on patch "
                   << m_lastCheckResult << ".\n";
    }

    gsMultiPatch<T> getMeshDisplacement() const
    {
        if (!m_isInitialized)
        {
            gsMultiPatch<T> zero;
            makeZeroLike(m_geometry, zero);
            return zero;
        }
        return m_meshDisplacement;
    }

    void resetToOriginal()
    {
        m_geometry = m_originalGeometry;
        if (m_isInitialized)
        {
            makeZeroLike(m_geometry, m_boundaryDisplacement);
            makeZeroLike(m_geometry, m_meshDisplacement);
        }
        m_lastCheckResult = -1;
    }

    gsOptionList& options() { return m_options; }
    const gsOptionList& options() const { return m_options; }

    void setVelocityBasis(const gsMultiBasis<T>* velBasis)
    {
        m_solverVelBasis = velBasis;
    }

    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("PoissonsRatio", "Poisson's ratio for linear elasticity", 0.3);
        opt.addSwitch("Check", "Check mesh geometry after deformation", true);
        return opt;
    }

    static bool isSupportedMethod(mesh_deform_method::method method)
    {
        return method == mesh_deform_method::HE  ||
               method == mesh_deform_method::IHE ||
               method == mesh_deform_method::LE  ||
               method == mesh_deform_method::BHE;
    }

    std::string methodName() const
    {
        switch (m_method)
        {
            case mesh_deform_method::HE:        return "Harmonic Extension";
            case mesh_deform_method::IHE:       return "Incremental Harmonic Extension";
            case mesh_deform_method::LE:        return "Linear Elasticity Extension";
            case mesh_deform_method::BHE:       return "Biharmonic Extension";
            case mesh_deform_method::ILE:
            case mesh_deform_method::TINE:
            case mesh_deform_method::TINE_StVK:
            case mesh_deform_method::IBHE:
            default:
                GISMO_NO_IMPLEMENTATION
        }
        return "Unknown";
    }

    bool isInitialized() const { return m_isInitialized; }
    index_t lastCheckResult() const { return m_lastCheckResult; }

    void saveState()
    {
        m_savedGeometry = m_geometry;
        m_savedMeshDisplacement = m_meshDisplacement;
        m_savedBoundaryDisplacement = m_boundaryDisplacement;
        m_savedCheckResult = m_lastCheckResult;
        m_hasSavedState = true;
    }

    void recoverState()
    {
        GISMO_ENSURE(m_hasSavedState, "No mesh-deformer state has been saved.");
        m_geometry = m_savedGeometry;
        m_meshDisplacement = m_savedMeshDisplacement;
        m_boundaryDisplacement = m_savedBoundaryDisplacement;
        m_lastCheckResult = m_savedCheckResult;
    }

protected:
    static void makeZeroLike(const gsMultiPatch<T>& source, gsMultiPatch<T>& target)
    {
        target.clear();
        for (size_t p = 0; p < source.nPatches(); ++p)
        {
            target.addPatch(source.patch(p).clone());
            target.patch(p).coefs().setZero();
        }
        if (target.nPatches() > 0)
            target.computeTopology();
    }

    void copyBoundaryDisplacement(const gsMultiPatch<T>& boundaryDisplacementValues)
    {
        makeZeroLike(m_geometry, m_boundaryDisplacement);

        const size_t nPatches = std::min(m_boundaryDisplacement.nPatches(),
                                         boundaryDisplacementValues.nPatches());
        for (size_t p = 0; p < nPatches; ++p)
        {
            gsMatrix<T>& target = m_boundaryDisplacement.patch(p).coefs();
            const gsMatrix<T>& source = boundaryDisplacementValues.patch(p).coefs();

            const index_t rows = std::min(target.rows(), source.rows());
            const index_t cols = std::min(target.cols(), source.cols());
            for (index_t i = 0; i < rows; ++i)
                for (index_t d = 0; d < cols; ++d)
                    target(i, d) = source(i, d);
        }
    }

    void solveAllPatches(bool incremental, gsMultiPatch<T>& result) const
    {
        GISMO_ENSURE(isSupportedMethod(m_method),
                     "Unsupported mesh deformation method: " << static_cast<int>(m_method));

        makeZeroLike(m_geometry, result);

        for (size_t p = 0; p < m_geometry.nPatches(); ++p)
        {
            gsMatrix<T> disp;
            if (m_method == mesh_deform_method::LE)
                solvePatchLinearElasticity(static_cast<index_t>(p), incremental, disp);
            else
                solvePatchScalarExtension(static_cast<index_t>(p), incremental, isBiharmonicMethod(), disp);
            result.patch(p).coefs() = disp;
        }
    }

    void solvePatchScalarExtension(index_t patch, bool incremental, bool biharmonic, gsMatrix<T>& disp) const
    {
        const gsMatrix<T>& orig = m_originalGeometry.patch(patch).coefs();
        const index_t nCoefs = orig.rows();
        const index_t dim = orig.cols();

        disp.setZero(nCoefs, dim);
        if (nCoefs == 0)
            return;

        std::vector<char> fixed;
        gsMatrix<T> values;
        collectConstraints(patch, incremental, fixed, values);

        index_t nFree = 0;
        std::vector<index_t> freeMap(nCoefs, -1);
        for (index_t i = 0; i < nCoefs; ++i)
        {
            if (!fixed[i])
                freeMap[i] = nFree++;
        }

        if (nFree == 0)
        {
            disp = values;
            return;
        }

        gsSparseMatrix<T> op;
        assembleScalarExtensionOperator(patch, biharmonic, op);

        std::vector<gsEigen::Triplet<T,index_t> > triplets;
        triplets.reserve(static_cast<size_t>(op.nonZeros()));
        gsMatrix<T> rhs(nFree, dim);
        rhs.setZero();

        for (index_t k = 0; k < op.outerSize(); ++k)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(op, k); it; ++it)
            {
                const index_t row = static_cast<index_t>(it.row());
                const index_t col = static_cast<index_t>(it.col());
                const T value = static_cast<T>(it.value());

                if (fixed[row])
                    continue;

                const index_t r = freeMap[row];
                if (fixed[col])
                {
                    for (index_t d = 0; d < dim; ++d)
                        rhs(r, d) -= value * values(col, d);
                }
                else
                {
                    triplets.push_back(gsEigen::Triplet<T,index_t>(r, freeMap[col], value));
                }
            }
        }

        gsSparseMatrix<T> matrix(nFree, nFree);
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        matrix.makeCompressed();

        gsSparseSolver<>::SimplicialLDLT solver(matrix);
        gsMatrix<T> freeSol = solver.solve(rhs);

        disp = values;
        for (index_t i = 0; i < nCoefs; ++i)
            if (!fixed[i])
                disp.row(i) = freeSol.row(freeMap[i]);
    }

    void solvePatchLinearElasticity(index_t patch, bool incremental, gsMatrix<T>& disp) const
    {
        const gsMatrix<T>& orig = m_originalGeometry.patch(patch).coefs();
        const index_t nCoefs = orig.rows();
        const index_t dim = orig.cols();

        disp.setZero(nCoefs, dim);
        if (nCoefs == 0)
            return;

        std::vector<char> fixed;
        gsMatrix<T> values;
        collectConstraints(patch, incremental, fixed, values);

        index_t nFree = 0;
        std::vector<index_t> freeMap(static_cast<size_t>(nCoefs * dim), -1);
        for (index_t i = 0; i < nCoefs; ++i)
            for (index_t d = 0; d < dim; ++d)
                if (!fixed[i])
                    freeMap[static_cast<size_t>(d * nCoefs + i)] = nFree++;

        if (nFree == 0)
        {
            disp = values;
            return;
        }

        gsSparseMatrix<T> op;
        assembleLinearElasticityOperator(patch, op);

        std::vector<gsEigen::Triplet<T,index_t> > triplets;
        triplets.reserve(static_cast<size_t>(op.nonZeros()));
        gsMatrix<T> rhs(nFree, 1);
        rhs.setZero();

        for (index_t k = 0; k < op.outerSize(); ++k)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(op, k); it; ++it)
            {
                const index_t row = static_cast<index_t>(it.row());
                const index_t col = static_cast<index_t>(it.col());
                const T value = static_cast<T>(it.value());
                const index_t rowCoef = row % nCoefs;
                const index_t rowComp = row / nCoefs;
                const index_t colCoef = col % nCoefs;
                const index_t colComp = col / nCoefs;

                if (fixed[rowCoef])
                    continue;

                const index_t r = freeMap[static_cast<size_t>(row)];
                if (fixed[colCoef])
                    rhs(r, 0) -= value * values(colCoef, colComp);
                else
                    triplets.push_back(gsEigen::Triplet<T,index_t>(
                        r, freeMap[static_cast<size_t>(col)], value));
            }
        }

        gsSparseMatrix<T> matrix(nFree, nFree);
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        matrix.makeCompressed();

        gsSparseSolver<>::SimplicialLDLT solver(matrix);
        gsMatrix<T> freeSol = solver.solve(rhs);

        disp = values;
        for (index_t i = 0; i < nCoefs; ++i)
            for (index_t d = 0; d < dim; ++d)
                if (!fixed[i])
                    disp(i, d) = freeSol(freeMap[static_cast<size_t>(d * nCoefs + i)], 0);
    }

    void collectConstraints(index_t patch, bool incremental,
                            std::vector<char>& fixed, gsMatrix<T>& values) const
    {
        const gsMatrix<T>& orig = m_originalGeometry.patch(patch).coefs();
        const index_t nCoefs = orig.rows();
        const index_t dim = orig.cols();

        fixed.assign(static_cast<size_t>(nCoefs), 0);
        values.setZero(nCoefs, dim);

        gsMatrix<index_t> allBoundary = m_originalGeometry.patch(patch).basis().allBoundary();
        for (index_t i = 0; i < allBoundary.size(); ++i)
        {
            const index_t id = allBoundary(i);
            if (id >= 0 && id < nCoefs)
                fixed[id] = 1;
        }

        for (size_t k = 0; k < m_movingBoundaries.size(); ++k)
        {
            if (m_movingBoundaries[k].patch != patch)
                continue;

            gsMatrix<index_t> ids =
                m_originalGeometry.patch(patch).basis().boundary(m_movingBoundaries[k].side());
            for (index_t i = 0; i < ids.size(); ++i)
            {
                const index_t id = ids(i);
                if (id < 0 || id >= nCoefs)
                    continue;

                fixed[id] = 1;
                for (index_t d = 0; d < dim && d < m_boundaryDisplacement.patch(patch).coefs().cols(); ++d)
                    values(id, d) = m_boundaryDisplacement.patch(patch).coefs()(id, d);

                if (incremental && patch < static_cast<index_t>(m_meshDisplacement.nPatches()))
                    for (index_t d = 0; d < dim && d < m_meshDisplacement.patch(patch).coefs().cols(); ++d)
                        values(id, d) -= m_meshDisplacement.patch(patch).coefs()(id, d);
            }
        }

        if (m_movingBoundaries.empty())
            makeBoundaryConstraintValues(patch, incremental, values);
    }

    void makeBoundaryConstraintValues(index_t patch, bool incremental, gsMatrix<T>& values) const
    {
        const index_t rows = std::min(values.rows(), m_boundaryDisplacement.patch(patch).coefs().rows());
        const index_t cols = std::min(values.cols(), m_boundaryDisplacement.patch(patch).coefs().cols());
        for (index_t i = 0; i < rows; ++i)
            for (index_t d = 0; d < cols; ++d)
                values(i, d) = m_boundaryDisplacement.patch(patch).coefs()(i, d);

        if (incremental && patch < static_cast<index_t>(m_meshDisplacement.nPatches()))
        {
            const index_t r = std::min(values.rows(), m_meshDisplacement.patch(patch).coefs().rows());
            const index_t c = std::min(values.cols(), m_meshDisplacement.patch(patch).coefs().cols());
            for (index_t i = 0; i < r; ++i)
                for (index_t d = 0; d < c; ++d)
                    values(i, d) -= m_meshDisplacement.patch(patch).coefs()(i, d);
        }
    }

    void assembleScalarExtensionOperator(index_t patch,
                                         bool biharmonic,
                                         gsSparseMatrix<T>& op) const
    {
        const gsGeometry<T>& geo = m_originalGeometry.patch(patch);
        const gsBasis<T>& basis = geo.basis();
        const index_t nCoefs = geo.coefs().rows();

        if (biharmonic)
            for (index_t d = 0; d < basis.dim(); ++d)
                GISMO_ENSURE(basis.degree(d) >= 2,
                             "Biharmonic mesh deformation requires degree at least 2 in every direction.");

        std::vector<gsEigen::Triplet<T,index_t> > triplets;
        gsQuadRule<T> rule;
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        gsMapData<T> md(biharmonic ? NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER
                                   : NEED_MEASURE | NEED_GRAD_TRANSFORM);
        gsMatrix<index_t> actives;
        std::vector<gsMatrix<T> > basisData;
        gsMatrix<T> physGrad, physLaplace;
        gsOptionList quadOptions = gsAssembler<T>::defaultOptions();

        rule = gsQuadrature::get(*basis.domain(), quadOptions);
        for (typename gsBasis<T>::domainIter domIt = basis.domain()->beginAll();
             domIt < basis.domain()->endAll(); ++domIt)
        {
            rule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);

            md.points = quNodes;
            basis.active_into(md.points.col(0), actives);
            basis.evalAllDers_into(md.points, biharmonic ? 2 : 1, basisData, true);
            geo.computeMap(md);

            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                const T weight = quWeights[q] * md.measure(q);
                if (biharmonic)
                {
                    transformLaplaceHgrad(md, q, basisData[1], basisData[2], physLaplace);
                    for (index_t i = 0; i < actives.rows(); ++i)
                        for (index_t j = 0; j < actives.rows(); ++j)
                            triplets.push_back(gsEigen::Triplet<T,index_t>(
                                actives(i), actives(j), weight * physLaplace(0, i) * physLaplace(0, j)));
                }
                else
                {
                    transformGradients(md, q, basisData[1], physGrad);
                    for (index_t i = 0; i < actives.rows(); ++i)
                        for (index_t j = 0; j < actives.rows(); ++j)
                            triplets.push_back(gsEigen::Triplet<T,index_t>(
                                actives(i), actives(j), weight * physGrad.col(i).dot(physGrad.col(j))));
                }
            }
        }

        op.resize(nCoefs, nCoefs);
        op.setFromTriplets(triplets.begin(), triplets.end());
        op.makeCompressed();
    }

    void assembleLinearElasticityOperator(index_t patch, gsSparseMatrix<T>& op) const
    {
        const gsGeometry<T>& geo = m_originalGeometry.patch(patch);
        const gsBasis<T>& basis = geo.basis();
        const index_t nCoefs = geo.coefs().rows();
        const index_t dim = geo.coefs().cols();
        const T nu = static_cast<T>(m_options.getReal("PoissonsRatio"));
        const T lambda = nu / ((T(1) + nu) * (T(1) - T(2) * nu));
        const T mu = T(1) / (T(2) * (T(1) + nu));

        GISMO_ENSURE(dim == basis.dim(), "Linear elasticity mesh deformation requires matching geometry and parametric dimensions.");
        GISMO_ENSURE(nu > T(-1) && nu < T(0.5), "Poisson's ratio must be in (-1, 0.5) for linear elasticity.");

        std::vector<gsEigen::Triplet<T,index_t> > triplets;
        gsQuadRule<T> rule;
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        gsMapData<T> md(NEED_MEASURE | NEED_GRAD_TRANSFORM);
        gsMatrix<index_t> actives;
        std::vector<gsMatrix<T> > basisData;
        gsMatrix<T> physGrad;
        gsOptionList quadOptions = gsAssembler<T>::defaultOptions();

        rule = gsQuadrature::get(*basis.domain(), quadOptions);
        for (typename gsBasis<T>::domainIter domIt = basis.domain()->beginAll();
             domIt < basis.domain()->endAll(); ++domIt)
        {
            rule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);

            md.points = quNodes;
            basis.active_into(md.points.col(0), actives);
            basis.evalAllDers_into(md.points, 1, basisData, true);
            geo.computeMap(md);

            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                const T weight = quWeights[q] * md.measure(q);
                transformGradients(md, q, basisData[1], physGrad);

                for (index_t i = 0; i < actives.rows(); ++i)
                {
                    for (index_t j = 0; j < actives.rows(); ++j)
                    {
                        const T gradDot = physGrad.col(i).dot(physGrad.col(j));
                        for (index_t a = 0; a < dim; ++a)
                        {
                            for (index_t b = 0; b < dim; ++b)
                            {
                                const T value =
                                    (a == b ? mu * gradDot : T(0)) +
                                    mu * physGrad(b, i) * physGrad(a, j) +
                                    lambda * physGrad(a, i) * physGrad(b, j);
                                triplets.push_back(gsEigen::Triplet<T,index_t>(
                                    a * nCoefs + actives(i),
                                    b * nCoefs + actives(j),
                                    weight * value));
                            }
                        }
                    }
                }
            }
        }

        op.resize(dim * nCoefs, dim * nCoefs);
        op.setFromTriplets(triplets.begin(), triplets.end());
        op.makeCompressed();
    }

    void addToDisplacement(const gsMultiPatch<T>& update)
    {
        const size_t nPatches = std::min(m_meshDisplacement.nPatches(), update.nPatches());
        for (size_t p = 0; p < nPatches; ++p)
            m_meshDisplacement.patch(p).coefs() += update.patch(p).coefs();
    }

    void applyDisplacementToGeometry()
    {
        const size_t nPatches = std::min(m_geometry.nPatches(), m_meshDisplacement.nPatches());
        for (size_t p = 0; p < nPatches; ++p)
        {
            if (p >= m_originalGeometry.nPatches())
                continue;
            if (m_originalGeometry.patch(p).coefs().rows() != m_meshDisplacement.patch(p).coefs().rows())
                continue;
            m_geometry.patch(p).coefs() =
                m_originalGeometry.patch(p).coefs() + m_meshDisplacement.patch(p).coefs();
        }
    }

    index_t checkCurrentGeometry() const
    {
        for (size_t p = 0; p < m_geometry.nPatches(); ++p)
        {
            gsMatrix<T> anchors = m_geometry.patch(p).basis().anchors();
            if (anchors.cols() == 0)
                continue;

            gsMapData<T> md(NEED_DERIV);
            md.points = anchors;
            m_geometry.patch(p).computeMap(md);

            for (index_t q = 0; q < anchors.cols(); ++q)
            {
                const gsMatrix<T> J = md.jacobian(q);
                if (J.rows() == J.cols() && J.determinant() <= T(0))
                    return static_cast<index_t>(p);
            }
        }
        return -1;
    }

    bool isIncrementalMethod() const
    {
        return m_method == mesh_deform_method::IHE ||
               m_method == mesh_deform_method::ILE ||
               m_method == mesh_deform_method::IBHE;
    }

    bool isBiharmonicMethod() const
    {
        return m_method == mesh_deform_method::BHE ||
               m_method == mesh_deform_method::IBHE;
    }

protected:
    gsMultiPatch<T>& m_geometry;
    gsMultiPatch<T> m_originalGeometry;
    gsMultiPatch<T> m_boundaryDisplacement;
    gsMultiPatch<T> m_meshDisplacement;
    gsMultiPatch<T> m_savedGeometry;
    gsMultiPatch<T> m_savedMeshDisplacement;
    gsMultiPatch<T> m_savedBoundaryDisplacement;

    std::vector<patchSide> m_movingBoundaries;
    mesh_deform_method::method m_method;
    gsOptionList m_options;
    const gsMultiBasis<T>* m_solverVelBasis;
    bool m_isInitialized;
    bool m_hasSavedState;
    index_t m_lastCheckResult;
    index_t m_savedCheckResult;
};

} // namespace gismo
