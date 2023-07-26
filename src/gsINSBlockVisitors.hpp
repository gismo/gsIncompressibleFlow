/** @file gsINSBlockVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sourek, H. Hornikova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSBlockVisitors.h>

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

// === PARENT === //

template <class T>
void gsINSBlockVisitor<T>::initialize(const gsBasisRefs<T>& basisRefs,
    const index_t patchIndex,
    const gsAssemblerOptions& options,
    gsQuadRule<T>& rule,
    unsigned& evFlags)
{
    const gsBasis<T>& basis = basisRefs.front();

    m_dim = basis.dim();
    m_patchIndex = patchIndex;

    // Setup Quadrature
    rule = gsGaussRule<T>(basis, options.quA, options.quB);

    initializeSpecific(evFlags);

    m_diffusionCoeff.setZero(1, rule.numNodes());

    for (int i = 0; i < rule.numNodes(); i++)
        m_diffusionCoeff(0, i) = m_viscosity;
}

template <class T>
void gsINSBlockVisitor<T>::initialize(const gsBasisRefs<T>& basisRefs,
    const index_t patchIndex,
    gsQuadRule<T>& rule,
    unsigned& evFlags)
{
    const gsBasis<T>& basis = basisRefs.front();

    m_dim = basis.dim();
    m_patchIndex = patchIndex;

    gsVector<index_t> numQuadNodes(m_dim);

    for (int i = 0; i < m_dim; ++i) // to do: improve
        numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

    // Setup Quadrature
    rule = gsGaussRule<T>(numQuadNodes);

    initializeSpecific(evFlags);

    m_diffusionCoeff.setZero(1, rule.numNodes());

    for (int i = 0; i < rule.numNodes(); i++)
        m_diffusionCoeff(0, i) = m_viscosity;
}

// === BLOCK A === //

template <class T>
void gsINSBlockAVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    // Evaluate basis functions on element  nodes
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.front().deriv_into(quNodes, m_basisGradsU);
}

template <class T>
void gsINSBlockAVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();

    m_localMat.setZero(numActU, numActU);//local_A

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, m_basisGradsU, m_physGradU);

        // Local block A
        m_localMat.template triangularView<gsEigen::Upper>() += weight * m_viscosity * (m_physGradU.transpose() * m_physGradU);
    }
}

template <class T>
void gsINSBlockAVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    const index_t usz = m_Umap.freeSize();

    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActU; ++j) // Build A-part of the matrix
            {
                const int jj = m_activesU(j);
                if (m_Umap.is_free_index(jj))
                {
                    if (j >= i) {
                        sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                    }
                }
                else // m_Umap.is_boundary_index(jj)
                {
                    const int bb = m_Umap.global_to_bindex(jj);
                    for (index_t s = 0; s != m_dim; ++s)
                        rhs(ii + s*usz, 0) -= // assuming single rhs
                        m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[0](bb, s);
                }
            }
        }
    }
}


// === BLOCKS B === //

template <class T>
void gsINSBlocksBVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    // Evaluate basis functions on element  nodes
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.back().active_into(quNodes.col(0), m_activesP);

    basisRefs.front().deriv_into(quNodes, m_basisGradsU);
    basisRefs.back().eval_into(quNodes, m_basisValsP);
}

template <class T>
void gsINSBlocksBVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();
    const index_t numActP = m_activesP.rows();

    m_localMat.resize(m_dim);
    for (index_t i = 0; i != m_dim; ++i)
        m_localMat[i].setZero(numActP, numActU);//local_B_i

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, m_basisGradsU, m_physGradU);

        // Local blocks B_i
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].noalias() += weight * (m_basisValsP.col(k) * m_physGradU.row(i));
    }
}

template <class T>
void gsINSBlocksBVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    const index_t usz = m_Umap.freeSize();
    const index_t ps = m_dim*usz;

    // Local Dofs to global dofs
    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActP; ++j) // Build B-part of the matrix
            {
                const int jj = m_activesP(j);
                if (m_Pmap.is_free_index(jj))
                {
                    for (index_t s = 0; s != m_dim; ++s)
                        sysBlock.coeffRef(jj, ii + s * usz) += m_localMat[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj)
                {
                    const int bb = m_Pmap.global_to_bindex(jj);
                    for (index_t s = 0; s < m_dim; ++s)
                        rhs(s * usz + ii, 0) += m_localMat[s](j, i) * eliminatedDofs[1](bb, 0);
                }
            }

        }
        else // m_Umap.is_boundary_index(ii)
        {
            const int bb = m_Umap.global_to_bindex(ii);
            for (index_t k = 0; k < numActP; ++k)
            {
                const int kk = m_activesP(k);
                if (m_Pmap.is_free_index(kk))
                {
                    T tmp = m_localMat[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += m_localMat[s](k, i) * eliminatedDofs[0](bb, s);
                    rhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
        }
    }
}


// === BLOCKS C === //

template <class T>
void gsINSBlocksCVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    // Evaluate basis functions on element  nodes
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.back().active_into(quNodes.col(0), m_activesP);

    basisRefs.front().eval_into(quNodes, m_basisValsU);
    basisRefs.back().deriv_into(quNodes, m_basisGradsP);
}

template <class T>
void gsINSBlocksCVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();
    const index_t numActP = m_activesP.rows();

    m_localMat.resize(m_dim);
    for (index_t i = 0; i != m_dim; ++i)
        m_localMat[i].setZero(numActP, numActU);//local_B_i

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, m_basisGradsP, m_physGradP);

        // Local blocks B_i
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].noalias() += weight * (m_physGradP.row(i).transpose() * m_basisValsU.col(k).transpose());
    }
}

template <class T>
void gsINSBlocksCVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    const index_t usz = m_Umap.freeSize();
    const index_t ps = m_dim*usz;

    // Local Dofs to global dofs
    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {

            for (index_t j = 0; j < numActP; ++j) // Build B-part of the matrix
            {
                const int jj = m_activesP(j);
                if (m_Pmap.is_free_index(jj))
                {
                    for (index_t s = 0; s != m_dim; ++s)
                        sysBlock.coeffRef(jj, ii + s * usz) += m_localMat[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj)
                {
                    const int bb = m_Pmap.global_to_bindex(jj);
                    for (index_t s = 0; s < m_dim; ++s)
                        rhs(s * usz + ii, 0) += m_localMat[s](j, i) * eliminatedDofs[1](bb, 0);
                }
            }

        }
        else // m_Umap.is_boundary_index(ii)
        {
            const int bb = m_Umap.global_to_bindex(ii);
            for (index_t k = 0; k < numActP; ++k)
            {
                const int kk = m_activesP(k);
                if (m_Pmap.is_free_index(kk)) {
                    T tmp = m_localMat[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += m_localMat[s](k, i) * eliminatedDofs[0](bb, s);
                    rhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
        }
    }
}


// === BLOCK N === //

template <class T>
void gsINSBlockNVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    // Evaluate basis functions on element nodes
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);

    // Evaluate solution on element nodes
    m_solUVals = m_solU.value(quNodes, m_patchIndex);
}

template <class T>
void gsINSBlockNVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();

    m_localMat.setZero(numActU, numActU);

    const gsMatrix<T> & basisValsU = m_basisDataU[0];
    const gsMatrix<T> & basisGradsU = m_basisDataU[1];

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, basisGradsU, m_physGradU);

        // Local block
        m_localMat.noalias() += weight * (basisValsU.col(k) * (m_solUVals.col(k).transpose() * m_physGradU));
    }
}

template <class T>
void gsINSBlockNVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    const index_t usz = m_Umap.freeSize();

    // Local Dofs to global dofs
    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActU; ++j)
            {
                const int jj = m_activesU(j);
                if (m_Umap.is_free_index(jj))
                {
                    sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                }
                else // m_Umap.is_boundary_index(jj)
                {
                    const int bb = m_Umap.global_to_bindex(jj);
                    for (index_t s = 0; s != m_dim; ++s)
                        rhs(ii + s*usz, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                }
            }
        }
    }
}


// === BLOCK M === //

template <class T>
void gsINSBlockMVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    m_basisDataU.resize(1);

    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.front().eval_into(quNodes, m_basisDataU.front());
}

template <class T>
void gsINSBlockMVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();

    m_localMat.setZero(numActU, numActU);//local_A

    const gsMatrix<T> & basisValsU = m_basisDataU.front();

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        const T weight = quWeights(k) * mapData.measure(k);

        m_localMat.template triangularView<gsEigen::Upper>() += weight * (basisValsU.col(k) * basisValsU.col(k).transpose());
    }
}

template <class T>
void gsINSBlockMVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    const index_t usz = m_Umap.freeSize();

    // Local Dofs to global dofs
    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActU; ++j) // Build A-part of the matrix
            {
                const int jj = m_activesU(j);
                if (m_Umap.is_free_index(jj))
                {
                    if (j >= i) {
                        sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                    }
                }
                else // m_Umap.is_boundary_index(jj)
                {
                    const int bb = m_Umap.global_to_bindex(jj);
                    for (index_t s = 0; s != m_dim; ++s)
                        rhs(ii + s*usz, 0) -= // assuming single rhs
                        m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[0](bb, s);
                }
            }
        }
    }
}


// === BLOCK Ap === //

template <class T>
void gsINSBlockApVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    basisRefs.back().active_into(quNodes.col(0), m_activesP);
    basisRefs.back().deriv_into(quNodes, m_basisGradsP);
}

template <class T>
void gsINSBlockApVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActP = m_activesP.rows();

    m_localMat.setZero(numActP, numActP);//local_A

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        const T weight = quWeights(k) * mapData.measure(k);

        transformGradients(mapData, k, m_basisGradsP, m_physGradP);

        m_localMat.template triangularView<gsEigen::Upper>() += weight * (m_physGradP.transpose() * m_physGradP);
    }
}

template <class T>
void gsINSBlockApVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    for (index_t i = 0; i < numActP; ++i)
    {
        const int ii = m_activesP(i);
        if (m_Pmap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActP; ++j) // Build A-part of the matrix
            {
                const int jj = m_activesP(j);
                if (m_Pmap.is_free_index(jj))
                {
                    if (j >= i) {
                        sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                    }
                }
                else // m_Pmap.is_boundary_index(jj)
                {
                    const int bb = m_Pmap.global_to_bindex(jj);
                    rhs(ii, 0) -= m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }
}


// === BLOCK Np === //

template <class T>
void gsINSBlockNpVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    basisRefs.back().active_into(quNodes.col(0), m_activesP);
    basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataP);

    // Evaluate solution on element nodes
    m_solUVals = m_solU.value(quNodes, m_patchIndex);
}

template <class T>
void gsINSBlockNpVisitor<T>::assemble(const gsMapData<T>& mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActP = m_activesP.rows();

    m_localMat.setZero(numActP, numActP);

    const gsMatrix<T> & basisValsP = m_basisDataP[0];
    const gsMatrix<T> & basisGradsP = m_basisDataP[1];

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, basisGradsP, m_physGradP);

        // Local block
        m_localMat.noalias() += weight * (basisValsP.col(k) * (m_solUVals.col(k).transpose() * m_physGradP));
    }
}

template <class T>
void gsINSBlockNpVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    for (index_t i = 0; i < numActP; ++i)
    {
        const int ii = m_activesP(i);
        if (m_Pmap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActP; ++j)
            {
                const int jj = m_activesP(j);
                if (m_Pmap.is_free_index(jj))
                {
                    sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                }
                else // m_Pmap.is_boundary_index(jj)
                {
                    const int bb = m_Pmap.global_to_bindex(jj);
                    rhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }
}

// === BLOCK Mp === //

template <class T>
void gsINSBlockMpVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, gsMatrix<T>& quNodes)
{
    basisRefs.back().active_into(quNodes.col(0), m_activesP);
    basisRefs.back().eval_into(quNodes, bVals_p);
}

template <class T>
void gsINSBlockMpVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActP = m_activesP.rows();

    m_localMat.setZero(numActP, numActP);//local_A

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        const T weight = quWeights(k) * mapData.measure(k);

        m_localMat.template triangularView<gsEigen::Upper>() += weight * (bVals_p.col(k) * bVals_p.col(k).transpose());
    }
}

template <class T>
void gsINSBlockMpVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & sysBlock,
    gsMatrix<T>           & rhs)
{
    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    for (index_t i = 0; i < numActP; ++i)
    {
        const int ii = m_activesP(i);
        if (m_Pmap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActP; ++j) // Build A-part of the matrix
            {
                const int jj = m_activesP(j);
                if (m_Pmap.is_free_index(jj))
                {
                    if (j >= i) {
                        sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                    }
                }
                else // m_Pmap.is_boundary_index(jj)
                {
                    const int bb = m_Pmap.global_to_bindex(jj);
                    rhs(ii, 0) -= m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }
}


// === RHS === //

template <class T>
void gsINSRhsVisitor<T>::evaluate(const gsBasisRefs<T>& basisRefs, const gsMapData<T> & mapData, gsMatrix<T>& quNodes)
{
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.front().eval_into(quNodes, m_basisValsU);

    // Evaluate right-hand side at the geometry points
    m_pRhsFcn->eval_into( mapData.values[0], m_rhsVals );
}

template <class T>
void gsINSRhsVisitor<T>::assemble(const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    const index_t numActU = m_activesU.rows();

    m_localRhsU.setZero(numActU, m_rhsVals.rows());

    const index_t nQuPoints = quWeights.rows();
    for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
    {
        // weight * abs(det J), where J is geometry Jacobian.
        const T weight = quWeights(k) * mapData.measure(k);

        // Right-hand side
        m_localRhsU.noalias() += weight * (m_basisValsU.col(k) *  m_rhsVals.col(k).transpose());
    }
}

template <class T>
void gsINSRhsVisitor<T>::localToGlobal(gsMatrix<T> & rhs)
{
    const index_t usz = m_Umap.freeSize();

    // Local Dofs to global dofs
    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    for (index_t i = 0; i < numActU; ++i)
    {
        const int ii = m_activesU(i);
        if (m_Umap.is_free_index(ii))
        {
            for (index_t s = 0; s != m_dim; ++s)
                rhs(ii + s*usz, 0) += m_localRhsU(i, s);
        }
    }
}


} // namespace gismo

