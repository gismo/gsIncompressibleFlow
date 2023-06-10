/** @file gsINSBlockVisitorsBnd.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s) : H.Hornikova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSBlockVisitorsBnd.h>

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo 
{

// === PARENT === //

template <class T>
void gsINSBlockVisitorBnd<T>::initialize(gsBasisRefs<T> const& basisRefs,
    const index_t patchIndex,
    gsQuadRule<T>& rule,
    unsigned& evFlags,
    boxSide side)
{
    const gsBasis<T>& basis = basisRefs.front();

    m_dim = basis.dim();
    m_patchIndex = patchIndex;

    m_side = side;

    const int dir = m_side.direction();
    gsVector<int> numQuadNodes(basis.dim());
    for (int i = 0; i < basis.dim(); ++i)
        numQuadNodes[i] = 2 * basis.degree(i) + 1;
    numQuadNodes[dir] = 1;

    // Setup Quadrature
    rule = gsGaussRule<T>(numQuadNodes);

    this->initializeSpecific(evFlags);
}


// === PCD Robin === //

template <class T>
void gsINSBlockVisitorRobinPCD<T>::evaluate(const gsBasisRefs<T>& basisRefs, const gsMapData<T> & mapData, gsMatrix<T>& quNodes)
{
    // Compute the active basis functions
    // Assumes actives are the same for all quadrature points on the current element
    basisRefs.back().active_into(quNodes.col(0), m_actives);
    const index_t numActive = m_actives.rows();

    basisRefs.back().eval_into(quNodes, m_bVals);

    // Evaluate solution on element nodes
    m_solUVals = m_solU.value(quNodes, m_patchIndex);

    // Initialize local matrix/rhs
    m_localMat.setZero(numActive, numActive);
}

template <class T>
void gsINSBlockVisitorRobinPCD<T>::assemble(gsDomainIterator<T>& element,
    const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {
        // Compute the outer normal vector on the side
        outerNormal(mapData, k, m_side, m_unormal);

        // Multiply quadrature weight by the geometry measure
        const T weight = quWeights[k] * m_unormal.norm();

        // Compute the unit normal vector 
        m_unormal.normalize();

        T coef = - m_solUVals.col(k).transpose() * m_unormal;

        m_localMat.noalias() += weight * coef * (m_bVals.col(k) * m_bVals.col(k).transpose());
    }
}

template <class T>
void gsINSBlockVisitorRobinPCD<T>::localToGlobal(gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{

    // Local Dofs to global dofs
    m_Pmap.localToGlobal(m_actives, m_patchIndex, m_actives);
    const index_t numActP = m_actives.rows();

    for (index_t i = 0; i < numActP; ++i)
    {
        const int ii = m_actives(i);
        if (m_Pmap.is_free_index(ii))
        {
            for (index_t j = 0; j < numActP; ++j)
            {
                const int jj = m_actives(j);
                if (m_Pmap.is_free_index(jj))
                    blockMatrix.coeffRef(ii, jj) += m_localMat(i, j);
            }
        }
    }
}


// === NITSCHE PARENT === //

template <class T>
void gsINSBlockVisitorNitsche<T>::initialize(gsBasisRefs<T> const& basisRefs,
    const index_t patchIndex,
    gsQuadRule<T>& rule,
    unsigned& evFlags,
    const boundary_condition<T> & dirCond)
{
    Base::initialize(basisRefs, patchIndex, rule, evFlags, dirCond.side());

    const gsBasis<T>& basis = basisRefs.front();

    const int deg = basis.maxDegree();
    m_penalty = (deg + basis.dim()) * (deg + 1) * T(2.5);
    m_pDirFcn = dirCond.function().get();
}


// === NITSCHE BLOCK A === //

template <class T>
void gsINSBlockAVisitorNitsche<T>::evaluate(const gsBasisRefs<T>& basisRefs, const gsMapData<T> & mapData, gsMatrix<T>& quNodes)
{
    // Compute the active basis functions
    // Assumes actives are the same for all quadrature points on the current element
    basisRefs.front().active_into(quNodes.col(0), m_actives);
    const index_t numActive = m_actives.rows();

    // Evaluate basis values and derivatives on element
    basisRefs.front().evalAllDers_into(quNodes, 1, m_basisData);

    // Evaluate the Dirichlet data
    m_pDirFcn->eval_into(mapData.values[0], m_dirData);

    // Initialize local matrix/rhs
    m_localMat.setZero(numActive, numActive);
    m_localRhs.setZero(numActive, m_pDirFcn->targetDim());

}

template <class T>
void gsINSBlockAVisitorNitsche<T>::assemble(gsDomainIterator<T>& element,
    const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    gsMatrix<T> & bGrads = m_basisData[1];
    

    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        const gsMatrix<T> bVals = m_basisData[0].col(k);

        // Compute the outer normal vector on the side
        outerNormal(mapData, k, m_side, m_unormal);

        // Multiply quadrature weight by the geometry measure
        const T weight = quWeights[k] * m_unormal.norm();

        // Compute the unit normal vector 
        m_unormal.normalize();

        // Compute physical gradients at k as a Dim x NumActive matrix
        transformGradients(mapData, k, bGrads, m_physGrad);

        // Get penalty parameter
        const T h = element.getCellSize();
        const T gamma = m_penalty / (0 != h ? h : 1);

        // Sum up quadrature point evaluations
        m_localRhs.noalias() -= weight * m_viscosity * ((m_physGrad.transpose() * m_unormal - gamma * bVals)
            * m_dirData.col(k).transpose());

        m_localMat.noalias() -= weight * m_viscosity * (bVals * m_unormal.transpose() * m_physGrad
            + (bVals * m_unormal.transpose() * m_physGrad).transpose()
            - gamma * bVals * bVals.transpose());
    }
}

template <class T>
void gsINSBlockAVisitorNitsche<T>::localToGlobal( gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    const index_t usz = m_Umap.freeSize();
    m_Umap.localToGlobal(m_actives, m_patchIndex, m_actives);
    const index_t numActive = m_actives.rows();

    // Push element contribution to the global matrix and load vector
    for (index_t i = 0; i != numActive; ++i)
    {
        const unsigned ii = m_actives(i);

        for (index_t s = 0; s != m_dim; s++)
            blockRhs(ii + s*usz, 0) += m_localRhs(i,s);

        for (index_t j = 0; j != numActive; ++j)
        {
            const unsigned jj = m_actives(j);
            blockMatrix.coeffRef(ii, jj) += m_localMat(i, j);
        }
    }

}


// === NITSCHE BLOCK B === //

template <class T>
void gsINSBlockBVisitorNitsche<T>::evaluate(const gsBasisRefs<T>& basisRefs, const gsMapData<T> & mapData, gsMatrix<T>& quNodes)
{
    // Compute the active basis functions
    // Assumes actives are the same for all quadrature points on the current element
    basisRefs.front().active_into(quNodes.col(0), m_activesU);
    basisRefs.back().active_into(quNodes.col(0), m_activesP);

    const index_t numActU = m_activesU.rows();
    const index_t numActP = m_activesP.rows();

    basisRefs.front().eval_into(quNodes, m_bValsU);
    basisRefs.back().eval_into(quNodes, m_bValsP);

    // Evaluate the Dirichlet data
    m_pDirFcn->eval_into(mapData.values[0], m_dirData);

    // Initialize local matrix/rhs
    m_localMatVec.resize(m_dim);
    m_localRhs.setZero(numActP, 1);

    for (index_t i = 0; i != m_dim; ++i)
        m_localMatVec[i].setZero(numActP, numActU); //local_B_i

}

template <class T>
void gsINSBlockBVisitorNitsche<T>::assemble(gsDomainIterator<T>& element,
    const gsMapData<T> & mapData,
    const gsVector<T>& quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {
        // Compute the outer normal vector on the side
        outerNormal(mapData, k, m_side, m_unormal);

        // Multiply quadrature weight by the geometry measure
        const T weight = quWeights[k] * m_unormal.norm();

        // Compute the unit normal vector 
        m_unormal.normalize();

        m_localRhs.noalias() -= weight * m_bValsP.col(k) * (m_dirData.col(k).transpose() * m_unormal);

        // Sum up quadrature point evaluations
        for (index_t i = 0; i != m_dim; ++i)
        {
            const T w = weight * m_unormal(i);

            m_localMatVec[i].noalias() -= w * (m_bValsP.col(k) * m_bValsU.col(k).transpose());
        }
    }
}

template <class T>
void gsINSBlockBVisitorNitsche<T>::localToGlobal(gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    const index_t usz = m_Umap.freeSize();
    const index_t pshift = m_dim * usz;

    m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
    const index_t numActU = m_activesU.rows();

    m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
    const index_t numActP = m_activesP.rows();

    // Push element contribution to the global matrix and load vector
    for (index_t i = 0; i != numActP; ++i)
    {
        const unsigned ii = m_activesP(i);

        if (m_Pmap.is_free_index(ii))
        {

            blockRhs(pshift + ii, 0) += m_localRhs(i, 0);

            for (index_t j = 0; j != numActU; ++j)
            {
                const unsigned jj = m_activesU(j);

                for (index_t s = 0; s != m_dim; s++)
                    blockMatrix.coeffRef(ii, jj + s * usz) += m_localMatVec[s](i, j);
            }
        }
    }
}

} // end namespace gismo
