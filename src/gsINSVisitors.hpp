/** @file gsINSVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSVisitors.h>

namespace gismo
{


template <class T, int MatOrder>
void gsINSVisitorUU<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);


    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < nComponents; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_localMat(i, j) * eliminatedDofs[m_shapeUnkID](bb, d);
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPU<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();


    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPU_withUPrhs<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                }
            }
        }
        else // part arising from block B (assuming that the offdiag. blocks are symmetric)
        {
            const int bb = m_dofMappers[m_testUnkID].global_to_bindex(ii);
            for (index_t k = 0; k < numActShape; k++)
            {
                const int kk = m_shapeFunActives(k);

                if (m_dofMappers[m_shapeUnkID].is_free_index(kk))
                {
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp += m_locMatVec[d](i, k) * eliminatedDofs[m_testUnkID](bb, d);

                    globalRhs(dim*uCompSize + kk, 0) += tmp;
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorUP<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_shapeUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii, jj + d*uCompSize) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);
                    
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp -= m_locMatVec[d](i, j) * eliminatedDofs[m_shapeUnkID](bb, d);

                    globalRhs(ii, 0) += tmp;
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPP<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const int ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const int jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorRhsU<T, MatOrder>::assemble()
{
    m_localMat.setZero(m_testFunActives.rows(), m_paramsPtr->getPde().dim());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}


template <class T, int MatOrder>
void gsINSVisitorRhsU<T, MatOrder>::localToGlobal(gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[0].freeSize(); // number of dofs for one velocity component

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

    index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t d = 0; d != dim; d++)
                globalRhs(ii + d*uCompSize, 0) += m_localMat(i, d);
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorRhsP<T, MatOrder>::assemble()
{
    m_localMat.setZero(m_testFunActives.rows(), 1);

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}


template <class T, int MatOrder>
void gsINSVisitorRhsP<T, MatOrder>::localToGlobal(gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[0].freeSize(); // number of dofs for one velocity component

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

   index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
            globalRhs(dim*uCompSize + ii, 0) += m_localMat(i);
    }
} 

} // namespace gismo