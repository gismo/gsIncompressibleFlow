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
void gsINSVisitorUU<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal_nonper.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);


    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    for (index_t d = 0; d < nComponents; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, d);
                }
            }
        }
    }
} 


template <class T, int MatOrder>
void gsINSVisitorUU<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal_per, matrix has to contain all components.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            bool iiElim = m_paramsPtr->getPerHelperPtr()->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr()->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    bool jjElim = m_paramsPtr->getPerHelperPtr()->isEliminated(jj);
                    const index_t jjMapped = m_paramsPtr->getPerHelperPtr()->map(jj);

                    // ii and jj are not eliminated periodic dofs:
                    if (!iiElim && !jjElim) 
                    {
                        for (index_t d = 0; d < nComponents; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped + d*uCompSize) += m_localMat(i, j);
                    }
                    // only ii or jj is eliminated periodic dof:
                    else if ( (!iiElim && jjElim) || (iiElim && !jjElim) )
                    {
                        for (int r = 0; r < nComponents; r++)
                        {
                            for (int s = 0; s < nComponents; s++)
                            {
                                T tmp = 0;

                                if (jjElim)
                                    tmp = m_periodicTransformMat(s, r);
                                else if (iiElim)
                                    tmp = m_periodicTransformMat(r, s);

                                tmp *= m_localMat(i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + r*uCompSize, jjMapped + s*uCompSize) += tmp;
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int q = 0; q < nComponents; q++)
                        {
                            for (int r = 0; r < nComponents; r++)
                            {
                                for (int s = 0; s < nComponents; s++)
                                {
                                    T tmp = m_periodicTransformMat(q, s) * m_periodicTransformMat(r, s) * m_localMat(i, j);
                                    
                                    if (tmp != 0)
                                        globalMat.coeffRef(iiMapped + q*uCompSize, jjMapped + r*uCompSize) += tmp;
                                }
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalRhs(iiMapped + d*uCompSize, 0) -= m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, d);
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            for (int t = 0; t < nComponents; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, s);

                                if (tmp != 0)
                                    globalRhs(iiMapped + t*uCompSize, 0) -= tmp;
                            }
                        }
                    }
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorUUrotation<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsINSVisitorUUrotation::localToGlobal_nonper.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);


    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj + uCompSize) -= m_omega * m_localMat(i, j);
                    globalMat.coeffRef(ii + uCompSize, jj) += m_omega * m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);
                    globalRhs(ii, 0) += m_omega * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 1);
                    globalRhs(ii + uCompSize, 0) -= m_omega * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }
} 


template <class T, int MatOrder>
void gsINSVisitorUUrotation<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsINSVisitorUUrotation::localToGlobal_per, matrix has to contain all components.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            bool iiElim = m_paramsPtr->getPerHelperPtr()->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr()->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    bool jjElim = m_paramsPtr->getPerHelperPtr()->isEliminated(jj);
                    const index_t jjMapped = m_paramsPtr->getPerHelperPtr()->map(jj);

                    // ii and jj are not eliminated periodic dofs:
                    if (!iiElim && !jjElim) 
                    {
                        globalMat.coeffRef(iiMapped, jjMapped + uCompSize) -= m_omega * m_localMat(i, j);
                        globalMat.coeffRef(iiMapped + uCompSize, jjMapped) += m_omega * m_localMat(i, j);
                    }
                    // only ii is eliminated periodic dof:
                    else if ( iiElim && !jjElim )
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            T tmp0 = -1.0 * m_periodicTransformMat(s, 0) * m_omega * m_localMat(i, j);
                            T tmp1 = m_periodicTransformMat(s, 1) * m_omega * m_localMat(i, j);

                            if (tmp0 != 0)
                                globalMat.coeffRef(iiMapped + s*uCompSize, jjMapped + uCompSize) += tmp0;

                            if (tmp1 != 0)
                                globalMat.coeffRef(iiMapped + s*uCompSize, jjMapped) += tmp1;

                        }
                    }
                    // only jj is eliminated periodic dof:
                    else if ( !iiElim && jjElim )
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            T tmp0 = m_periodicTransformMat(s, 0) * m_omega * m_localMat(i, j);
                            T tmp1 = -1.0 * m_periodicTransformMat(s, 1) * m_omega * m_localMat(i, j);

                            if (tmp0 != 0)
                                globalMat.coeffRef(iiMapped + uCompSize, jjMapped + s*uCompSize) += tmp0;

                            if (tmp1 != 0)
                                globalMat.coeffRef(iiMapped, jjMapped + s*uCompSize) += tmp1;

                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int q = 0; q < nComponents; q++)
                        {
                            for (int r = 0; r < nComponents; r++)
                            {
                                T tmp0 = -1.0 * m_periodicTransformMat(q, 0) * m_periodicTransformMat(r, 1) * m_omega * m_localMat(i, j);
                                T tmp1 = m_periodicTransformMat(q, 1) * m_periodicTransformMat(r, 0) * m_omega * m_localMat(i, j);

                                if (tmp0 != 0)
                                    globalMat.coeffRef(iiMapped + q*uCompSize, jjMapped + r*uCompSize) += tmp0;

                                if (tmp1 != 0)
                                    globalMat.coeffRef(iiMapped + q*uCompSize, jjMapped + r*uCompSize) += tmp1;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalRhs(iiMapped + d*uCompSize, 0) -= m_omega * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, d);
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            T tmp0 = -1.0 * m_periodicTransformMat(s, 0) * m_omega * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 1);
                            T tmp1 = m_periodicTransformMat(s, 1) * m_omega * m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 0);

                            if (tmp0 != 0)
                                globalRhs(iiMapped + s*uCompSize) -= tmp0;

                            if (tmp1 != 0)
                                globalRhs(iiMapped + s*uCompSize) -= tmp1;

                        }
                    }
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPU<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_nonper.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();


    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }
} 


template <class T, int MatOrder>
void gsINSVisitorPU<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_per.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            bool iiElim = m_paramsPtr->getPerHelperPtr()->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr()->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jj) += m_locMatVec[d](i, j);
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(s, t) * m_locMatVec[t](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + s*uCompSize, jj) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalRhs(iiMapped + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(s, t) * m_locMatVec[t](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                                
                                if (tmp != 0)
                                    globalRhs(iiMapped + s*uCompSize, 0) -= tmp;
                            }
                        }
                    }
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPU_withUPrhs<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_nonper.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
        else // part arising from block B (assuming that the offdiag. blocks are symmetric)
        {
            const int bb = m_paramsPtr->getMapper(m_testUnkID).global_to_bindex(ii);
            for (index_t k = 0; k < numActTrial; k++)
            {
                const int kk = m_trialFunActives(k);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(kk))
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


template <class T, int MatOrder>
void gsINSVisitorPU_withUPrhs<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU_withUPrhs::localToGlobal_per.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            bool iiElim = m_paramsPtr->getPerHelperPtr()->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr()->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jj) += m_locMatVec[d](i, j);
                    }
                    // ii is eliminated periodic dof:
                    else
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(s, t) * m_locMatVec[t](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + s*uCompSize, jj) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    // ii is not eliminated periodic dof:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalRhs(iiMapped + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(s, t) * m_locMatVec[t](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                                
                                if (tmp != 0)
                                    globalRhs(iiMapped + s*uCompSize, 0) -= tmp;
                            }
                        }
                    }
                }
            }
        }
        else // part arising from block B (assuming that the offdiag. blocks are symmetric)
        {
            const int bb = m_paramsPtr->getMapper(m_testUnkID).global_to_bindex(ii);

            for (index_t k = 0; k < numActTrial; k++)
            {
                const int kk = m_trialFunActives(k);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(kk))
                {
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp += m_locMatVec[d](i, k) * eliminatedDofs[m_testUnkID](bb, d);

                    globalRhs(kk + dim*uCompSize, 0) += tmp;
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorUP<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_trialUnkID).freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal_nonper.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii, jj + d*uCompSize) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);
                    
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, d);

                    globalRhs(ii, 0) += tmp;
                }
            }
        }
    }
} 


template <class T, int MatOrder>
void gsINSVisitorUP<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal_per.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    bool jjElim = m_paramsPtr->getPerHelperPtr()->isEliminated(jj);
                    const index_t jjMapped = m_paramsPtr->getPerHelperPtr()->map(jj);

                    // jj is not eliminated periodic dof:
                    if (!jjElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(ii, jjMapped + d*uCompSize) += m_locMatVec[d](i, j);
                    }
                    // jj is eliminated periodic dof:
                    else
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(s, t) * m_locMatVec[t](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(ii, jjMapped + s*uCompSize) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);
                    
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, d);

                    globalRhs(ii, 0) += tmp;
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPP<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const int ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const int jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
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
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_localMat);
}


template <class T, int MatOrder>
void gsINSVisitorRhsU<T, MatOrder>::localToGlobal_nonper(gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(0).freeSize(); // number of dofs for one velocity component

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

    index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t d = 0; d != dim; d++)
                globalRhs(ii + d*uCompSize, 0) += m_localMat(i, d);
        }
    }
} 


template <class T, int MatOrder>
void gsINSVisitorRhsU<T, MatOrder>::localToGlobal_per(gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

    index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            bool iiElim = m_paramsPtr->getPerHelperPtr()->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr()->map(ii);

            // ii is not eliminated periodic dof
            if (!iiElim) 
            {
                for (index_t d = 0; d < dim; d++)
                    globalRhs(iiMapped + d*uCompSize, 0) += m_localMat(i, d);
            }
            // ii is eliminated periodic dof:
            else
            {
                for (int s = 0; s < dim; s++)
                {
                    for (int t = 0; t < dim; t++)
                    {
                        T tmp = m_periodicTransformMat(s, t) * m_localMat(i, t);

                        if (tmp != 0)
                            globalRhs.coeffRef(iiMapped + s*uCompSize, 0) += tmp;
                    }
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorRhsP<T, MatOrder>::assemble()
{
    m_localMat.setZero(m_testFunActives.rows(), 1);

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_localMat);
}

template <class T, int MatOrder>
void gsINSVisitorRhsP<T, MatOrder>::localToGlobal_nonper(gsMatrix<T>& globalRhs)
{
    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

   index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
            globalRhs(ii, 0) += m_localMat(i);
    }
} 

} // namespace gismo