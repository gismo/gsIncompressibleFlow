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
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal_nonper.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);


    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < nComponents; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

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
    const index_t uCompSize = m_testPeriodicHelperPtr->numFreeDofs(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal_per, matrix has to contain all components.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            bool iiElim = m_testPeriodicHelperPtr->isEliminated(ii);
            const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    bool jjElim = m_trialPeriodicHelperPtr->isEliminated(jj);
                    const index_t jjMapped = m_trialPeriodicHelperPtr->map(jj);

                    // ii and jj are not eliminated periodic dofs:
                    if (!iiElim && !jjElim) 
                    {
                        for (index_t d = 0; d < nComponents; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped + d*uCompSize) += m_localMat(i, j);
                    }
                    // only ii or jj is eliminated periodic dof:
                    else if ( (!iiElim && jjElim) || (iiElim && !jjElim) )
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            for (int t = 0; t < nComponents; t++)
                            {
                                T tmp = 0;

                                if (jjElim)
                                    tmp = m_periodicTransformMat(t, s);
                                else if (iiElim)
                                    tmp = m_periodicTransformMat(s, t);

                                tmp *= m_localMat(i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + s*uCompSize, jjMapped + t*uCompSize) += tmp;
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int s = 0; s < nComponents; s++)
                        {
                            for (int t = 0; t < nComponents; t++)
                            {
                                for (int u = 0; u < nComponents; u++)
                                {
                                    T tmp = m_periodicTransformMat(u, s) * m_periodicTransformMat(t, s) * m_localMat(i, j);
                                    
                                    if (tmp != 0)
                                        globalMat.coeffRef(iiMapped + u*uCompSize, jjMapped + t*uCompSize) += tmp;
                                }
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

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
void gsINSVisitorPU<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_nonper.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();


    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

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
    const index_t uCompSize = m_testPeriodicHelperPtr->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_per.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            bool iiElim = m_testPeriodicHelperPtr->isEliminated(ii);
            const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    bool jjElim = m_trialPeriodicHelperPtr->isEliminated(jj);
                    const index_t jjMapped = m_trialPeriodicHelperPtr->map(jj);

                    // ii and jj are not eliminated periodic dofs or only jj is eliminated:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped) += m_locMatVec[d](i, j);
                    }
                    // only ii is eliminated periodic dof:
                    else if ( (iiElim && !jjElim) )
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + t*uCompSize, jjMapped) += tmp;
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);
                                
                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + t*uCompSize, jjMapped) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

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
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                                
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
void gsINSVisitorPU_withUPrhs<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal_nonper.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++)
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
        else // part arising from block B (assuming that the offdiag. blocks are symmetric)
        {
            const int bb = m_dofMappers[m_testUnkID].global_to_bindex(ii);
            for (index_t k = 0; k < numActTrial; k++)
            {
                const int kk = m_trialFunActives(k);

                if (m_dofMappers[m_trialUnkID].is_free_index(kk))
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
    const index_t uCompSize = m_testPeriodicHelperPtr->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU_withUPrhs::localToGlobal_per.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            bool iiElim = m_testPeriodicHelperPtr->isEliminated(ii);
            const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    bool jjElim = m_trialPeriodicHelperPtr->isEliminated(jj);
                    const index_t jjMapped = m_trialPeriodicHelperPtr->map(jj);

                    // ii and jj are not eliminated periodic dofs or only jj is eliminated:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped) += m_locMatVec[d](i, j);
                    }
                    // only ii is eliminated periodic dof:
                    else if ( iiElim && !jjElim )
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + t*uCompSize, jjMapped) += tmp;
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);
                                
                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped + t*uCompSize, jjMapped) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

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
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j)* eliminatedDofs[m_trialUnkID](bb, 0);
  
                                if (tmp != 0)
                                    globalRhs(iiMapped + t*uCompSize, 0) -= tmp;
                            }
                        }
                    }
                }
            }
        }
        else // part arising from block B (assuming that the offdiag. blocks are symmetric)
        {
            const int bb = m_dofMappers[m_testUnkID].global_to_bindex(ii);

            for (index_t k = 0; k < numActTrial; k++)
            {
                const int kk = m_trialFunActives(k);

                if (m_dofMappers[m_trialUnkID].is_free_index(kk))
                {
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp += m_locMatVec[d](i, k) * eliminatedDofs[m_testUnkID](bb, d);

                    globalRhs(m_trialPeriodicHelperPtr->map(kk) + dim*uCompSize, 0) += tmp;
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
    const index_t uCompSize = m_dofMappers[m_trialUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal_nonper.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii, jj + d*uCompSize) += m_locMatVec[d](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);
                    
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
    const index_t uCompSize = m_trialPeriodicHelperPtr->numFreeDofs(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal_per.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            bool iiElim = m_testPeriodicHelperPtr->isEliminated(ii);
            const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    bool jjElim = m_trialPeriodicHelperPtr->isEliminated(jj);
                    const index_t jjMapped = m_trialPeriodicHelperPtr->map(jj);

                    // ii and jj are not eliminated periodic dofs or only ii is eliminated:
                    if (!iiElim) 
                    {
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped, jjMapped + d*uCompSize) += m_locMatVec[d](i, j);
                    }
                    // only jj is eliminated periodic dof:
                    else if ( !iiElim && jjElim )
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);

                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped, jjMapped + t*uCompSize) += tmp;
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int s = 0; s < dim; s++)
                        {
                            for (int t = 0; t < dim; t++)
                            {
                                T tmp = m_periodicTransformMat(t, s) * m_locMatVec[s](i, j);
                                
                                if (tmp != 0)
                                    globalMat.coeffRef(iiMapped, jjMapped + t*uCompSize) += tmp;
                            }
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);
                    
                    T tmp = 0;

                    for (index_t d = 0; d < dim; d++)
                        tmp -= m_locMatVec[d](i, j) * eliminatedDofs[m_trialUnkID](bb, d);

                    globalRhs(iiMapped, 0) += tmp;
                }
            }
        }
    }
} 

// ===================================================================================================================

template <class T, int MatOrder>
void gsINSVisitorPP<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const int ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const int jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }
} 

template <class T, int MatOrder>
void gsINSVisitorPP<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);
    
    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const int ii = m_testFunActives(i);
        const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActTrial; ++j)
            {
                const int jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    globalMat.coeffRef(iiMapped, m_trialPeriodicHelperPtr->map(jj)) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    globalRhs(iiMapped, 0) -= m_localMat(i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
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


template <class T, int MatOrder>
void gsINSVisitorRhsU<T, MatOrder>::localToGlobal_per(gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_testPeriodicHelperPtr->numFreeDofs(); // number of dofs for one velocity component

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

    index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            bool iiElim = m_testPeriodicHelperPtr->isEliminated(ii);
            const index_t iiMapped = m_testPeriodicHelperPtr->map(ii);

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
                        T tmp = m_periodicTransformMat(t, s) * m_localMat(i, t);

                        if (tmp != 0)
                            globalRhs.coeffRef(iiMapped + t*uCompSize, 0) += tmp; // t or s???
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
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

   index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
            globalRhs(ii, 0) += m_localMat(i);
    }
} 


template <class T, int MatOrder>
void gsINSVisitorRhsP<T, MatOrder>::localToGlobal_per(gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);

   index_t numActTest = m_testFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
            globalRhs(m_testPeriodicHelperPtr->map(ii), 0) += m_localMat(i);
    }
} 

} // namespace gismo