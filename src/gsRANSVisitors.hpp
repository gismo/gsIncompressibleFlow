/** @file gsRANSVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsRANSVisitors.h>

namespace gismo
{

template <class T, int MatOrder>
void gsRANSVisitorUU<T, MatOrder>::initMembers()
{
    m_viscosity = m_paramsPtr->getPde().viscosity();
}

template <class T, int MatOrder>
void gsRANSVisitorUU<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    } 
}

template <class T, int MatOrder>
void gsRANSVisitorUU<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    }  
}

template <class T, int MatOrder>
void gsRANSVisitorUU<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    // diagonal blocks
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j); // block A + block Eii

                    // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                    for (index_t s = 0; s < dim; s++)
                    {
                        for (index_t t = 0; t < dim; t++)
                        {
                            if (s == t) // diagonal block 
                                continue;

                            T matCoef = (s < t) ? m_locMatVec[dim+s+t](i, j) : m_locMatVec[dim+s+t](j, i); // Ets = Est^T
                            globalRhs(ii + s*uCompSize, 0) -= matCoef * m_solution(jj + t*uCompSize, 0);
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++) 
                        globalRhs(ii + d*uCompSize, 0) -= (m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j)) * eliminatedDofs[m_trialUnkID](bb, d);    // block A + Eii

                    for (index_t s = 0; s < dim; s++)
                    {
                        for (index_t t = 0; t < dim; t++)
                        {
                            if (s == t) // diagonal block 
                                continue;

                            T matCoef = (s < t) ? m_locMatVec[dim+s+t](i, j) : m_locMatVec[dim+s+t](j, i); 
                            globalRhs(ii + s*uCompSize, 0) -= matCoef * eliminatedDofs[m_trialUnkID](bb, t);
                        }
                    }
                }
            }
        }
    }
}


template <class T, int MatOrder>
void gsRANSVisitorUU<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
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
                        // diagonal blocks
                        for (index_t d = 0; d < dim; d++)
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j); // block A + block Eii

                        // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                        for (index_t s = 0; s < dim; s++)
                        {
                            for (index_t t = 0; t < dim; t++)
                            {
                                if (s == t) // diagonal block 
                                    continue;

                                T matCoef = (s < t) ? m_locMatVec[dim+s+t](i, j) : m_locMatVec[dim+s+t](j, i); // Ets = Est^T
                                globalRhs(iiMapped + s*uCompSize, 0) -= matCoef * m_solution(jjMapped + t*uCompSize, 0);
                            }
                        }
                    }
                    // only jj is eliminated periodic dof:
                    else if ( (!iiElim && jjElim) )
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        if (t == r)     // diagonal blocks
                                            globalMat.coeffRef(iiMapped + r*uCompSize, jjMapped + s*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[r+1](i, j) );
                                        else if (r < t) // upper off-diagonal blocks
                                            globalRhs(iiMapped + r*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](i, j) * m_solution(jjMapped + s*uCompSize, 0);
                                        else            // lower off-diagonal blocks
                                            globalRhs(iiMapped + r*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](j, i) * m_solution(jjMapped + s*uCompSize, 0);
                                    }
                                }
                            }
                        }
                    }
                    // only ii is eliminated periodic dof:
                    else if ( (iiElim && !jjElim) )
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        if (t == r)     // diagonal blocks
                                            globalMat.coeffRef(iiMapped + s*uCompSize, jjMapped + r*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[r+1](i, j) );
                                        else if (t < r) // upper off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](i, j) * m_solution(jjMapped + r*uCompSize, 0);
                                        else            // lower off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](j, i) * m_solution(jjMapped + r*uCompSize, 0);
                                    }
                                }
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int q = 0; q < dim; q++)
                        {
                            for (int r = 0; r < dim; r++)
                            {
                                for (int s = 0; s < dim; s++)
                                {
                                    for (int t = 0; t < dim; t++)
                                    {
                                        T rotCoef = m_periodicTransformMat(q, s) * m_periodicTransformMat(r, t);

                                        if (rotCoef != 0)
                                        {
                                            if (s == t)     // diagonal blocks
                                                globalMat.coeffRef(iiMapped + q*uCompSize, jjMapped + r*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[s+1](i, j) );
                                            else if (s < t) // upper off-diagonal blocks
                                                globalRhs(iiMapped + q*uCompSize, 0) -= rotCoef * m_locMatVec[dim+s+t](i, j) * m_solution(jjMapped + r*uCompSize, 0);
                                            else            // lower off-diagonal blocks
                                                globalRhs(iiMapped + q*uCompSize, 0) -= rotCoef * m_locMatVec[dim+s+t](j, i) * m_solution(jjMapped + r*uCompSize, 0);
                                        }
                                    }
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
                            globalRhs(iiMapped + d*uCompSize, 0) -= (m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j)) * eliminatedDofs[m_trialUnkID](bb, d);   // block A + Eii

                        for (index_t s = 0; s < dim; s++)
                        {
                            for (index_t t = 0; t < dim; t++)
                            {
                                if (s == t) // diagonal block 
                                    continue;

                                T matCoef = (s < t) ? m_locMatVec[dim+s+t](i, j) : m_locMatVec[dim+s+t](j, i); 
                                globalRhs(iiMapped + s*uCompSize, 0) -= matCoef * eliminatedDofs[m_trialUnkID](bb, t);
                            }
                        }
                        
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        if (t == r)     // diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[r+1](i, j) ) * eliminatedDofs[m_trialUnkID](bb, r);
                                        else if (t < r) // upper off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](i, j) * eliminatedDofs[m_trialUnkID](bb, r);
                                        else            // lower off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[dim+r+t](j, i) * eliminatedDofs[m_trialUnkID](bb, r);
                                    }
                                }
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
void gsRANSVisitorUU_full<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getMapper(m_testUnkID).freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsRANSVisitorUU_full::localToGlobal.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    index_t blockIndex = 0; // index of current block in m_locMatVec

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
                    // diagonal blocks
                    for (index_t d = 0; d < dim; d++)
                    {
                        blockIndex = d * dim + d + 1; // index of block Edd in m_locMatVec
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j); // block A + block Edd
                    }

                    // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                    for (index_t s = 0; s < dim; s++)
                    {
                        for (index_t t = 0; t < dim; t++)
                        {
                            if (s == t) // diagonal block 
                                continue;

                            blockIndex = s * dim + t + 1; // index of block Est in m_locMatVec
                            globalRhs(ii + s*uCompSize, 0) -= m_locMatVec[blockIndex](i, j) * m_solution(jj + t*uCompSize, 0);
                        }
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++) 
                    {
                        blockIndex = d * dim + d + 1;
                        globalRhs(ii + d*uCompSize, 0) -= (m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j)) * eliminatedDofs[m_trialUnkID](bb, d);    // block A + Edd
                    }

                    for (index_t s = 0; s < dim; s++)
                    {
                        for (index_t t = 0; t < dim; t++)
                        {
                            if (s == t) // diagonal block 
                                continue;

                            blockIndex = s * dim + t + 1; 
                            globalRhs(ii + s*uCompSize, 0) -= m_locMatVec[blockIndex](i, j) * eliminatedDofs[m_trialUnkID](bb, t);
                        }
                    }
                }
            }
        }
    }
}


template <class T, int MatOrder>
void gsRANSVisitorUU_full<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr()->numFreeDofs(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsINSVisitorUU_full::localToGlobal_per, matrix has to contain all components.");

    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActTrial = m_trialFunActives.rows();

    index_t blockIndex = 0; // index of current block in m_locMatVec

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
                        // diagonal blocks
                        for (index_t d = 0; d < dim; d++)
                        {
                            blockIndex = d * dim + d + 1;
                            globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j); // block A + block Edd
                        }

                        // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                        for (index_t s = 0; s < dim; s++)
                        {
                            for (index_t t = 0; t < dim; t++)
                            {
                                if (s == t) // diagonal block 
                                    continue;

                                blockIndex = s * dim + t + 1; 
                                globalRhs(iiMapped + s*uCompSize, 0) -= m_locMatVec[blockIndex](i, j) * m_solution(jjMapped + t*uCompSize, 0);
                            }
                        }
                    }
                    // only jj is eliminated periodic dof:
                    else if ( (!iiElim && jjElim) )
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        blockIndex = r * dim + t + 1; // index of block Ert in m_locMatVec

                                        if (t == r)     // diagonal blocks
                                            globalMat.coeffRef(iiMapped + r*uCompSize, jjMapped + s*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j) );
                                        else // off-diagonal blocks
                                            globalRhs(iiMapped + r*uCompSize, 0) -= rotCoef * m_locMatVec[blockIndex](i, j) * m_solution(jjMapped + s*uCompSize, 0);
                                    }
                                }
                            }
                        }
                    }
                    // only ii is eliminated periodic dof:
                    else if ( (iiElim && !jjElim) )
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        blockIndex = r * dim + t + 1; // index of block Ert in m_locMatVec

                                        if (t == r)     // diagonal blocks
                                            globalMat.coeffRef(iiMapped + s*uCompSize, jjMapped + r*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j) );
                                        else // off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[blockIndex](i, j) * m_solution(jjMapped + r*uCompSize, 0);
                                    }
                                }
                            }
                        }
                    }
                    // both ii and jj are eliminated periodic dofs:
                    else 
                    {
                        for (int q = 0; q < dim; q++)
                        {
                            for (int r = 0; r < dim; r++)
                            {
                                for (int s = 0; s < dim; s++)
                                {
                                    for (int t = 0; t < dim; t++)
                                    {
                                        T rotCoef = m_periodicTransformMat(q, s) * m_periodicTransformMat(r, t);

                                        if (rotCoef != 0)
                                        {
                                            blockIndex = s * dim + t + 1; // index of block Est in m_locMatVec

                                            if (s == t)     // diagonal blocks
                                                globalMat.coeffRef(iiMapped + q*uCompSize, jjMapped + r*uCompSize) += rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j) );
                                            else // off-diagonal blocks
                                                globalRhs(iiMapped + q*uCompSize, 0) -= rotCoef * m_locMatVec[blockIndex](i, j) * m_solution(jjMapped + r*uCompSize, 0);
                                        }
                                    }
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
                            globalRhs(iiMapped + d*uCompSize, 0) -= (m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j)) * eliminatedDofs[m_trialUnkID](bb, d);   // block A + Eii

                        for (index_t s = 0; s < dim; s++)
                        {
                            for (index_t t = 0; t < dim; t++)
                            {
                                if (s == t) // diagonal block 
                                    continue;

                                blockIndex = s * dim + t + 1;
                                globalRhs(iiMapped + s*uCompSize, 0) -= m_locMatVec[blockIndex](i, j) * eliminatedDofs[m_trialUnkID](bb, t);
                            }
                        }
                        
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        for (int r = 0; r < dim; r++)
                        {
                            for (int s = 0; s < dim; s++)
                            {
                                for (int t = 0; t < dim; t++)
                                {
                                    T rotCoef = m_periodicTransformMat(s, t);

                                    if (rotCoef != 0)
                                    {
                                        blockIndex = r * dim + t + 1; // index of block Ert in m_locMatVec

                                        if (t == r)     // diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * ( m_locMatVec[0](i, j) + m_locMatVec[blockIndex](i, j) ) * eliminatedDofs[m_trialUnkID](bb, r);
                                        else // off-diagonal blocks
                                            globalRhs(iiMapped + s*uCompSize, 0) -= rotCoef * m_locMatVec[blockIndex](i, j) * eliminatedDofs[m_trialUnkID](bb, r);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


// ================================================================================================================
// For T-CSD stabilization
//

template <class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_time<T, MatOrder>::initMembers()
{
    m_viscosity = m_paramsPtr->getPde().viscosity();
}

template<class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_time<T, MatOrder>::initialize()
{
    Base::initialize();  
}

template <class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_time<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2));

    gsFlowTerm_TCSDStabilization_time<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_time<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(velocities);
    }
}

template <class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_time<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2));

    gsFlowTerm_TCSDStabilization_time<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_time<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(velocities);
    }
}

// ================================================================================================================

template <class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_advection<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2));

    gsFlowTerm_TCSDStabilization_advection<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_advection<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(velocities);
    }
}

template <class T, int MatOrder>
void gsRANSVisitorTCSDStabilization_advection<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->evalTurbulentViscosity(m_quNodes, m_patchID);
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2));

    gsFlowTerm_TCSDStabilization_advection<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_advection<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(velocities);
    }
}

}