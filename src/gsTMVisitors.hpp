/** @file gsTMVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once

#include <gsIncompressibleFlow/src/gsTMVisitors.h>

namespace gismo
{

template <class T, int MatOrder>
void gsTMVisitorLinearSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_trialFunActives.rows();
    
    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);
  
        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);
    
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

// ===================================================================================================================

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setCurrentSolution(USolField);
    } 

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setCurrentSolution(USolField);
    }   

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::assemble()
{
    m_locMatVec.resize(2);

    m_locMatVec[0].setZero(m_testFunActives.rows(), m_trialFunActives.rows());
    m_terms[0]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[0]);

    m_locMatVec[1].setZero(m_testFunActives.rows(), 1);
    m_terms[1]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[1]);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            globalRhs(ii, 0) += m_locMatVec[1](i, 0); 
            
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj) += m_locMatVec[0](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_locMatVec[0](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }

}

// ===================================================================================================================

template<class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::initialize()
{
    Base::initialize();
}


template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);

}


template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);

}

template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::assemble()
{
    GISMO_ASSERT((size_t) (m_numLhsTerms + m_numRhsTerms) == m_terms.size(), "Incorrect number of nonlinear terms for turbulent model!");
    
    m_locMatVec.resize(m_numLhsTerms + m_numRhsTerms);

    for (index_t i = 0; i < m_numLhsTerms; i++)
    {
        m_locMatVec[i].setZero(m_testFunActives.rows(), m_trialFunActives.rows());
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[i]);
    }
        
    for (index_t i = 0; i < m_numRhsTerms; i++)
    {
        m_locMatVec[i + m_numLhsTerms].setZero(m_testFunActives.rows(), 1);
        m_terms[i + m_numLhsTerms]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[i + m_numLhsTerms]);
    }

}

template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t m_dofshift = 0;
    for (short_t i = 2; i < m_unknown; i++)
        m_dofshift += m_dofMappers[i].freeSize();
    
    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_trialUnkID].localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        // nonlinear terms going to lhs
        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            // blended coeff going directly to rhs
            for (index_t k = m_numLhsTerms; k < (m_numLhsTerms + m_numRhsTerms); k++)
                globalRhs(ii, 0) += m_locMatVec[k](i, 0); 
            
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalMat.coeffRef(ii, jj) += m_locMatVec[k](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalRhs(ii, 0) -= m_locMatVec[k](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }
}

}