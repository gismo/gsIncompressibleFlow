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
    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_trialFunActives.rows();
    
    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);
  
        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);
    
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
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::assemble()
{
    m_locMatVec.resize(1);

    //m_locMatVec[0].setZero(m_testFunActives.rows(), m_trialFunActives.rows());
    //m_terms[0]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[0]);

    m_locMatVec[0].setZero(m_testFunActives.rows(), 1);
    m_terms[0]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec[0]);
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    // index_t numActShape = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            globalRhs(ii, 0) += m_locMatVec[0](i, 0); 
            
            /*for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    globalMat.coeffRef(ii, jj) += m_locMatVec[0](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_locMatVec[0](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }*/
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

    gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms[i]);
        if (termPtr)
        {
            termPtr->setCurrentSolution(USolField);
        } 
    }

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);

}


template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms[i]);
        if (termPtr)
        {
            termPtr->setCurrentSolution(USolField);
        } 
    }

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);

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
        m_dofshift += m_paramsPtr->getMapper(i).freeSize();
    
    m_paramsPtr->getMapper(m_testUnkID).localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_paramsPtr->getMapper(m_trialUnkID).localToGlobal(m_trialFunActives, m_patchID, m_trialFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_trialFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        // nonlinear terms going to lhs
        if (m_paramsPtr->getMapper(m_testUnkID).is_free_index(ii))
        {
            // blended coeff going directly to rhs
            for (index_t k = m_numLhsTerms; k < (m_numLhsTerms + m_numRhsTerms); k++)
                globalRhs(ii, 0) += m_locMatVec[k](i, 0); 
            
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_paramsPtr->getMapper(m_trialUnkID).is_free_index(jj))
                {
                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalMat.coeffRef(ii, jj) += m_locMatVec[k](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_paramsPtr->getMapper(m_trialUnkID).global_to_bindex(jj);

                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalRhs(ii, 0) -= m_locMatVec[k](i, j) * eliminatedDofs[m_trialUnkID](bb, 0);
                }
            }
        }
    }
}


// ================================================================================================================
// For T-CSD stabilization
//

template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_time<T, MatOrder>::initMembers()
{
    m_viscosity = m_paramsPtr->getPde().viscosity();
}

template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_time<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->updateModel(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    gsMatrix<T> USolVals = m_TMModelPtr->getUSolVals();
    gsMatrix<T> OSolVals = m_TMModelPtr->getOSolVals();
    gsVector<T> F1 = m_TMModelPtr->getF1Vals();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t sigmaK1 = m_TMModelPtr->get_sigmaK1();
    real_t sigmaK2 = m_TMModelPtr->get_sigmaK2();
    real_t sigmaO1 = m_TMModelPtr->get_sigmaO1();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    real_t eps = m_TMModelPtr->get_eps();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        if (m_unknown == 2)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaK1 + (1 - F1(i)) * sigmaK2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max(betaStar * OSolVals(0, i), eps));
        else if (m_unknown == 3)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaO1 + (1 - F1(i)) * sigmaO2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max((F1(i) * beta1 + (1 - F1(i)) * beta2) * OSolVals(0, i), eps));
        else
            gsWarn << "Wrong unknown identification!\n";

    gsFlowTerm_TCSDStabilization_time<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_time<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(USolVals);
    }

}


template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_time<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);

    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    gsMatrix<T> USolVals = m_TMModelPtr->getUSolVals();
    gsMatrix<T> OSolVals = m_TMModelPtr->getOSolVals();
    gsVector<T> F1 = m_TMModelPtr->getF1Vals();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t sigmaK1 = m_TMModelPtr->get_sigmaK1();
    real_t sigmaK2 = m_TMModelPtr->get_sigmaK2();
    real_t sigmaO1 = m_TMModelPtr->get_sigmaO1();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    real_t eps = m_TMModelPtr->get_eps();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        if (m_unknown == 2)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaK1 + (1 - F1(i)) * sigmaK2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max(betaStar * OSolVals(0, i), eps));
        else if (m_unknown == 3)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaO1 + (1 - F1(i)) * sigmaO2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max((F1(i) * beta1 + (1 - F1(i)) * beta2) * OSolVals(0, i), eps));
        else
            gsWarn << "Wrong unknown identification!\n";

    gsFlowTerm_TCSDStabilization_time<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_time<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(USolVals);
    }

}

template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_advection<T, MatOrder>::initMembers()
{
    m_viscosity = m_paramsPtr->getPde().viscosity();
}

template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_advection<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->updateModel(m_quNodes, this->m_quRule.numNodes(), m_patchID);
    
    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    gsMatrix<T> USolVals = m_TMModelPtr->getUSolVals();
    gsMatrix<T> OSolVals = m_TMModelPtr->getOSolVals();
    gsVector<T> F1 = m_TMModelPtr->getF1Vals();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t sigmaK1 = m_TMModelPtr->get_sigmaK1();
    real_t sigmaK2 = m_TMModelPtr->get_sigmaK2();
    real_t sigmaO1 = m_TMModelPtr->get_sigmaO1();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    real_t eps = m_TMModelPtr->get_eps();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        if (m_unknown == 2)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaK1 + (1 - F1(i)) * sigmaK2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max(betaStar * OSolVals(0, i), eps));
        else if (m_unknown == 3)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaO1 + (1 - F1(i)) * sigmaO2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max((F1(i) * beta1 + (1 - F1(i)) * beta2) * OSolVals(0, i), eps));
        else
            gsWarn << "Wrong unknown identification!\n";

    gsFlowTerm_TCSDStabilization_advection<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_advection<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(USolVals);
    }

}


template <class T, int MatOrder>
void gsTMVisitorSSTTCSDStabilization_advection<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->updateModel(m_mapData.points, this->m_quRule.numNodes(), m_mapData.patchId);

    m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    gsMatrix<T> USolVals = m_TMModelPtr->getUSolVals();
    gsMatrix<T> OSolVals = m_TMModelPtr->getOSolVals();
    gsVector<T> F1 = m_TMModelPtr->getF1Vals();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t sigmaK1 = m_TMModelPtr->get_sigmaK1();
    real_t sigmaK2 = m_TMModelPtr->get_sigmaK2();
    real_t sigmaO1 = m_TMModelPtr->get_sigmaO1();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    real_t eps = m_TMModelPtr->get_eps();

    const index_t nQuPoints = m_quNodes.cols();
    gsMatrix<T> physPoints = m_mapData.values[0];
    real_t h = (physPoints.col(physPoints.cols()-1) - physPoints.col(0)).norm();
    
    gsMatrix<T> velocities = m_solution.value(m_quNodes);
    m_tauS.resize(1, nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        if (m_unknown == 2)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaK1 + (1 - F1(i)) * sigmaK2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max(betaStar * OSolVals(0, i), eps));
        else if (m_unknown == 3)
            m_tauS(0, i) = 1 / math::sqrt(math::pow(2 * velocities.col(i).norm() / h, 2) + 9 * math::pow(4 * (m_viscosity + (F1(i) * sigmaO1 + (1 - F1(i)) * sigmaO2) * m_TurbulentViscosityVals(i)) / math::pow(h, 2), 2) + math::max((F1(i) * beta1 + (1 - F1(i)) * beta2) * OSolVals(0, i), eps));
        else
            gsWarn << "Wrong unknown identification!\n";

    gsFlowTerm_TCSDStabilization_advection<T>* termPtr = dynamic_cast< gsFlowTerm_TCSDStabilization_advection<T>* > (m_terms.back());
    if (termPtr)
    {
        termPtr->setTauS(m_tauS);
        termPtr->setUSolVals(USolVals);
    }

}

}