/** @file gsRANSVisitors.hpp

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

/*template<class T, int MatOrder>
void gsTMVisitorLinearSST<T, MatOrder>::initialize()
{
    defineTestShapeUnknowns();  
    
    Base::deleteTerms();
    defineTerms();
    Base::gatherEvalFlags();
    m_mapData.flags = m_geoFlags;
}
*/

template <class T, int MatOrder>
void gsTMVisitorLinearSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    //index_t dim = m_paramsPtr->getPde().dim();
    //const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    //index_t nComponents = globalMat.rows() / uCompSize;
    //index_t m_dofshift = 0;
    //for (short_t i = 2; i < m_unknown; i++)
    //    m_dofshift += m_dofMappers[i].freeSize();

    //GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal.");

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
                    globalMat.coeffRef(ii/* + m_dofshift*/, jj/*+ m_dofshift*/) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);
   
                    globalRhs(ii/* + m_dofshift*/, 0) -= m_localMat(i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                }
            }
        }
    }
    } 

// ===================================================================================================================

/*template<class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::initialize()
{
    defineTestShapeUnknowns();  
    
    Base::deleteTerms();
    defineTerms();
    Base::gatherEvalFlags();
    m_mapData.flags = m_geoFlags;
}
    */

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    m_USolField = m_paramsPtr->getVelocitySolution();

    gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_USolField);
    } 
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_USolField = m_paramsPtr->getVelocitySolution();

    gsTMTerm_VecCoeffGradVal<T>* termPtr = dynamic_cast< gsTMTerm_VecCoeffGradVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_USolField);
    }   
}

template <class T, int MatOrder>
void gsTMVisitorTimeIterationSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    //index_t dim = m_paramsPtr->getPde().dim();
    //const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    //index_t m_dofshift = 0;
    //for (short_t i = 2; i < m_unknown; i++)
    //    m_dofshift += m_dofMappers[i].freeSize();
    //index_t nComponents = globalMat.rows() / uCompSize;

    //GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    globalMat.coeffRef(ii/* + m_dofshift*/, jj/*+ m_dofshift*/) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    globalRhs(ii/* + m_dofshift*/, 0) -= m_localMat(i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                }
            }
        }
    }

}

// ===================================================================================================================

template<class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::initialize()
{
    
    //getAssembler()->setSSTModelEvaluator(m_SSTPtr);

    Base::initialize();
}


template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);

    // gsTMTerm_CoeffGradGrad<T>* termPtr1 = dynamic_cast< gsTMTerm_CoeffGradGrad<T>* > (m_terms[1]);
    // termPtr1->setDistanceField();
    // gsTMTerm_CoeffValVal<T>* termPtr2 = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms[2]);
    // termPtr2->setDistanceField();
    // gsTMTerm_BlendCoeffRhs<T>* termPtr3 = dynamic_cast< gsTMTerm_BlendCoeffRhs<T>* > (m_terms[3]);
    // termPtr3->setDistanceField();
    // gsTMTerm_ProductionRhs<T>* termPtr4 = dynamic_cast< gsTMTerm_ProductionRhs<T>* > (m_terms[4]);
    // termPtr4->setDistanceField();

    
    //gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    //if (termPtr)
    //{
    //    termPtr->setCurrentSolution(m_currentSol);
    //} 
}


template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMModelPtr->updateModel(m_mapData.points, m_mapData.patchId);

    // gsTMTerm_CoeffGradGrad<T>* termPtr1 = dynamic_cast< gsTMTerm_CoeffGradGrad<T>* > (m_terms[1]);
    // termPtr1->setDistanceField();
    // gsTMTerm_CoeffValVal<T>* termPtr2 = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms[2]);
    // termPtr2->setDistanceField();
    // gsTMTerm_BlendCoeffRhs<T>* termPtr3 = dynamic_cast< gsTMTerm_BlendCoeffRhs<T>* > (m_terms[3]);
    // termPtr3->setDistanceField();
    // gsTMTerm_ProductionRhs<T>* termPtr4 = dynamic_cast< gsTMTerm_ProductionRhs<T>* > (m_terms[4]);
    // termPtr4->setDistanceField();

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    //gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    //if (termPtr)
    //{
    //    termPtr->setCurrentSolution(m_currentSol);
    //}  
}

template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::assemble()
{
    GISMO_ASSERT((size_t) (m_numLhsTerms + m_numRhsTerms) == m_terms.size(), "Incorrect number of nonlinear terms for turbulent model!");
    
    m_locMatVec.resize(m_numLhsTerms + m_numRhsTerms);

    for (index_t i = 0; i < m_numLhsTerms; i++)
    {
        m_locMatVec[i].setZero(m_testFunActives.rows(), m_shapeFunActives.rows());
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_locMatVec[i]);
    }
        
    for (index_t i = 0; i < m_numRhsTerms; i++)
    {
        m_locMatVec[i + m_numLhsTerms].setZero(m_testFunActives.rows(), 1);
        m_terms[i + m_numLhsTerms]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_locMatVec[i + m_numLhsTerms]);
    }

    //gsInfo << m_locMatVec[2].sum() << std::endl;
        
    //m_SSTPtr->setNotCurrent();
}

template <class T, int MatOrder>
void gsTMVisitorNonlinearSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    //index_t dim = m_paramsPtr->getPde().dim();
    //const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    index_t m_dofshift = 0;
    for (short_t i = 2; i < m_unknown; i++)
        m_dofshift += m_dofMappers[i].freeSize();
    //index_t nComponents = globalMat.rows() / uCompSize;

    //GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

    m_dofMappers[m_testUnkID].localToGlobal(m_testFunActives, m_patchID, m_testFunActives);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);

    index_t numActTest = m_testFunActives.rows();
    index_t numActShape = m_shapeFunActives.rows();

    for (index_t i = 0; i < numActTest; ++i)
    {
        const index_t ii = m_testFunActives(i);

        // production term and blended coeff going directly to rhs
        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            //globalRhs(ii, 0) += m_locMatVec[2](i);
            for (index_t k = m_numLhsTerms; k < (m_numLhsTerms + m_numRhsTerms); k++)
                globalRhs(ii, 0) += m_locMatVec[k](i);           
        }
        
        // nonlinear terms going to lhs
        if (m_dofMappers[m_testUnkID].is_free_index(ii))
        {
            for (index_t j = 0; j < numActShape; ++j)
            {
                const index_t jj = m_shapeFunActives(j);

                if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
                {
                    // globalMat.coeffRef(ii, jj) += m_locMatVec[0](i, j);
                    // globalMat.coeffRef(ii, jj) += m_locMatVec[1](i, j);
                    // globalMat.coeffRef(ii, jj) += m_locMatVec[2](i, j);
                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalMat.coeffRef(ii, jj) += m_locMatVec[k](i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    // globalRhs(ii, 0) -= m_locMatVec[0](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                    // globalRhs(ii, 0) -= m_locMatVec[1](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                    // globalRhs(ii, 0) -= m_locMatVec[2](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                    for (index_t k = 0; k < m_numLhsTerms; k++)
                        globalRhs(ii, 0) -= m_locMatVec[k](i, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
                    // production term and blended coeff going directly to rhs
                    //globalRhs(ii, 0) -= m_locMatVec[1](j, i) * m_solution(jj + m_dofshift, 0); 
                }
            }
        }
    }
}

// ==================================================================================================================



// ===========================================================================================================================

/*
template<class T, int MatOrder>
void gsTMVisitorProductuionRhsSST<T, MatOrder>::initialize()
{
    Base::initialize();  
}

template <class T, int MatOrder>
void gsTMVisitorProductuionRhsSST<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    

    // evaluatuion of a turbulent viscosity
    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    // evaluation of velocity gradients

    // evaluation of a strainrate magnitude

    // evaluation of k and omega solutions

    // evaluation of a production

    gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_currentSol);
    } 
}

template <class T, int MatOrder>
void gsTMVisitorProductuionRhsSST<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsBasisRefs<T> basisRefsRANS = m_paramsPtr->getBases();
    gsBasisRefs<T> basisRefsTM = m_paramsPtr->getBasesTM();

    gsMatrix<T> activesU, activesK, activesO;
    basisRefsRANS.front().active_into(m_quNodes, activesU);
    basisRefsTM.front().active_into(m_quNodes, activesK);
    basisRefsTM.back().active_into(m_quNodes, activesO);

    gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_currentSol);
    }  
}

template <class T, int MatOrder>
void gsTMVisitorProductuionRhsSST<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    //const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    //index_t m_dofshift = 0;
    //for (short_t i = 0; i < m_unknown; i++)
    //    m_dofshift += m_dofMappers[i].freeSize();
    //index_t nComponents = globalMat.rows() / uCompSize;

    

    //GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    globalMat.coeffRef(ii + m_dofshift, jj + m_dofshift) += m_localMat(i, j);
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    globalRhs(ii + m_dofshift, 0) -= m_localMat(i, j) * eliminatedDofs[m_shapeUnkID](bb, d);
                }
            }
        }
    }
}
*/

}