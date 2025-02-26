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
/*
template <class T, int MatOrder>
void gsTMVisitorVelocityAdvection<T, MatOrder>::initMembers()
{
    
}
*/
template<class T, int MatOrder>
void gsTMVisitorVelocityAdvection<T, MatOrder>::initialize()
{
    Base::initialize();  
}

/*
template <class T, int MatOrder>
void gsTMVisitorVelocityAdvection<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    } 
}
*/

/*
template <class T, int MatOrder>
void gsTMVisitorVelocityAdvection<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    }  
}
*/

template <class T, int MatOrder>
void gsTMVisitorVelocityAdvection<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    index_t m_dofshift = 0;
    for (short_t i = 0; i < m_unknown; i++)
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

// ===================================================================================================================

template<class T, int MatOrder>
void gsTMVisitorDiffusion<T, MatOrder>::initialize()
{
    Base::initialize();  
}

template <class T, int MatOrder>
void gsTMVisitorDiffusion<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsTMTerm_CoeffGradGrad<T>* termPtr = dynamic_cast< gsTMTerm_CoeffGradGrad<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    } 
}

template <class T, int MatOrder>
void gsTMVisitorDiffusion<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsTMTerm_CoeffGradGrad<T>* termPtr = dynamic_cast< gsTMTerm_CoeffGradGrad<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    }  
}

template <class T, int MatOrder>
void gsTMVisitorDiffusion<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    index_t m_dofshift = 0;
    for (short_t i = 0; i < m_unknown; i++)
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

// ==================================================================================================================

template<class T, int MatOrder>
void gsTMVisitorReaction<T, MatOrder>::initialize()
{
    Base::initialize();  
}

template <class T, int MatOrder>
void gsTMVisitorReaction<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_currentSol);
    } 
}

template <class T, int MatOrder>
void gsTMVisitorReaction<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    //m_TMsolverPtr->evalTurbulentViscosity(m_quNodes);
    //m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsTMTerm_CoeffValVal<T>* termPtr = dynamic_cast< gsTMTerm_CoeffValVal<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setCurrentSolution(m_currentSol);
    }  
}

template <class T, int MatOrder>
void gsTMVisitorReaction<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for given turbulent variable
    index_t m_dofshift = 0;
    for (short_t i = 0; i < m_unknown; i++)
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

// ===========================================================================================================================

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

}