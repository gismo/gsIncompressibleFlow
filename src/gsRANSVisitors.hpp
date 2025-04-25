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
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::initMembers()
{
    m_viscosity = m_paramsPtr->getPde().viscosity();
}

template<class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::initialize()
{
    Base::initialize();  
}

template <class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::evaluate(index_t testFunID)
{
    Base::evaluate(testFunID);

    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes, m_patchID);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    } 
}

template <class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    Base::evaluate(domIt);

    m_TMsolverPtr->evalTurbulentViscosity(m_quNodes, m_patchID);
    m_TurbulentViscosityVals = m_TMsolverPtr->getTurbulentViscosity();

    gsRANSTerm_SymmetricGradient<T>* termPtr = dynamic_cast< gsRANSTerm_SymmetricGradient<T>* > (m_terms.back());

    if (termPtr)
    {
        termPtr->setViscosity(m_viscosity);
        termPtr->setTurbulentViscosityVals(m_TurbulentViscosityVals);
    }  
}

template <class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    // diagonal blocks
                    for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j); // block A + block Eii

                    // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                    globalRhs(ii, 0) -= m_locMatVec[dim+1](i, j) * m_solution(jj + uCompSize, 0); // - E12*u_2^*
                    globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+1](j, i) * m_solution(jj, 0);     // - E21*u_1^*
                    
                    if (dim == 3)
                    {
                        globalRhs(ii, 0) -= m_locMatVec[dim+2](i, j) * m_solution(jj + 2*uCompSize, 0); // - E13*u_3^*
                        globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+3](i, j) * m_solution(jj + 2*uCompSize, 0);     // - E23*u_3^*
                        globalRhs(ii + 2*uCompSize, 0) -= m_locMatVec[dim+2](j, i) * m_solution(jj, 0)                  // - E31*u_1^*
                                                        + m_locMatVec[dim+3](j, i) * m_solution(jj + uCompSize, 0);   // - E32*u_2^* 
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_trialUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++) 
                    {
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[0](i, j) * eliminatedDofs[m_trialUnkID](bb, d);   // block A
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d+1](i, j) * eliminatedDofs[m_trialUnkID](bb, d);   // block Eii
                    }

                    globalRhs(ii, 0) -= m_locMatVec[dim+1](i, j) * eliminatedDofs[m_trialUnkID](bb, 1);                 // - E12*u_2^*
                    globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+1](j, i) * eliminatedDofs[m_trialUnkID](bb, 0);     // - E21*u_1^*
                    
                    if (dim == 3)
                    {
                        globalRhs(ii, 0) -= m_locMatVec[dim+2](i, j) * eliminatedDofs[m_trialUnkID](bb, 2);                 // - E13*u_3^*
                        globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+3](i, j) * eliminatedDofs[m_trialUnkID](bb, 2);     // - E23*u_3^*
                        globalRhs(ii + 2*uCompSize, 0) -= m_locMatVec[dim+2](j, i) * eliminatedDofs[m_trialUnkID](bb, 0);   // - E31*u_1^*
                        globalRhs(ii + 2*uCompSize, 0) -= m_locMatVec[dim+3](j, i) * eliminatedDofs[m_trialUnkID](bb, 1);   // - E32*u_2^*
                    }
                }
            }
        }
    }
}


template <class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradient<T, MatOrder>::localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_paramsPtr->getPerHelperPtr(m_testUnkID)->numFreeDofs(); // number of dofs for one velocity component
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
            bool iiElim = m_paramsPtr->getPerHelperPtr(m_testUnkID)->isEliminated(ii);
            const index_t iiMapped = m_paramsPtr->getPerHelperPtr(m_testUnkID)->map(ii);

            for (index_t j = 0; j < numActTrial; ++j)
            {
                const index_t jj = m_trialFunActives(j);

                if (m_dofMappers[m_trialUnkID].is_free_index(jj))
                {
                    bool jjElim = m_paramsPtr->getPerHelperPtr(m_trialUnkID)->isEliminated(jj);
                    const index_t jjMapped = m_paramsPtr->getPerHelperPtr(m_trialUnkID)->map(jj);

                    // ii and jj are not eliminated periodic dofs:
                    if (!iiElim && !jjElim) 
                    {
                        // diagonal blocks
                        for (index_t d = 0; d < dim; d++)
                        globalMat.coeffRef(iiMapped + d*uCompSize, jjMapped + d*uCompSize) += m_locMatVec[0](i, j) + m_locMatVec[d+1](i, j); // block A + block Eii

                        // off-diagonal blocks asrising from symmetric gradient are put directly to right-hand side of the system
                        globalRhs(iiMapped, 0) -= m_locMatVec[dim+1](i, j) * m_solution(jjMapped + uCompSize, 0); // - E12*u_2^*
                        globalRhs(iiMapped + uCompSize, 0) -= m_locMatVec[dim+1](j, i) * m_solution(jjMapped, 0);     // - E21*u_1^*
                        
                        if (dim == 3)
                        {
                            globalRhs(iiMapped, 0) -= m_locMatVec[dim+2](i, j) * m_solution(jjMapped + 2*uCompSize, 0); // - E13*u_3^*
                            globalRhs(iiMapped + uCompSize, 0) -= m_locMatVec[dim+3](i, j) * m_solution(jjMapped + 2*uCompSize, 0);     // - E23*u_3^*
                            globalRhs(iiMapped + 2*uCompSize, 0) -= m_locMatVec[dim+2](j, i) * m_solution(jjMapped, 0)                  // - E31*u_1^*
                                                            + m_locMatVec[dim+3](j, i) * m_solution(jjMapped + uCompSize, 0);   // - E32*u_2^* 
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

                                    }
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
                        
                    }
                    // ii is eliminated periodic dof:
                    else 
                    {
                        
                    }
                }
            }
        }
    }
}

}