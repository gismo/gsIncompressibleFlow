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
void gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    {
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_locMatVec[0](i, j);     // block A
                        globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_locMatVec[d+1](i, j);   // block Eii
                    }
                }
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    for (index_t d = 0; d < dim; d++) 
                    {
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[0](i, j) * eliminatedDofs[m_shapeUnkID](bb, d);   // block A
                        globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d+1](i, j) * eliminatedDofs[m_shapeUnkID](bb, d);   // block Eii
                    }
                }
            }
        }
    }

}

// ===========================================================================================================

template <class T, int MatOrder>
void gsRANSVisitorUUSymmetricGradientOffdiag<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_paramsPtr->getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsRANSVisitorUU::localToGlobal.");

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
                    globalMat.coeffRef(ii, jj + uCompSize) += m_locMatVec[dim+1](i, j);   // block E12
                    globalMat.coeffRef(ii + uCompSize, jj) += m_locMatVec[dim+1](i, j);   // block E21=sym(E12)
                    if (dim == 3)
                    {
                        globalMat.coeffRef(ii, jj + 2*uCompSize) += m_locMatVec[dim+2](i, j);   // block E13
                        globalMat.coeffRef(ii + 2*uCompSize, jj) += m_locMatVec[dim+2](i, j);   // block E31=sym(E13)
                        globalMat.coeffRef(ii + uCompSize, jj + 2*uCompSize) += m_locMatVec[dim+3](i, j);   // block E23
                        globalMat.coeffRef(ii + 2*uCompSize, jj + uCompSize) += m_locMatVec[dim+3](i, j);   // block E32=sym(E23)
                    }
                }
                /* as off-diagonal blocks Ekl go directly to rhs, nothing originating from boundary conditions for these block is added to rhs
                else // is_boundary_index(jj)
                {
                    const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                    globalRhs(ii, 0) -= m_locMatVec[dim+1](i, j) * eliminatedDofs[m_shapeUnkID](bb, 1); // - E12*u_2^*
                    globalRhs(ii, 0) -= m_locMatVec[dim+2](i, j) * eliminatedDofs[m_shapeUnkID](bb, 2); // - E13*u_3^*
                    if (dim == 3)
                    {
                        globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+1](j, i) * eliminatedDofs[m_shapeUnkID](bb, 0);     // - E21*u_1^*
                        globalRhs(ii + uCompSize, 0) -= m_locMatVec[dim+3](i, j) * eliminatedDofs[m_shapeUnkID](bb, 2);     // - E23*u_3^*
                        globalRhs(ii + 2*uCompSize, 0) -= m_locMatVec[dim+2](j, i) * eliminatedDofs[m_shapeUnkID](bb, 0);   // - E31*u_1^*
                        globalRhs(ii + 2*uCompSize, 0) -= m_locMatVec[dim+3](j, i) * eliminatedDofs[m_shapeUnkID](bb, 1);   // - E32*u_2^*
                    }
                }
                */
            }
        }
    }

}

}