/** @file gsINSAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova), B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMAssemblerBase.h>

namespace gismo
{

template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::initMembers()
{ 
    Base::initMembers();

    m_viscosity = m_paramsPtr->getPde().viscosity();

    m_bases = m_paramsPtr->getBasesTM();

    m_nnzPerRowTM = 1;
    index_t maxDeg = 0;
    for (short_t i = 0; i < m_tarDim; i++)
        for (short_t j = 0; j < numTMvars; j++)
            if (m_bases[j].maxDegree(i) > maxDeg)
                maxDeg = m_bases[j].maxDegree(i);
    m_nnzPerRowTM = m_tarDim * (2 * maxDeg + 1);

    m_bInitialized = false;

    updateSizes();
}

template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::updateSizes()
{
    m_dofs = 0;
    for (short_t i = 0; i < numTMvars; i++)
    {
        m_kdofs[i] = m_dofMappers[i].freeSize();
        m_dofs += m_kdofs[i];

        m_ddof[i].setZero(m_dofMappers[i].boundarySize(), 1);

        if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            this->computeDirichletDofs(i, i, m_ddof[i]);
        }
    }
    
    m_solution.setZero(m_dofs, 1);

    //m_currentSolField = constructSolution(m_solution, 0);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);
}


template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
        m_solution = solVector;

    m_currentSolField = constructSolution(solVector, 0);
    //m_visitorUUnonlin.setCurrentSolution(m_currentVelField);
}


template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat, const index_t unk)
{
    if (sourceMat.rows() == m_dofs) // sourceMat is of the size of the whole globalMat matrix
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row(), it.col()) += it.value();
    }
    else // sourceMat is a block for one turbulent unknown
    {
        index_t dofs = 0;
        for (index_t i = 0; i < unk; i++)
            dofs += m_kdofs[i];
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                    globalMat.coeffRef(it.row() + dofs, it.col() + dofs) += it.value();
    }
}

template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
{
    m_dofMappers[unk] = gsDofMapper(getBases().at(unk), getBCs(), unk);

    if (getAssemblerOptions().intStrategy == iFace::conforming)
        for (gsBoxTopology::const_iiterator it = getPatches().iBegin(); it != getPatches().iEnd(); ++it)
        {
            getBases().at(unk).matchInterface(*it, m_dofMappers[unk]);
        }

    for (size_t i = 0; i < boundaryDofs.size(); i++)
        m_dofMappers[unk].markBoundary(i, boundaryDofs[i]);

    m_dofMappers[unk].finalize();

    this->updateDofMappers();
    this->updateSizes();
}


template<class T, int MatOrder>
gsField<T> gsTMAssemblerBase<T, MatOrder>::constructSolution(const gsMatrix<T>& solVector, index_t unk) const
{
    GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

    gsMultiPatch<T>* result = new gsMultiPatch<T>;

    const gsDofMapper& mapper = m_dofMappers[unk];

    //const index_t dim = (unk == 0 ? m_tarDim : 1);
    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    index_t dofs = 0;
    for (index_t i = 0; i < unk; i++)
        dofs += m_kdofs[i];
    gsAsConstMatrix<T> sol(solVector.data() + dofs, m_kdofs[unk], 1);
    
    for (size_t p = 0; p < this->getPatches().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const index_t sz = this->getBases().at(unk).piece(p).size();
        coeffs.resize(sz, m_tarDim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, p)) // DoF value is in the solVector
            {
                coeffs.row(i) = sol.row(mapper.index(i, p));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof[unk].row(mapper.bindex(i, p));
            }
        }

        result->addPatch(this->getBases().at(unk).piece(p).makeGeometry(coeffs));
    }

    return gsField<T>(this->getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
}

template<class T, int MatOrder>
index_t gsTMAssemblerBase<T, MatOrder>::numDofsUnk(size_t i)
{
    if ((i >= 0) && (i < m_kdofs.size()))
        return m_kdofs[i];
    else
        GISMO_ERROR("numDofsUnk(i): i must be between 0 and the total number of turbulent variables.");
}

} // namespace gismo