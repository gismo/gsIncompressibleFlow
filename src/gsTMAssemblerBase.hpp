/** @file gsTMAssemblerBase.hpp

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

    GISMO_ASSERT(m_numVars > 0, "Number of TM variables not set yet.");

    m_viscosity = m_paramsPtr->getPde().viscosity();

    //m_bases = m_paramsPtr->getBases();
    //m_numVars = m_bases.size() - 2;
    m_kdofs.resize(m_numVars);

    m_dofMappers.resize(m_numVars+2);
    m_ddof.resize(m_numVars+2);
    gsMatrix<T> ddof;
    m_ddof[0] = ddof.setZero(1, 1);
    m_ddof[1] = ddof.setZero(1, 1);
    
    m_nnzPerRowTM = 1;
    index_t maxDeg;
    // for (short_t i = 0; i < m_tarDim; i++)
    //     for (short_t j = 0; j < m_numVars; j++)
    //         if (m_bases[j].maxDegree(i) > maxDeg)
    //             maxDeg = m_bases[j].maxDegree(i);
    // m_nnzPerRowTM = 2 * maxDeg + 1;
    // m_nnzPerRowTM = 2 * m_nnzPerRowTM;
    for (short_t i = 0; i < m_tarDim; i++)
    {
        maxDeg = 0;
        for (short_t j = 0; j < m_numVars; j++)
            if (getBasis(2+j).maxDegree(i) > maxDeg)
                maxDeg = getBasis(2+j).maxDegree(i);

        m_nnzPerRowTM *= 2 * maxDeg + 1;
    }
    
    m_isInitialized = false;

}

template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::updateSizes()
{
    m_dofs = 0;
    gsMatrix<T> ddof;
    for (short_t i = 0; i < m_numVars; i++)
    {
        m_kdofs[i] = m_dofMappers[i+2].freeSize();
        m_dofs += m_kdofs[i];

        ddof.setZero(m_dofMappers[i+2].boundarySize(), 1);
        m_ddof[i+2] = ddof;

        if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            this->computeDirichletDofs(i+2, i+2, m_ddof[i+2]);
        }
    }
    
    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);
}


template<class T, int MatOrder>
void gsTMAssemblerBase<T, MatOrder>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
{
    m_dofMappers[unk] = gsDofMapper(getBasis(unk), getBCs(), unk);

    if (getAssemblerOptions().intStrategy == iFace::conforming)
        for (gsBoxTopology::const_iiterator it = getPatches().iBegin(); it != getPatches().iEnd(); ++it)
        {
            getBasis(unk).matchInterface(*it, m_dofMappers[unk]);
        }

    for (size_t i = 0; i < boundaryDofs.size(); i++)
        m_dofMappers[unk].markBoundary(i, boundaryDofs[i]);

    m_dofMappers[unk].finalize();

    this->updateDofMappers();
    this->updateSizes();
}

// TODO: update with respect to the new version in gsINSAssembler.hpp
template<class T, int MatOrder>
gsField<T> gsTMAssemblerBase<T, MatOrder>::constructSolution(const gsMatrix<T>& solVector, index_t unk, bool customSwitch) const
{
    GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

    gsMultiPatch<T>* result = new gsMultiPatch<T>;

    const gsDofMapper& mapper = m_dofMappers[unk];

    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    index_t dofs = 0;
    for (index_t i = 2; i < unk; i++)
        dofs += m_kdofs[i-2];
    
    for (size_t p = 0; p < this->getPatches().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const index_t sz = this->getBasis(unk).piece(p).size();
        coeffs.resize(sz, 1);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, p)) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row(mapper.index(i, p) + dofs);
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof[unk].row(mapper.bindex(i, p));
            }
        }

        result->addPatch(this->getBasis(unk).piece(p).makeGeometry(coeffs));
    }

    return gsField<T>(this->getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
}

template<class T, int MatOrder>
index_t gsTMAssemblerBase<T, MatOrder>::numDofsUnk(size_t i)
{
    if ((i >= 0) && (i < m_kdofs.size()))
        return m_kdofs[i];
    else
        GISMO_ERROR("numDofsUnk(i): i must be greater than or equal to 0 and less than the total number of turbulent variables.");
}

} // namespace gismo