/** @file gsINSAssemblerBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSAssemblerBase.h>

namespace gismo
{

template<class T>
void gsINSAssemblerBase<T>::addPressureOutletCondition(int patch, boxSide side)
{
    gsMatrix<unsigned> boundaryIndicesOutlet = m_blockAssembler.getBases().back().basis(patch).boundary(side);

    std::vector< gsMatrix< index_t > > boundaryDofsToEliminate;
    boundaryDofsToEliminate.resize(m_blockAssembler.getBases().back().nBases());
    for (size_t i = 0; i < boundaryDofsToEliminate.size(); ++i)
    {
        boundaryDofsToEliminate[i].setZero(0, 0);
        if ( i == static_cast<size_t>(patch) )
            boundaryDofsToEliminate[i] = boundaryIndicesOutlet;
    }

    m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofsToEliminate, 1);

    reinitMembers();
}


template<class T>
void gsINSAssemblerBase<T>::fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
{
    Ap = m_blockAssembler.getPressurePoissonMatrix(assembAp, lumping);

    gsSparseMatrix<T> ApforFp;
    if (assembAp == assembFp)
        ApforFp = Ap;
    else
        ApforFp = m_blockAssembler.getPressurePoissonMatrix(assembFp, lumping);

    if(Fp.nonZeros()) // the matrix is not empty
        Fp += getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();
    else 
        Fp = getViscosity() * ApforFp + m_blockAssembler.getPressureConvectionMatrix();

    m_blockAssembler.applyPCDboundaryConditions(Ap, Fp, bcType);
}


} //namespace gismo
