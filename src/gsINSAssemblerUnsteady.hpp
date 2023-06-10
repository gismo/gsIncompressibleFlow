/** @file gsINSAssemblerUnsteady.hpp
  
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSAssemblerUnsteady.h>

namespace gismo
{

template<class T>
void gsINSAssemblerUnsteady<T>::initMembers()
{
    int dofs = Base::numDofs();

    m_baseMatrix.resize(dofs, dofs);
    m_baseRhs.setZero(dofs, 1);

    m_matrix.resize(dofs, dofs);
    m_rhs.setZero(dofs, 1);

    m_solution.setZero(dofs, 1);

    m_bInitialized = false;
    m_bMatrixReady = false;
    m_bRhsReady = false;

    m_blockAssembler.setUnsteady(true);
}


template<class T>
void gsINSAssemblerUnsteady<T>::fillBase()
{
    const T invTimeStep = 1. / m_timeStepSize;
    int numDofs = Base::numDofs();
    int uDofs = Base::getUdofs();
    int tarDim = Base::getTarDim();

    gsSparseMatrix<T> stokesMatrix(numDofs, numDofs);

    m_blockAssembler.fillStokesSystem_into(stokesMatrix, m_baseRhs);

    gsVector<int> nonZerosPerColumnVector;
    nonZerosPerColumnVector.setZero(numDofs);
    for (int s = 0; s < tarDim; ++s)
        for (int i = 0; i < uDofs; i++)
            nonZerosPerColumnVector(i + s * uDofs) = m_blockAssembler.getBlockNpattern().col(i).nonZeros();

    gsSparseMatrix<T> blockNpatternMatrix(numDofs, numDofs);
    blockNpatternMatrix.reserve(nonZerosPerColumnVector);

    #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
    for (index_t col = 0; col < uDofs; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern(), col); it; ++it)
            for (index_t s = 0; s < tarDim; ++s)
                blockNpatternMatrix.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;

    m_baseMatrix = stokesMatrix + blockNpatternMatrix;

    if (m_blockAssembler.isUnsteady())
    {
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += invTimeStep * it.value();
    }

    m_bMatrixReady = false;
    m_bRhsReady = false;
}


template<class T>
void gsINSAssemblerUnsteady<T>::fillMatrix()
{
    int uDofs = Base::getUdofs();

    m_matrix = m_baseMatrix;
    #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
    for (index_t col = 0; col < uDofs; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
            for (index_t s = 0; s < Base::getTarDim(); ++s)
                m_matrix.coeffRef(it.row() + s*uDofs, it.col() + s*uDofs) += it.value();

    m_bMatrixReady = true;

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();
}


template<class T>
void gsINSAssemblerUnsteady<T>::fillRhs()
{
    const T invTimeStep = 1. / m_timeStepSize;
    int uDofs = Base::getUdofs();

    m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();

    if (m_blockAssembler.isUnsteady())
    {
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < Base::getTarDim(); ++s)
            m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
    }

    m_bRhsReady = true;
}


template<class T>
void gsINSAssemblerUnsteady<T>::changeTimeStep(const T timeStepSize)
{
    GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

    m_timeStepSize = timeStepSize;

    m_matrix.resize(Base::numDofs(), Base::numDofs());
    m_rhs.setZero();

    fillBase();
}


} // namespace gismo
