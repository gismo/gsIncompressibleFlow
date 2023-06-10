/** @file gsINSAssemblerSteady.hpp
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSAssemblerSteady.h>

namespace gismo
{

template<class T>
void gsINSAssemblerSteady<T>::initMembers()
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
}


template<class T>
void gsINSAssemblerSteady<T>::fillBase()
{
    int uDofs = Base::getUdofs();
    int numDofs = Base::numDofs();
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

    m_bMatrixReady = false;
    m_bRhsReady = false;
}


template<class T>
void gsINSAssemblerSteady<T>::fillMatrix()
{
    int uDofs = Base::getUdofs();

    m_matrix = m_baseMatrix;

    #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
    for (index_t col = 0; col < uDofs; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
            for (index_t s = 0; s < Base::getTarDim(); ++s)
                m_matrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_bMatrixReady = true;
}

} 