/** @file uwbTMCrosswindAssemblerKOmegaLinSteady.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMAssemblerKOmegaLinSteady.h"
#include "uwbTMSUPGAssemblerKOmegaLinSteady.h"

namespace gismo
{

template<class T>
class uwbTMCrosswindAssemblerKOmegaLinSteady : public uwbTMSUPGAssemblerKOmegaLinSteady<T>
{

public:
    typedef uwbTMSUPGAssemblerKOmegaLinSteady<T> Base;

public:
    uwbTMCrosswindAssemblerKOmegaLinSteady(uwbINSSolverParams<T>& params) :
        Base(params)
    { }

    virtual ~uwbTMCrosswindAssemblerKOmegaLinSteady()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        Base::initAssembly(uSolField);
        m_blockAssembler.assembleCrosswindPatternBlock_kOmega();
    }

    virtual void updateAssembly()
    {
        Base::updateAssembly();
        m_blockAssembler.assembleCrosswindBlock_kOmega();
    }

    virtual void updatePicardAssembly()
    { 
        Base::updatePicardAssembly();
        m_blockAssembler.assembleCrosswindBlock_kOmega();
    }

public:
    //===================================================== fillBase ================================================================

    virtual void fillBase()  //+pattern
    {
        Base::fillBase();
        
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlockCrosswindPattern().col(i).nonZeros();

        gsSparseMatrix<T> matrixCrosswind(varDofs, dofs);
        matrixCrosswind.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswindPattern(), col); it; ++it)
                    matrixCrosswind.insert(it.row(), it.col()) = 0.;

        m_baseMatrix += matrixCrosswind;

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        Base::fillSystem();

        int dofs = this->numDofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswind(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += it.value();

        m_rhs += m_blockAssembler.getRhsCrosswind();

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        m_bSystemReady = true;
    } //end fillSystem

protected:
    //using Base::m_pBlockAssembler;
    using Base::getBlockAssembler;

    // members from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_baseRhs;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;

}; //uwbTMCrosswindAssemblerKOmegaLinSteady

} //namespace gismo
