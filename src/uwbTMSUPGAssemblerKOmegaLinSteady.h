/** @file uwbTMSUPGAssemblerKOmegaLinSteady.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMAssemblerKOmegaLinSteady.h"

namespace gismo
{

template<class T>
class uwbTMSUPGAssemblerKOmegaLinSteady : public uwbTMAssemblerKOmegaLinSteady<T>
{

public:
    typedef uwbTMAssemblerKOmegaLinSteady<T> Base;

public:
    uwbTMSUPGAssemblerKOmegaLinSteady(uwbINSSolverParams<T>& params) :
        Base(params)
    {
    }

    virtual ~uwbTMSUPGAssemblerKOmegaLinSteady()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        Base::initAssembly(uSolField);
        m_blockAssembler.assembleSUPGpatternBlocks_kOmega();
    }

    virtual void updateAssembly()
    {
        Base::updateAssembly();
        m_blockAssembler.assembleNonlinearPartSUPG_kOmega();
    }

    virtual void updatePicardAssembly()
    { 
        Base::updatePicardAssembly();
        m_blockAssembler.assembleNonlinearPartSUPG_kOmega();
    }

public:
    void updateTurbCoeffSolField(const gsMatrix<T> & solVector, const gsField<T> & uSolField)
    {
        //GISMO_ASSERT(m_blockAssembler.isTMsupg(), "SUPG stabilization must be set.");
        std::vector<gsField<T> > kOmegaDiffusionCoeffSol;
        kOmegaDiffusionCoeffSol.push_back(this->constructTurbCoeffSol(solVector, uSolField, "kDiffusionCoeff"));
        kOmegaDiffusionCoeffSol.push_back(this->constructTurbCoeffSol(solVector, uSolField, "oDiffusionCoeff"));
        m_blockAssembler.setKOmegaDiffusionCoeffSolField(kOmegaDiffusionCoeffSol);
    }
    //===================================================== fillBase ================================================================

    virtual void fillBase()  //+pattern
    {
        Base::fillBase();
        
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;

        //blockN
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlockN_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> matrixN_SUPG(varDofs, dofs);
        matrixN_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN_SUPG(), col); it; ++it)
                matrixN_SUPG.insert(it.row(), it.col()) = it.value();

        //block Apattern_SUPG
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksAnonlinPattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> matrixAnonlinPattern_SUPG(varDofs, dofs);
        matrixAnonlinPattern_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlinPattern_SUPG(), col); it; ++it)
                matrixAnonlinPattern_SUPG.insert(it.row(), it.col()) = 0.;

        m_baseMatrix += matrixN_SUPG + matrixAnonlinPattern_SUPG;

        m_baseRhs += m_blockAssembler.getRhsN_SUPG();

        if (!m_baseMatrix.isCompressed() && !m_blockAssembler.isTMcrosswind())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        Base::fillSystem();

        int dofs = this->numDofs();

        //SUPG block A nonlin.
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlin_SUPG(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += it.value();

        m_rhs += m_blockAssembler.getRhsAnonlin_SUPG();

        if (!m_matrix.isCompressed() && !m_blockAssembler.isTMcrosswind())
            m_matrix.makeCompressed();

        if (!m_blockAssembler.isTMcrosswind())
            m_bSystemReady = true;
    } //end fillSystem

protected:
    // members from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_baseRhs;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;

}; //uwbTMSUPGAssemblerKOmegaLinSteady

} //namespace gismo
