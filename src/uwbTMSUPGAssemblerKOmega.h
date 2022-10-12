/** @file uwbTMSUPGAssemblerKOmega.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMAssemblerKOmega.h"
//#include "uwbTMStabilizationBlockAssembler.h"

namespace gismo
{

template<class T>
class uwbTMSUPGAssemblerKOmega : public uwbTMAssemblerKOmega<T>
{

public:
    typedef uwbTMAssemblerKOmega<T> Base;

public:
    uwbTMSUPGAssemblerKOmega(uwbINSSolverParams<T>& params) : Base(params)
    {
    }

    virtual ~uwbTMSUPGAssemblerKOmega()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        Base::initAssembly(uSolField);
        m_blockAssembler.assembleSUPGpatternBlocks_kOmega();
    }

    virtual void updateAssembly(const gsField<T> & uSolField)
    {
        Base::updateAssembly(uSolField);
        //m_blockAssembler.assembleNewTimestepPartSUPG_kOmega();
        //m_blockAssembler.assembleNonlinearPartSUPG_kOmega();

        m_blockAssembler.assembleSUPGblocks_kOmega();
    }

    virtual void updatePicardAssembly()
    { 
        Base::updatePicardAssembly();
        //m_blockAssembler.assembleNonlinearPartSUPG_kOmega();

        m_blockAssembler.assembleSUPGblocks_kOmega();
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

        //blocks pattern_SUPG
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksPattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> matrixPattern_SUPG(varDofs, dofs);
        matrixPattern_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksPattern_SUPG(), col); it; ++it)
                matrixPattern_SUPG.insert(it.row(), it.col()) = 0.;

        m_baseMatrix += matrixPattern_SUPG;

        if (!m_baseMatrix.isCompressed() && !m_blockAssembler.isTMcrosswind())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

public:
    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        Base::fillSystem();

        int dofs = this->numDofs();
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        const T invTimeStep = 1. / m_timeStepSize;

        //SUPG blocks
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocks_SUPG(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += it.value();

        m_rhs += m_blockAssembler.getRhs_SUPG();
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < numVar; ++s)
            m_rhs.col(s).noalias() += invTimeStep * m_blockAssembler.getBlockM_SUPG().middleCols(s * varDofs, varDofs) * m_solution.col(s);

        if (!m_matrix.isCompressed() && !m_blockAssembler.isTMcrosswind())
            m_matrix.makeCompressed();

        if (!m_blockAssembler.isTMcrosswind())
            m_bSystemReady = true;

    } //end fillSystem

protected:
    //uwbTMStabilizationBlockAssembler<T>* m_pBlockAssembler;

    // members from uwbTMAssemblerBaseUnsteady
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseMatrix;
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseRhs;
    using uwbTMAssemblerBaseUnsteady<T>::m_timeStepSize;

    // members from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;

}; //uwbTMSUPGAssemblerKOmega

} //namespace gismo
