/** @file uwbADRAssemblerNonlinear.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRAssemblerBase.h"
#include "uwbADRBlockVisitors.h"

namespace gismo
{

template<class T>
class uwbADRAssemblerNonlinear : public uwbADRAssemblerBaseUnsteady<T>
{

public:
    typedef uwbADRAssemblerBaseUnsteady<T> Base;

public:
    uwbADRAssemblerNonlinear(uwbADRSolverParams<T>& params) : Base(params, 1) { }

    uwbADRAssemblerNonlinear(uwbADRSolverParams<T>& params, gsField<T>& advectionField,
                             gsField<T>& diffusionField) : Base(params, 1)
    {
        m_blockAssembler.setAdvectionDiffusionCoeffFields(advectionField, diffusionField);
        //m_blockAssembler.setCoeffGeometry(patchesCoeffs, basesCoeffs);
    }

    uwbADRAssemblerNonlinear(uwbADRSolverParams<T>& params, gsField<T>& advectionField,
                             gsField<T>& diffusionField, gsField<T>& reactionField) : Base(params, 1)
    {
        m_blockAssembler.setADRcoeffFields(advectionField, diffusionField, reactionField);
        //m_blockAssembler.setCoeffGeometry(patchesCoeffs, basesCoeffs);
    }

    uwbADRAssemblerNonlinear(uwbADRSolverParams<T>& params, T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff) : Base(params, 1)
    {
        m_blockAssembler.setLinearADRCoefficients(diffusionCoeff, advectionCoeff, reactionCoeff);
    }

    virtual ~uwbADRAssemblerNonlinear() { }

protected:

    virtual void initSteadyAssembly()
    {
        m_blockAssembler.assemblePatternBlocks();
    }

    virtual void initAssembly()
    {
        if (m_blockAssembler.isUnsteady())
            m_blockAssembler.assembleMassMatrix();
        m_blockAssembler.assemblePatternBlocks();
    }

    virtual void updateAssembly()
    {
        m_blockAssembler.assembleNonlinearBlocks(m_solution);
        if (m_blockAssembler.isSUPG() && m_blockAssembler.isUnsteady())
            m_blockAssembler.assembleMassMatrixSUPG();
    }

    virtual void updatePicardAssembly(const gsMatrix<T> & solVector)
    { 
        m_blockAssembler.assembleNonlinearBlocks(solVector, false);
        if (m_blockAssembler.isSUPG() && m_blockAssembler.isUnsteady())
            m_blockAssembler.assembleMassMatrixSUPG();
    }

//==========================================================================================================================
public:

    virtual void fillBase()
    {
        fillBaseLinear();
        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = true;
    }


    virtual void fillBaseLinear()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        const T invTimeStep = 1. / m_timeStepSize;

        reserveSpace();

        if (m_blockAssembler.isUnsteady())
        {
            //block mass matrix
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < varDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrix(), col); it; ++it)
                    for (index_t s = 0; s < numVar; ++s)
                        m_baseMatrix.coeffRef(it.row(), it.col() + s*varDofs) = invTimeStep * it.value();
        }

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

    } //end fillBaseLinear

    virtual void fillSystem()
    {
        m_matrix = m_baseMatrix;
        fillBlocks();
        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        m_bSystemReady = true;
    }

    //only unsteady case
    virtual void fillRhs()
    {
        int numVar = this->getNumVar();
        const T invTimeStep = 1. / m_timeStepSize;

        m_rhs = m_baseRhs;

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < numVar; ++s)
            m_rhs.col(s).noalias() += invTimeStep * m_blockAssembler.getMassMatrix() * m_solution.col(s);
        //================= SUPG part ===========================
        if (m_blockAssembler.isSUPG())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < numVar; ++s)
                m_rhs.col(s).noalias() += invTimeStep * m_blockAssembler.getMassMatrix_SUPG() * m_solution.col(s);
        }
    }

protected:
    virtual void fillBlocks()
    {
        int dofs = this->numDofs();
        int varDofs = this->numVarDofs();
        int numVar = this->getNumVar();
        const T invTimeStep = 1. / m_timeStepSize;

        //ADR matrix
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrix(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += it.value();

        //================= SUPG part ===========================
        if (m_blockAssembler.isSUPG())
        {
            //SUPG block A
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrix_SUPG(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();

            if (m_blockAssembler.isUnsteady())
            {
                //SUPG mass matrix
                #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
                for (index_t col = 0; col < varDofs; ++col)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrix_SUPG(), col); it; ++it)
                        for (index_t s = 0; s < numVar; ++s)
                            m_matrix.coeffRef(it.row(), it.col() + s*varDofs) += invTimeStep * it.value();
            }
        }

        //================= CROSSWIND part ===========================
        if (m_blockAssembler.isCrosswind())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getCrosswindMatrix(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();
        }

        //=========== isotropic artificial diffusion part =================
        if (m_blockAssembler.isIsoArtificialDiff())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getIsoArtificialDiffMatrix(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();
        }

        //=========== artificial diffusion part =================
        if (m_blockAssembler.isArtificialDiff())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getArtificialDiffMatrix(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();
        }

        //---------------------------
        //--------- fillRhs ---------
        //---------------------------
        m_baseRhs = m_blockAssembler.getRhsADR();

        //================= SUPG part ===========================
        if (m_blockAssembler.isSUPG())
            m_baseRhs += m_blockAssembler.getRhsADR_SUPG();
        if (m_blockAssembler.isCrosswind())
            m_baseRhs += m_blockAssembler.getRhsCrosswind();
        if (m_blockAssembler.isIsoArtificialDiff())
            m_baseRhs += m_blockAssembler.getRhsIsoArtificialDiff();
        if (m_blockAssembler.isArtificialDiff())
            m_baseRhs += m_blockAssembler.getRhsArtificialDiff();


    } //end

    void reserveSpace()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t s = 0; s < numVar; ++s)
            for (index_t i = 0; i < varDofs; ++i)
                nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getADRMatrixPattern().col(i).nonZeros();

        if (m_blockAssembler.isUnsteady())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getMassMatrix().col(i).nonZeros();
        }

        if (m_blockAssembler.isSUPG())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getADRMatrixPattern_SUPG().col(i).nonZeros();
            if (m_blockAssembler.isUnsteady())
            {
                for (index_t s = 0; s < numVar; ++s)
                    for (index_t i = 0; i < varDofs; ++i)
                        nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getMassMatrixPattern_SUPG().col(i).nonZeros();
            }
        }

        if (m_blockAssembler.isCrosswind())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getCrosswindPatternMatrix().col(i).nonZeros();
        }

        if (m_blockAssembler.isIsoArtificialDiff())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getIsoArtificialDiffPatternMatrix().col(i).nonZeros();
        }

        if (m_blockAssembler.isArtificialDiff())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getArtificialDiffPatternMatrix().col(i).nonZeros();
        }

        m_baseMatrix.reserve(nonZerosPerColumnVector);

        // ADR matrix
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrixPattern(), col); it; ++it)
                    m_baseMatrix.insert(it.row(), it.col()) = 0.;

        if (m_blockAssembler.isUnsteady())
        {
            //Mass matrix
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < varDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrix(), col); it; ++it)
                    for (index_t s = 0; s < numVar; ++s)
                        m_baseMatrix.coeffRef(it.row(), it.col() + s*varDofs) = 0.;
        }

        //================= SUPG part ===========================
        if (m_blockAssembler.isSUPG())
        {
            if (m_blockAssembler.isUnsteady())
            {
                //SUPG mass matrix
                #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
                for (index_t col = 0; col < varDofs; ++col)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrixPattern_SUPG(), col); it; ++it)
                        for (index_t s = 0; s < numVar; ++s)
                            m_baseMatrix.coeffRef(it.row(), it.col() + s*varDofs) = 0.;
            }

            //SUPG ADR matrix
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrixPattern_SUPG(), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row(), it.col()) = 0.;
        }

        //================= CROSSWIND part ===========================
        if (m_blockAssembler.isCrosswind())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getCrosswindPatternMatrix(), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row(), it.col()) = 0.;
        }

        //========= isotropic artificial diffusion part ================
        if (m_blockAssembler.isIsoArtificialDiff())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getIsoArtificialDiffPatternMatrix(), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row(), it.col()) = 0.;
        }

        //========= artificial diffusion part ================
        if (m_blockAssembler.isArtificialDiff())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getArtificialDiffPatternMatrix(), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row(), it.col()) = 0.;
        }
    }


protected:
    // members from uwbADRAssemblerBaseUnsteady
    using Base::m_timeStepSize;

    // members from uwbADRAssemblerSteady
    using uwbADRAssemblerSteady<T>::m_blockAssembler;
    using uwbADRAssemblerSteady<T>::m_baseMatrix;
    using uwbADRAssemblerSteady<T>::m_baseRhs;
    using uwbADRAssemblerSteady<T>::m_matrix;
    using uwbADRAssemblerSteady<T>::m_rhs;
    using uwbADRAssemblerSteady<T>::m_solution;
    using uwbADRAssemblerSteady<T>::m_bSystemReady;

}; //uwbADRAssemblerNonlinear

} //namespace gismo
