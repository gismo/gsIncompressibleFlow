/** @file uwbADRAssemblerLinear.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRAssemblerBase.h"
#include "uwbADRBlockVisitors.h"

namespace gismo
{

template<class T>
class uwbADRAssemblerLinear : public uwbADRAssemblerBaseUnsteady<T>
{

public:
    typedef uwbADRAssemblerBaseUnsteady<T> Base;

public:
    uwbADRAssemblerLinear(uwbADRSolverParams<T>& params, T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff) :
        Base(params, 1)
    {
        m_blockAssembler.setLinearADRCoefficients(diffusionCoeff, advectionCoeff, reactionCoeff);
    }

    uwbADRAssemblerLinear(uwbADRSolverParams<T>& params) : Base(params, 1)
    { }

    virtual ~uwbADRAssemblerLinear() { }

protected:

    virtual void initSteadyAssembly()
    {
        m_blockAssembler.assembleADR();
    }

    virtual void initAssembly()
    {
        m_blockAssembler.assembleMassMatrix();
        m_blockAssembler.assembleADR();
    }

public:
    //===================================================== fillSystem ==============================================================
    virtual void fillBase()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        const T invTimeStep = 1. / m_timeStepSize;

        reserveSpace();

        fillADRMatrix();

        if (m_blockAssembler.isUnsteady())
        {
            //mass matrix
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < varDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrix(), col); it; ++it)
                    for (index_t s = 0; s < numVar; ++s)
                        m_matrix.coeffRef(it.row(), it.col() + s*varDofs) += invTimeStep * it.value();

            //================= SUPG part ===========================
            if (m_blockAssembler.isSUPG())
            {
                //SUPG mass matrix
                #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
                for (index_t col = 0; col < varDofs; ++col)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getMassMatrix_SUPG(), col); it; ++it)
                        for (index_t s = 0; s < numVar; ++s)
                            m_matrix.coeffRef(it.row(), it.col() + s*varDofs) += invTimeStep * it.value();
            }
        }

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        m_bSystemReady = true;
    }

    virtual void fillADRMatrix()
    {
        int dofs = this->numDofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrix(), col); it; ++it)
                m_matrix.insert(it.row(), it.col()) = it.value();

        //================= SUPG part ===========================
        if (m_blockAssembler.isSUPG())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getADRMatrix_SUPG(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();
        }

        //========= artificial diffusion part =====================
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
        //========= artificial diffusion part ===================
        if (m_blockAssembler.isArtificialDiff())
            m_baseRhs += m_blockAssembler.getRhsArtificialDiff();

    } //end

    //only in unsteady case
    void fillRhs()
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


    void reserveSpace()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        if (m_blockAssembler.isUnsteady())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getMassMatrix().col(i).nonZeros() +
                                                             m_blockAssembler.getADRMatrix().col(i).nonZeros();
            if (m_blockAssembler.isSUPG())
            {
                for (index_t s = 0; s < numVar; ++s)
                    for (index_t i = 0; i < varDofs; ++i)
                        nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getMassMatrix_SUPG().col(i).nonZeros() +
                                                                 m_blockAssembler.getADRMatrix_SUPG().col(i).nonZeros();
            }
        }
        else
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getADRMatrix().col(i).nonZeros();
            if (m_blockAssembler.isSUPG())
            {
                for (index_t s = 0; s < numVar; ++s)
                    for (index_t i = 0; i < varDofs; ++i)
                        nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getADRMatrix_SUPG().col(i).nonZeros();
            }
        }

        if (m_blockAssembler.isArtificialDiff())
        {
            for (index_t s = 0; s < numVar; ++s)
                for (index_t i = 0; i < varDofs; ++i)
                    nonZerosPerColumnVector(i + s*varDofs) += m_blockAssembler.getArtificialDiffMatrix().col(i).nonZeros();
        }

        m_matrix.reserve(nonZerosPerColumnVector);
    }

protected:
    // members from uwbADRAssemblerBaseUnsteady
    using Base::m_baseRhs;
    using Base::m_timeStepSize;

    // members from uwbADRAssemblerSteady
    using uwbADRAssemblerSteady<T>::m_blockAssembler;
    using uwbADRAssemblerSteady<T>::m_matrix;
    using uwbADRAssemblerSteady<T>::m_rhs;
    using uwbADRAssemblerSteady<T>::m_solution;
    using uwbADRAssemblerSteady<T>::m_bSystemReady;

}; //uwbADRAssemblerLinear

} //namespace gismo
