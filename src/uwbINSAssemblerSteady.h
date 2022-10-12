/** @file uwbINSAssemblerSteady.h
    
    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include "uwbINSAssemblerBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerSteady : public uwbINSAssemblerBase<T>
{

public:
    typedef uwbINSAssemblerBase<T> Base;

public:
    uwbINSAssemblerSteady(uwbINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSAssemblerSteady()
    {
    }

protected:

    void initMembers()
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

    virtual void reinitMembers() { initMembers(); }

public:

    virtual void initialize()
    {
        m_blockAssembler.assembleLinearStokesPart();
        m_blockAssembler.assembleBlockNpattern();

        fillBase();

        if (m_blockAssembler.isSUPG())
        {
            m_blockAssembler.assembleBlockNonlinearSUPGpattern();

            fillSUPGBase();
        }

        m_bInitialized = true;
    }

    virtual void updateAssembly()
    {
        Base::updateAssembly();

        if (m_blockAssembler.isSUPG())
            m_blockAssembler.assembleNonlinearSUPGPart();
        //if (m_blockAssembler.isPSPG())
        //    m_blockAssembler.assemblePSPGPart();

    }

protected:

    virtual void fillBase()
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

protected:

    virtual void fillMatrix()
    {
        int uDofs = Base::getUdofs();

        m_matrix = m_baseMatrix;

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < Base::getTarDim(); ++s)
                    m_matrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();

        if (m_blockAssembler.isSUPG())
            fillSUPGMatrix();

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        m_bMatrixReady = true;
    }

protected:

    virtual void fillRhs()
    {
        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();

        if (m_blockAssembler.isSUPG())
            fillSUPGRhs();

        //if (m_blockAssembler.isPSPG())
        //    m_rhs.noalias() += m_blockAssembler.getRhsPSPG();

        m_bRhsReady = true;
    }

    //----------------------- SUPG part -----------------------------------
    void fillSUPGBase()
    {
        int uDofs = Base::getUdofs();
        int pDofs = Base::getPdofs();
        int tarDim = Base::getTarDim();
        int numDofs = Base::numDofs();
        int pShift = Base::getPshift();

        gsVector<int> nonZerosPerColumnVector;

        //------ A_SUPG -----------
        nonZerosPerColumnVector.setZero(numDofs);
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(s * uDofs + i) = m_blockAssembler.getBlockApattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> blockApatternMatrix_SUPG(numDofs, numDofs);
        blockApatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockApattern_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockApatternMatrix_SUPG.insert(s * uDofs + it.row(), s * uDofs + it.col()) = 0.;


        //------ N_SUPG -----------
        nonZerosPerColumnVector.setZero(numDofs);
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(i + s * uDofs) = m_blockAssembler.getBlockNpattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> blockNpatternMatrix_SUPG(numDofs, numDofs);
        blockNpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockNpatternMatrix_SUPG.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;

        //------ B_SUPG -----------
        nonZerosPerColumnVector.setZero(numDofs);
        for (int i = 0; i < pDofs; i++)
            for (int s = 0; s < tarDim; ++s)
                nonZerosPerColumnVector(i + pShift) += m_blockAssembler.getBlockBpattern_SUPG(s).col(i).nonZeros();

        gsSparseMatrix<T> blockBpatternMatrix_SUPG(numDofs, numDofs);
        blockBpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < pDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockBpattern_SUPG(s), col); it; ++it)
                    blockBpatternMatrix_SUPG.insert(it.row() + s * uDofs, it.col() + pShift) = 0.;


        m_baseMatrix += blockApatternMatrix_SUPG + blockNpatternMatrix_SUPG + blockBpatternMatrix_SUPG;

        if (m_blockAssembler.isRotation())
        {
            //------ M_SUPG -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(i + s * uDofs) = m_blockAssembler.getBlockMpattern_SUPG().col(i).nonZeros();

            gsSparseMatrix<T> blockMpatternMatrix_SUPG(numDofs, numDofs);
            blockMpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockMpattern_SUPG(), col); it; ++it)
                {
                    blockMpatternMatrix_SUPG.insert(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) = 0.;
                    blockMpatternMatrix_SUPG.insert(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) = 0.;
                }

            m_baseMatrix += blockMpatternMatrix_SUPG;
        }

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    void fillSUPGMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        //------ A_SUPG -----------
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_matrix.coeffRef(s * uDofs + it.row(), s * uDofs + it.col()) += it.value();

        //------ N_SUPG -----------
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_matrix.coeffRef(it.row() + s* uDofs, it.col() + s * uDofs) += it.value();

        //------ B_SUPG -----------
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < Base::getPdofs(); ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockB_SUPG(s), col); it; ++it)
                {
                    m_matrix.coeffRef(it.row() + s*uDofs, it.col() + Base::getPshift()) += it.value();
                }

        if (m_blockAssembler.isRotation())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM_SUPG(), col); it; ++it)
                {
                    m_matrix.coeffRef(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) += -m_blockAssembler.getOmega() * it.value();
                    m_matrix.coeffRef(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) += m_blockAssembler.getOmega() * it.value();
                }
        }
    }

    void fillSUPGRhs()
    {
        int tarDim = Base::getTarDim();
        int uDofs = Base::getUdofs();

        m_rhs.noalias() += m_blockAssembler.getRhsN_SUPG() + m_blockAssembler.getRhsA_SUPG() + m_blockAssembler.getRhsB_SUPG();

        // BC for rotation term
        if (m_blockAssembler.isRotation())
        {
            m_rhs.middleRows((tarDim - 2) * uDofs, uDofs) -= m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 1) * uDofs, uDofs);
            m_rhs.middleRows((tarDim - 1) * uDofs, uDofs) += m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 2) * uDofs, uDofs);
        }
    }
    //---------------------------------------------------------------------

public:

    virtual const gsSparseMatrix<T> & matrix() const
    {
        GISMO_ASSERT(m_bMatrixReady, "Matrix not ready, update() must be called first");
        return m_matrix;
    }

    virtual const gsMatrix<T> & rhs() const
    {
        GISMO_ASSERT(m_bRhsReady, "Rhs not ready, update() must be called first");
        return m_rhs;
    }

    void evalResiduum(const gsMatrix<T> & solVector, std::vector<T> & residuum)
    {
        m_blockAssembler.evaluateINSsteadyResiduumL2norm(solVector, residuum);
    }


protected:

    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;

    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;

    bool m_bMatrixReady;
    bool m_bRhsReady;

    // members from uwbINSAssemblerBase
    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;

}; // class uwbINSAssemblerSteady

} // namespace gismo


