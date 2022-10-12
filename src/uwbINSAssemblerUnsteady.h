/** @file uwbINSAssemblerUnsteady.h

Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include "uwbINSAssemblerBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerUnsteady : public uwbINSAssemblerBase<T>
{

public:
    typedef uwbINSAssemblerBase<T> Base;

public:
    uwbINSAssemblerUnsteady(uwbINSSolverParams<T>& params) : Base(params)
    {
        m_timeStepSize = params.settings().get(constantsINS::timeStep);
        initMembers();
    }

    virtual ~uwbINSAssemblerUnsteady()
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

        m_blockAssembler.setUnsteady(true);
    }

    virtual void reinitMembers() { initMembers(); }

public:

    virtual void initialize()
    {
        m_blockAssembler.assembleLinearStokesPart();
        m_blockAssembler.assembleBlockNpattern();

        fillBase();

        if (m_blockAssembler.isCROSSWIND())
        {
            m_blockAssembler.assembleBlockCROSSWINDpattern();
            fillCROSSWINDBase();
        }

        if (m_blockAssembler.isSUPG())
        {
            m_blockAssembler.assembleBlockNonlinearSUPGpattern();
            fillSupgCsdBase();
        }

        if (m_blockAssembler.isTCSD())
        {
            m_blockAssembler.assembleBlockNonlinearCSDpattern();
            fillSupgCsdBase();
        }

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bInitialized = true;
    }

    virtual void initialize(gsMatrix<T> & initialConditionCoeffVector)
    {
        initialize();
        Base::update(initialConditionCoeffVector);
    }

    virtual void updateAssembly()
    {
        Base::updateAssembly();

        if (m_blockAssembler.isSUPG())
            m_blockAssembler.assembleNonlinearSUPGPart();

        if (m_blockAssembler.isTCSD())
            m_blockAssembler.assembleNonlinearCSDPart();

        if (m_blockAssembler.isCROSSWIND())
            m_blockAssembler.assembleCrosswindPart();
    }

    virtual void updatePicardAssembly()
    {
        m_blockAssembler.assembleNonlinearPart();

        if (m_blockAssembler.isSUPG())
            m_blockAssembler.assembleNonlinearSUPGPart();

        if (m_blockAssembler.isTCSD())
            m_blockAssembler.assembleNonlinearCSDPart();

        if (m_blockAssembler.isCROSSWIND())
            m_blockAssembler.assembleCrosswindPart();
    }

    virtual void changeTimeStep(const T timeStepSize)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_timeStepSize = timeStepSize;

        m_matrix.resize(Base::numDofs(), Base::numDofs());
        m_rhs.setZero();

        fillBase();
    }

protected:

    virtual void fillBase()
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

        if (m_blockAssembler.isTimeDerTerm())
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

    virtual void fillMatrix()
    {
        int uDofs = Base::getUdofs();

        m_matrix = m_baseMatrix;
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < Base::getTarDim(); ++s)
                    m_matrix.coeffRef(it.row() + s*uDofs, it.col() + s*uDofs) += it.value();

        if (m_blockAssembler.isSUPG() || m_blockAssembler.isTCSD())
            fillSupgCsdMatrix();
        if (m_blockAssembler.isCROSSWIND())
            fillCROSSWINDMatrix();

        m_bMatrixReady = true;

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();
    }

    void fillRhs()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();

        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();

        if (m_blockAssembler.isTimeDerTerm())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < Base::getTarDim(); ++s)
                m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
        }

        if (m_blockAssembler.isSUPG() || m_blockAssembler.isTCSD())
            fillSupgCsdRhs();
        if (m_blockAssembler.isCROSSWIND())
            m_rhs.noalias() += m_blockAssembler.getRhsCrosswind();

        m_bRhsReady = true;
    }

    //--------------------------------------------------------------------------------------------

    void fillSupgCsdBase()
    {
        int numDofs = Base::numDofs();
        int uDofs = Base::getUdofs();
        int pDofs = Base::getPdofs();
        int tarDim = Base::getTarDim();
        int pShift = Base::getPshift();

        gsVector<int> nonZerosPerColumnVector;
        gsSparseMatrix<T> blockApatternMatrix_SUPG(numDofs, numDofs);
        gsSparseMatrix<T> blockBpatternMatrix_SUPG(numDofs, numDofs);
        if (m_blockAssembler.isSUPG())
        {
            //------ A_SUPG -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(s * uDofs + i) = m_blockAssembler.getBlockApattern_SUPG().col(i).nonZeros();

            blockApatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockApattern_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        blockApatternMatrix_SUPG.insert(s * uDofs + it.row(), s * uDofs + it.col()) = 0.;

            //------ B_SUPG -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int i = 0; i < pDofs; i++)
                for (int s = 0; s < tarDim; ++s)
                    nonZerosPerColumnVector(i + pShift) += m_blockAssembler.getBlockBpattern_SUPG(s).col(i).nonZeros();

            blockBpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < pDofs; ++col)
                for (index_t s = 0; s < tarDim; ++s)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockBpattern_SUPG(s), col); it; ++it)
                        blockBpatternMatrix_SUPG.insert(it.row() + s*uDofs, it.col() + pShift) = 0.;
        }

        //------ N_SUPG -----------
        nonZerosPerColumnVector.setZero(numDofs);
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < uDofs; i++)
                nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockNpattern_SUPG().col(i).nonZeros();

        gsSparseMatrix<T> blockNpatternMatrix_SUPG(numDofs, numDofs);
        blockNpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    blockNpatternMatrix_SUPG.insert(it.row() + s* uDofs, it.col() + s* uDofs) = 0.;

        if (m_blockAssembler.isTimeDerTerm())
        {
            //------ M_SUPG unsteady-----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockMpattern_SUPG().col(i).nonZeros();

            gsSparseMatrix<T> blockMpatternMatrix_SUPG(numDofs, numDofs);
            blockMpatternMatrix_SUPG.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockMpattern_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        blockMpatternMatrix_SUPG.insert(it.row() + s* uDofs, it.col() + s* uDofs) = 0.;

            m_baseMatrix += blockMpatternMatrix_SUPG;
        }

        if (m_blockAssembler.isTCSD())
            m_baseMatrix += blockNpatternMatrix_SUPG;
        else
            m_baseMatrix += blockNpatternMatrix_SUPG + blockBpatternMatrix_SUPG + blockApatternMatrix_SUPG;

        if (m_blockAssembler.isRotation())
        {
            //------ M_SUPG rotation -----------
            nonZerosPerColumnVector.setZero(numDofs);
            for (int s = 0; s < tarDim; ++s)
                for (int i = 0; i < uDofs; i++)
                    nonZerosPerColumnVector(i + s *  uDofs) = m_blockAssembler.getBlockMpattern_SUPG().col(i).nonZeros();

            gsSparseMatrix<T> blockMpatternMatrix_SUPG_rot(numDofs, numDofs);
            blockMpatternMatrix_SUPG_rot.reserve(nonZerosPerColumnVector);

            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockMpattern_SUPG(), col); it; ++it)
                {
                    blockMpatternMatrix_SUPG_rot.insert(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) = 0.;
                    blockMpatternMatrix_SUPG_rot.insert(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) = 0.;
                }

            m_baseMatrix += blockMpatternMatrix_SUPG_rot;
        }

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    void fillSupgCsdMatrix()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        if (m_blockAssembler.isSUPG())
        {
            //------ A_SUPG -----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        m_matrix.coeffRef(s * uDofs + it.row(), s * uDofs + it.col()) += it.value();

            //------ B_SUPG -----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < Base::getPdofs(); ++col)
                for (index_t s = 0; s < tarDim; ++s)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockB_SUPG(s), col); it; ++it)
                        m_matrix.coeffRef(it.row() + s*uDofs, it.col() + Base::getPshift()) += it.value();
        }

        //------ N_SUPG -----------
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN_SUPG(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_matrix.coeffRef(it.row() + s* uDofs, it.col() + s* uDofs) += it.value();

        if (m_blockAssembler.isTimeDerTerm())
        {
            //------ M_SUPG unsteady-----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < tarDim; ++s)
                        m_matrix.coeffRef(it.row() + s* uDofs, it.col() + s* uDofs) += invTimeStep * it.value();
        }

        if (m_blockAssembler.isRotation())
        {
            //------ M_SUPG rotation-----------
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < uDofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM_SUPG(), col); it; ++it)
                {
                    m_matrix.coeffRef(it.row() + ((tarDim - 2) * uDofs), it.col() + ((tarDim - 1) * uDofs)) = -m_blockAssembler.getOmega() * it.value();
                    m_matrix.coeffRef(it.row() + ((tarDim - 1) * uDofs), it.col() + ((tarDim - 2) * uDofs)) = m_blockAssembler.getOmega() * it.value();
                }

        }
    }

    virtual void fillSupgCsdRhs()
    {
        const T invTimeStep = 1. / m_timeStepSize;
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        if (m_blockAssembler.isSUPG())
            m_rhs.noalias() += m_blockAssembler.getRhsN_SUPG() + m_blockAssembler.getRhsB_SUPG() + m_blockAssembler.getRhsA_SUPG();
        else
            m_rhs.noalias() += m_blockAssembler.getRhsN_SUPG();

        if (m_blockAssembler.isTimeDerTerm())
        {
            //m_rhs.noalias() += invTimeStep * m_blockAssembler.getRhsM_SUPG();

            // unsteady term on the RHS
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < tarDim; ++s)
                m_rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * m_blockAssembler.getBlockM_SUPG() * m_blockAssembler.getSolution().middleRows(s * uDofs, uDofs);
        }

        // BC for rotation term
        if (m_blockAssembler.isRotation())
        {
            m_rhs.middleRows((tarDim - 2) * uDofs, uDofs) -= m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 1) * uDofs, uDofs);
            m_rhs.middleRows((tarDim - 1) * uDofs, uDofs) += m_blockAssembler.getOmega() * m_blockAssembler.getRhsM_SUPG().middleRows((tarDim - 2) * uDofs, uDofs);
        }
    }

    void fillCROSSWINDBase()
    {
        int numDofs = Base::numDofs();
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);
        for (int i = 0; i < uDofs; i++)
            for (int s = 0; s < tarDim; ++s)
                nonZerosPerColumnVector(i + s * uDofs) += m_blockAssembler.getBlockCrosswindPattern(s).col(i).nonZeros();

        gsSparseMatrix<T> blockCROSSWINDMatrix(numDofs, numDofs);
        blockCROSSWINDMatrix.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswindPattern(s), col); it; ++it)
                    blockCROSSWINDMatrix.insert(it.row() + s * uDofs, it.col() + s * uDofs) = 0.;


        m_baseMatrix += blockCROSSWINDMatrix;

        m_bMatrixReady = false;
        m_bRhsReady = false;
    }

    void fillCROSSWINDMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswind(s), col); it; ++it)
                    m_matrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    //---------------------------------------------------------------------------------------------------

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

    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        Fp = (1. / m_timeStepSize) * m_blockAssembler.getBlockMp();

     
        Base::fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
    }

    const T getTimeStepSize() const { return m_timeStepSize; }

protected:

    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;

    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;

    bool m_bMatrixReady;
    bool m_bRhsReady;

    T m_timeStepSize;

    // members from uwbINSAssemblerBase
    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;

}; // class uwbINSAssemblerUnsteady

} // namespace gismo
