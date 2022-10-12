/** @file uwbINSAssemblerDecoupled.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSAssemblerBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerDecoupled1 : public uwbINSAssemblerBase<T>
{

public:
    typedef uwbINSAssemblerBase<T> Base;

public:
    uwbINSAssemblerDecoupled1(uwbINSSolverParams<T>& params) : Base(params)
    {
        m_alpha_u = params.settings().get(constantsINS::alpha_u);
        m_alpha_p = params.settings().get(constantsINS::alpha_p);

        initMembers();
    }

    virtual ~uwbINSAssemblerDecoupled1()
    {
    }

protected:

    void initMembers()
    {
        short_t tarDim = Base::getTarDim();
        int udofs = Base::getUdofs();
        int pdofs = Base::getPdofs();

        m_baseBlock1.resize(tarDim * udofs, tarDim * udofs);

        m_matrix1.resize(tarDim * udofs, tarDim * udofs);
        m_matrix2.resize(pdofs, pdofs);
        m_matrix3.resize(tarDim * udofs, tarDim * udofs);
        m_rhs1.setZero(tarDim * udofs, 1);
        m_rhs2.setZero(pdofs, 1);
        m_rhs3.setZero(tarDim * udofs, 1);

        m_solution.setZero(Base::numDofs(), 1);

        m_bInitialized = false;
        m_bMatrix1Ready = false;
        m_bRhs1Ready = false;
        m_bRhs2Ready = false;
        m_bRhs3Ready = false;
    }

    virtual void reinitMembers() { initMembers(); }

public:

    virtual void initialize()
    {
        m_blockAssembler.assembleAllLinearBlocks();

        m_relax = (1 - m_alpha_u) / m_alpha_u;

        fillBase();

        m_matrix2 = - m_blockAssembler.getBlockAp();

        fillBlockOnDiagonal_into(m_blockAssembler.getBlockM(), m_matrix3);

        m_bInitialized = true;
    }

    virtual void update(const gsMatrix<T> & solVector)
    {
        Base::update(solVector);

        m_bRhs2Ready = false;
        m_bRhs3Ready = false;
    }

protected:

    void fillBase()
    {
        m_baseBlock1 = m_blockAssembler.getBlockA() + m_relax * m_blockAssembler.getBlockM();

        m_bMatrix1Ready = false;
        m_bRhs1Ready = false;
        m_bRhs2Ready = false;
        m_bRhs3Ready = false;
    }

    void fillMatrix()
    {
        fillBlockOnDiagonal_into(m_baseBlock1 + m_blockAssembler.getBlockN(), m_matrix1, m_blockAssembler.isRotation());

        m_bMatrix1Ready = true;
    }

    virtual void fillRhs()
    {
        #pragma omp parallel for
        for (index_t s = 0; s < Base::getTarDim(); ++s) {
            m_rhs1.middleRows(s * Base::getUdofs(), Base::getUdofs()) =
                m_relax * m_blockAssembler.getBlockM() * m_solution.middleRows(s * Base::getUdofs(), Base::getUdofs())
                + m_blockAssembler.getRhsA().middleRows(s * Base::getUdofs(), Base::getUdofs())
                + m_blockAssembler.getRhsN().middleRows(s * Base::getUdofs(), Base::getUdofs())
                + m_blockAssembler.getRhsF().middleRows(s * Base::getUdofs(), Base::getUdofs())
                + m_blockAssembler.getRhsB().middleRows(s * Base::getUdofs(), Base::getUdofs())
                - m_blockAssembler.getBlockMinusBT(s) * m_solution.middleRows(Base::getPshift(), Base::getPdofs());
        }

        if (m_blockAssembler.isRotation())
        {
            m_rhs1.middleRows((Base::getTarDim() - 2) * Base::getUdofs(), Base::getUdofs()) -= m_blockAssembler.getOmega() * m_blockAssembler.getRhsM().middleRows((Base::getTarDim() - 1) * Base::getUdofs(), Base::getUdofs());
            m_rhs1.middleRows((Base::getTarDim() - 1) * Base::getUdofs(), Base::getUdofs()) += m_blockAssembler.getOmega() * m_blockAssembler.getRhsM().middleRows((Base::getTarDim() - 2) * Base::getUdofs(), Base::getUdofs());
        }

        m_bRhs1Ready = true;
    }

public:
    virtual void updateRhs2(const gsMatrix<T> & solVector1)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_rhs2.setZero(Base::getPdofs(), 1);
        for (index_t s = 0; s < Base::getTarDim(); ++s) {
            m_rhs2.noalias() += m_relax * m_blockAssembler.getBlockB(s) * solVector1.middleRows(s * Base::getUdofs(), Base::getUdofs());
        }
        m_rhs2.noalias() -=  m_relax *m_blockAssembler.getRhsB().middleRows(Base::getPshift(), Base::getPdofs());

        m_bRhs2Ready = true;
    }

    virtual void updateRhs3(const gsMatrix<T> & solVector2)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        #pragma omp parallel for
        for (index_t s = 0; s < Base::getTarDim(); ++s) {
            //m_rhs3.col(s).noalias() = - m_blockAssembler.getBlockMinusBT(s) * solVector2;
            m_rhs3.middleRows(s * Base::getUdofs(), Base::getUdofs()) = - (1 / m_relax) * m_blockAssembler.getBlockCT(s) * solVector2;
        }

        m_bRhs3Ready = true;
    }

    virtual void createSolVector_into(const gsMatrix<T> & solVector1, const gsMatrix<T> & solVector2, const gsMatrix<T> & solVector3, gsMatrix<T> & result) const
    {
        result.setZero(Base::numDofs(), 1);

        result.middleRows(0, Base::getTarDim() * Base::getUdofs()) = solVector1 + solVector3;

        result.middleRows(Base::getPshift(), Base::getPdofs()) = m_solution.middleRows(Base::getPshift(), Base::getPdofs()) + m_alpha_p * solVector2;
    }

    void changeRelaxParametres(const T alpha_u, const T alpha_p)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_alpha_u = alpha_u;
        m_alpha_p = alpha_p;
        m_relax = (1 - alpha_u) / alpha_u;

        fillBase();
    }

protected:

    void fillBlockOnDiagonal_into(const gsSparseMatrix<T> & block, gsSparseMatrix<T> & matrix, const bool addRotation = false) const
    {
        short_t tarDim = Base::getTarDim();
        int udofs = Base::getUdofs();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(tarDim * udofs);
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < udofs; i++)
                nonZerosPerColumnVector(i + s * udofs) = block.col(i).nonZeros();

        if (addRotation)
            for (int s = tarDim - 2; s < tarDim; ++s)
                for (int i = 0; i < udofs; i++)
                    nonZerosPerColumnVector(i + s * udofs) += m_blockAssembler.getBlockM().col(i).nonZeros();

        matrix.resize(tarDim * udofs, tarDim * udofs);
        matrix.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(block, col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    matrix.insert(it.row() + s * udofs, it.col() + s * udofs) = it.value();

        if (addRotation)
        {
            GISMO_ENSURE(tarDim == 2 || tarDim == 3, "Rotation is implemented only for 2D and 3D.");

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < udofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM(), col); it; ++it) {
                    matrix.insert(it.row() + ((tarDim - 2) * udofs), it.col() + ((tarDim - 1) * udofs)) = - m_blockAssembler.getOmega() * it.value();
                    matrix.insert(it.row() + ((tarDim - 1) * udofs), it.col() + ((tarDim - 2) * udofs)) = m_blockAssembler.getOmega() * it.value();
                }
        }

        matrix.makeCompressed();
    }

public:

    virtual const gsSparseMatrix<T> & matrix1() const
    {
        GISMO_ASSERT(m_bMatrix1Ready, "Matrix not ready, update() must be called first");
        return m_matrix1;
    }

    virtual const gsMatrix<T> & rhs1() const
    {
        GISMO_ASSERT(m_bRhs1Ready, "Rhs not ready, update() must be called first");
        return m_rhs1;
    }

    virtual const gsSparseMatrix<T> & matrix2() const
    {
        GISMO_ASSERT(m_bInitialized, "Matrix not ready, initialize() must be called first");
        return m_matrix2;
    }

    virtual const gsMatrix<T> & rhs2() const
    {
        GISMO_ASSERT(m_bRhs2Ready, "Rhs not ready, updateRhs2() must be called first");

        return m_rhs2;
    }

    virtual const gsSparseMatrix<T> & matrix3() const
    {
        GISMO_ASSERT(m_bInitialized, "Matrix not ready, initialize() must be called first");
        return m_matrix3;
    }

    virtual const gsMatrix<T> & rhs3() const
    {
        GISMO_ASSERT(m_bRhs3Ready, "Rhs not ready, updateRhs3() must be called first");

        return m_rhs3;
    }

protected:

    gsSparseMatrix<T> m_baseBlock1;

    gsSparseMatrix<T> m_matrix1;
    gsSparseMatrix<T> m_matrix2;
    gsSparseMatrix<T> m_matrix3;

    gsMatrix<T> m_rhs1;
    gsMatrix<T> m_rhs2;
    gsMatrix<T> m_rhs3;

    bool m_bMatrix1Ready;
    bool m_bRhs1Ready;
    bool m_bRhs2Ready;
    bool m_bRhs3Ready;

    T m_alpha_u;
    T m_alpha_p;
    T m_relax;

    // members from uwbINSAssemblerBase
    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;
};

} // namespace gismo


