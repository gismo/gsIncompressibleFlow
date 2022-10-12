/** @file uwbINSAssemblerDecoupled2.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSAssemblerDecoupled1.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerDecoupled2 : public uwbINSAssemblerDecoupled1<T>
{

public:
    typedef uwbINSAssemblerDecoupled1<T> Base;

public:
    uwbINSAssemblerDecoupled2(uwbINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSAssemblerDecoupled2()
    {
    }

protected:

    void initMembers()
    {
        m_matrix3.resize(Base::getPdofs(), Base::getPdofs());
        m_rhs3.setZero(Base::getPdofs(), 1);

        m_solution2.setZero(Base::getPdofs(), 1);

        m_bRhs3Ready = false;
    }

    virtual void reinitMembers() 
    {
        Base::reinitMembers();
        initMembers();
    }

public:

    virtual void initialize()
    {
        Base::initialize();

        m_blockAssembler.assemblePressureMassMatrix();

        m_matrix3 = m_blockAssembler.getBlockMp();
    }

    virtual void update(const gsMatrix<T> & solVector, const gsMatrix<T> & solVector2)
    {
        m_solution2 = solVector2;

        Base::update(solVector);

        m_bRhs2Ready = false;
        m_bRhs3Ready = false;
    }

protected:

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
                - m_blockAssembler.getBlockMinusBT(s) * m_solution.middleRows(Base::getPshift(), Base::getPdofs())
                - m_blockAssembler.getBlockMinusBT(s) * m_solution2;
        }

        if (m_blockAssembler.isRotation())
        {
            m_rhs1.middleRows((Base::getTarDim() - 2) * Base::getUdofs(), Base::getUdofs()) -= m_blockAssembler.getOmega() * m_blockAssembler.getRhsM().middleRows((Base::getTarDim() - 1) * Base::getUdofs(), Base::getUdofs());
            m_rhs1.middleRows((Base::getTarDim() - 1) * Base::getUdofs(), Base::getUdofs()) += m_blockAssembler.getOmega() * m_blockAssembler.getRhsM().middleRows((Base::getTarDim() - 2) * Base::getUdofs(), Base::getUdofs());
        }

        m_bRhs1Ready = true;
    }

public:

    virtual void updateRhs3(const gsMatrix<T> & solVector1, const gsMatrix<T> & solVector2)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_rhs3.setZero(Base::getPdofs(), 1);
        for (index_t s = 0; s < Base::getTarDim(); ++s) {
            m_rhs3.noalias() -= m_blockAssembler.getViscosity() * m_blockAssembler.getBlockB(s) * solVector1.middleRows(s * Base::getUdofs(), Base::getUdofs());
        }
        m_rhs3.noalias() += m_blockAssembler.getBlockMp() * (m_solution.middleRows(Base::getPshift(), Base::getPdofs()) + solVector2)
            + m_blockAssembler.getViscosity() * m_blockAssembler.getRhsB().middleRows(Base::getPshift(), Base::getPdofs());

        m_bRhs3Ready = true;
    }

    virtual void createSolVector_into(const gsMatrix<T> & solVector1, const gsMatrix<T> & solVector2, const gsMatrix<T> & solVector3, gsMatrix<T> & result) const
    {
        result.setZero(Base::numDofs(), 1);

        result.middleRows(0, Base::getTarDim() * Base::getUdofs()) = solVector1;
        result.middleRows(Base::getPshift(), Base::getPdofs()) = solVector3;
    }

public:

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

    gsSparseMatrix<T> m_matrix3;

    gsMatrix<T> m_rhs3;

    gsMatrix<T> m_solution2;

    bool m_bRhs3Ready;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_solution;
    using uwbINSAssemblerBase<T>::m_bInitialized;

    // members from uwbINSAssemblerDecoupled1
    using Base::m_rhs1;
    using Base::m_bRhs1Ready;
    using Base::m_bRhs2Ready;
    using Base::m_relax;
   
};

} // namespace gismo


