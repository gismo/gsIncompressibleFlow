/** @file uwbINSAssemblerDecoupledPeriodic.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSAssemblerDecoupled1.h"
#include "uwbINSAssemblerPeriodicDecoupledBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerDecoupled1Periodic : public uwbINSAssemblerDecoupled1<T>, public uwbINSAssemblerPeriodicDecoupledBase<T>
{

public:
    typedef uwbINSAssemblerDecoupled1<T> Base;
    typedef uwbINSAssemblerPeriodicDecoupledBase<T> BasePeriodic;

public:

    uwbINSAssemblerDecoupled1Periodic(uwbINSSolverParams<T>& params) : Base(params), BasePeriodic(params, this->m_blockAssembler)
    {
        initMembers();
    }

    virtual ~uwbINSAssemblerDecoupled1Periodic()
    {
    }

protected:

    void initMembers()
    {
        int pdofs = this->getPDofsPer();
        int pshift = this->getPShiftPer();

        m_perMatrix1.resize(pshift, pshift);
        m_nonZerosPerColumnVector.setZero(pshift);

        m_perMatrix2.resize(pdofs, pdofs);
        m_perMatrix3.resize(pshift, pshift);

        m_perRhs1.setZero(pshift, 1);
        m_perRhs2.setZero(pdofs, 1);
        m_perRhs3.setZero(pshift, 1);
    }

    virtual void reinitMembers()
    {
        Base::reinitMembers();
        BasePeriodic::reinitMembers();
        initMembers();
    }

public:

    virtual void initialize()
    {
        Base::initialize();

        this->applyPeriodicConditionsOnPressure_into(m_matrix2, m_rhs2, m_perMatrix2, m_perRhs2, true);
        this->applyPeriodicConditionsOnVelocity_into(m_matrix3, m_rhs3, m_perMatrix3, m_perRhs3, true);
    }

public:
    virtual void update(const gsMatrix<T> & solVector) {
        m_perSolution = solVector;
        gsMatrix<T> nonPerSolution;
        this->per2nonper_into(m_perSolution, nonPerSolution);

        Base::update(nonPerSolution);

        if (!m_bNonZerosExact) {
            m_nonZerosPerColumnVector.setZero(this->getPShiftPer());
            for (index_t col = 0; col < Base::getPshift(); ++col) {
                m_nonZerosPerColumnVector(this->map(col)) += m_matrix1.col(col).nonZeros() + this->getPerUDofs();
            }
        }

        m_perMatrix1.resize(this->getPShiftPer(), this->getPShiftPer());
        m_perMatrix1.reserve(m_nonZerosPerColumnVector);
        this->applyPeriodicConditionsOnVelocity_into(m_matrix1, m_rhs1, m_perMatrix1, m_perRhs1, false);

        if (!m_bNonZerosExact) {
            for (index_t col = 0; col < this->getPShiftPer(); ++col)
                m_nonZerosPerColumnVector(col) = m_perMatrix1.col(col).nonZeros();
            m_bNonZerosExact = true;
        }
    }

    virtual void updateRhs2(const gsMatrix<T> & solVector1)
    {
        gsMatrix<T> nonPerSolVector1;
        this->per2nonperVelocity_into(solVector1, nonPerSolVector1);

        Base::updateRhs2(nonPerSolVector1);

        this->nonper2perPressure_into(m_rhs2, m_perRhs2);
    }

    virtual void updateRhs3(const gsMatrix<T> & solVector2)
    {
        gsMatrix<T> nonPerSolVector2;
        this->per2nonperPressure_into(solVector2, nonPerSolVector2);

        Base::updateRhs3(nonPerSolVector2);

        this->nonper2perVelocity_into(m_rhs3, m_perRhs3);
    }

    virtual void createSolVector_into(const gsMatrix<T> & solVector1, const gsMatrix<T> & solVector2, const gsMatrix<T> & solVector3, gsMatrix<T> & result) const
    {
        result.setZero(this->numDofsPer(), 1);

        result.middleRows(0, this->getPShiftPer()) = solVector1 + solVector3;

        result.middleRows(this->getPShiftPer(), this->getPDofsPer()) = m_perSolution.middleRows(this->getPShiftPer(), this->getPDofsPer()) + this->m_alpha_p * solVector2;
    }

    virtual const gsSparseMatrix<T> & matrix1() const
    {
        GISMO_ASSERT(this->m_bMatrix1Ready, "Matrix not ready, update() must be called first");
        return m_perMatrix1;
    }

    virtual const gsMatrix<T> & rhs1() const
    {
        GISMO_ASSERT(this->m_bRhs1Ready, "Rhs not ready, update() must be called first");
        return m_perRhs1;
    }

    virtual const gsSparseMatrix<T> & matrix2() const
    {
        GISMO_ASSERT(m_bInitialized, "Matrix not ready, initialize() must be called first");
        return m_perMatrix2;
    }

    virtual const gsMatrix<T> & rhs2() const
    {
        GISMO_ASSERT(this->m_bRhs2Ready, "Rhs not ready, updateRhs2() must be called first");

        return m_perRhs2;
    }

    virtual const gsSparseMatrix<T> & matrix3() const
    {
        GISMO_ASSERT(m_bInitialized, "Matrix not ready, initialize() must be called first");
        return m_perMatrix3;
    }

    virtual const gsMatrix<T> & rhs3() const
    {
        GISMO_ASSERT(this->m_bRhs3Ready, "Rhs not ready, updateRhs3() must be called first");

        return m_perRhs3;
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& perSolVector, int unk, bool relative = false) const
    {
        return BasePeriodic::constructSolution(perSolVector, unk, relative);
    }

    virtual gsField<T> constructSolutionCombined(const gsMatrix<T>& perSolVector, gsVector<size_t> relPatches) const
    {

        return BasePeriodic::constructSolutionCombined(perSolVector, relPatches);
    }

    // computes flow rate through a side of a patch
    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> perSolVector) const
    {
        return BasePeriodic::computeFlowRate(patch, side, perSolVector);
    }

    virtual T computeDimensionlessWallDistance(gsMatrix<T> perSolVector, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return BasePeriodic::computeDimensionlessWallDistance(perSolVector, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }


    virtual const gsMatrix<T> & getSolution() const { return m_perSolution; }
    virtual int numDofs() const { return this->numDofsPer(); }
    virtual int getUdofs() const { return this->getUDofsPer(); }
    virtual int getPdofs() const { return this->getPDofsPer(); }
    virtual int getPshift() const { return this->getPShiftPer(); }

protected:
    gsSparseMatrix<T> m_perMatrix1;
    gsSparseMatrix<T> m_perMatrix2;
    gsSparseMatrix<T> m_perMatrix3;

    gsMatrix<T> m_perRhs1;
    gsMatrix<T> m_perRhs2;
    gsMatrix<T> m_perRhs3;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_bInitialized;

    // members from uwbINSAssemblerDecoupled
    using Base::m_matrix1;
    using Base::m_matrix2;
    using Base::m_matrix3;
    using Base::m_rhs1;
    using Base::m_rhs2;
    using Base::m_rhs3;

    // members from uwbINSAssemblerPeriodicBase
    using BasePeriodic::Base::m_transformMatrix;
    using BasePeriodic::Base::m_invTransformMatrix;
    using BasePeriodic::Base::m_perDofs2_u;
    using BasePeriodic::Base::m_perDofs2_p;
    using BasePeriodic::Base::m_uMap;
    using BasePeriodic::Base::m_pMap;
    using BasePeriodic::Base::m_uInvMap;
    using BasePeriodic::Base::m_pInvMap;
    using BasePeriodic::Base::m_bNonZerosExact;
    using BasePeriodic::Base::m_nonZerosPerColumnVector;
    using BasePeriodic::Base::m_perSolution;

}; //class uwbINSAssemblerDecoupledPeriodic

} // namespace gismo
