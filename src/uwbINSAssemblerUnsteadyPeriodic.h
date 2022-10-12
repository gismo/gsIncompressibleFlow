/** @file uwbINSAssemblerUnsteadyPeriodic.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSAssemblerUnsteady.h"
#include "uwbINSAssemblerPeriodicCoupledBase.h"

// this class is almost identical to uwbINSAssemblerSteady (to be improved later)

namespace gismo
{

template<class T>
class uwbINSAssemblerUnsteadyPeriodic : public uwbINSAssemblerUnsteady<T>, public uwbINSAssemblerPeriodicCoupledBase<T>
{

public:
    typedef uwbINSAssemblerUnsteady<T> Base;
    typedef uwbINSAssemblerPeriodicCoupledBase<T> BasePeriodic;

public:

    uwbINSAssemblerUnsteadyPeriodic(uwbINSSolverParams<T>& params) : Base(params), BasePeriodic(params, this->m_blockAssembler)
    {
    }

    virtual ~uwbINSAssemblerUnsteadyPeriodic()
    {
    }

protected:

    virtual void reinitMembers()
    {
        Base::reinitMembers();
        BasePeriodic::reinitMembers();
    }

public:

    virtual void update(const gsMatrix<T> & solVector)
    {
        m_perSolution = solVector;
        gsMatrix<T> nonPerSolution;
        this->per2nonper_into(m_perSolution, nonPerSolution);

        Base::update(nonPerSolution);

        this->applyPeriodicConditions(m_matrix, m_rhs);
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        gsMatrix<T> nonPerSolution;
        this->per2nonper_into(solVector, nonPerSolution);

        Base::updatePicard(nonPerSolution);

        this->applyPeriodicConditions(m_matrix, m_rhs);
    }

    virtual const gsSparseMatrix<T> & matrix() const
    {
        GISMO_ASSERT(this->m_bMatrixReady, "Matrix not ready, update() must be called first");
        return m_perMatrix;
    }

    virtual const gsMatrix<T> & rhs() const
    {
        GISMO_ASSERT(this->m_bRhsReady, "Rhs not ready, update() must be called first");
        return m_perRhs;
    }

    virtual const gsSparseMatrix<T>& getVelocityMassMatrix()
    {
        return BasePeriodic::getVelocityMassMatrix();
    }

    virtual const gsSparseMatrix<T>& getPressureMassMatrix()
    {
        return BasePeriodic::getPressureMassMatrix();
    }

    virtual void fillPCDblocks(gsSparseMatrix<T>& perAp, gsSparseMatrix<T>& perFp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        gsSparseMatrix<T> Ap, Fp;
        Base::fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);

        this->applyPeriodicConditionsOnPressure_into(Ap, perAp);
        this->applyPeriodicConditionsOnPressure_into(Fp, perFp);
    }

    virtual void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        uwbINSAssemblerPeriodicCoupledBase<T>::fillStokesSystem_into(stokesMatrix, stokesRhs);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& perSolVector, int unk, bool relative = false) const
    {
        return uwbINSAssemblerPeriodicBase<T>::constructSolution(perSolVector, unk, relative);
    }

    virtual gsField<T> constructSolution(const gsMatrix<T>& perSolVector, int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        return uwbINSAssemblerPeriodicBase<T>::constructSolution(perSolVector, unk, patchNumbers, sides, relative);
    }

    virtual gsField<T> constructSolutionCombined(const gsMatrix<T>& perSolVector, gsVector<size_t> relPatches) const
    {

        return uwbINSAssemblerPeriodicBase<T>::constructSolutionCombined(perSolVector, relPatches);
    }

    // computes flow rate through a side of a patch
    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> perSolVector) const
    {
        return uwbINSAssemblerPeriodicBase<T>::computeFlowRate(patch, side, perSolVector);
    }

    virtual T computeDimensionlessWallDistance(gsMatrix<T> perSolVector, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return uwbINSAssemblerPeriodicBase<T>::computeDimensionlessWallDistance(perSolVector, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

public:
    virtual const gsMatrix<T> & getSolution() const { return m_perSolution; }
    virtual int numDofs() const { return this->numDofsPer(); }
    virtual int getUdofs() const { return this->getUDofsPer(); }
    virtual int getPdofs() const { return this->getPDofsPer(); }
    virtual int getPshift() const { return this->getPShiftPer(); }

protected:

    // members from uwbINSAssemblerPeriodicCoupledBase
    using BasePeriodic::m_perMatrix;
    using BasePeriodic::m_perRhs;

    // members from uwbINSAssemblerPeriodicBase
    using BasePeriodic::Base::m_perSolution;

    // members from uwbINSAssemblerSteady
    using Base::m_matrix;
    using Base::m_rhs;
};

} // namespace gismo


