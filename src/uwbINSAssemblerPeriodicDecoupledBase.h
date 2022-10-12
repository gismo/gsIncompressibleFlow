/** @file uwbINSAssemblerPeriodicDecoupledBase.h

Author(s): H. Hornikova
*/

#pragma once
#include "uwbINSAssemblerPeriodicBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerPeriodicDecoupledBase : public uwbINSAssemblerPeriodicBase<T>
{
public:
    typedef uwbINSAssemblerPeriodicBase<T> Base;

public:
    uwbINSAssemblerPeriodicDecoupledBase(uwbINSSolverParams<T>& params, uwbINSBlockAssembler<T>& blockAssembler) :
        Base(params, blockAssembler)
    {  }

protected:

    void nonper2perPressure_into(const gsMatrix<T> & nonperVector, gsMatrix<T> & perVector) const
    {
        perVector.setZero(this->getPDofsPer(), 1);

        for (index_t row = 0; row < nonperVector.rows(); ++row)
            perVector(m_pMap(row), 0) += nonperVector(row, 0);
    }

    void per2nonperPressure_into(const gsMatrix<T> & perVector, gsMatrix<T> & nonperVector) const
    {
        nonperVector.setZero(m_blockAssemblerRef.getPdofs(), 1);

        for (index_t row = 0; row < this->getPDofsPer(); ++row)
            nonperVector(m_pInvMap(row), 0) += perVector(row, 0);

        for (int i = 0; i < this->getPerPDofs(); i++)
            nonperVector(m_perDofs2_p[i], 0) = perVector(m_pMap(m_perDofs2_p[i]), 0);
    }

    void nonper2perVelocity_into(const gsMatrix<T> & nonperVector, gsMatrix<T> & perVector) const
    {
        perVector.setZero(this->getPShiftPer(), 1);

        short_t tarDim = m_blockAssemblerRef.getTarDim();
        int udofs = m_blockAssemblerRef.getUdofs();

        for (int s = 0; s < tarDim; s++) {
            for (index_t row = 0; row < udofs; ++row)
                if (!(isEliminated(row)))
                    perVector(map(row) + s * getUDofsPer(), 0) += nonperVector(row + s * udofs, 0);
                else {
                    for (int t = 0; t < tarDim; t++)
                        perVector(m_uMap(row) + t * getUDofsPer(), 0) += m_transformMatrix(t, s) * nonperVector(row + s * udofs, 0);
                }
        }
    }

    void per2nonperVelocity_into(const gsMatrix<T> & perVector, gsMatrix<T> & nonperVector) const
    {
        nonperVector.setZero(m_blockAssemblerRef.getPshift(), 1);

        short_t tarDim = m_blockAssemblerRef.getTarDim();
        int udofs = m_blockAssemblerRef.getUdofs();

        for (int s = 0; s < tarDim; s++)
            for (index_t row = 0; row < getUDofsPer(); ++row)
                nonperVector(m_uInvMap(row) + s * udofs, 0) += perVector(row + s * getUDofsPer(), 0);

        for (int i = 0; i < this->getPerUDofs(); i++) {
            gsMatrix<T> u1(tarDim, 1);
            for (int s = 0; s < tarDim; s++) {
                u1(s, 0) = perVector(m_uMap(m_perDofs2_u[i]) + s * getUDofsPer(), 0);
            }
            gsMatrix<T> u2(tarDim, 1);
            u2 = m_invTransformMatrix * u1;
            for (int s = 0; s < tarDim; s++) {
                nonperVector(m_perDofs2_u[i] + s * udofs, 0) = u2(s, 0);
            }
        }
    }

    void applyPeriodicConditionsOnVelocity_into(const gsSparseMatrix<T> & matrix, const gsMatrix<T> & rhs, gsSparseMatrix<T> & resultMatrix, gsMatrix<T> & resultRhs, bool reserveMemory = true) {

        nonper2perVelocity_into(rhs, resultRhs);

        Base::applyPeriodicConditionsOnVelocity_into(matrix, resultMatrix, reserveMemory);
    }

    void applyPeriodicConditionsOnPressure_into(const gsSparseMatrix<T> & matrix, const gsMatrix<T> & rhs, gsSparseMatrix<T> & resultMatrix, gsMatrix<T> & resultRhs, bool reserveMemory = true) {

        nonper2perPressure_into(rhs, resultRhs);

        Base::applyPeriodicConditionsOnPressure_into(matrix, resultMatrix, reserveMemory);
    }

protected:

    // members from uwbINSAssemblerPeriodicBase
    using Base::m_transformMatrix;
    using Base::m_invTransformMatrix;
    using Base::m_perDofs2_u;
    using Base::m_perDofs2_p;
    using Base::m_uMap;
    using Base::m_pMap;
    using Base::m_uInvMap;
    using Base::m_pInvMap;
    using Base::m_bNonZerosExact;
    using Base::m_nonZerosPerColumnVector;
    using Base::m_blockAssemblerRef;

    // member functions from uwbINSAssemblerPeriodicBase
    using Base::map;
    using Base::isEliminated;
    using Base::getUDofsPer;

}; // class uwbINSAssemblerPeriodicDecoupledBase

} // namespace gismo