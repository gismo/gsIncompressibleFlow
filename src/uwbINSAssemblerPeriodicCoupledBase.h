/** @file uwbINSAssemblerPeriodicCoupledBase.h

Author(s): H. Hornikova
*/

#pragma once
#include "uwbINSAssemblerPeriodicBase.h"

namespace gismo
{

template<class T>
class uwbINSAssemblerPeriodicCoupledBase : public uwbINSAssemblerPeriodicBase<T>
{
public:
    typedef uwbINSAssemblerPeriodicBase<T> Base;

public:
    uwbINSAssemblerPeriodicCoupledBase(uwbINSSolverParams<T>& params, uwbINSBlockAssembler<T>& blockAssembler) :
        Base(params, blockAssembler)
    { 
        initMembers();
    }

protected:
    void initMembers()
    {
        int dofsPer = numDofsPer();

        m_perMatrix.resize(dofsPer, dofsPer);
        m_perRhs.setZero(dofsPer, 1);

        m_bMassPerReady = false;

        m_nonZerosPerColumnVector.setZero(dofsPer);
    }

    virtual void reinitMembers() 
    { 
        Base::reinitMembers();
        initMembers(); 
    }

    void applyPeriodicConditions_into(const gsSparseMatrix<T>& matrix, gsSparseMatrix<T>& resultMatrix, bool reserveMemory = true) const
    {
        int numDofs = m_blockAssemblerRef.numDofs();
        int udofs = m_blockAssemblerRef.getUdofs();
        int pshift = m_blockAssemblerRef.getPshift();
        short_t tarDim = m_blockAssemblerRef.getTarDim();

        if (reserveMemory)
        {
            resultMatrix.resize(numDofsPer(), numDofsPer());

            gsVector<int> nonZerosPerColumnVector;
            nonZerosPerColumnVector.setZero(numDofsPer());

            //compute maximum non zero element number
            //from original matrix
            for (int i = 0; i < numDofs; i++)
                if (!isEliminated(i))
                    nonZerosPerColumnVector(map(i)) = matrix.col(i).nonZeros();

            //from col additions
            //U
            for (int j = 0; j < getPerUDofs(); j++)
                for (int s = 0; s < tarDim; s++)
                    for (int t = 0; t < tarDim; t++)
                        nonZerosPerColumnVector(map(m_perDofs2_u[j] + t * udofs)) += matrix.col(m_perDofs2_u[j] + s * udofs).nonZeros();
            //P
            for (int k = 0; k < getPerPDofs(); k++)
                nonZerosPerColumnVector(map(m_perDofs2_p[k] + pshift)) += matrix.col(m_perDofs2_p[k] + pshift).nonZeros();

            //from row additions
            nonZerosPerColumnVector.topRows(getPShiftPer()) += gsVector<int>::Constant(getPShiftPer(), 1, tarDim * getPerUDofs() + getPerPDofs());
            nonZerosPerColumnVector.bottomRows(this->getPDofsPer()) += gsVector<int>::Constant(this->getPDofsPer(), 1, tarDim * getPerPDofs());
        
            resultMatrix.reserve(nonZerosPerColumnVector);
        }

        // MATRIX TRANSFORMATION
        // map non-eliminated elements
        for (index_t col = 0; col < numDofs; ++col)
            if (!(isEliminated(col)))
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, col); it; ++it)
                    if (!(isEliminated(it.row())))
                        resultMatrix.insert(map(it.row()), map(it.col())) = it.value();

        //col additions
        //U
        for (int i = 0; i < getPerUDofs(); i++)
        {
            for (int s = 0; s < tarDim; s++)
            {
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, m_perDofs2_u[i] + s * udofs); it; ++it)
                {
                    if (!(isEliminated(it.row())))
                    {
                        for (int t = 0; t < tarDim; t++)
                            resultMatrix.coeffRef(map(it.row()), m_uMap(m_perDofs2_u[i]) + t * getUDofsPer()) += m_transformMatrix(t, s) * it.value();
                    }
                    //cross elements
                    else
                    {
                        if (it.row() < pshift)
                        {
                            for (int t = 0; t < tarDim; t++)
                                for (int u = 0; u < tarDim; u++)
                                    resultMatrix.coeffRef(m_uMap(it.row() % udofs) + u * getUDofsPer(), m_uMap(m_perDofs2_u[i]) + t * getUDofsPer()) += m_transformMatrix(u, it.row() / udofs) * m_transformMatrix(t, s) * it.value();
                        }
                        else
                        {
                            for (int t = 0; t < tarDim; t++)
                                resultMatrix.coeffRef(map(it.row()), m_uMap(m_perDofs2_u[i]) + t * getUDofsPer()) += m_transformMatrix(t, s) * it.value();
                        }
                    }
                }
            }
        }

        //P
        for (int i = 0; i < getPerPDofs(); i++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(matrix, m_perDofs2_p[i] + pshift); it; ++it)
            {
                if (!(isEliminated(it.row())))
                {
                    resultMatrix.coeffRef(map(it.row()), m_pMap(m_perDofs2_p[i]) + getPShiftPer()) += it.value();
                }
                //cross elements
                else
                {
                    if (it.row() < pshift)
                    {
                        for (int u = 0; u < tarDim; u++)
                            resultMatrix.coeffRef(m_uMap(it.row() % udofs) + u * getUDofsPer(), m_pMap(m_perDofs2_p[i]) + getPShiftPer()) += m_transformMatrix(u, it.row() / udofs) * it.value();
                    }
                    else
                        resultMatrix.coeffRef(map(it.row()), m_pMap(m_perDofs2_p[i]) + getPShiftPer()) += it.value();
                }
            }
        }

        //row additions
        for (index_t col = 0; col < numDofs; ++col)
        {
            if (!(isEliminated(col)))
            {
                for (typename gsSparseMatrix<T>::InnerIterator it(matrix, col); it; ++it)
                {
                    if (isEliminated(it.row()))
                    {
                        if (it.row() < pshift)
                        {
                            for (int t = 0; t < tarDim; t++)
                                resultMatrix.coeffRef(m_uMap(it.row() % udofs) + t * getUDofsPer(), map(col)) += m_transformMatrix(t, it.row() / udofs) * it.value();
                        }
                        else
                            resultMatrix.coeffRef(map(it.row()), map(col)) += it.value();
                    }
                }
            }
        }

        resultMatrix.makeCompressed();
    }

    void applyPeriodicConditions(const gsSparseMatrix<T>& matrix, const gsMatrix<T>& rhs)
    {
        int numDofs = m_blockAssemblerRef.numDofs();
        int udofs = m_blockAssemblerRef.getUdofs();
        int pshift = m_blockAssemblerRef.getPshift();
        short_t tarDim = m_blockAssemblerRef.getTarDim();

        m_perMatrix.resize(numDofsPer(), numDofsPer());
        m_perRhs.setZero(numDofsPer(), 1);

        if (!m_bNonZerosExact) {
            //compute maximum non zero element number
            //from original matrix
            for (int i = 0; i < numDofs; i++)
                if (!isEliminated(i))
                    m_nonZerosPerColumnVector(map(i)) = matrix.col(i).nonZeros();

            //from col additions
            //U
            for (int j = 0; j < getPerUDofs(); j++)
                for (int s = 0; s < tarDim; s++)
                    for (int t = 0; t < tarDim; t++)
                        m_nonZerosPerColumnVector(map(m_perDofs2_u[j] + t * udofs)) += matrix.col(m_perDofs2_u[j] + s * udofs).nonZeros();
            //P
            for (int k = 0; k < getPerPDofs(); k++)
                m_nonZerosPerColumnVector(map(m_perDofs2_p[k] + pshift)) += matrix.col(m_perDofs2_p[k] + pshift).nonZeros();

            //from row additions
            m_nonZerosPerColumnVector.topRows(getPShiftPer()) += gsVector<int>::Constant(getPShiftPer(), 1, tarDim * getPerUDofs() + getPerPDofs());
            m_nonZerosPerColumnVector.bottomRows(this->getPDofsPer()) += gsVector<int>::Constant(this->getPDofsPer(), 1, tarDim * getPerPDofs());
        }

        m_perMatrix.reserve(m_nonZerosPerColumnVector);

        applyPeriodicConditions_into(matrix, m_perMatrix, false);

        if (!m_bNonZerosExact)
        {
            for (index_t col = 0; col < numDofsPer(); ++col)
            {
                m_nonZerosPerColumnVector(col) = m_perMatrix.col(col).nonZeros();
            }
            m_bNonZerosExact = true;
        }

        this->nonper2per_into(rhs, m_perRhs);
    }

public:

    void applyPeriodicConditionsOnMass()
    {
        const gsSparseMatrix<T>& velM = m_blockAssemblerRef.getBlockM();
        const gsSparseMatrix<T>& presM = m_blockAssemblerRef.getBlockMp();

        GISMO_ASSERT((velM.nonZeros() > 0 && presM.nonZeros() > 0), "One or both mass matrices were not assembled yet.");

        int numDofs = m_blockAssemblerRef.numDofs();
        int udofs = m_blockAssemblerRef.getUdofs();
        int pdofs = m_blockAssemblerRef.getPdofs();
        int pshift = m_blockAssemblerRef.getPshift();
        short_t tarDim = m_blockAssemblerRef.getTarDim();

        // non-zeros in velocity mass matrix
        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (int s = 0; s < tarDim; ++s)
            for (int i = 0; i < udofs; i++)
                nonZerosPerColumnVector(i + s * udofs) = velM.col(i).nonZeros();

        // non-zeros in pressure mass matrix
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (int i = 0; i < pdofs; i++)
            nonZerosPerColumnVector(pshift + i) = presM.col(i).nonZeros();

        // matrix with tarDim * velocity mass matrix and pressure mass matrix on its block-diagonal
        gsSparseMatrix<T> massMatrices(numDofs, numDofs);
        massMatrices.reserve(nonZerosPerColumnVector);

        // fill velocity mass matrix block on the diagonal
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(velM, col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    massMatrices.insert(it.row() + s * udofs, it.col() + s * udofs) = it.value();

        // fill pressure mass matrix block on the diagonal
        #pragma omp parallel for num_threads(m_blockAssemblerRef.getNumThreads())
        for (index_t col = 0; col < pdofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(presM, col); it; ++it)
                massMatrices.insert(pshift + it.row(), pshift + it.col()) = it.value();

        gsSparseMatrix<T> perMassMatrices;
        applyPeriodicConditions_into(massMatrices, perMassMatrices);

        m_perVelMass = perMassMatrices.block(0, 0, tarDim*getUDofsPer(), tarDim*getUDofsPer());
        m_perVelMass_xComp = perMassMatrices.block(0, 0, getUDofsPer(), getUDofsPer());
        m_perPresMass = perMassMatrices.block(getPShiftPer(), getPShiftPer(), getPDofsPer(), getPDofsPer());

        m_bMassPerReady = true;
    }

    virtual void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        gsSparseMatrix<T> stokesMatrixNonper;
        gsMatrix<T> stokesRhsNonper;

        m_blockAssemblerRef.fillStokesSystem_into(stokesMatrixNonper, stokesRhsNonper);

        applyPeriodicConditions_into(stokesMatrixNonper, stokesMatrix);
        this->nonper2per_into(stokesRhsNonper, stokesRhs);
    }

    virtual const gsSparseMatrix<T>& getVelocityMassMatrix()
    {
        if (!m_bMassPerReady)
            applyPeriodicConditionsOnMass();

        return m_perVelMass; 
    }

    const gsSparseMatrix<T>& getVelocityMassMatrix_xComp()
    {
        if (!m_bMassPerReady)
            applyPeriodicConditionsOnMass();

        return m_perVelMass_xComp;
    }

    virtual const gsSparseMatrix<T>& getPressureMassMatrix()
    {
        if (!m_bMassPerReady)
            applyPeriodicConditionsOnMass();

        return m_perPresMass;
    }

protected:
    gsSparseMatrix<T> m_perMatrix;
    gsMatrix<T> m_perRhs;

    gsSparseMatrix<T> m_perVelMass, m_perVelMass_xComp, m_perPresMass;
    bool m_bMassPerReady;

    // members from uwbINSAssemblerPeriodicBase
    using Base::m_transformMatrix;
    using Base::m_perDofs2_u;
    using Base::m_perDofs2_p;
    using Base::m_uMap;
    using Base::m_pMap;
    using Base::m_bNonZerosExact;
    using Base::m_nonZerosPerColumnVector;
    using Base::m_blockAssemblerRef;

    // member functions from uwbINSAssemblerPeriodicBase
    using Base::map;
    using Base::isEliminated;
    using Base::numDofsPer;
    using Base::getUDofsPer;
    using Base::getPDofsPer;
    using Base::getPShiftPer;
    using Base::getPerUDofs;
    using Base::getPerPDofs;


}; // class uwbINSAssemblerPeriodicCoupledBase

} // namespace gismo