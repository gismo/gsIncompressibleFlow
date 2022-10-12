/** @file uwbRANSAssembler.h

Author(s): E. Turnerova, H. Hornikova
*/

#pragma once

#include "uwbINSAssemblerUnsteady.h"
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

template<class T>
class uwbRANSAssembler : public uwbINSAssemblerUnsteady<T>
{

public:
    typedef uwbINSAssemblerUnsteady<T> Base;

public:
    uwbRANSAssembler(uwbINSSolverParams<T>& params, uwbTMSolverBase<T>* pTMsolver) : 
        Base(params), m_pTurbulenceSolver(pTMsolver)
    {
        m_blockAssembler.setRANS();
    }

    virtual ~uwbRANSAssembler()
    {
    }

public:

    virtual void updateAssembly()
    {
        Base::updateAssembly();

        m_blockAssembler.assembleRANSPart();

        fillNMatrix();

        if (m_blockAssembler.isRANScrosswind())
        {
            m_blockAssembler.assembleRANScrosswindPart();

            fillRANScrosswindMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsCrosswind();
        }
        if (m_blockAssembler.isRANStanhCSD())
        {
            m_blockAssembler.assembleRANSTanhCSDpart();
            fillRANStanhCSDmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsTanhCSD();
        }
        if (m_blockAssembler.isRANSisoAD())
        {
            m_blockAssembler.assembleRANSisoADpart();
            fillRANSisoADmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsIsoAD();
        }
        if (m_blockAssembler.isRANSad())
        {
            m_blockAssembler.assembleBlockRANSad();

            fillRANSadMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhs_RANS_AD();
        }
    }

    virtual void updatePicardAssembly()
    {
        Base::updatePicardAssembly();

        m_baseMatrix = m_NMatrix;
        m_baseRhs = m_NRhs;

        if (m_blockAssembler.isRANScrosswind())
        {
            m_blockAssembler.assembleRANScrosswindPart();

            fillRANScrosswindMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsCrosswind();
        }
        if (m_blockAssembler.isRANStanhCSD())
        {
            m_blockAssembler.assembleRANSTanhCSDpart();
            fillRANStanhCSDmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsTanhCSD();
        }
        if (m_blockAssembler.isRANSisoAD())
        {
            m_blockAssembler.assembleRANSisoADpart();
            fillRANSisoADmatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhsIsoAD();
        }
        if (m_blockAssembler.isRANSad())
        {
            m_blockAssembler.assembleBlockRANSad();

            fillRANSadMatrix();
            m_baseRhs.noalias() += m_blockAssembler.getRhs_RANS_AD();
        }
    }

protected:

    virtual void fillBase()
    {
        Base::fillBase();

        m_StokesMatrix = m_baseMatrix;
        m_StokesRhs = m_baseRhs;
    }
    
    virtual void fillNMatrix()
    {
        m_baseMatrix = m_StokesMatrix;
        m_baseRhs = m_StokesRhs;

        int udofs = m_blockAssembler.getUdofs();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_RANS(), col); it; ++it)
                for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                    m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < udofs; ++col)
            for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockEdiag_RANS(s), col); it; ++it)
                    m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

        //m_baseRhs += m_blockAssembler.getRhsA_RANS();
        m_baseRhs += m_blockAssembler.getRhsA_RANS() + m_blockAssembler.getRhsEdiag_RANS() + m_blockAssembler.getRhsEnondiag_RANS();

/*        T kKoef = -2 / 3;
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
            m_baseRhs.middleRows(s * udofs, udofs).noalias() += kKoef * m_blockAssembler.getBlockMinusBT(s) * m_pTurbulenceSolver->getSolutionK_full(m_blockAssembler.getMappers().back());
*/

        if (m_blockAssembler.isSUPG())
        {
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < udofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockA_RANS_SUPG(), col); it; ++it)
                    for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
                        m_baseMatrix.coeffRef(s*udofs + it.row(), s*udofs + it.col()) += it.value();

            m_baseRhs += m_blockAssembler.getRhsA_RANS_SUPG();

            //#pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            //for (index_t s = 0; s < m_blockAssembler.getTarDim(); ++s)
            //    m_baseRhs.middleRows(s * udofs, udofs).noalias() += kKoef * m_blockAssembler.getBlockB_SUPG(s) * m_pTurbulenceSolver->getSolutionK_full(m_blockAssembler.getMappers().back());

        }

        m_NMatrix = m_baseMatrix;
        m_NRhs = m_baseRhs;
    }

    void fillRANScrosswindMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockCrosswind(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANStanhCSDmatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockTanhCSD(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANSisoADmatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (index_t s = 0; s < tarDim; ++s)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockIsoAD(s), col); it; ++it)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

    void fillRANSadMatrix()
    {
        int uDofs = Base::getUdofs();
        int tarDim = Base::getTarDim();

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < uDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlock_RANS_AD(), col); it; ++it)
                for (index_t s = 0; s < tarDim; ++s)
                    m_baseMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += it.value();
    }

protected:
    uwbTMSolverBase<T>* m_pTurbulenceSolver;

    gsSparseMatrix<T> m_StokesMatrix;
    gsMatrix<T> m_StokesRhs;

    gsSparseMatrix<T> m_NMatrix;
    gsMatrix<T> m_NRhs;

    // members from uwbINSAssemblerBase
    using uwbINSAssemblerBase<T>::m_blockAssembler;
    using uwbINSAssemblerBase<T>::m_solution;
    using uwbINSAssemblerBase<T>::m_bInitialized;

    // members from uwbINSAssemblerUnsteady
    using Base::m_baseMatrix;
    using Base::m_baseRhs;
    

}; // class uwbRANSAssembler

} // namespace gismo
