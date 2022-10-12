/** @file uwbTMAssemblerKOmega_AFC.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMAssemblerKOmega.h"

namespace gismo
{

template<class T>
class uwbTMAssemblerKOmega_AFC : public uwbTMAssemblerKOmega<T>
{

public:
    typedef uwbTMAssemblerKOmega<T> Base;

public:
    uwbTMAssemblerKOmega_AFC(uwbINSSolverParams<T>& params) : // assuming dirichlet::none from solver
        Base(params)
    {
        GISMO_ASSERT(params.settings().get(constantsINS::timeDerTerm) == true, "AFC method implemented only for unsteady case.");

        m_bAFC_HO = params.settings().get(constantsINS::TMafcHO);

        uwbINSSolverParams<T> params1(params);
        params1.getAssemblerOptions().dirStrategy = dirichlet::elimination;
        m_pBlockAssemblerBC = new uwbTMBlockAssembler<T>(params1, 2);

    }

    virtual ~uwbTMAssemblerKOmega_AFC()
    { 
        if (m_pBlockAssemblerBC)
        {
            delete m_pBlockAssemblerBC;
            m_pBlockAssemblerBC = NULL;
        }
    }

    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        m_solPicard = solVector;
        uwbTMAssemblerBase<T>::updatePicard(m_solPicard);
    }

public:
    //===================================================== fillExplicit ================================================================
    virtual void fillExplicitPartMatrix()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        m_nExplicitPartMatrix.setZero();
        
        const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN();

        // block N (AFC modified)
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
        {
            for (index_t s = 0; s < numVar; ++s)
            {
                m_nExplicitPartMatrix.coeffRef(col, col + s*varDofs) += blockNref.at(col, col);
            }

            for (typename gsSparseMatrix<T>::InnerIterator it(blockNref, col); it; ++it)
            {
                int row = it.row();

                if (col > row)
                {
                    T temp1 = it.value(); // N(row, col)
                    T temp2 = blockNref.at(col, row);
                    T delta = math::max(0.0, math::max(temp1, temp2));

                    for (index_t s = 0; s < numVar; ++s)
                    {
                        m_nExplicitPartMatrix.coeffRef(row, col + s*varDofs) += temp1 - delta;
                        m_nExplicitPartMatrix.coeffRef(col, row + s*varDofs) += temp2 - delta;
                        m_nExplicitPartMatrix.coeffRef(row, row + s*varDofs) += delta;
                        m_nExplicitPartMatrix.coeffRef(col, col + s*varDofs) += delta;
                    }
                }
            }
        }

        //block A nonlin.
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlin(), col); it; ++it)
                m_nExplicitPartMatrix.coeffRef(it.row(), it.col()) += it.value();

        //blockReaction
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksReaction(), col); it; ++it)
                m_nExplicitPartMatrix.coeffRef(it.row(), it.col()) += it.value();
    }

    //===================================================== fillBase ================================================================

    virtual void fillBase()  //+pattern
    {

        //---------------------------
        //--------- fillBase --------
        //---------------------------
        
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        const T invTimeStep = 1. / m_timeStepSize;

        // lumped mass matrix M
        m_lumpedM.resize(varDofs, dofs);
        fillLumpedMassMatrix();

        // pattern blocks
        gsSparseMatrix<T> matrixNpattern(varDofs, dofs);
        gsSparseMatrix<T> matrixAnonlinPattern(varDofs, dofs);
        gsSparseMatrix<T> matrixReactionPattern(varDofs, dofs);
        this->fillPatternBlocks(matrixNpattern, matrixAnonlinPattern, matrixReactionPattern);

        m_baseMatrix.resize(varDofs, dofs);
        m_baseMatrix = invTimeStep * m_lumpedM + matrixNpattern + matrixAnonlinPattern + matrixReactionPattern;

        m_nExplicitPartMatrix = matrixNpattern + matrixAnonlinPattern + matrixReactionPattern;

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

protected:
    void fillLumpedMassMatrix()
    {
        int varDofs = this->numVarDofs();
        int numVar = this->getNumVar();

        const gsSparseMatrix<T>& blockMref = m_blockAssembler.getBlockM();

        m_lumpedM.reserve(gsVector<int>::Constant(m_lumpedM.cols(), 1));

        for (int j = 0; j < varDofs; j++)
            for (int s = 0; s < numVar; s++)
                m_lumpedM.insert(j, j + s*varDofs) = blockMref.at(j, j);  

        for (int j = 0; j < varDofs; j++)
            for (typename gsSparseMatrix<T>::InnerIterator it(blockMref, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                if ((blockMref.at(k, l) != 0) && k != l)
                {
                    for (int s = 0; s < numVar; s++)
                        m_lumpedM.coeffRef(k, k + s*varDofs) += blockMref.at(k, l);
                }
            }
    }

    //==================================================== fillNBase ===============================================================

public:
    virtual void fillNBase() 
    {
        int numVar = this->getNumVar();
        const T invTimeStep = 1. / m_timeStepSize;
        m_solPicard = m_solution;

        m_nBaseMatrix = m_baseMatrix;

        int varDofs = this->numVarDofs();

        const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN();

        // block N (AFC modified)
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
        {
            for (index_t s = 0; s < numVar; ++s)
            {
                m_nBaseMatrix.coeffRef(col, col + s*varDofs) += m_theta * blockNref.at(col, col);
            }

            for (typename gsSparseMatrix<T>::InnerIterator it(blockNref, col); it; ++it)
            {
                int row = it.row();

                if (col > row)
                {
                    T temp1 = it.value(); // N(row, col)
                    T temp2 = blockNref.at(col, row);
                    T delta = math::max(0.0, math::max(temp1, temp2));

                    for (index_t s = 0; s < numVar; ++s)
                    {
                        m_nBaseMatrix.coeffRef(row, col + s*varDofs) += m_theta * (temp1 - delta);
                        m_nBaseMatrix.coeffRef(col, row + s*varDofs) += m_theta * (temp2 - delta);
                        m_nBaseMatrix.coeffRef(row, row + s*varDofs) += m_theta * delta;
                        m_nBaseMatrix.coeffRef(col, col + s*varDofs) += m_theta * delta;
                    }
                }
            }
        }

        m_nBaseRhs = m_blockAssembler.getRhsN() + (1 - m_theta) * m_nExplicitPartRhs;

        // block M
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t s = 0; s < numVar; ++s)
            m_nBaseRhs.col(s).noalias() += invTimeStep * m_lumpedM.leftCols(varDofs) * m_solution.col(s);

    } //end fillNBase

    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        int dofs = this->numDofs();

        //---------------------------
        //--------- fillMatrix ------
        //---------------------------
        m_matrix = m_nBaseMatrix;

        //block A nonlin.
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlin(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += m_theta * it.value();

        //blockReaction
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksReaction(), col); it; ++it)
                m_matrix.coeffRef(it.row() , it.col()) += m_theta * it.value();

        if (!m_matrix.isCompressed())
            m_matrix.makeCompressed();

        //---------------------------
        //--------- fillRhs ---------
        //---------------------------

        m_rhs = m_nBaseRhs + m_blockAssembler.getRhsAnonlin() + m_blockAssembler.getRhsReaction()
              + m_theta * m_blockAssembler.getRhsF();

        if (m_bAFC_HO)
        {
            gsMatrix<T> afcFlux(this->numVarDofs(), this->getNumVar());
            fillAfcFluxes(afcFlux);
            m_rhs += afcFlux;
        }

        //blending term on the right hand side
        m_rhs.col(1).noalias() += m_blockAssembler.getBlockBlend() * m_solution.col(1)
                                + m_blockAssembler.getRhsBlend();

        applyDirichletConditions();

        m_bSystemReady = true;

    } //end fillSystem

protected:
    void fillAfcFluxes(gsMatrix<T>& ff)
    {
        int varDofs = this->numVarDofs();
        int numVar = this->getNumVar();
        int dofs = this->numDofs();

        const T invTimeStep = 1. / m_timeStepSize;

        const gsSparseMatrix<T>& blockMref = m_blockAssembler.getBlockM();
        const gsSparseMatrix<T>& blockNref = m_blockAssembler.getBlockN();

        gsSparseMatrix<T> flux;
        flux.resize(varDofs, dofs);
        flux.leftCols(varDofs) = blockMref + blockNref;
        flux.rightCols(varDofs) = blockMref + blockNref;

        gsMatrix<T> Pp, Pm, Qp, Qm;
        Pp.setZero(varDofs, numVar);
        Pm.setZero(varDofs, numVar);
        Qp.setZero(varDofs, numVar);
        Qm.setZero(varDofs, numVar);

        for (int j = 0; j < varDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(flux, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                T temp1 = blockNref.at(k, l);
                T temp2 = blockNref.at(k, l);
                T delta = math::max(0.0, math::max(temp1, temp2));

                for (int s = 0; s < numVar; s++)
                {
                    flux.coeffRef(k, l + s*varDofs) =
                        invTimeStep * (blockMref.at(k, l) * (m_solPicard(k, s) - m_solPicard(l, s)) - blockMref.at(k, l) * (m_solution(k, s) - m_solution(l, s))) +
                        m_theta * delta * (m_solPicard(k, s) - m_solPicard(l, s)) +
                        (1 - m_theta) * delta * (m_solution(k, s) - m_solution(l, s));

                    if (k != l)
                    {
                        Pp.coeffRef(k, s) += math::max(0.0, flux.at(k, l + s*varDofs));
                        Pm.coeffRef(k, s) += math::min(0.0, flux.at(k, l + s*varDofs));

                        Qp.coeffRef(k, s) = math::max(Qp.coeffRef(k, s), m_solPicard(l, s) - m_solPicard(k, s));
                        Qm.coeffRef(k, s) = math::min(Qm.coeffRef(k, s), m_solPicard(l, s) - m_solPicard(k, s));
                    }

                }
            }
        }

        gsMatrix<T> Rp, Rm;
        Rp.setZero(varDofs, numVar);
        Rm.setZero(varDofs, numVar);

        for (int k = 0; k < varDofs; k++)
        {
            for (int s = 0; s < numVar; s++)
            {
                if (Pp(k, s) == 0)
                    Rp.coeffRef(k, s) = 1.0;
                else
                    Rp.coeffRef(k, s) = math::min(1.0, m_lumpedM.at(k, k) * Qp(k, s) / Pp(k, s));
                if (Pm(k, s) == 0)
                    Rm.coeffRef(k, s) = 1.0;
                else
                    Rm.coeffRef(k, s) = math::min(1.0, m_lumpedM.at(k, k) * Qm(k, s) / Pm(k, s));
            }
        }

        ff.setZero(varDofs, numVar);
        for (int j = 0; j < varDofs; j++)
        {
            for (typename gsSparseMatrix<T>::InnerIterator it(flux, j); it; ++it)
            {
                int k = it.row();
                int l = it.col();

                for (int s = 0; s < numVar; s++)
                {
                    if (flux.at(k, l + s*varDofs) > 0.0)
                    {
                        ff.coeffRef(k, s) += math::min(Rp(k, s), Rm(l, s))*flux.at(k, l + s*varDofs);
                    }
                    else if (flux.at(k, l + s*varDofs) < 0.0)
                    {
                        ff.coeffRef(k, s) += math::min(Rm(k, s), Rp(l, s))*flux.at(k, l + s*varDofs);
                    }
                }
            }
        }
    }

    void applyDirichletConditions()
    {
        int varDofs = this->numVarDofs();
        int numVar = this->getNumVar();

        const std::vector<gsMatrix<T> > & ddofs = m_pBlockAssemblerBC->getDirichletDofs();
        const gsMultiBasis<T>& basis = m_pBlockAssemblerBC->getSolBasis();

        for (typename gsBoundaryConditions<>::const_iterator it = m_pBlockAssemblerBC->getBCs().dirichletBegin(); it != m_pBlockAssemblerBC->getBCs().dirichletEnd(); ++it)
        {
            gsMatrix<unsigned> boundary = basis.piece(it->patch()).boundary(it->side());

            for (index_t k = 0; k != boundary.size(); ++k)
            {
                // dof number with dirichlet::none
                index_t ii = m_blockAssembler.getMapper().index(boundary.at(k), it->patch());

                for (int col = 0; col < varDofs; ++col)
                    for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix, col); it; ++it)
                    {
                        if (it.row() == ii)
                        {
                            for (index_t s = 0; s < numVar; ++s)
                                m_matrix.coeffRef(ii, col + s * varDofs) = 0.0;
                        }
                    }

                for (index_t s = 0; s < numVar; ++s)
                {
                    m_matrix.coeffRef(ii, ii + s*varDofs) = 1.0;

                    index_t bi = m_pBlockAssemblerBC->getMapper().bindex(boundary.at(k), it->patch());
                    m_rhs(ii, s) = ddofs[s](bi, 0);
                }
            }
        }
    }

protected:
    uwbTMBlockAssembler<T>* m_pBlockAssemblerBC;
    gsSparseMatrix<T> m_lumpedM;
    gsMatrix<T> m_solPicard;
    bool m_bAFC_HO; //high order

    // members from uwbTMAssemblerBaseUnsteady
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseMatrix;
    using uwbTMAssemblerBaseUnsteady<T>::m_nBaseRhs;
    using uwbTMAssemblerBaseUnsteady<T>::m_timeStepSize;
    using uwbTMAssemblerBaseUnsteady<T>::m_theta;
    using uwbTMAssemblerBaseUnsteady<T>::m_nExplicitPartRhs;
    using uwbTMAssemblerBaseUnsteady<T>::m_nExplicitPartMatrix;


    // members from uwbTMAssemblerBase
    using Base::getBlockAssembler;
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;

}; //uwbTMAssemblerKOmega_AFC

} //namespace gismo
