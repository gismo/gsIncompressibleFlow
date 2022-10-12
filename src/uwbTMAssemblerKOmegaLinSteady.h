/** @file uwbTMAssemblerKOmegaLinSteady.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMBlockVisitorsKOmega.h"

namespace gismo
{

template<class T>
class uwbTMAssemblerKOmegaLinSteady : public uwbTMAssemblerBase<T>
{

public:
    typedef uwbTMAssemblerBase<T> Base;

public:
    uwbTMAssemblerKOmegaLinSteady(uwbINSSolverParams<T>& params) :
        Base(params, 2)
    { }

    virtual ~uwbTMAssemblerKOmegaLinSteady()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        m_blockAssembler.updateVelocitySolution(uSolField);
        m_blockAssembler.assembleNewTimestepPart_kOmega();
        m_blockAssembler.assemblePatternBlocks_kOmega();
    }

    virtual void updateAssembly()
    {
        m_blockAssembler.assembleNonlinearPart_kOmega();
    }

    virtual void updatePicardAssembly()
    { 
        m_blockAssembler.assembleNonlinearPart_kOmega();
    }

public:
    //===================================================== fillBase ================================================================

    virtual void fillBase()  //+pattern
    {
        //---------------------------
        //--------- fillBase --------
        //---------------------------
        
        int varDofs = this->numVarDofs();
        int numVar = this->getNumVar();
        int dofs = this->numDofs();

        // blockN
        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t s = 0; s < numVar; ++s)
            for (index_t i = 0; i < varDofs; ++i)
                nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getBlockN().col(i).nonZeros();

        gsSparseMatrix<T> matrixN(varDofs, dofs);
        matrixN.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    matrixN.insert(it.row(), it.col() + s*varDofs) = it.value();

        // block Apattern
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksAnonlinPattern().col(i).nonZeros();

        gsSparseMatrix<T> matrixAnonlinPattern(varDofs, dofs);
        matrixAnonlinPattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlinPattern(), col); it; ++it)
                matrixAnonlinPattern.insert(it.row(), it.col()) = 0.;

        m_baseMatrix.resize(varDofs, dofs);
        m_baseMatrix = matrixN + matrixAnonlinPattern;

        // base rhs
        m_baseRhs = m_blockAssembler.getRhsN();

        if (!m_baseMatrix.isCompressed() && !m_blockAssembler.isStabilization())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        int dofs = this->numDofs();

        //---------------------------
        //--------- fillMatrix ------
        //---------------------------
        m_matrix = m_baseMatrix;

        // block A nonlin.
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlin(), col); it; ++it)
                m_matrix.coeffRef(it.row(), it.col()) += it.value();

        //---------------------------
        //--------- fillRhs ---------
        //---------------------------

        m_rhs = m_baseRhs + m_blockAssembler.getRhsF() + m_blockAssembler.getRhsAnonlin();

        if (!m_matrix.isCompressed() && !m_blockAssembler.isStabilization())
            m_matrix.makeCompressed();

        if (!m_blockAssembler.isStabilization())
            m_bSystemReady = true;
    } //end fillSystem

    virtual gsMatrix<T> getSolutionK_full(const gsDofMapper& pMapper, const gsMatrix<T>& solVector) const
    {
        const gsDofMapper& mapper = m_blockAssembler.getMapper();
        const gsMatrix<T>& ddofs = m_blockAssembler.getDirichletDofs().at(0); // Dirichlet dofs coefs for k

        index_t totalDofs = mapper.size();
        gsMatrix<T> solKFull;
        solKFull.setZero(totalDofs, 1);

        for (unsigned int p = 0; p < m_blockAssembler.getPatches().nPatches(); ++p)
        {
            const int sz = m_blockAssembler.getSolBasis().piece(p).size();
            for (index_t i = 0; i < sz; ++i)
            {
                index_t ii = mapper.index(i, p);
                index_t jj = pMapper.index(i, p);

                if (mapper.is_free_index(ii))
                {
                    solKFull(jj, 0) = solVector(ii, 0);
                }
                else
                {
                    solKFull(jj, 0) = ddofs(mapper.global_to_bindex(ii), 0);
                }
            }
        }

        return solKFull;
    }

    virtual void evalTurbCoeff_into(index_t patchIndex, gsMatrix<T>& points, gsMatrix<T> solution, std::vector<gsMatrix<T> >& solUGrads, gsVector<T>& turbCoeffVals, std::string coeffType = "turbViscosity")
    {
        getTMEvaluator()->initialize(getViscosity(), points.cols());
        getTMEvaluator()->setAveraging(false);

        gsField<T> solTM = constructSolution(solution);
        gsMatrix<T> solKOVals = solTM.value(points, patchIndex);

        index_t nQuPoints = points.cols();

        //--- KOGrads
        gsGeometry<T>* patchKO = &m_blockAssembler.getPatches().patch(patchIndex);
        const gsMultiBasis<T> basisKO = m_blockAssembler.getSolBasis();
        unsigned evFlagsKO = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlagsKO, *patchKO));
        gsMatrix<T> solActKOCoeffs, physGradKO, bGradsKO;
        gsMatrix<index_t> activesKO;

        basisKO.basis(patchIndex).active_into(points.col(0), activesKO);
        basisKO.basis(patchIndex).deriv_into(points, bGradsKO);

        geoEval->evaluateAt(points);

        const index_t numActKO = activesKO.rows();

        solActKOCoeffs.setZero(2, numActKO);
        for (int j = 0; j < numActKO; j++)
            solActKOCoeffs.col(j) = solTM.coefficientVector(patchIndex).row(activesKO(j)).transpose();

        std::vector<gsMatrix<T> > solKOGrads(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval->transformGradients(k, bGradsKO, physGradKO);
            solKOGrads[k].noalias() = solActKOCoeffs * physGradKO.transpose();
        }
        //---

        getTMEvaluator()->initAtElement(solUGrads, solKOVals, solKOGrads);

        std::string tmEvaluator = m_blockAssembler.getTMEvaluator();
        if (m_blockAssembler.checkWallDistanceBasedTM())
        {
            getTMEvaluator()->setKOmegaVariant(tmEvaluator);

            const gsField<T> solPoisson = m_blockAssembler.getPoissonSolution();
            gsGeometry<T>* patch = &m_blockAssembler.getPoissonPatches().patch(patchIndex);
            const gsMultiBasis<T> basisPoisson = m_blockAssembler.getPoissonBasis();
            unsigned evFlags = NEED_VALUE | NEED_GRAD_TRANSFORM | NEED_DERIV2;
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlags, *patch));
            gsMatrix<T> solPoissonVals, solActPoissonCoeffs, physGradPoisson, bGradsPoisson;
            gsMatrix<index_t> activesPoisson;

            basisPoisson.basis(patchIndex).active_into(points.col(0), activesPoisson);
            basisPoisson.basis(patchIndex).deriv_into(points, bGradsPoisson);

            geoEvalPoisson->evaluateAt(points);

            const index_t numActPoisson = activesPoisson.rows();

            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(patchIndex).row(activesPoisson(j)).transpose();

            solPoissonVals = solPoisson.value(points, patchIndex);
            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, bGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            getTMEvaluator()->evalWallDistance(solPoissonVals, solPoissonGrads);
       }

        if (tmEvaluator == "koSAS" || tmEvaluator == "koSAS_SS" || tmEvaluator == "koSAS_SO" || tmEvaluator == "koSAS_OO")
        {
            std::vector<gsMatrix<T> > solULaplaces(nQuPoints);
            gsMatrix<T> physLaplacianU;
            std::vector<gsMatrix<T> > basisDataU;
            m_blockAssembler.getBases().front().basis(patchIndex).evalAllDers_into(points, 2, basisDataU);

            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval->transformLaplaceHgrad(k, basisDataU[1], basisDataU[2], physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                solULaplaces[k].noalias() = m_actCoeffsU * physLaplacianU.transpose();
            }
            getTMEvaluator()->setULaplacian(solULaplaces);
        }

        if (coeffType == "kDiffusionCoeff")
        {
            getTMEvaluator()->evalQuantities_diffusionCoeff();
            turbCoeffVals = getTMEvaluator()->getKDiffusionCoefficient();
        }
        else if (coeffType == "oDiffusionCoeff")
        {
            getTMEvaluator()->evalQuantities_diffusionCoeff();
            turbCoeffVals = getTMEvaluator()->getOmegaDiffusionCoefficient();
        }
        else if (coeffType == "kReactionCoeff")
        {
            getTMEvaluator()->evalQuantities_reactionCoeff();
            turbCoeffVals = getTMEvaluator()->getBetaStar();
        }
        else if (coeffType == "oReactionCoeff")
        {
            getTMEvaluator()->evalQuantities_reactionCoeff();
            turbCoeffVals = getTMEvaluator()->getBeta();
        }
        else if (coeffType == "blendCoeff")
        {
            getTMEvaluator()->evalQuantities_blendCoeff();
            turbCoeffVals = getTMEvaluator()->getBlendCoeff();
        }
        else if (coeffType == "kSourceCoeff")
        {
            getTMEvaluator()->evalQuantities_rhsPart();
            turbCoeffVals = getTMEvaluator()->getRhsK();
        }
        else if (coeffType == "oSourceCoeff")
        {
            getTMEvaluator()->evalQuantities_rhsPart();
            turbCoeffVals = getTMEvaluator()->getRhsOmega();
        }
        else if (coeffType == "F1")
        {
            getTMEvaluator()->evalTurbVariables();
            turbCoeffVals = getTMEvaluator()->getF1();
        }
        else if (coeffType == "vorticityMagnitude")
        {
            getTMEvaluator()->evalTurbVariables();
            turbCoeffVals = getTMEvaluator()->getVorticityMagnitude();
        }
        else if (coeffType == "strainRateMagnitude")
        {
            getTMEvaluator()->evalTurbVariables();
            turbCoeffVals = getTMEvaluator()->getStrainRateMagnitude();
        }
        else if (coeffType == "strainRateUGrad")
        {
            getTMEvaluator()->evalTurbVariables();
            turbCoeffVals = getTMEvaluator()->getStrainRateUGrad();
        }
        else if (coeffType == "sourceQSAS")
        {
            getTMEvaluator()->evalTurbVariables();
            turbCoeffVals = getTMEvaluator()->getSourceQSAS();
        }
        else if (coeffType == "viscostyRatio")
        {
            getTMEvaluator()->evalQuantities_turbViscosity();
            turbCoeffVals = getTMEvaluator()->getTurbViscosity();
            turbCoeffVals = turbCoeffVals / m_blockAssembler.getViscosity();
        }
        else if (coeffType == "turbViscosity")
        {
            getTMEvaluator()->evalQuantities_turbViscosity();
            turbCoeffVals = getTMEvaluator()->getTurbViscosity();
        }
        else
        {
            gsWarn << "Wrong 'coeffType' chosen. Default variant selected (i.e. 'turbViscosity')";
            getTMEvaluator()->evalQuantities_turbViscosity();
            turbCoeffVals = getTMEvaluator()->getTurbViscosity();
        }
    }

public:
    virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        std::string tmEvaluator = m_blockAssembler.getTMEvaluator();
        GISMO_ASSERT(tmEvaluator != "", "No evaluator type set in the visitor.");

        if (tmEvaluator == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else //(tmEvaluator == "koSST")
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }

protected:
    // members from uwbTMAssemblerBase
    using Base::m_blockAssembler;
    using Base::m_baseMatrix;
    using Base::m_baseRhs;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_solution;
    using Base::m_bSystemReady;

    using uwbTMAssemblerBase<T>::m_pEvaluator;
    using uwbTMAssemblerBase<T>::m_actCoeffsU;

    using uwbTMAssemblerBase<T>::getViscosity;
    using uwbTMAssemblerBase<T>::constructSolution;
    using uwbTMAssemblerBase<T>::getBlockAssembler;

}; //uwbTMAssemblerKOmegaLinSteady

} //namespace gismo
