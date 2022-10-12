/** @file uwbTMAssemblerKOmega.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMAssemblerBase.h"
#include "uwbTMBlockVisitorsKOmega.h"

namespace gismo
{

template<class T>
class uwbTMAssemblerKOmega : public uwbTMAssemblerBaseUnsteady<T>
{

public:
    typedef uwbTMAssemblerBaseUnsteady<T> Base;

public:
    uwbTMAssemblerKOmega(uwbINSSolverParams<T>& params) :
        Base(params, 2)
    { }

    virtual ~uwbTMAssemblerKOmega()
    { }

protected:

    virtual void initAssembly(const gsField<T> & uSolField)
    {
        m_blockAssembler.updateVelocitySolution(uSolField);
        if (m_blockAssembler.isTimeDerTerm())
            m_blockAssembler.assembleLinearPart_kOmega();
        m_blockAssembler.assemblePatternBlocks_kOmega(m_bAssembleAllBlocks);
    }

    virtual void updateAssembly(const gsField<T> & uSolField)
    {
        m_blockAssembler.updateVelocitySolution(uSolField);
        m_blockAssembler.assembleNewTimestepPart_kOmega();
        m_blockAssembler.assembleNonlinearPart_kOmega(m_bAssembleAllBlocks);
    }

    virtual void updatePicardAssembly()
    { 
        m_blockAssembler.assembleNonlinearPart_kOmega(m_bAssembleAllBlocks);
    }

    //===================================================== fillExplicit ================================================================
    virtual void fillExplicitPartMatrix()
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        m_nExplicitPartMatrix.setZero();

        // block N
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    m_nExplicitPartMatrix.coeffRef(it.row(), it.col() + s*varDofs) += it.value();

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

        if (m_bAssembleAllBlocks)
        {
            // pattern blocks
            gsSparseMatrix<T> matrixNpattern(varDofs, dofs);
            gsSparseMatrix<T> matrixPattern(varDofs, dofs);
            fillPatternBlocks(matrixNpattern, matrixPattern);

            m_baseMatrix = matrixNpattern + matrixPattern;

            if (m_blockAssembler.isTimeDerTerm())
            {
                // block M
                gsSparseMatrix<T> matrixM(varDofs, dofs);
                fillBlockM(matrixM);
                m_baseMatrix += matrixM;
            }
        }
        else
        {
            // pattern blocks
            gsSparseMatrix<T> matrixNpattern(varDofs, dofs);
            gsSparseMatrix<T> matrixAnonlinPattern(varDofs, dofs);
            gsSparseMatrix<T> matrixReactionPattern(varDofs, dofs);
            fillPatternBlocks(matrixNpattern, matrixAnonlinPattern, matrixReactionPattern);

            m_baseMatrix = matrixNpattern + matrixAnonlinPattern + matrixReactionPattern;

            if (m_blockAssembler.isTimeDerTerm())
            {
                // block M
                gsSparseMatrix<T> matrixM(varDofs, dofs);
                fillBlockM(matrixM);
                m_baseMatrix += matrixM;
            }

            m_nExplicitPartMatrix = matrixNpattern + matrixAnonlinPattern + matrixReactionPattern;
        }

        if (!m_baseMatrix.isCompressed()  && !m_blockAssembler.isStabilization())
            m_baseMatrix.makeCompressed();

        m_bSystemReady = false;

    } //end fillBase

protected:
    void fillBlockM(gsSparseMatrix<T>& matrixM)
    {
        const T invTimeStep = 1. / m_timeStepSize;

        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t s = 0; s < numVar; ++s)
            for (index_t i = 0; i < varDofs; ++i)
                nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getBlockM().col(i).nonZeros();

        matrixM.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockM(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    matrixM.insert(it.row(), it.col() + s*varDofs) = invTimeStep * it.value();
    }

    void fillPatternBlocks(gsSparseMatrix<T>& matrixNpattern, gsSparseMatrix<T>& matrixAnonlinPattern, gsSparseMatrix<T>& matrixReactionPattern)
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;

        //blockNpattern
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t s = 0; s < numVar; ++s)
            for (index_t i = 0; i < varDofs; ++i)
                nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getBlockNpattern().col(i).nonZeros();

        matrixNpattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    matrixNpattern.insert(it.row(), it.col() + s*varDofs) = 0.;

        //block Apattern
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksAnonlinPattern().col(i).nonZeros();

        matrixAnonlinPattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksAnonlinPattern(), col); it; ++it)
                matrixAnonlinPattern.insert(it.row(), it.col()) = 0.;

        //block reactionPattern
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksReactionPattern().col(i).nonZeros();

        matrixReactionPattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksReactionPattern(), col); it; ++it)
                matrixReactionPattern.insert(it.row(), it.col()) = 0.;
    }

    void fillPatternBlocks(gsSparseMatrix<T>& matrixNpattern, gsSparseMatrix<T>& matrixPattern)
    {
        int numVar = this->getNumVar();
        int varDofs = this->numVarDofs();
        int dofs = this->numDofs();

        gsVector<int> nonZerosPerColumnVector;

        //blockNpattern
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t s = 0; s < numVar; ++s)
            for (index_t i = 0; i < varDofs; ++i)
                nonZerosPerColumnVector(i + s*varDofs) = m_blockAssembler.getBlockNpattern().col(i).nonZeros();

        matrixNpattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockNpattern(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    matrixNpattern.insert(it.row(), it.col() + s*varDofs) = 0.;

        //blocks Anonlin + reaction + blend
        nonZerosPerColumnVector.setZero(dofs);
        for (index_t i = 0; i < dofs; ++i)
            nonZerosPerColumnVector(i) = m_blockAssembler.getBlocksPattern().col(i).nonZeros();

        matrixPattern.reserve(nonZerosPerColumnVector);

        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < dofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocksPattern(), col); it; ++it)
                matrixPattern.insert(it.row(), it.col()) = 0.;
    }

    //==================================================== fillNBase ===============================================================

public:
    virtual void fillNBase() 
    {
        int numVar = this->getNumVar();
        const T invTimeStep = 1. / m_timeStepSize;

        m_nBaseMatrix = m_baseMatrix;

        int varDofs = this->numVarDofs();

        // block N
        #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
        for (index_t col = 0; col < varDofs; ++col)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlockN(), col); it; ++it)
                for (index_t s = 0; s < numVar; ++s)
                    m_nBaseMatrix.coeffRef(it.row(), it.col() + s*varDofs) += m_theta * it.value();

        m_nBaseRhs = m_blockAssembler.getRhsN() + (1 - m_theta) * m_nExplicitPartRhs;

        if(m_blockAssembler.isTimeDerTerm())
        {
            // block M
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t s = 0; s < numVar; ++s)
                m_nBaseRhs.col(s).noalias() += invTimeStep * m_blockAssembler.getBlockM() * m_solution.col(s);
        }

    } //end fillNBase

    //===================================================== fillSystem ==============================================================

    virtual void fillSystem()
    {
        int dofs = this->numDofs();

        //---------------------------
        //--------- fillMatrix ------
        //---------------------------
        m_matrix = m_nBaseMatrix;

        if(m_bAssembleAllBlocks)
        {
            //block A nonlin. + reaction + blend
            #pragma omp parallel for num_threads(m_blockAssembler.getNumThreads())
            for (index_t col = 0; col < dofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockAssembler.getBlocks(), col); it; ++it)
                    m_matrix.coeffRef(it.row(), it.col()) += it.value();

            //---------------------------
            //--------- fillRhs ---------
            //---------------------------
            m_rhs = m_nBaseRhs + m_blockAssembler.getRhs();
        }
        else
        {
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

            //---------------------------
            //--------- fillRhs ---------
            //---------------------------

            m_rhs = m_nBaseRhs + m_theta * m_blockAssembler.getRhsF() + m_blockAssembler.getRhsAnonlin()
                  + m_blockAssembler.getRhsReaction();
            //blending term on the right hand side
            m_rhs.col(1).noalias() -= m_blockAssembler.getBlockBlend() * m_solution.col(1);
            m_rhs.col(1).noalias() += m_blockAssembler.getRhsBlend();
        }

        if (!m_matrix.isCompressed() && !m_blockAssembler.isStabilization())
            m_matrix.makeCompressed();

        if (!m_blockAssembler.isStabilization())
            m_bSystemReady = true;

    } //end fillSystem


    virtual gsMatrix<T> getSolutionK_full(const gsDofMapper& pMapper, const gsMatrix<T>& solVector) const
    {
        const gsDofMapper& mapper = m_blockAssembler.getMapper();
        const gsMatrix<T>& ddofs = m_blockAssembler.getDirichletDofs().at(0); // Dirichlet dofs coefs for k

        index_t totalDofs = pMapper.size();
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
        unsigned evFlagsKO = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlagsKO, *patchKO));
        gsMatrix<T> solActKOCoeffs, physGradKO;  //, bGradsKO;
        std::vector<gsMatrix<T> > basisDataKO;
        gsMatrix<index_t> activesKO;

        basisKO.basis(patchIndex).active_into(points.col(0), activesKO);
        basisKO.basis(patchIndex).evalAllDers_into(points, 2, basisDataKO);
        //basisKO.basis(patchIndex).deriv_into(points, bGradsKO);

        geoEval->evaluateAt(points);

        const index_t numActKO = activesKO.rows();

        solActKOCoeffs.setZero(2, numActKO);
        for (int j = 0; j < numActKO; j++)
            solActKOCoeffs.col(j) = solTM.coefficientVector(patchIndex).row(activesKO(j)).transpose();

        std::vector<gsMatrix<T> > solKOGrads(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval->transformGradients(k, basisDataKO[1], physGradKO);
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

        turbCoeffVals.setZero(nQuPoints);
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
            gsVector<T> betaStar = getTMEvaluator()->getBetaStar();
            for (int ii = 0; ii < nQuPoints; ii++)
                turbCoeffVals(ii) = betaStar(ii)*solKOVals.coeff(1, ii);
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
            turbCoeffVals = turbCoeffVals / getViscosity();
        }
        else if (coeffType == "turbViscosity")
        {
            getTMEvaluator()->evalQuantities_turbViscosity();
            turbCoeffVals = getTMEvaluator()->getTurbViscosity();
        }
        else if (coeffType == "cwParamK" || coeffType == "cwParamO")
        {
            getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
            getTMEvaluator()->evalQuantities_rhsPart();

            gsMatrix<T> solUVals = this->m_uSol.value(points, patchIndex);

            gsField<T> oldSolTM = constructSolution(this->m_oldSolution);
            gsMatrix<T> oldSolKOVals = oldSolTM.value(points, patchIndex);

            gsMatrix<T> physLaplacian;
            gsVector<T> residual;
            residual.setZero(nQuPoints);
            gsVector<T> cwStabParam;
            cwStabParam.setZero(nQuPoints);

            T elemDiam = (geoEval->value(0) - geoEval->value(nQuPoints-1)).norm();
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval->transformLaplaceHgrad(k, basisDataKO[1], basisDataKO[2], physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActKOCoeffs * physLaplacian.transpose();

                T normAdvection = (solUVals.col(k)).norm();
                T denom;

                T tau_s = elemDiam / (2. * normAdvection);

                if (coeffType == "cwParamK")
                {
                    residual(k) = solKOGrads[k].row(0) * solUVals.col(k) //advection
                                - getTMEvaluator()->getKDiffusionCoefficient(k) * solLaplacian(0, 0) //diffusion
                                + getTMEvaluator()->getBetaStar(k) * solKOVals(1, k) * solKOVals(0, k) //reaction
                                - getTMEvaluator()->getRhsK(k); //source
                    if (m_blockAssembler.isUnsteady())
                        residual(k) += 1./m_timeStepSize * (solKOVals(0, k) - oldSolKOVals(0, k));

                    denom = normAdvection * (solKOGrads[k].row(0)).norm() + math::abs(residual(k));
                    cwStabParam(k) = tau_s * math::pow(normAdvection, 2) * math::abs(residual(k)) /
                                  denom;
                }
                else
                {
                    residual(k) = solKOGrads[k].row(1) * solUVals.col(k) //advection
                                - getTMEvaluator()->getOmegaDiffusionCoefficient(k) * solLaplacian(1, 0) //diffusion
                                + getTMEvaluator()->getBeta(k) * solKOVals(1, k) * solKOVals(1, k) // reaction
                                - getTMEvaluator()->getRhsOmega(k) // source
                                - getTMEvaluator()->getBlendCoeff(k) / math::max(solKOVals(1, k), math::pow(10, -15)) *
                                  solKOGrads[k].row(0).dot(solKOGrads[k].row(1)); //blend
                    if (m_blockAssembler.isUnsteady())
                        residual(k) += 1./m_timeStepSize * (solKOVals(1, k) - oldSolKOVals(1, k));

                    denom = normAdvection * (solKOGrads[k].row(1)).norm() + math::abs(residual(k));
                    cwStabParam(k) = tau_s * math::pow(normAdvection, 2) * math::abs(residual(k)) /
                                  denom;
                }
                if (denom == 0 || math::isnan(cwStabParam(k)))
                    cwStabParam(k) = 0.;
            }

            turbCoeffVals = cwStabParam;
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
    // members from uwbTMAssemblerBaseUnsteady
    using Base::m_nBaseMatrix;
    using Base::m_nBaseRhs;
    using Base::m_nExplicitPartRhs;
    using Base::m_nExplicitPartMatrix;
    using Base::m_timeStepSize;
    using Base::m_theta;

    using Base::m_bAssembleAllBlocks;

    // members from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::getBlockAssembler;
    using uwbTMAssemblerBase<T>::m_blockAssembler;
    using uwbTMAssemblerBase<T>::m_baseMatrix;
    using uwbTMAssemblerBase<T>::m_matrix;
    using uwbTMAssemblerBase<T>::m_rhs;
    using uwbTMAssemblerBase<T>::m_solution;
    using uwbTMAssemblerBase<T>::m_bSystemReady;
    using uwbTMAssemblerBase<T>::m_pEvaluator;
    using uwbTMAssemblerBase<T>::m_actCoeffsU;

    // functions from uwbTMAssemblerBase
    using uwbTMAssemblerBase<T>::getViscosity;
    using uwbTMAssemblerBase<T>::constructSolution;

}; //uwbTMAssemblerKOmega

} //namespace gismo
