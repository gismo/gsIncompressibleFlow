/** @file uwbTMCrosswindVisitors.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbTMBlockVisitorsKOmega.h"
#include "uwbStabilizationEvaluator.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbTMCrosswindVisitorKOmega : public uwbTMBlockVisitor<T>
{

public:
    typedef uwbTMBlockVisitor<T> Base;

public:
    uwbTMCrosswindVisitorKOmega(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_deg(0), m_tauStabType(0), m_crosswindType(0), m_cwResidualType(0),
        m_r(2), m_timeStep(0.),
        m_bOldSolSet(false)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
    }

    ~uwbTMCrosswindVisitorKOmega()
    {
        if (m_pStabEvaluator)
        {
            delete m_pStabEvaluator;
            m_pStabEvaluator = NULL;
        }
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        Base::initialize(basisRefs, patchIndex, options, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        Base::initialize(basisRefs, patchIndex, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    virtual inline void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }


    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);

        basisRefs.back().evalAllDers_into(quNodes, 2, m_basisDataKO);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        m_deg = basisRefs.back().maxDegree();

        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        const index_t numActKOmega = m_activesKO.rows();
        const index_t numActU = m_activesU.rows();

        // Evaluate basis functions on element  nodes
        basisRefs.front().deriv_into(quNodes, m_basisGradsU);

        // Evaluate solution on element nodes
        m_solKOmegaVals = m_solTM.value(quNodes, m_patchIndex);
        gsMatrix<T> solKOmegaValsOld = m_solKOmegaVals;
        if (m_bUnsteady)
        {
            GISMO_ASSERT(m_bOldSolSet, "Old time step solution not set in Crosswind visitor.");
            solKOmegaValsOld = m_oldSolTM.value(quNodes, m_patchIndex);
        }

        gsMatrix<T> solActUCoeffs;
        solActUCoeffs.setZero(m_dim, numActU);
        for (int j = 0; j < numActU; j++) {
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();
        }
        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        index_t nQuPoints = quNodes.cols();
        m_solUGrads.resize(nQuPoints);
        m_solKOmegaGrads.resize(nQuPoints);
        gsMatrix<T> physGrad_u;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, m_basisGradsU, physGrad_u);
            m_solUGrads[k].noalias() = solActUCoeffs * physGrad_u.transpose();
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);
            m_solKOmegaGrads[k].noalias() = solActKOmegaCoeffs * m_physGradKO.transpose();
        }

        getTMEvaluator()->initAtElement(m_solUGrads, m_solKOmegaVals, m_solKOmegaGrads);
        if (this->checkWallDistanceBasedTM())
        {
            getTMEvaluator()->setKOmegaVariant(m_evaluatorType);

            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = m_patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            m_basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            m_basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            m_solPoissonVals = m_solPoisson.value(quNodes, m_patchIndex);
            gsMatrix<T> solActPoissonCoeffs;
            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = m_solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            getTMEvaluator()->evalWallDistance(m_solPoissonVals, solPoissonGrads);

        }

        if (m_evaluatorType == "koSAS" || m_evaluatorType == "koSAS_SS" || m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO")
        {
            std::vector<gsMatrix<T> > solULaplaces(nQuPoints);
            gsMatrix<T> basisHessian, physLaplacianU;
            basisRefs.front().deriv2_into(quNodes, basisHessian);

            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, m_basisGradsU, basisHessian, physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                solULaplaces[k].noalias() = solActUCoeffs * physLaplacianU.transpose();
            }
            getTMEvaluator()->setULaplacian(solULaplaces);
        }

        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
        getTMEvaluator()->evalQuantities_rhsPart();

        gsMatrix<T> physLaplacian;
        std::vector<gsVector<T> > residual(m_numVar);
        for (int var = 0; var < m_numVar; var++)
            residual[var].setZero(nQuPoints);
        switch (m_cwResidualType)
        {
        case 0:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, m_basisDataKO[1], m_basisDataKO[2], physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActKOmegaCoeffs * physLaplacian.transpose();

                residual[0](k) = m_solKOmegaGrads[k].row(0) * m_solUVals.col(k) //advection
                               - getTMEvaluator()->getKDiffusionCoefficient(k) * solLaplacian(0, 0) //diffusion
                               + getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k) * m_solKOmegaVals(0, k) //reaction
                               - getTMEvaluator()->getRhsK(k); //source

                residual[1](k) = m_solKOmegaGrads[k].row(1) * m_solUVals.col(k) //advection
                               - getTMEvaluator()->getOmegaDiffusionCoefficient(k) * solLaplacian(1, 0) //diffusion
                               + getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k) * m_solKOmegaVals(1, k) // reaction
                               - getTMEvaluator()->getRhsOmega(k) // source
                               - getTMEvaluator()->getBlendCoeff(k) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15)) *
                                 m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1)); //blend

                if (m_bUnsteady)
                {
                    residual[0](k) += 1./m_timeStep * (m_solKOmegaVals(0, k) - solKOmegaValsOld(0, k));
                    residual[1](k) += 1./m_timeStep * (m_solKOmegaVals(1, k) - solKOmegaValsOld(1, k));
                }
            }
            break;
        case 1:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                residual[0](k) = m_solKOmegaGrads[k].row(0) * m_solUVals.col(k); //advection
                residual[1](k) = m_solKOmegaGrads[k].row(1) * m_solUVals.col(k); //advection
            }
            break;
        case 2:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                residual[0](k) = m_solKOmegaGrads[k].row(0) * m_solUVals.col(k); //advection
                residual[1](k) = m_solKOmegaGrads[k].row(1) * m_solUVals.col(k); //advection
                if (m_bUnsteady)
                {
                    residual[0](k) += 1./m_timeStep * (m_solKOmegaVals(0, k) - solKOmegaValsOld(0, k));
                    residual[1](k) += 1./m_timeStep * (m_solKOmegaVals(1, k) - solKOmegaValsOld(1, k));
                }
            }
            break;
        default:
            gsWarn << "Wrong residualType choosen for TM crosswind stabilization parameter. Choosen default case 1.\n";
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                residual[0](k) = m_solKOmegaGrads[k].row(0) * m_solUVals.col(k); //advection
                residual[1](k) = m_solKOmegaGrads[k].row(1) * m_solUVals.col(k); //advection
            }
        }


        //getTMEvaluator()->computeAverage(residual[0]);
        //getTMEvaluator()->computeAverage(residual[1]);

        m_pStabEvaluator->initAtElement(m_solUVals);
        m_pStabEvaluator->setCrosswindVars(m_crosswindType, m_deg, m_timeStep, m_solKOmegaGrads, residual, m_elementLength);
        m_pStabEvaluator->setTauType(m_tauStabType, m_r);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);

        //const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & basisGradsKO = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsKO, m_physGradKO);

            T cwStabParamK;
            T cwStabParamO;
            if (m_tauStabType >= 12)
            {
                /*cwStabParamK = m_pStabEvaluator->getCrosswindStabParam(k, 0, getTMEvaluator()->getKDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBetaStar(k));
                cwStabParamO = m_pStabEvaluator->getCrosswindStabParam(k, 1, getTMEvaluator()->getOmegaDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBeta(k));*/

                cwStabParamK = m_pStabEvaluator->getCrosswindStabParam(k, 0, getTMEvaluator()->getKDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k));
                cwStabParamO = m_pStabEvaluator->getCrosswindStabParam(k, 1, getTMEvaluator()->getOmegaDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k));
            }
            else {
                cwStabParamK = m_pStabEvaluator->getCrosswindStabParam(k, 0, getTMEvaluator()->getKDiffusionCoefficient(k));
                cwStabParamO = m_pStabEvaluator->getCrosswindStabParam(k, 1, getTMEvaluator()->getOmegaDiffusionCoefficient(k));
            }
            gsMatrix<T> proj = m_pStabEvaluator->getCrosswindProjection(k);

            m_localMat[0].noalias() += weight * cwStabParamK * ((proj * m_physGradKO).transpose() * m_physGradKO);
            m_localMat[1].noalias() += weight * cwStabParamO * ((proj * m_physGradKO).transpose() * m_physGradKO);
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        const index_t KOmegaSize = m_mapperTM.freeSize();

        // Local Dofs to global dofs
        m_mapperTM.localToGlobal(m_activesKO, m_patchIndex, m_activesKO);
        const index_t numActKOmega = m_activesKO.rows();

        for (index_t i = 0; i < numActKOmega; ++i)
        {
            const int ii = m_activesKO(i);
            if (m_mapperTM.is_free_index(ii))
            {
                for (index_t j = 0; j < numActKOmega; ++j)
                {
                    const int jj = m_activesKO(j);
                    if (m_mapperTM.is_free_index(jj))
                    {
                        for (index_t s = 0; s != m_numVar; ++s)
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //CW = [CW_k, CW_omega]
                    }
                    else // m_mapperTM.is_boundary_index(jj)
                    {
                        const int bb = m_mapperTM.global_to_bindex(jj);
                        for (index_t s = 0; s != m_numVar; ++s)
                            rhs(ii, s) -= m_localMat[s](i, j) * eliminatedDofs[s](bb, 0); //rhs = [rhs_k, rhs_omega]
                    }
                }
            }
        }
    }

    void setCrosswind(const int crosswindType, const int cwResidualType, bool unsteady = false, const T timeStep = 0.)
    {
        m_crosswindType = crosswindType;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setCrosswind(const int crosswindType, const int cwResidualType, const int tauType = 1, bool unsteady = false,
                      const T timeStep = 0., const int r = 2)
    {
        m_crosswindType = crosswindType;
        m_tauStabType = tauType;
        m_r = r;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setOldSolutionField(gsField<T>& solution)
    {
        m_oldSolTM = solution;
        m_bOldSolSet = true;
    }

    virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");

        if (m_evaluatorType == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else // (m_evaluatorType == "koSST") or its variants
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }

protected:
    int m_deg;
    int m_tauStabType;
    int m_crosswindType;
    int m_cwResidualType;
    int m_r;
    T m_timeStep;
    bool m_bUnsteady;
    bool m_bOldSolSet;

    // Basis values
    gsMatrix<index_t> m_activesKO;
    std::vector<gsMatrix<T> > m_basisDataKO;
    gsMatrix<T> m_solUVals;
    gsMatrix<T> m_physGradKO;
    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_solKOmegaVals;
    gsMatrix<T> m_solPoissonVals;
    gsMatrix<T> m_basisGradsU;

    std::vector<gsMatrix<T> > m_solUGrads;
    std::vector<gsMatrix<T> > m_solKOmegaGrads;
    std::vector<gsMatrix<T> > m_localMat;

    gsField<T> m_oldSolTM;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    // members from uwbTMBlockVisitor
    using Base::m_patchIndex;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_solTM;
    using Base::m_evaluatorType;
    using Base::m_pEvaluator;
    using Base::m_numVar;
    using Base::m_solPoisson;
    using Base::m_patchesPoisson;
    using Base::m_basesPoisson;
    using Base::m_mapperTM;

    using uwbVisitorBase<T>::m_elementLength;
};

} // namespace gismo

