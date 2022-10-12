/** @file uwbTMSUPGBlockVisitorsKOmega.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once
#include "uwbTMBlockVisitorsKOmega.h"
#include "uwbStabilizationEvaluator.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbTMSUPGBlockVisitorKOmega : public uwbTMBlockVisitor<T>
{

public:
    typedef uwbTMBlockVisitor<T> Base;

public:
    uwbTMSUPGBlockVisitorKOmega(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_deg(0), m_tauStabType(0), m_timeStep(0.), m_r(2)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
        m_bAverageTau = false;
    }

    ~uwbTMSUPGBlockVisitorKOmega()
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

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);

        basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataKO);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        m_deg = basisRefs.back().maxDegree();

        m_pStabEvaluator->initAtElement(m_solUVals);
        m_pStabEvaluator->setSUPGvars(m_tauStabType, m_deg, m_timeStep, m_elementLength, m_r);
    }

    void setKOmegaDiffusionCoeffSolField(const std::vector<gsField<T> >& kOmegaDiffusionCoeffSolField)
    { m_kDiffusionCoeffSolField = kOmegaDiffusionCoeffSolField[0]; m_omegaDiffusionCoeffSolField = kOmegaDiffusionCoeffSolField[1]; }

    void setTauStabType(const int tauStabType, bool unsteady, const T timeStep = 0., const int r = 2)
    {
        m_tauStabType = tauStabType;
        m_timeStep = timeStep;
        if (!unsteady)
            m_timeStep = 0.;
        m_r = r;
    }

protected:
    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

protected:
    int m_deg;
    int m_tauStabType;
    T m_timeStep;
    int m_r;

    bool m_bAverageTau;

    // Basis values
    gsMatrix<index_t> m_activesKO;
    std::vector<gsMatrix<T> > m_basisDataKO;
    gsMatrix<T> m_solUVals;
    gsMatrix<T> m_physGradKO;

    gsVector<T> m_kTau, m_oTau;

    gsField<T> m_kDiffusionCoeffSolField;
    gsField<T> m_omegaDiffusionCoeffSolField;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    using Base::m_viscosity;
    using Base::m_patchIndex;
    using Base::m_solU;
    using Base::m_dim;

    using uwbVisitorBase<T>::m_elementLength;
};

// ============================================================= PARENT SUPG nonlinear k-o ============================================================= //
template <class T>
class uwbTMSUPGBlockNonlinVisitorKOmega : public uwbTMSUPGBlockVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockVisitorKOmega<T> Base;

public:
    uwbTMSUPGBlockNonlinVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) : Base(dofMappers, viscosity){ }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        const index_t numActKOmega = m_activesKO.rows();
        const index_t numActU = m_activesU.rows();

        // Evaluate basis functions on element  nodes
        basisRefs.front().deriv_into(quNodes, m_basisGradsU);

        // Evaluate solution on element nodes
        m_solKOmegaVals = m_solTM.value(quNodes, m_patchIndex);

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
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t nQuPoints = quWeights.rows();
        m_kTau.setZero(nQuPoints);
        m_oTau.setZero(nQuPoints);

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            if (this->m_tauStabType >= 12)
            {
                /*m_kTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getKDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBetaStar(k));
                m_oTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getOmegaDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBeta(k));*/

                m_kTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getKDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBetaStar(k) * this->m_solKOmegaVals(1, k));
                m_oTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getOmegaDiffusionCoefficient(k),
                                                                       getTMEvaluator()->getBeta(k) * this->m_solKOmegaVals(1, k));
            }
            else {
                m_kTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getKDiffusionCoefficient(k));
                m_oTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getOmegaDiffusionCoefficient(k));
            }

            //T m_kTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getKDiffusionCoefficient(k));
            //T m_oTau(k) = m_pStabEvaluator->getTauS(k, getTMEvaluator()->getOmegaDiffusionCoefficient(k));
        }

        if (Base::m_bAverageTau)
        {
            T sumTauK = 0.;
            T sumTauO = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
            {
                sumTauK += m_kTau(k);
                sumTauO += m_oTau(k);
            }
            m_kTau.setConstant(nQuPoints, sumTauK / nQuPoints);
            m_oTau.setConstant(nQuPoints, sumTauO / nQuPoints);
        }
    }

    virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");

        if (m_evaluatorType == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else // (m_evaluatorType == "koSST")
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }

protected:
    // Basis values
    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_solKOmegaVals;
    gsMatrix<T> m_solPoissonVals;

    gsMatrix<T> m_basisGradsU;
    gsMatrix<T> m_physGradKO;

    std::vector<gsMatrix<T> > m_solUGrads;
    std::vector<gsMatrix<T> > m_solKOmegaGrads;

    // memberes from uwbTMSUPGBlockVisitorKOmega
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_solU;
    using Base::m_pStabEvaluator;

    using Base::m_kTau;
    using Base::m_oTau;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_pEvaluator;
    using uwbTMBlockVisitor<T>::m_solTM;
    using uwbTMBlockVisitor<T>::m_dim; // number of unknown variables of turbulent model
    using uwbTMBlockVisitor<T>::m_patchIndex;
    using uwbTMBlockVisitor<T>::m_evaluatorType;
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_solPoisson;
    using uwbTMBlockVisitor<T>::m_patchesPoisson;
    using uwbTMBlockVisitor<T>::m_basesPoisson;
};

// ============================================================= BLOCK M ============================================================= //
template <class T>
class uwbTMSUPGBlockMVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlockMVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & basisGradsKO = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsKO, m_physGradKO);

            m_localMat[0].noalias() += weight * m_kTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
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
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //M = [NM_k, NM_omega]
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

protected:

    std::vector<gsMatrix<T> > m_localMat;

    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::getTMEvaluator;

    using NonlinBase::m_kTau;
    using NonlinBase::m_oTau;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_patchIndex;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;
};

// ============================================================= BLOCK N_TMkomega ============================================================= //
template <class T>
class uwbTMSUPGBlocksNVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlocksNVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega); // local [N_k(u^n), N_omega(u^n)]

        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            m_localMat[0].noalias() += weight * m_kTau(k) * ((m_solUVals.col(k).transpose() * m_physGradKO).transpose() * (m_solUVals.col(k).transpose() * m_physGradKO));
            m_localMat[1].noalias() += weight * m_oTau(k) * ((m_solUVals.col(k).transpose() * m_physGradKO).transpose() * (m_solUVals.col(k).transpose() * m_physGradKO));
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
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //N = [N_k, N_omega]
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

protected:

    std::vector<gsMatrix<T> > m_localMat;

    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::getTMEvaluator;

    using NonlinBase::m_kTau;
    using NonlinBase::m_oTau;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_patchIndex;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;

};

// ============================================================= BLOCK NonlinearM ============================================================= //
template <class T>
class uwbTMSUPGBlocksReactionVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlocksReactionVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        NonlinBase::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate solution on element nodes
        m_solUVals = this->m_solU.value(quNodes, m_patchIndex);

        this->m_deg = basisRefs.back().maxDegree(); // front or back???
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local NM_komega = [local NMk, local NMomega] + TMBlend(k^k, omega^k)

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            //--- nonlinearM block
            m_localMat[0].noalias() += weight * m_kTau(k) * (getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k)) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * m_oTau(k) * (getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k)) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
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
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //NM_komega = [NM_k, NM_omega]
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

protected:
    std::vector<gsMatrix<T> > m_localMat;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // members from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::m_solKOmegaVals;
    using NonlinBase::m_patchIndex;
    using NonlinBase::m_solKOmegaGrads;

    using NonlinBase::m_kTau;
    using NonlinBase::m_oTau;
    
    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;
};

// ============================================================= BLOCK Blend ============================================================= //
template <class T>
class uwbTMSUPGBlockBlendVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlockBlendVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        NonlinBase::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate solution on element nodes
        m_solUVals = this->m_solU.value(quNodes, m_patchIndex);

        this->m_deg = basisRefs.back().maxDegree(); // front or back???
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.setZero(numActKOmega, numActKOmega);//local NM_komega = TMBlend(k^k, omega^k)

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            //--- Blend block
            if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
            {
                m_localMat.noalias() -= weight * m_oTau(k) *
                        math::max((getTMEvaluator()->getBlendCoeff(k) / math::max(math::pow(m_solKOmegaVals(1, k), 2), math::pow(10, -15))) *
                        (m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1))), 0.) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                        basisValsKO.col(k).transpose());
            }
            else
            {
                m_localMat.noalias() -= weight * m_oTau(k) *
                        (getTMEvaluator()->getBlendCoeff(k) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15))) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                         (m_solKOmegaGrads[k].row(0) * m_physGradKO));
            }
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
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
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j); //NM_komega = [blend_omega]
                    else // m_mapperTM.is_boundary_index(jj)
                    {
                        const int bb = m_mapperTM.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[1](bb, 0); //rhs = [rhs_omega]
                    }
                }
            }
        }
    }

protected:
    gsMatrix<T> m_localMat;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // members from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::m_solKOmegaVals;
    using NonlinBase::m_patchIndex;
    using NonlinBase::m_solKOmegaGrads;

    using NonlinBase::m_oTau;

    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;
};

// ============================================================= BLOCK NonlinearA ============================================================= //
template <class T>
class uwbTMSUPGBlockANonlinearVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlockANonlinearVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        NonlinBase::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate solution on element nodes
        m_solUVals = this->m_solU.value(quNodes, m_patchIndex);

        this->m_deg = basisRefs.back().maxDegree();

        gsMatrix<T> deriv2KO;
        basisRefs.back().deriv2_into(quNodes, deriv2KO);

        m_basisDataKO.push_back(deriv2KO);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        gsMatrix<T> solActKDiffusionCoeffs, solActOmegaDiffusionCoeffs;
        solActKDiffusionCoeffs.setZero(1, numActKOmega);
        solActOmegaDiffusionCoeffs.setZero(1, numActKOmega);
        /*for (int j = 0; j < numActKOmega; j++)
        {
            solActKDiffusionCoeffs.col(j) = m_kDiffusionCoeffSolField.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();
            solActOmegaDiffusionCoeffs.col(j) = m_omegaDiffusionCoeffSolField.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();
        }*/

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local A_komega

        const gsMatrix<T> & basisGradsKO = m_basisDataKO[1];
        const gsMatrix<T> & basisHessianKO = m_basisDataKO[2];

        const index_t nQuPoints = quWeights.rows();
        //gsMatrix<T> solKDiffusionCoeffGrad;
        //gsMatrix<T> solOmegaDiffusionCoeffGrad;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsKO, m_physGradKO);
            geoEval.transformLaplaceHgrad(k, basisGradsKO, basisHessianKO, m_physLaplacianKO); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n

            //solKDiffusionCoeffGrad.noalias() = solActKDiffusionCoeffs * m_physGradKO.transpose();
            //solOmegaDiffusionCoeffGrad.noalias() = solActOmegaDiffusionCoeffs * m_physGradKO.transpose();

            //--- A block ----
            //laplacian part
            m_localMat[0].noalias() -= weight * m_kTau(k) * getTMEvaluator()->getKDiffusionCoefficient(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * m_physLaplacianKO); //local A_komega = [local A_k(nu_T^k), ...
            m_localMat[1].noalias() -= weight * m_oTau(k) * getTMEvaluator()->getOmegaDiffusionCoefficient(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * m_physLaplacianKO); // ... local A_omega(nu_T^k)]
            // part with the derivative of the diffusion coefficient
            // nuT averaged over elements => derivative is zero
            //m_localMat[0].noalias() -= weight * m_kTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * (solKDiffusionCoeffGrad * m_physGradKO));
            //m_localMat[1].noalias() -= weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * (solOmegaDiffusionCoeffGrad * m_physGradKO));
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
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //A_komega = [A_k, A_omega]
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

protected:
    std::vector<gsMatrix<T> > m_localMat;

    std::vector<gsMatrix<T> > m_kGradDiffusionCoeff;
    std::vector<gsMatrix<T> > m_omegaGradDiffusionCoeff;
    
    gsMatrix<T> m_physLaplacianKO;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_kDiffusionCoeffSolField;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_omegaDiffusionCoeffSolField;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::getTMEvaluator;
    //using NonlinBase::m_solUGrads;
    using NonlinBase::m_patchIndex;

    using NonlinBase::m_kTau;
    using NonlinBase::m_oTau;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;
};

// ============================================================= RHS: [f_k; f_omega] ============================================================= //
template <class T>
class uwbTMSUPGRhsVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{

public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> Base;

    uwbTMSUPGRhsVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) : 
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);
        getTMEvaluator()->evalQuantities_rhsPart();
        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
        //getTMEvaluator()->evalKOmegaDiffusionCoefficient();
    }

protected:
    gsVector<T> m_localRhsK;
    gsVector<T> m_localRhsO;

protected:
    using Base::getTMEvaluator;

};


template <class T>
class uwbTMSUPGRhsFVisitorKOmega : public uwbTMSUPGRhsVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGRhsVisitorKOmega<T> RhsBase;

    uwbTMSUPGRhsFVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        RhsBase(dofMappers, viscosity)
    { }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        RhsBase::evaluate(basisRefs, geoEval, quNodes);

        basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataKO);

        // Evaluate solution on element nodes
        m_solUVals = this->m_solU.value(quNodes, m_patchIndex);

        this->m_deg = basisRefs.back().maxDegree(); // front or back???
    }

    void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        uwbTMSUPGBlockNonlinVisitorKOmega<T>::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        m_localRhsK.setZero(numActKOmega, 1);
        m_localRhsO.setZero(numActKOmega, 1);

        const gsMatrix<T> & m_basisGradsKO = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, m_basisGradsKO, m_physGradKO);

            // Right-hand side
            gsVector<T> physQuadPoint = geoEval.value(k);
            if (m_limitTMProduction && (physQuadPoint(0) < m_productionXPoint)) //-0.16
            {
                m_localRhsK.noalias() += 0. * (m_physGradKO.transpose() * m_solUVals.col(k));
                m_localRhsO.noalias() += 0. * (m_physGradKO.transpose() * m_solUVals.col(k));
            }
            else
            {
                m_localRhsK.noalias() += weight * m_kTau(k) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *  getTMEvaluator()->getRhsK(k));

                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                {
                    m_localRhsO.noalias() += weight * m_oTau(k) *
                                            ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                                             (getTMEvaluator()->getRhsOmega(k) + getTMEvaluator()->getReactionAtRhs(k)));
                }
                else
                {
                    m_localRhsO.noalias() += weight * m_oTau(k) *
                            ((m_physGradKO.transpose() * m_solUVals.col(k)) *  getTMEvaluator()->getRhsOmega(k));
                }
                if (checkSASTypeTM())
                    m_localRhsO.noalias() += weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * getTMEvaluator()->getSourceQSAS(k));
            }
        }
    }

    void localToGlobal(gsMatrix<T> & rhs)
    {
        // Local Dofs to global dofs
        m_mapperTM.localToGlobal(m_activesKO, m_patchIndex, m_activesKO);
        const index_t numActKOmega = m_activesKO.rows();

        for (index_t i = 0; i < numActKOmega; ++i)
        {
            const int ii = m_activesKO(i);
            if (m_mapperTM.is_free_index(ii))
            {
                rhs(ii, 0) += m_localRhsK(i);
                rhs(ii, 1) += m_localRhsO(i);
            }
        }
    }

    bool checkSASTypeTM()
    {
        return (m_evaluatorType == "koSAS" || m_evaluatorType == "koSAS_SS" ||
                m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO");
    }

protected:

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    using uwbTMSUPGBlockNonlinVisitorKOmega<T>::m_kTau;
    using uwbTMSUPGBlockNonlinVisitorKOmega<T>::m_oTau;

    // members from uwbTMSUPGRhsVisitorKOmega
    using RhsBase::m_localRhsK;
    using RhsBase::m_localRhsO;

    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using uwbTMSUPGBlockNonlinVisitorKOmega<T>::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_patchIndex;
    using uwbTMBlockVisitor<T>::m_mapperTM;
    using uwbTMBlockVisitor<T>::m_evaluatorType;
    using uwbTMBlockVisitor<T>::m_limitTMProduction;
    using uwbTMBlockVisitor<T>::m_productionXPoint;

};

// ================================================================================================================================================= //
// ============================================================= SUPG BLOCK TM ============================================================= //
// ================================================================================================================================================= //
template <class T>
class uwbTMSUPGBlocksVisitorKOmega : public uwbTMSUPGBlockNonlinVisitorKOmega<T>
{
public:
    typedef uwbTMSUPGBlockNonlinVisitorKOmega<T> NonlinBase;

    uwbTMSUPGBlocksVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        NonlinBase(dofMappers, viscosity)
    { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        NonlinBase::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate solution on element nodes
        m_solUVals = this->m_solU.value(quNodes, m_patchIndex);

        this->m_deg = basisRefs.back().maxDegree();

        gsMatrix<T> deriv2KO;
        basisRefs.back().deriv2_into(quNodes, deriv2KO);

        m_basisDataKO.push_back(deriv2KO);

        getTMEvaluator()->evalQuantities_rhsPart();
        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        NonlinBase::assemble(geoEval, quWeights);

        const index_t numActKOmega = m_activesKO.rows();

        gsMatrix<T> solActKDiffusionCoeffs, solActOmegaDiffusionCoeffs;
        solActKDiffusionCoeffs.setZero(1, numActKOmega);
        solActOmegaDiffusionCoeffs.setZero(1, numActKOmega);
        /*for (int j = 0; j < numActKOmega; j++)
        {
            solActKDiffusionCoeffs.col(j) = m_kDiffusionCoeffSolField.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();
            solActOmegaDiffusionCoeffs.col(j) = m_omegaDiffusionCoeffSolField.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();
        }*/

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local A_komega

        m_localRhsK.setZero(numActKOmega, 1);
        m_localRhsO.setZero(numActKOmega, 1);

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & basisGradsKO = m_basisDataKO[1];
        const gsMatrix<T> & basisHessianKO = m_basisDataKO[2];

        const index_t nQuPoints = quWeights.rows();
        //gsMatrix<T> solKDiffusionCoeffGrad;
        //gsMatrix<T> solOmegaDiffusionCoeffGrad;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsKO, m_physGradKO);
            geoEval.transformLaplaceHgrad(k, basisGradsKO, basisHessianKO, m_physLaplacianKO); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n

            //--- M SUPG block ----
            m_localMat[0].noalias() += weight * m_kTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());


            //--- A SUPG block ----
            //laplacian part
            m_localMat[0].noalias() -= weight * m_kTau(k) *
                    getTMEvaluator()->getKDiffusionCoefficient(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * m_physLaplacianKO); //local A_komega = [local A_k(nu_T^k), ...
            m_localMat[1].noalias() -= weight * m_oTau(k) *
                    getTMEvaluator()->getOmegaDiffusionCoefficient(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * m_physLaplacianKO); // ... local A_omega(nu_T^k)]
            // part with the derivative of the diffusion coefficient
            // nuT averaged over elements => derivative is zero
            //solKDiffusionCoeffGrad.noalias() = solActKDiffusionCoeffs * m_physGradKO.transpose();
            //solOmegaDiffusionCoeffGrad.noalias() = solActOmegaDiffusionCoeffs * m_physGradKO.transpose();
            //m_localMat[0].noalias() -= weight * m_kTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * (solKDiffusionCoeffGrad * m_physGradKO));
            //m_localMat[1].noalias() -= weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * (solOmegaDiffusionCoeffGrad * m_physGradKO));


            //--- N SUPG block ----
            m_localMat[0].noalias() += weight * m_kTau(k) * ((m_solUVals.col(k).transpose() * m_physGradKO).transpose() * (m_solUVals.col(k).transpose() * m_physGradKO));
            m_localMat[1].noalias() += weight * m_oTau(k) * ((m_solUVals.col(k).transpose() * m_physGradKO).transpose() * (m_solUVals.col(k).transpose() * m_physGradKO));


            //--- reaction SUPG block ----
            m_localMat[0].noalias() += weight * m_kTau(k) * (getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k)) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * m_oTau(k) * (getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k)) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * basisValsKO.col(k).transpose());


            //--- blend SUPG block ----
            if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
            {
                m_localMat[1].noalias() += weight * m_oTau(k) *
                        math::max((-1.)*(getTMEvaluator()->getBlendCoeff(k) / math::max(math::pow(m_solKOmegaVals(1, k), 2), math::pow(10, -15))) *
                        (m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1))), 0.) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                        basisValsKO.col(k).transpose());
            }
            else
            {
                m_localMat[1].noalias() -= weight * m_oTau(k) *
                        (getTMEvaluator()->getBlendCoeff(k) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15))) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                         (m_solKOmegaGrads[k].row(0) * m_physGradKO));
            }


            // SUPG Right-hand side
            gsVector<T> physQuadPoint = geoEval.value(k);
            if (m_limitTMProduction && (physQuadPoint(0) < m_productionXPoint)) //-0.16
            {
                m_localRhsK.noalias() += 0. * (m_physGradKO.transpose() * m_solUVals.col(k));
                m_localRhsO.noalias() += 0. * (m_physGradKO.transpose() * m_solUVals.col(k));
                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                    m_localRhsO.noalias() += weight * m_oTau(k) *
                                            ((m_physGradKO.transpose() * m_solUVals.col(k)) * getTMEvaluator()->getReactionAtRhs(k));
            }
            else
            {
                m_localRhsK.noalias() += weight * m_kTau(k) *
                        ((m_physGradKO.transpose() * m_solUVals.col(k)) *  getTMEvaluator()->getRhsK(k));

                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                {
                    m_localRhsO.noalias() += weight * m_oTau(k) *
                                            ((m_physGradKO.transpose() * m_solUVals.col(k)) *
                                             (getTMEvaluator()->getRhsOmega(k) + getTMEvaluator()->getReactionAtRhs(k)));
                }
                else
                {
                    m_localRhsO.noalias() += weight * m_oTau(k) *
                            ((m_physGradKO.transpose() * m_solUVals.col(k)) *  getTMEvaluator()->getRhsOmega(k));
                }
                if (checkSASTypeTM())
                    m_localRhsO.noalias() += weight * m_oTau(k) * ((m_physGradKO.transpose() * m_solUVals.col(k)) * getTMEvaluator()->getSourceQSAS(k));
            }
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
                            sysBlock.coeffRef(ii, jj + s * KOmegaSize) += m_localMat[s](i, j); //A_komega = [A_k, A_omega]
                    }
                    else // m_mapperTM.is_boundary_index(jj)
                    {
                        const int bb = m_mapperTM.global_to_bindex(jj);
                        for (index_t s = 0; s != m_numVar; ++s)
                            rhs(ii, s) -= m_localMat[s](i, j) * eliminatedDofs[s](bb, 0); //rhs = [rhs_k, rhs_omega]
                    }
                }
                rhs(ii, 0) += m_localRhsK(i);
                rhs(ii, 1) += m_localRhsO(i);
            }
        }
    }

    bool checkSASTypeTM()
    {
        return (m_evaluatorType == "koSAS" || m_evaluatorType == "koSAS_SS" ||
                m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO");
    }

    /*virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");

        if (m_evaluatorType == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else // (m_evaluatorType == "koSST")
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }*/

protected:
    std::vector<gsMatrix<T> > m_localMat;
    gsVector<T> m_localRhsK;
    gsVector<T> m_localRhsO;

    std::vector<gsMatrix<T> > m_kGradDiffusionCoeff;
    std::vector<gsMatrix<T> > m_omegaGradDiffusionCoeff;

    gsMatrix<T> m_physLaplacianKO;

protected:
    // member functions from uwbTMSUPGBlockNonlinVisitorKOmega
    using uwbTMSUPGBlockNonlinVisitorKOmega<T>::getTMEvaluator;

    // members from uwbTMSUPGBlockNonlinVisitorKOmega
    using NonlinBase::m_kTau;
    using NonlinBase::m_oTau;
    using NonlinBase::m_solKOmegaVals;
    using NonlinBase::m_solKOmegaGrads;

    // members from uwbTMSUPGBlockVisitorKOmega
    using uwbTMSUPGBlockVisitorKOmega<T>::m_solUVals;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_activesKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_basisDataKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_physGradKO;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_kDiffusionCoeffSolField;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_omegaDiffusionCoeffSolField;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_patchIndex;
    using uwbTMSUPGBlockVisitorKOmega<T>::m_pStabEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_mapperTM;
    using uwbTMBlockVisitor<T>::m_evaluatorType;
    using uwbTMBlockVisitor<T>::m_limitTMProduction;
    using uwbTMBlockVisitor<T>::m_productionXPoint;
};


} // namespace gismo
