/** @file uwbTMBlockVisitorsKOmega.h

    Author(s): E. Turnerova

	k-omega Wilcox 2006 - LRN
	https://turbmodels.larc.nasa.gov/wilcox.html

*/

#pragma once
#include "uwbINSBlockVisitors.h"
#include "uwbTMEvaluators.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbTMBlockVisitor : public uwbVisitorBase<T>
{
protected:

    const T m_viscosity;
    const gsDofMapper & m_mapperTM;
    gsField<T> m_solU;
    gsField<T> m_solTM;
    index_t m_dim; // number of unknown variables of turbulent model
    index_t m_numVar; // number of unknown variables of turbulent model - i.e. [k, omega]
    int m_patchIndex;
    gsField<T> m_solPoisson;
    gsMultiPatch<T> m_patchesPoisson;
    gsMultiBasis<T> m_basesPoisson;

    std::string m_evaluatorType;
    uwbTMEvaluator<T>* m_pEvaluator;

    bool m_limitTMProduction;
    T m_productionXPoint;

    bool m_bOFreactionTreatment;

public:
	uwbTMBlockVisitor(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        m_viscosity(viscosity),
        m_mapperTM(dofMappers.back()),
        m_limitTMProduction(false),
        m_productionXPoint(0.)
    { 
        m_evaluatorType = "";
        m_pEvaluator = NULL;
        m_bOFreactionTreatment = false;
    }

    ~uwbTMBlockVisitor()
        {
            if (m_pEvaluator)
            {
                delete m_pEvaluator;
                m_pEvaluator = NULL;
            }
        }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions & options,
        gsQuadRule<T>    & rule,
        unsigned         & evFlags)
    {
        const gsBasis<T>& basis = basisRefs.back();

        m_dim = basisRefs.front().dim();
        m_numVar = 2;
        m_patchIndex = patchIndex;

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options.quA, options.quB);

        if (getTMEvaluator() != NULL)
            getTMEvaluator()->initialize(m_viscosity, rule.numNodes());

        GISMO_ASSERT(this->m_bSolutionSet, "No solution set in the visitor.");

        initializeSpecific(evFlags);
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>    & rule,
        unsigned         & evFlags)
    {
        const gsBasis<T>& basis = basisRefs.back();

        m_dim = basisRefs.front().dim();
        m_numVar = 2;
        m_patchIndex = patchIndex;

        gsVector<index_t> numQuadNodes(m_dim);
        for (int i = 0; i < m_dim; ++i)
            numQuadNodes[i] = basis.maxDegree() + 1; 

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);

        if (getTMEvaluator() != NULL)
            getTMEvaluator()->initialize(m_viscosity, rule.numNodes());

        GISMO_ASSERT(this->m_bSolutionSet, "No solution set in the visitor.");

        initializeSpecific(evFlags);
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    { GISMO_NO_IMPLEMENTATION }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    virtual inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_solU = solutions.front();
        m_solTM = solutions.back();
        this->m_bSolutionSet = true;
    }

    void createTMEvaluator(std::string evaluatorType)
    {
        m_evaluatorType = evaluatorType;
        if (m_evaluatorType == "koWilcoxLRN")
        {
            m_pEvaluator = new uwbTMEvaluatorKOmegaWilcoxLRN<T>();
            m_bOFreactionTreatment = false;
        }
        else
        {
            m_pEvaluator = new uwbTMEvaluatorKOmegaSST<T>();
            m_bOFreactionTreatment = true;
        }
    }

    bool checkWallDistanceBasedTM()
    {
        return (m_evaluatorType == "koSST" || m_evaluatorType == "koSSTMenter2009" || m_evaluatorType == "koSAS"
                || m_evaluatorType == "koSAS_SS" || m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO");
    }

    void setPoissonSolution(gsField<T>& solPoisson, const gsMultiPatch<T>& patchesPoisson, const gsMultiBasis<T>& basesPoisson)
    { m_solPoisson = solPoisson; m_patchesPoisson = patchesPoisson; m_basesPoisson = basesPoisson; }

    void setTMProductionLimiter(bool limitTMProduction, T productionXPoint)
    { m_limitTMProduction = limitTMProduction; m_productionXPoint = productionXPoint; }

    virtual uwbTMEvaluator<T>* getTMEvaluator() { return m_pEvaluator; }

    virtual inline void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }

};

// ============================================================= PARENT nonlinear k-o ============================================================= //
template <class T>
class uwbTMBlockNonlinVisitorKOmega : public uwbTMBlockVisitor<T>
{
public:
    typedef uwbTMBlockVisitor<T> Base;

public:
    uwbTMBlockNonlinVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) : Base(dofMappers, viscosity) { }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        const index_t numActKOmega = m_activesKO.rows(); // number of active basis functions
        const index_t numActU = m_activesU.rows();

        // Evaluate basis functions on element  nodes
        basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataKO);
        basisRefs.front().deriv_into(quNodes, m_basisGradsU);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solKOmegaVals = m_solTM.value(quNodes, m_patchIndex);

        m_solActUCoeffs.setZero(m_dim, numActU);
        for (int j = 0; j < numActU; j++) {
            m_solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();
        }
        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        index_t nQuPoints = quNodes.cols();
        std::vector<gsMatrix<T> > solUGrads(nQuPoints);
        std::vector<gsMatrix<T> > solKOmegaGrads(nQuPoints);
        gsMatrix<T> physGrad_u;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, m_basisGradsU, physGrad_u);
            solUGrads[k].noalias() = m_solActUCoeffs * physGrad_u.transpose();
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);
            solKOmegaGrads[k].noalias() = solActKOmegaCoeffs * m_physGradKO.transpose();
        }

        getTMEvaluator()->initAtElement(solUGrads, m_solKOmegaVals, solKOmegaGrads);
        getTMEvaluator()->setKOmegaVariant(m_evaluatorType);

        if (this->checkWallDistanceBasedTM())
        {
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
        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
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
    // Basis values
    gsMatrix<index_t> m_activesKO;
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataKO;
    gsMatrix<T> m_solKOmegaVals;
    gsMatrix<T> m_solPoissonVals;
    gsMatrix<T> m_solActUCoeffs;

    gsMatrix<T> m_basisGradsU;
    gsMatrix<T> m_physGradKO;

    // members from uwbTMBlockVisitor
    using Base::m_solU;
    using Base::m_solTM;
    using Base::m_dim; // number of unknown variables of turbulent model
    using Base::m_patchIndex;
    using Base::m_evaluatorType;
    using Base::m_pEvaluator;
    using Base::m_numVar;
    using Base::m_solPoisson;
    using Base::m_patchesPoisson;
    using Base::m_basesPoisson;
};


// ============================================================= BLOCK M ============================================================= //
template <class T>
class uwbTMBlockMVisitorKOmega : public uwbTMBlockVisitor<T>
{
public:
    typedef uwbTMBlockVisitor<T> Base;

public:

    uwbTMBlockMVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);
        // Evaluate basis functions on element  nodes

        basisRefs.back().eval_into(quNodes, m_basisValsKO);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);
    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.setZero(numActKOmega, numActKOmega);//local M_k = M_omega

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            m_localMat.noalias() += weight * (m_basisValsKO.col(k) * m_basisValsKO.col(k).transpose());
        }
    }

    virtual inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
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
                    {
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j); //M_k = M_omega
                    }
                    else // m_mapperTM.is_boundary_index(jj)
                    {
                        const int bb = m_mapperTM.global_to_bindex(jj);
                        for (index_t s = 0; s != m_numVar; ++s)
                            rhs(ii, s) -= m_localMat(i, j) * eliminatedDofs[s](bb, 0); //rhs = [rhs_k, rhs_omega]
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesKO;
    gsMatrix<T> m_basisValsKO;

    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_mapperTM;
    using Base::m_numVar; // number of unknown variables of turbulent model - i.e. [k, omega]
    using Base::m_patchIndex;
};


// ============================================================= BLOCK N_TMkomega ============================================================= //
template <class T>
class uwbTMBlocksNVisitorKOmega : public uwbTMBlockVisitor<T>
{

public:
    typedef uwbTMBlockVisitor<T> Base;

public:

    uwbTMBlocksNVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);

        // Evaluate basis functions on element  nodes
        basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataKO);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);
    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.setZero(numActKOmega, numActKOmega);//local N_k(u^n) = local N_omega(u^n)

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            m_localMat.noalias() += weight * (basisValsKO.col(k) * (m_solUVals.col(k).transpose() * m_physGradKO));
        }
    }

    virtual inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
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
                    {
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j); //N_k(u^n) = N_omega(u^n)
                    }
                    else // m_mapperTM.is_boundary_index(jj)
                    {
                        const int bb = m_mapperTM.global_to_bindex(jj);
                        for (index_t s = 0; s != m_numVar; ++s)
                            rhs(ii, s) -= m_localMat(i, j) * eliminatedDofs[s](bb, 0); //rhs = [rhs_k, rhs_omega]
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesKO;
    std::vector<gsMatrix<T> > m_basisDataKO;
    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_physGradKO;
    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_mapperTM;
    using Base::m_solU;
    using Base::m_numVar; // number of unknown variables of turbulent model - i.e. [k, omega]
    using Base::m_patchIndex;
};


// ============================================================= BLOCK NonlinearM ============================================================= //
template <class T>
class uwbTMBlocksReactionVisitorKOmega : public uwbTMBlockNonlinVisitorKOmega<T>
{

public:
    typedef uwbTMBlockNonlinVisitorKOmega<T> Base;

public:

    uwbTMBlocksReactionVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local NM_komega = [local NMk, local NMomega] + TMBlend(k^k, omega^k)

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        gsMatrix<T> solKOmegaGrads;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            solKOmegaGrads.noalias() = solActKOmegaCoeffs * m_physGradKO.transpose();

            //--- nonlinearM block
            m_localMat[0].noalias() += weight * (getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k)) * (basisValsKO.col(k) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * (getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k)) * (basisValsKO.col(k) * basisValsKO.col(k).transpose());
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

protected:
    // members from uwbTMBlockNonlinVisitorKOmega
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_solKOmegaVals;
    using Base::m_physGradKO;

    // member functions from uwbTMBlockNonlinVisitorKOmega
    using Base::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_mapperTM;
    using uwbTMBlockVisitor<T>::m_solTM;
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_patchIndex;
};

// ============================================================= BLOCK Blend ============================================================= //
template <class T>
class uwbTMBlockBlendVisitorKOmega : public uwbTMBlockNonlinVisitorKOmega<T>
{

public:
    typedef uwbTMBlockNonlinVisitorKOmega<T> Base;

public:

    uwbTMBlockBlendVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        m_localMat.setZero(numActKOmega, numActKOmega);//local TMBlend(k^k, omega^k)

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        gsMatrix<T> solKOmegaGrads;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            solKOmegaGrads.noalias() = solActKOmegaCoeffs * m_physGradKO.transpose();

            //--- Blend block
            if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
            {
                m_localMat.noalias() -= weight *
                        math::max((getTMEvaluator()->getBlendCoeff(k) / math::max(math::pow(m_solKOmegaVals(1, k), 2), math::pow(10, -15))) *
                        (solKOmegaGrads.row(0).dot(solKOmegaGrads.row(1))), 0.) *
                        (basisValsKO.col(k) * basisValsKO.col(k).transpose());
            }
            else
            {
                m_localMat.noalias() -= weight *
                        (getTMEvaluator()->getBlendCoeff(k) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15))) *
                        (basisValsKO.col(k) * (solKOmegaGrads.row(0) * m_physGradKO));
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
                            sysBlock.coeffRef(ii, jj) += m_localMat(i, j); // blending term in omega eq.
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

protected:
    // members from uwbTMBlockNonlinVisitorKOmega
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_solKOmegaVals;
    using Base::m_physGradKO;

    // member functions from uwbTMBlockNonlinVisitorKOmega
    using Base::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_mapperTM;
    using uwbTMBlockVisitor<T>::m_solTM;
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_patchIndex;
};

// ============================================================= BLOCK NonlinearA ============================================================= //
template <class T>
class uwbTMBlockANonlinearVisitorKOmega : public uwbTMBlockNonlinVisitorKOmega<T>
{

public:
    typedef uwbTMBlockNonlinVisitorKOmega<T> Base;

public:

    uwbTMBlockANonlinearVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity){ }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local A_komega

        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            //--- A block
            m_localMat[0].noalias() += weight * getTMEvaluator()->getKDiffusionCoefficient(k) * (m_physGradKO.transpose() * m_physGradKO); //local A_komega = [local A_k(nu_T^k), ...
            m_localMat[1].noalias() += weight * getTMEvaluator()->getOmegaDiffusionCoefficient(k) * (m_physGradKO.transpose() * m_physGradKO); // ... local A_omega(nu_T^k)]
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

protected:
    // members from uwbTMBlockNonlinVisitorKOmega
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_physGradKO;

    // member functions from uwbTMBlockNonlinVisitorKOmega
    using Base::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_mapperTM;
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_patchIndex;
};

// ============================================================= RHS: [f_k; f_omega] ============================================================= //
template <class T>
class uwbTMRhsVisitorKOmega : public uwbTMBlockVisitor<T>
{

public:
    typedef uwbTMBlockVisitor<T> Base;

    uwbTMRhsVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity){ }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), m_activesKO);

        const index_t numActKOmega = m_activesKO.rows(); // number of active basis functions
        const index_t numActU = m_activesU.rows();

        std::vector<gsMatrix<T> > basisDataKO;
        basisRefs.back().evalAllDers_into(quNodes, 1, basisDataKO);
        basisRefs.front().deriv_into(quNodes, m_basisGradsU);

        m_basisValsKO = basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = basisDataKO[1];

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            m_solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        m_solKOmegaVals.noalias() = m_solActKOmegaCoeffs * m_basisValsKO;

        m_solActUCoeffs.setZero(m_dim, numActU);
        for (int j = 0; j < numActU; j++) {
            m_solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();
        }

        index_t nQuPoints = quNodes.cols();
        gsMatrix<T> physGradU;
        m_solUGrads.resize(nQuPoints);
        m_solKOmegaGrads.resize(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, m_basisGradsU, physGradU);
            m_solUGrads[k].noalias() = m_solActUCoeffs * physGradU.transpose();
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);
            m_solKOmegaGrads[k].noalias() = m_solActKOmegaCoeffs * m_physGradKO.transpose();
        }

        getTMEvaluator()->initAtElement(m_solUGrads, m_solKOmegaVals, m_solKOmegaGrads);
        getTMEvaluator()->setKOmegaVariant(m_evaluatorType);

        if (this->checkWallDistanceBasedTM())
        {
            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = m_patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            m_basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            m_basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            m_solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                m_solActPoissonCoeffs.col(j) = m_solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            m_solPoissonVals.noalias() = m_solPoisson.value(quNodes, m_patchIndex);

            m_solPoissonGrads.resize(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                m_solPoissonGrads[k].noalias() = m_solActPoissonCoeffs * physGradPoisson.transpose();
            }

            getTMEvaluator()->evalWallDistance(m_solPoissonVals, m_solPoissonGrads);
        }

        if (checkSASTypeTM())
        {
            std::vector<gsMatrix<T> > solULaplaces(nQuPoints);
            gsMatrix<T> basisHessian, physLaplacianU;
            basisRefs.front().deriv2_into(quNodes, basisHessian);


            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, m_basisGradsU, basisHessian, physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                solULaplaces[k].noalias() = m_solActUCoeffs * physLaplacianU.transpose();
            }
            getTMEvaluator()->setULaplacian(solULaplaces);
        }

        getTMEvaluator()->evalQuantities_rhsPart();
        //-----
        getTMEvaluator()->evalQuantities_nonlinearBlocksPart();
        //-----
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localRhsK.setZero(numActKOmega, 1);
        m_localRhsO.setZero(numActKOmega, 1);

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // Right-hand side
            gsVector<T> physQuadPoint = geoEval.value(k);
            if (m_limitTMProduction && (physQuadPoint(0) < m_productionXPoint))
            {
                m_localRhsK.noalias() += 0. * m_basisValsKO.col(k);
                m_localRhsO.noalias() += 0. * m_basisValsKO.col(k);
            }
            else
            {
                m_localRhsK.noalias() += weight * (m_basisValsKO.col(k) *  getTMEvaluator()->getRhsK(k));
                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                    m_localRhsO.noalias() += weight * (m_basisValsKO.col(k) *  (getTMEvaluator()->getRhsOmega(k) + getTMEvaluator()->getReactionAtRhs(k)));
                else
                    m_localRhsO.noalias() += weight * (m_basisValsKO.col(k) *  getTMEvaluator()->getRhsOmega(k));
                if (checkSASTypeTM())
                    m_localRhsO.noalias() += weight * (m_basisValsKO.col(k) *  getTMEvaluator()->getSourceQSAS(k));
            }
        }
    }

    inline void localToGlobal(gsMatrix<T> & rhs)
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

    virtual uwbTMEvaluatorKOmega<T>* getTMEvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");

        if (m_evaluatorType == "koWilcoxLRN")
            return dynamic_cast<uwbTMEvaluatorKOmegaWilcoxLRN<T>*>(m_pEvaluator);
        else // (m_evaluatorType == "koSST") or its variants
            return dynamic_cast<uwbTMEvaluatorKOmegaSST<T>*>(m_pEvaluator);
    }

    bool checkSASTypeTM()
    {
        return (m_evaluatorType == "koSAS" || m_evaluatorType == "koSAS_SS" ||
                m_evaluatorType == "koSAS_SO" || m_evaluatorType == "koSAS_OO");
    }

protected:
        // Basis values
        gsMatrix<index_t> m_activesKO;
        gsMatrix<index_t> m_activesU;
        gsMatrix<T> m_basisValsKO;
        gsMatrix<T> m_solActKOmegaCoeffs;
        gsMatrix<T> m_solActPoissonCoeffs;
        gsMatrix<T> m_solActUCoeffs;
        gsMatrix<T> m_solKOmegaVals;
        gsMatrix<T> m_solPoissonVals;
        gsMatrix<T> m_basisGradsU;
        std::vector<gsMatrix<T> > m_solUGrads;
        std::vector<gsMatrix<T> > m_solKOmegaGrads;
        std::vector<gsMatrix<T> > m_solPoissonGrads;

        gsMatrix<T> m_physGradKO;
        gsVector<T> m_localRhsK;
        gsVector<T> m_localRhsO;

protected:
    using Base::m_mapperTM;
    using Base::m_solU;
    using Base::m_solTM;
    using Base::m_solPoisson;
    using Base::m_patchesPoisson;
    using Base::m_basesPoisson;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;
    using Base::m_evaluatorType;
    using Base::m_pEvaluator;
    using Base::m_numVar;
    using Base::m_limitTMProduction;
    using Base::m_productionXPoint;

};


// =============================================================================================================================================== //
// ============================================================= TM visitor ============================================================= //
// =============================================================================================================================================== //

template <class T>
class uwbTMBlockVisitorKOmega : public uwbTMBlockNonlinVisitorKOmega<T>
{

public:
    typedef uwbTMBlockNonlinVisitorKOmega<T> Base;

public:

    uwbTMBlockVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity){ }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate solution on element nodes
        //m_solUVals = m_solU.value(quNodes, m_patchIndex);

        if (checkSASTypeTM())
        {
            index_t nQuPoints = quNodes.cols();
            std::vector<gsMatrix<T> > solULaplaces(nQuPoints);
            gsMatrix<T> basisHessian, physLaplacianU;
            basisRefs.front().deriv2_into(quNodes, basisHessian);


            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, m_basisGradsU, basisHessian, physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                solULaplaces[k].noalias() = m_solActUCoeffs * physLaplacianU.transpose();
            }
            getTMEvaluator()->setULaplacian(solULaplaces);
        }

        getTMEvaluator()->evalQuantities_rhsPart();
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        gsMatrix<T> solActKOmegaCoeffs;
        solActKOmegaCoeffs.setZero(m_numVar, numActKOmega);
        for (int j = 0; j < numActKOmega; j++)
            solActKOmegaCoeffs.col(j) = m_solTM.coefficientVector(m_patchIndex).row(m_activesKO(j)).transpose();

        m_localMat.resize(m_numVar);
        for (index_t s = 0; s != m_numVar; ++s)
            m_localMat[s].setZero(numActKOmega, numActKOmega);//local A_komega//local TMBlend(k^k, omega^k)//local NM_komega = [local NMk, local NMomega]

        m_localRhsK.setZero(numActKOmega, 1);
        m_localRhsO.setZero(numActKOmega, 1);

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        gsMatrix<T> solKOmegaGrads;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);
            solKOmegaGrads.noalias() = solActKOmegaCoeffs * m_physGradKO.transpose();

            //--- A block
            m_localMat[0].noalias() += weight * getTMEvaluator()->getKDiffusionCoefficient(k) * (m_physGradKO.transpose() * m_physGradKO); //local A_komega = [local A_k(nu_T^k), ...
            m_localMat[1].noalias() += weight * getTMEvaluator()->getOmegaDiffusionCoefficient(k) * (m_physGradKO.transpose() * m_physGradKO); // ... local A_omega(nu_T^k)]

            //--- Blend block
            if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
            {
                m_localMat[1].noalias() += weight *
                        math::max((-1.)*(getTMEvaluator()->getBlendCoeff(k) / math::max(math::pow(m_solKOmegaVals(1, k), 2), math::pow(10, -15))) *
                        (solKOmegaGrads.row(0).dot(solKOmegaGrads.row(1))), 0.) *
                        (basisValsKO.col(k) * basisValsKO.col(k).transpose());
            }
            else
            {
                m_localMat[1].noalias() -= weight *
                        (getTMEvaluator()->getBlendCoeff(k) / math::max(m_solKOmegaVals(1, k), math::pow(10, -15))) *
                        (basisValsKO.col(k) * (solKOmegaGrads.row(0) * m_physGradKO));
            }

            //--- nonlinearM block
            m_localMat[0].noalias() += weight * (getTMEvaluator()->getBetaStar(k) * m_solKOmegaVals(1, k)) * (basisValsKO.col(k) * basisValsKO.col(k).transpose());
            m_localMat[1].noalias() += weight * (getTMEvaluator()->getBeta(k) * m_solKOmegaVals(1, k)) * (basisValsKO.col(k) * basisValsKO.col(k).transpose());

            //--- N_TM block - assembled only once per time step
            //m_localMat[0].noalias() += weight * (basisValsKO.col(k) * (m_solUVals.col(k).transpose() * m_physGradKO));
            //m_localMat[1].noalias() += weight * (basisValsKO.col(k) * (m_solUVals.col(k).transpose() * m_physGradKO));

            // Right-hand side
            gsVector<T> physQuadPoint = geoEval.value(k);
            if (m_limitTMProduction && (physQuadPoint(0) < m_productionXPoint))
            {
                m_localRhsK.noalias() += 0. * basisValsKO.col(k);
                m_localRhsO.noalias() += 0. * basisValsKO.col(k);
                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                    m_localRhsO.noalias() += weight * (basisValsKO.col(k) *  getTMEvaluator()->getReactionAtRhs(k));
            }
            else
            {
                m_localRhsK.noalias() += weight * (basisValsKO.col(k) *  getTMEvaluator()->getRhsK(k));
                if (uwbTMBlockVisitor<T>::m_bOFreactionTreatment) //OF treatment
                    m_localRhsO.noalias() += weight * (basisValsKO.col(k) *  (getTMEvaluator()->getRhsOmega(k) + getTMEvaluator()->getReactionAtRhs(k)));
                else
                    m_localRhsO.noalias() += weight * (basisValsKO.col(k) *  getTMEvaluator()->getRhsOmega(k));
                if (checkSASTypeTM())
                    m_localRhsO.noalias() += weight * (basisValsKO.col(k) *  getTMEvaluator()->getSourceQSAS(k));
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

protected:
    std::vector<gsMatrix<T> > m_localMat;
    gsVector<T> m_localRhsK;
    gsVector<T> m_localRhsO;

    //gsMatrix<T> m_solUVals;

protected:
    // members from uwbTMBlockNonlinVisitorKOmega
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_physGradKO;
    using Base::m_solKOmegaVals;
    using Base::m_solActUCoeffs;
    using Base::m_basisGradsU;

    // member functions from uwbTMBlockNonlinVisitorKOmega
    using Base::getTMEvaluator;

    // members from uwbTMBlockVisitor
    using uwbTMBlockVisitor<T>::m_mapperTM;
     using uwbTMBlockVisitor<T>::m_solTM;
    using uwbTMBlockVisitor<T>::m_numVar;
    using uwbTMBlockVisitor<T>::m_patchIndex;
    //using uwbTMBlockVisitor<T>::m_solU;
    using uwbTMBlockVisitor<T>::m_evaluatorType;
    using uwbTMBlockVisitor<T>::m_limitTMProduction;
    using uwbTMBlockVisitor<T>::m_productionXPoint;
};



} // namespace gismo

