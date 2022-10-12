/** @file uwbResiduumEvaluators.h

    Author(s): E. Turnerova
*/

#pragma once

namespace gismo
{
// ========================================================== SUPER CLASS ========================================================== //
template <class T>
class uwbEvaluatorBase
{
public:
    uwbEvaluatorBase()
    {
        m_bSolutionSet = false;
        m_elementLength = 0.;
    }

    virtual void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { GISMO_NO_IMPLEMENTATION }

    void setElementLength(T elementLength) { m_elementLength = elementLength; }

protected: 
    bool m_bSolutionSet;
    T m_elementLength;
};

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbINSResiduumEvaluator : public uwbEvaluatorBase<T>
{
public:
    uwbINSResiduumEvaluator(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        m_viscosity(viscosity)//,
        //m_Umap(dofMappers.front()),
        //m_Pmap(dofMappers.back())
    { }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basis.dim();
        m_patchIndex = patchIndex;

        // Setup Quadrature
        rule = gsGaussRule<T>(basis, options.quA, options.quB);// harmless slicing occurs here

        initializeSpecific(evFlags);

        m_diffusionCoeff.setZero(1, rule.numNodes());
        for (int i = 0; i < rule.numNodes(); i++)
            m_diffusionCoeff(0, i) = m_viscosity;
    }

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basis.dim();
        m_patchIndex = patchIndex;

        gsVector<index_t> numQuadNodes(m_dim);
        for (int i = 0; i < m_dim; ++i) // to do: improve
            numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        initializeSpecific(evFlags);

        m_diffusionCoeff.setZero(1, rule.numNodes());
        for (int i = 0; i < rule.numNodes(); i++)
            m_diffusionCoeff(0, i) = m_viscosity;
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        gsMatrix<index_t> activesP;
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), activesP);

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);
        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);

        gsMatrix<T> basisGradsP;
        basisRefs.back().deriv_into(quNodes, basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate velocity gradients
        index_t numActU = m_activesU.rows();
        const index_t numActP = activesP.rows();
        index_t nQuPoints = quNodes.cols();

        m_solActUCoeffs.resize(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            m_solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        gsMatrix<T> solActPCoeffs;
        solActPCoeffs.setZero(1, numActP);
        for (int j = 0; j < numActP; j++)
            solActPCoeffs.col(j) = m_solP.coefficientVector(m_patchIndex).row(activesP(j)).transpose();

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        //gsMatrix<T> solUValsOld = m_solUVals;
        /*if (m_bUnsteady)
        {
            GISMO_ASSERT(m_bOldSolSet, "Old time step solution not set in Crosswind visitor.");
            solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);
        }*/

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        m_solUGrads.resize(nQuPoints);
        m_solPGrads.resize(nQuPoints);
        gsMatrix<T> physGradU, physGradP;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsU, physGradU);
            m_solUGrads[k].noalias() = m_solActUCoeffs * physGradU.transpose();
            geoEval.transformGradients(k, basisGradsP, physGradP);
            m_solPGrads[k].noalias() = solActPCoeffs * physGradP.transpose();
        }
    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_solU = solutions.front();
        m_solP = solutions.back();
        this->m_bSolutionSet = true;
    }

    /*void setOldSolutionField(gsField<T>& solution)
    {
        m_oldSolU = solution;
        m_bOldSolSet = true;
    }*/

    T getElementResiduum() const { return m_elemResiduum; }

protected:
    virtual inline void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }

protected:

    T m_elemResiduum;
    const T m_viscosity;
    gsField<T> m_solU, m_solP;
    index_t m_dim; // Velocity vector dimension
    int m_patchIndex;
    gsMatrix<T> m_diffusionCoeff;
    //gsField<T> m_oldSolU
    //bool m_bOldSolSet;
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    std::vector<gsMatrix<T> > m_solUGrads, m_solPGrads;

    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_solActUCoeffs;
};

// ============================================================= Residuum steady NS ===================================================== //
template <class T>
class uwbINSsteadyResiduumEvaluatorL2norm : public uwbINSResiduumEvaluator<T>
{

public:
    typedef uwbINSResiduumEvaluator<T> Base;

public:

    uwbINSsteadyResiduumEvaluatorL2norm(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity and pressure solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> physLaplacian;

        m_residuumMomentumEq.resize(m_dim);
        index_t nQuPoints = quNodes.cols();
        for (int var = 0; var < m_dim; var++)
            m_residuumMomentumEq[var].setZero(nQuPoints);
        m_residuumContinuityEq.setZero(nQuPoints);

            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = m_solActUCoeffs * physLaplacian.transpose();

                for (int var = 0; var < m_dim; var++)
                {
                    m_residuumMomentumEq[var](k) = m_solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                                  - m_viscosity * solLaplacian(var, 0) //diffusion
                                                  + m_solPGrads[k](0, var); //pressure term
                                                  //- f; //source term

                    //if (m_bUnsteady)
                    //    m_residuumMomentumEq[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));

                    m_residuumContinuityEq(k) += m_solUGrads[k](var, var);
                }
            }
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t nQuPoints = quWeights.rows();

        T elResMomentum = 0.;
        T elResContinuity = 0.;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            for (index_t var = 0; var != m_dim; ++var)
                elResMomentum += weight * math::pow(m_residuumMomentumEq[var](k), 2);

            elResContinuity += weight * math::pow(m_residuumContinuityEq(k), 2);
        }
        m_elemResiduum = (math::sqrt(elResMomentum) + math::sqrt(elResContinuity)) * m_elementLength;
    }

protected:
    std::vector<gsVector<T> > m_residuumMomentumEq;
    gsVector<T> m_residuumContinuityEq;

protected:
    using Base::m_elemResiduum;
    using Base::m_viscosity;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_solActUCoeffs;
    using Base::m_solUGrads;
    using Base::m_solPGrads;
    using Base::m_basisDataU;
    using Base::m_solUVals;

    using uwbEvaluatorBase<T>::m_elementLength;

};

} // namespace gismo

