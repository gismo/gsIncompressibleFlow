/** @file uwbINSPSPGBlockVisitors.h

Author(s): E. Turnerova
*/

#pragma once
#include "uwbINSBlockVisitors.h"
#include "uwbStabilizationEvaluator.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbINSPSPGBlockVisitors : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:
    uwbINSPSPGBlockVisitors(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_degU(0), m_pspgResidualType(0), m_tauStabType(0),
        m_timeStep(0.), m_r(2), m_bOldSolSet(false), m_bPSolSet(false)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
    }

    ~uwbINSPSPGBlockVisitors()
    {
        if (m_pStabEvaluator)
        {
            delete m_pStabEvaluator;
            m_pStabEvaluator = NULL;
        }
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
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
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), m_activesP);

        const index_t numActU = m_activesU.rows();
        const index_t numActP = m_activesP.rows();

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);
        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);
        m_basisDataU.push_back(deriv2_u);

        basisRefs.back().deriv_into(quNodes, m_basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);
        gsMatrix<T> solUValsOld = m_solUVals;
        if (m_bUnsteady)
        {
            GISMO_ASSERT(m_bOldSolSet, "Old time step solution not set in Crosswind visitor.");
            solUValsOld = m_oldSolU.value(quNodes, m_patchIndex);
        }

        m_degU = basisRefs.front().maxDegree();

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        gsMatrix<T> solActUCoeffs;
        solActUCoeffs.setZero(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        gsMatrix<T> solActPCoeffs;
        solActPCoeffs.setZero(1, numActP);
        GISMO_ASSERT(m_bPSolSet, "Pressure solution not set in Crosswind visitor.");
        for (int j = 0; j < numActP; j++)
            solActPCoeffs.col(j) = m_solP.coefficientVector(m_patchIndex).row(m_activesP(j)).transpose();

        index_t nQuPoints = quNodes.cols();

        std::vector<gsMatrix<T> > solUGrads, solPGrads;
        solUGrads.resize(nQuPoints);
        solPGrads.resize(nQuPoints);
        gsMatrix<T> physGradU, physGradP;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsU, physGradU);
            solUGrads[k].noalias() = solActUCoeffs * physGradU.transpose();
            geoEval.transformGradients(k, m_basisGradsP, physGradP);
            solPGrads[k].noalias() = solActPCoeffs * physGradP.transpose();
        }

        gsMatrix<T> physLaplacian;
        m_residual.resize(m_dim);
        for (int var = 0; var < m_dim; var++)
            m_residual[var].setZero(nQuPoints);

        switch (m_pspgResidualType)
        {
        case 0:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActUCoeffs * physLaplacian.transpose();

                for (int var = 0; var < m_dim; var++)
                {
                    m_residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                     - m_viscosity * solLaplacian(var, 0) //diffusion
                                     + solPGrads[k](0, var); //pressure term

                    if (m_bUnsteady)
                        m_residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));
                }
            }
            break;
        case 1:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                for (int var = 0; var < m_dim; var++)
                    m_residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
            break;
        case 2:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                for (int var = 0; var < m_dim; var++)
                {
                    m_residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
                    if (m_bUnsteady)
                        m_residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));
                }
            }
            break;
        default:
            gsWarn << "Wrong residualType choosen for crosswind stabilization parameter. Choosen default case 1.\n";
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                for (int var = 0; var < m_dim; var++)
                    m_residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
        }

        /*gsInfo << "m_residual[0] = \n" << m_residual[0] << "\n";
        gsInfo << "m_residual[1] = \n" << m_residual[1] << "\n";
        gsInfo << "m_bUnsteady = " << m_bUnsteady << "\n";
        getchar();*/

        m_pStabEvaluator->initAtElement(m_solUVals, this->m_diffusionCoeff);
        m_pStabEvaluator->setSUPGvars(m_tauStabType, m_degU, m_timeStep, m_elementLength, m_r);
        m_pStabEvaluator->setTauType(m_tauStabType, m_r);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActP = m_activesP.rows();

        m_localMat.setZero(numActP, 1);

        const index_t nQuPoints = quWeights.rows();

        gsMatrix<T> physGradP;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, m_basisGradsP, physGradP);

            T pspgStabParam;
            for (index_t var = 0; var != m_dim; ++var)
            {
                pspgStabParam = m_pStabEvaluator->getTauS(k);
                m_localMat.noalias() += weight * pspgStabParam * ((physGradP.row(var)).transpose() * m_residual[var](k));

            }
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        const index_t usz = m_Umap.freeSize();

        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        for (index_t i = 0; i < numActP; ++i)
        {
            const int ii = m_activesP(i);
            if (m_Pmap.is_free_index(ii))
                rhs(ii + (m_dim)*usz, 0) += m_localMat(i, 0);
        }
    }

    void setPSPG(const int tauType = 1, bool unsteady = false, const T timeStep = 0., const int residualType = 0, const int r = 2)
    {
        m_pspgResidualType = residualType;
        m_tauStabType = tauType;
        m_r = r;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setOldSolutionField(gsField<T>& solution)
    {
        m_oldSolU = solution;
        m_bOldSolSet = true;
    }

    void setPressureSolutionField(gsField<T>& solution)
    {
        m_solP = solution;
        m_bPSolSet = true;
    }

protected:
    int m_degU;
    int m_pspgResidualType;
    int m_tauStabType;
    T m_timeStep;
    int m_r;
    bool m_bUnsteady;
    bool m_bOldSolSet;
    bool m_bPSolSet;

    gsMatrix<T> m_localMat;
    std::vector<gsVector<T> > m_residual;

    // Basis values
    gsMatrix<index_t> m_activesU, m_activesP;
    std::vector<gsMatrix<T> > m_basisDataU;
    gsMatrix<T> m_basisGradsP;
    gsMatrix<T> m_solUVals;

    gsField<T> m_oldSolU;
    gsField<T> m_solP;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    using Base::m_viscosity;
    using Base::m_patchIndex;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_Umap;
    using Base::m_Pmap;

    using uwbVisitorBase<T>::m_elementLength;
};

} // namespace gismo

