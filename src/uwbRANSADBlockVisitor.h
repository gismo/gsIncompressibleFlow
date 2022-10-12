/** @file uwbRANSADBlockVisitors.h

Author(s): E. Turnerova
*/

#pragma once
#include "uwbRANSBlockVisitors.h"
#include "uwbStabilizationEvaluator.h"
//#include "uwbINSSUPGBlockVisitors.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbRANSADBlockVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;

public:
    uwbRANSADBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity),
        m_degU(0), m_tauStabType(0), m_r(2), m_timeStep(0.)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
    }

    ~uwbRANSADBlockVisitor()
    {
        if (m_pStabEvaluator)
        {
            delete m_pStabEvaluator;
            m_pStabEvaluator = NULL;
        }
    }

    /*virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }*/
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

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);
        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);
        m_basisDataU.push_back(deriv2_u);

        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        m_degU = basisRefs.front().maxDegree();

        /*gsMatrix<T> diffusionCoeff(1, quNodes.cols());
        for (int i = 0; i < quNodes.cols(); i++)
            diffusionCoeff(0, i) = m_viscosity + this->getTurbViscosityVal(i);*/
        m_pStabEvaluator->initAtElement(m_solUVals, this->m_diffusionCoeff);
        m_pStabEvaluator->setADVars(m_tauStabType, m_degU, m_timeStep, m_elementLength, m_r);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2];

        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        gsMatrix<T> physGradU, solUGrad, physDeriv2U;
        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, physGradU);

            //geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, m_physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            geoEval.transformDeriv2Hgrad(k, basisGradsU, bHessian_u, physDeriv2U);

            solUGrad.noalias() = solActUCoeffs * physGradU.transpose();

            T tau_s = m_pStabEvaluator->getTauS(k);

            m_localMat.noalias() += weight * tau_s *
                                    ((m_solUVals.col(k).transpose() * physGradU).transpose()) *
                                    (m_solUVals.col(k).transpose() * physGradU);
            //m_localMat.noalias() += weight * tau_s * (physGradU.transpose() * physGradU);

            /*m_localMat.noalias() -= weight * tau_s *
                                    basisValsU.col(k) *
                                    ((m_solUVals(0, k)*m_solUVals(0, k)) * physDeriv2U.col(0).transpose() +
                                     (2.*m_solUVals(0, k)*m_solUVals(1, k)) * physDeriv2U.col(2).transpose() +
                                     (m_solUVals(1, k)*m_solUVals(1, k)) * physDeriv2U.col(1).transpose());*/

            m_localMat.noalias() -= weight * tau_s *
                                    basisValsU.col(k) *
                                    ((m_solUVals.col(k).transpose() * solUGrad.transpose()) * physGradU);
        }
    }//end assemble

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {

        const index_t usz = m_Umap.freeSize();

        // Local Dofs to global dofs
        m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
        const index_t numActU = m_activesU.rows();

        for (index_t i = 0; i < numActU; ++i)
        {
            const int ii = m_activesU(i);
            if (m_Umap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActU; ++j)
                {
                    const int jj = m_activesU(j);
                    if (m_Umap.is_free_index(jj))
                    {
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(s*usz + ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                    }
                }//end for j
            }//end if free_index
        }//end for i
    }//end local2Global

    void setADType(const int tauStabType, const T timeStep = 0., const int r = 2)
    { m_tauStabType = tauStabType; m_timeStep = timeStep; m_r = r; }

protected:
    int m_degU;
    int m_tauStabType;
    int m_r;
    T m_timeStep;

    gsMatrix<T> m_solUVals;
    gsMatrix<T> m_localMat; // Local matrix

    std::vector<gsMatrix<T> > m_basisDataU;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    //gsField<T> m_turbViscosityField;

    // Members from Base (uwbRANSBlockVisitor)
    using Base::m_parGradsU;
    using Base::m_activesU;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_dim;
    using uwbINSBlockVisitor<T>::m_Umap;

    // Members from uwbVisitorBase
    using uwbVisitorBase<T>::m_elementLength;
};

} // namespace gismo
