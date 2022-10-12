/** @file uwbRANSisoADVisitor.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbRANSBlockVisitors.h"
#include "uwbStabilizationEvaluator.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbRANSisoADVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;

public:
    uwbRANSisoADVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_degU(0), m_isoType(0), m_isoResidualType(0), //m_tauStabType(2),
        m_timeStep(0.), m_alpha(1.), m_r(2), m_bOldSolSet(false), m_bPSolSet(false)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
    }

    ~uwbRANSisoADVisitor()
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
        uwbINSBlockVisitor<T>::initialize(basisRefs, patchIndex, options, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags)
    {
        uwbINSBlockVisitor<T>::initialize(basisRefs, patchIndex, rule, evFlags);
        if (m_pStabEvaluator != NULL)
            m_pStabEvaluator->initialize(rule.numNodes(), m_dim);
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<index_t> activesP;
        // Evaluate basis functions on element nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), activesP);

        const index_t numActU = m_activesU.rows();
        const index_t numActP = activesP.rows();

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);
        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);

        gsMatrix<T> basisGradsP;
        basisRefs.back().deriv_into(quNodes, basisGradsP);

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
            solActPCoeffs.col(j) = m_solP.coefficientVector(m_patchIndex).row(activesP(j)).transpose();

        index_t nQuPoints = quNodes.cols();

        std::vector<gsMatrix<T> > solUGrads, solPGrads;
        solUGrads.resize(nQuPoints);
        solPGrads.resize(nQuPoints);
        gsMatrix<T> physGradU, physGradP;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsU, physGradU);
            solUGrads[k].noalias() = solActUCoeffs * physGradU.transpose();
            geoEval.transformGradients(k, basisGradsP, physGradP);
            solPGrads[k].noalias() = solActPCoeffs * physGradP.transpose();
        }

        gsMatrix<T> physLaplacian, physDeriv2U;
        T diffCoeff;

        std::vector<gsVector<T> > residual(m_dim);
        std::vector<gsVector<T> > residual_adv(m_dim);
        for (int var = 0; var < m_dim; var++)
        {
            residual[var].setZero(nQuPoints);
            residual_adv[var].setZero(nQuPoints);
        }

        switch (m_isoResidualType)
        {
        case 0:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActUCoeffs * physLaplacian.transpose();

                for (int var = 0; var < m_dim; var++)
                {
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                     - (m_viscosity + this->getTurbViscosityVal(k)) * solLaplacian(var, 0) //diffusion
                                     + solPGrads[k](0, var); //pressure term

                    if (m_bUnsteady)
                        residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));
                }
            }
            break;
        case 1:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                for (int var = 0; var < m_dim; var++)
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
            break;
        case 2:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                for (int var = 0; var < m_dim; var++)
                {
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
                    if (m_bUnsteady)
                        residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));
                }
            }
            break;
        case 3:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActUCoeffs * physLaplacian.transpose();
                geoEval.transformDeriv2Hgrad(k, basisGradsU, bHessian_u, physDeriv2U);

                diffCoeff = m_viscosity + this->getTurbViscosityVal(k);

                for (int var = 0; var < m_dim; var++)
                {
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                     - diffCoeff * solLaplacian(var, 0) //diffusion
                                     + solPGrads[k](0, var) //pressure term
                                     - diffCoeff * solActUCoeffs.row(var) * physDeriv2U.col(var); //E-diag term

                    if (m_bUnsteady)
                        residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));
                }
                residual[0](k) -= diffCoeff * solActUCoeffs.row(1) * physDeriv2U.col(m_dim);
                residual[1](k) -= diffCoeff * solActUCoeffs.row(0) * physDeriv2U.col(m_dim);
                if (m_dim == 3)
                {
                    residual[0](k) -= diffCoeff * solActUCoeffs.row(2) * physDeriv2U.col(4);
                    residual[2](k) -= diffCoeff * solActUCoeffs.row(0) * physDeriv2U.col(4);
                    residual[1](k) -= diffCoeff * solActUCoeffs.row(2) * physDeriv2U.col(5);
                    residual[2](k) -= diffCoeff * solActUCoeffs.row(1) * physDeriv2U.col(5);
                }
            }
            break;
        /*case 4:
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
                gsMatrix<T> solLaplacian = solActUCoeffs * physLaplacian.transpose();

                for (int var = 0; var < m_dim; var++)
                {
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k) //advection
                                     - (m_viscosity + this->getTurbViscosityVal(k)) * solLaplacian(var, 0) //diffusion
                                     + solPGrads[k](0, var); //pressure term

                    if (m_bUnsteady)
                        residual[var](k) += 1./m_timeStep * (m_solUVals(var, k) - solUValsOld(var, k));

                    residual_adv[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
                }
            }
            break;*/
        default:
            gsWarn << "Wrong residualType choosen for crosswind stabilization parameter. Choosen default case 1.\n";
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
                for (int var = 0; var < m_dim; var++)
                    residual[var](k) = solUGrads[k].row(var) * m_solUVals.col(k); //advection
        }

        //gsMatrix<T> diffusionCoeff(1, quNodes.cols());
        //for (int i = 0; i < quNodes.cols(); i++)
        //    diffusionCoeff(0, i) = m_viscosity + this->getTurbViscosityVal(i);
        m_pStabEvaluator->initAtElement(m_solUVals, this->m_diffusionCoeff);
        //if (m_isoResidualType == 4)
        //    m_pStabEvaluator->setIsoADvars(m_isoType, m_degU, m_timeStep, solUGrads, residual, residual_adv, m_elementLength, m_alpha);
        //else
            m_pStabEvaluator->setIsoADvars(m_isoType, m_degU, m_timeStep, solUGrads, residual, m_elementLength, m_alpha);
        //m_pStabEvaluator->setTauType(m_tauStabType, m_r);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.resize(m_dim);
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].setZero(numActU, numActU);

        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        const index_t nQuPoints = quWeights.rows();

        gsMatrix<T> physGradU;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, physGradU);

            T isoStabParam;
            for (index_t var = 0; var != m_dim; ++var)
            {
                isoStabParam = m_pStabEvaluator->getRansIsoADStabParam(k, var);
                m_localMat[var].noalias() += weight * isoStabParam * physGradU.transpose() * physGradU;
            }
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        const index_t usz = m_Umap.freeSize();

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
                            for (index_t s = 0; s != m_dim; ++s)
                                sysBlock.coeffRef(ii, jj + s * usz) += m_localMat[s](i, j);
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(ii + s*usz, 0) -= m_localMat[s](i, j) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

    void setIsoAD(const int isoType, const int isoResidualType, bool unsteady = false, const T timeStep = 0.)
    {
        m_isoType = isoType;
        m_isoResidualType = isoResidualType;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setIsoAD(const int isoType, const int isoResidualType, bool unsteady = false,
                      const T timeStep = 0., const T alpha = 1., const int r = 2)
    {
        m_isoType = isoType;
        m_isoResidualType = isoResidualType;
        m_r = r;
        m_alpha = alpha;
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
    int m_isoType;
    int m_isoResidualType;
    //int m_tauStabType;
    T m_timeStep;
    T m_alpha;
    int m_r;
    bool m_bUnsteady;
    bool m_bOldSolSet;
    bool m_bPSolSet;

    std::vector<gsMatrix<T> > m_localMat;

    // Basis values
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    gsMatrix<T> m_solUVals;

    gsField<T> m_oldSolU;
    gsField<T> m_solP;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_dim;
    using uwbINSBlockVisitor<T>::m_Umap;

    using uwbVisitorBase<T>::m_elementLength;
};

} // namespace gismo

