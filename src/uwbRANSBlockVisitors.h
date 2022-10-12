/** @file uwbRANSBlockVisitors.h

Author(s): E. Turnerova
*/

#pragma once
#include "uwbINSBlockVisitors.h"
#include "uwbTMEvaluators.h"
#include "uwbTMSolverBase.h"
#include "uwbStabilizationEvaluator.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbRANSBlockVisitor : public uwbINSBlockVisitor<T>
{
public:
    typedef memory::unique_ptr<uwbRANSBlockVisitor> uPtr;

    typedef uwbINSBlockVisitor<T> Base;

public:
    uwbRANSBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_bEffectiveViscSet(false)
    {
        m_pTMsolver = NULL;
        //m_bIncreaseNuT = false;
        //m_yPlusLimit = 5;
        m_sTMEvaluatorType = "";
    }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        basisRefs.front().deriv_into(quNodes, m_parGradsU);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate velocity gradients
        index_t numActU = m_activesU.rows();
        index_t nQuPoints = quNodes.cols();
        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        std::vector<gsMatrix<T> > solUGrads(nQuPoints);
        m_physGradsU.resize(nQuPoints);

        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, m_parGradsU, m_physGradsU[k]);
            solUGrads[k].noalias() = solActUCoeffs * m_physGradsU[k].transpose();
        }

        // Evaluate turbulent viscosity at element nodes
        evalTurbulentViscosity(solUGrads, quNodes, basisRefs, geoEval);


    }

    void evalTurbulentViscosity(std::vector<gsMatrix<T> >& solUGrads, gsMatrix<T>& quNodes, gsBasisRefs<T> const & basisRefs,
                                gsGeometryEvaluator<T> & geoEval)
    {
        GISMO_ENSURE(m_pTMsolver != NULL, "uwbRANSBlockVisitor: No turbulent model set!");
        m_sTMEvaluatorType = m_pTMsolver->getTMEvaluator();
        GISMO_ASSERT(m_sTMEvaluatorType != "", "No evaluator type set.");

        typename uwbTMEvaluator<T>::uPtr evaluatorTM = uwbTMEvaluator<T>::make(m_sTMEvaluatorType);

        evaluatorTM->initialize(m_viscosity, quNodes.cols());
        evaluatorTM->setKOmegaVariant(m_sTMEvaluatorType);

        gsField<T> solTMfield = m_pTMsolver->constructSolution();
        gsMatrix<T> solKOVals = solTMfield.value(quNodes, m_patchIndex);

        //--- KOGrads
        gsMatrix<T> solActKOCoeffs, physGradKO, bGradsKO;
        gsMatrix<index_t> activesKO;

        basisRefs.back().active_into(quNodes.col(0), activesKO);
        basisRefs.back().deriv_into(quNodes, bGradsKO);

        const index_t numActKO = activesKO.rows();

        solActKOCoeffs.setZero(2, numActKO);
        for (int j = 0; j < numActKO; j++)
            solActKOCoeffs.col(j) = solTMfield.coefficientVector(m_patchIndex).row(activesKO(j)).transpose();

        index_t nQuPoints = quNodes.cols();
        std::vector<gsMatrix<T> > solKOGrads(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, bGradsKO, physGradKO);
            solKOGrads[k].noalias() = solActKOCoeffs * physGradKO.transpose();
        }

        if (this->checkWallDistanceBasedTM())
        {
            gsField<T> solPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonSolution();
            gsMultiPatch<T> patchesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonPatches();
            gsMultiBasis<T> basesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonBasis();

            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            gsMatrix<T> solPoissonVals = solPoisson.value(quNodes, m_patchIndex);
            gsMatrix<T> solActPoissonCoeffs;
            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            evaluatorTM->evalWallDistance(solPoissonVals, solPoissonGrads);
        }

        evaluatorTM->initAtElement(solUGrads, solKOVals, solKOGrads);
        evaluatorTM->evalQuantities_turbViscosity();
        m_turbViscosityVals = evaluatorTM->getTurbViscosity();

        /*for (int ii = 0; ii < (evaluatorTM->getWallDistance()).rows(); ii++)
        {
            if (m_turbViscosityVals(ii) > 1e-4)
            {
                gsInfo << "m_turbViscosityVals(ii) = \n" << m_turbViscosityVals(ii) << "\n";
                gsInfo << "evaluatorTM->getWallDistance(ii) = \n" << evaluatorTM->getWallDistance(ii) << "\n";
                getchar();
            }
        }*/

        /*if (m_bIncreaseNuT)
        {
            gsMatrix<T> dirGrads;
            gsVector<T> unormal;
            geoEval.outerNormal(k, side, unormal);
            unormal = (-1) * unormal/unormal.norm();

            dirGrads = uGrads * unormal;

            T shearStress = m_blockAssembler.getViscosity() * (std::signbit(solUVals(k)) == true ? -1 : 1) * dirGrads.norm();

        }*/

        T uFreeStream = 44.2;
        T inletWidth = 8 * 0.0127;
        T Re, skinFrCoeff, tau_w, u_tau, yPlus;
        for (int i = 0; i < quNodes.cols(); i++)
        {
            m_diffusionCoeff(0, i) = m_viscosity + m_turbViscosityVals(i);

            /*if (m_bIncreaseNuT)
            {
                Re = uFreeStream * inletWidth / m_diffusionCoeff(0, i);
                skinFrCoeff = math::pow(2.0*math::log10(Re) - 0.65, -2.3);
                tau_w = 0.5 * skinFrCoeff * math::pow(uFreeStream, 2);
                u_tau = math::sqrt(tau_w);
                yPlus = evaluatorTM->getWallDistance(i) * u_tau / m_diffusionCoeff(0, i);

                if (math::isnan(yPlus))
                {
                    yPlus = 10000000.;
                }

                //if (yPlus < m_yPlusLimit)
                //{
                    //m_diffusionCoeff(0, i) += 0.01 * math::exp(-(1/15)*yPlus*yPlus);
                    //m_turbViscosityVals(i) += 0.01 * math::exp(-(1/15)*yPlus*yPlus);
                    //m_diffusionCoeff(0, i) += 0.001 * math::exp(-(1/5)*yPlus*yPlus);
                    //m_turbViscosityVals(i) += 0.001 * math::exp(-(1/5)*yPlus*yPlus);
                    //m_diffusionCoeff(0, i) += 0.0001 * math::exp(-(1/5)*yPlus*yPlus);
                    //m_turbViscosityVals(i) += 0.0001 * math::exp(-(1/5)*yPlus*yPlus);
                    //m_diffusionCoeff(0, i) += 0.0001 * math::exp(-2*yPlus*yPlus);
                    //m_turbViscosityVals(i) += 0.0001 * math::exp(-2*yPlus*yPlus);
                    m_diffusionCoeff(0, i) += 0.001 * math::exp(-2*yPlus*yPlus);
                    m_turbViscosityVals(i) += 0.001 * math::exp(-2*yPlus*yPlus);
                //}


            }*/



        }

        //for (int i = 0; i < quNodes.cols(); i++)
        //    m_diffusionCoeff(0, i) = m_viscosity + m_turbViscosityVals(i);

        m_bEffectiveViscSet = true;
    }

    void setTurbulenceSolver(uwbTMSolverBase<T> * pTMsolver) { m_pTMsolver = pTMsolver; }

    gsMatrix<T> getTurbViscosityVals()
    { 
        GISMO_ASSERT(m_turbViscosityVals.rows() > 0, "Turbulent viscosity not evaluated yet.");
        gsMatrix<T> turbViscosity = m_turbViscosityVals;
        return turbViscosity;
    }

    T getTurbViscosityVal(int k)
    {
        GISMO_ASSERT(m_turbViscosityVals.rows() > 0, "Turbulent viscosity not evaluated yet.");
        return m_turbViscosityVals(k);
    }

    bool checkWallDistanceBasedTM()
    {
        return (m_sTMEvaluatorType == "koSST" || m_sTMEvaluatorType == "koSSTMenter2009" || m_sTMEvaluatorType == "koSAS"
                || m_sTMEvaluatorType == "koSAS_SS" || m_sTMEvaluatorType == "koSAS_SO" || m_sTMEvaluatorType == "koSAS_OO");
    }

protected:
    gsVector<T> m_turbViscosityVals;
    uwbTMSolverBase<T> * m_pTMsolver;

    // Basis values
    gsMatrix<index_t>  m_activesU;
    gsMatrix<T>         m_parGradsU;
    std::vector<gsMatrix<T> > m_physGradsU;

    bool m_bEffectiveViscSet;

    //bool m_bIncreaseNuT;
    //T m_yPlusLimit;

    std::string m_sTMEvaluatorType;

    // Members from uwbINSBlockVisitor
    using Base::m_solU;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;
    using Base::m_viscosity;
    using Base::m_diffusionCoeff;
};


// ============================================================= BLOCK A_RANS ============================================================= //

template <class T>
class uwbRANSBlockAsymVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;

public:

    uwbRANSBlockAsymVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        GISMO_ENSURE(this->m_bEffectiveViscSet, "uwbRANSBlockVisitor: effective viscosity not set!");

        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);//local_A

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // Local block A
            m_localMat.template triangularView<Eigen::Upper>() += weight * m_turbViscosityVals(k) * (m_physGradsU[k].transpose() * m_physGradsU[k]);
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
                for (index_t j = 0; j < numActU; ++j) // Build A-part of the matrix
                {
                    const int jj = m_activesU(j);
                    if (m_Umap.is_free_index(jj))
                    {
                        if (j >= i) {
                            sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                        }
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(ii + s*usz, 0) -= m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

protected:

    gsMatrix<T> m_localMat; // Local matrix

    // Members from uwbRANSBlockVisitor
    using Base::m_activesU;
    using Base::m_physGradsU;
    using Base::m_turbViscosityVals;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_dim; // Velocity vector dimension
    using uwbINSBlockVisitor<T>::m_patchIndex;
};


// ============================================================= BLOCK EsymDiag_RANS ============================================================= //

template <class T>
class uwbRANSBlockEdiagsymVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;


public:

    uwbRANSBlockEdiagsymVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        GISMO_ENSURE(this->m_bEffectiveViscSet, "uwbRANSBlockVisitor: effective viscosity not set!");

        const index_t numActU = m_activesU.rows();

        m_localMat.resize(m_dim);
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].setZero(numActU, numActU);

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // Local block Ediag
            for (index_t s = 0; s != m_dim; ++s)
                m_localMat[s].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (m_physGradsU[k].row(s).transpose() * m_physGradsU[k].row(s));
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
                            {
                                sysBlock.coeffRef(ii, jj + s * usz) += m_localMat[s](i, j);
                            }
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

protected:

    std::vector<gsMatrix<T> > m_localMat;

    // Members from uwbRANSBlockVisitor
    using Base::m_activesU;
    using Base::m_physGradsU;
    using Base::m_turbViscosityVals;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_dim; // Velocity vector dimension
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
};


// ============================================================= BLOCK E12 E13 E23 RANS ============================================================= //

template <class T>
class uwbRANSRhsEnondiagVisitor : public uwbRANSBlockVisitor<T>
{

public:
    typedef uwbRANSBlockVisitor<T> Base;


public:

    uwbRANSRhsEnondiagVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        GISMO_ENSURE(this->m_bEffectiveViscSet, "uwbRANSBlockVisitor: effective viscosity not set!");

        const index_t numActU = m_activesU.rows();

        switch (m_dim)
        {
        case 2:
            m_numBlocks = m_dim - 1;
            break;

        case 3:
            m_numBlocks = m_dim;
            break;

        default:
            GISMO_ERROR("RANS are implemented only for 2D and 3D.");
            break;
        }

        m_localMat.resize(m_numBlocks);
        for (index_t i = 0; i != m_numBlocks; ++i)
            m_localMat[i].setZero(numActU, numActU);

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // Local blocks
            m_localMat[0].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (m_physGradsU[k].row(1).transpose() * m_physGradsU[k].row(0)); //dv1/dy * du2/dx

            if (m_numBlocks == m_dim)
            {
                m_localMat[1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (m_physGradsU[k].row(2).transpose() * m_physGradsU[k].row(0)); //dv1/dz * du3/dx
                m_localMat[2].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (m_physGradsU[k].row(2).transpose() * m_physGradsU[k].row(1)); //dv2/dz * du3/dy
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
                        sysBlock.coeffRef(ii, jj) += m_localMat[0](i, j);
                        if (m_numBlocks == m_dim)
                        {
                            sysBlock.coeffRef(ii, usz + jj) += m_localMat[1](i, j);
                            sysBlock.coeffRef(ii, 2 * usz + jj) += m_localMat[2](i, j);
                        }
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat[0](i, j) * eliminatedDofs[0](bb, 1);
                        rhs(usz + ii, 0) -= m_localMat[0](j, i) * eliminatedDofs[0](bb, 0);
                        if (m_numBlocks == m_dim)
                        {
                            rhs(ii, 0) -= m_localMat[1](i, j) * eliminatedDofs[0](bb, 2);
                            rhs(2 * usz + ii, 0) -= m_localMat[1](j, i) * eliminatedDofs[0](bb, 0);
                            rhs(usz + ii, 0) -= m_localMat[2](i, j) * eliminatedDofs[0](bb, 2);
                            rhs(2 * usz + ii, 0) -= m_localMat[2](j, i) * eliminatedDofs[0](bb, 1);
                        }
                    }
                }
            }
        }
    }

protected:

    std::vector<gsMatrix<T> > m_localMat;
    index_t m_numBlocks;

    // Members from uwbRANSBlockVisitor
    using Base::m_activesU;
    using Base::m_physGradsU;
    using Base::m_turbViscosityVals;

    // Members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_dim; // Velocity vector dimension
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
};

} // namespace gismo
