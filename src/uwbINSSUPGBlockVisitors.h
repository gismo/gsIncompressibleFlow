/** @file uwbINSSUPGBlockVisitors.h

Author(s): E. Turnerova
*/

#pragma once
#include "uwbINSBlockVisitors.h"
#include "uwbStabilizationEvaluator.h"
#include "uwbTMEvaluators.h"
#include "uwbTMSolverBase.h"

namespace gismo
{

// ============================================================= PARENT ============================================================= //
template <class T>
class uwbINSSUPGBlockVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:
    uwbINSSUPGBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_degU(0), m_tauStabType(0), m_timeStep(0.), m_r(2)
    {
        m_pStabEvaluator = new uwbStabilizationEvaluator<T>();
        m_pTMsolver = NULL;
        m_bRANSsupg = false;

        m_bAverageTau = false;
        m_dirElemLength = false;
    }

    ~uwbINSSUPGBlockVisitor()
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
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        const index_t numActU = m_activesU.rows();

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);
        //bLaplace_u = laplacian(quNodes);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        m_degU = basisRefs.front().maxDegree();

        index_t nQuPoints = quNodes.cols();
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        //-------------------------------------------------
        if (m_bRANSsupg)
        {
            std::vector<gsMatrix<T> > solUGrads(nQuPoints);
            std::vector< gsMatrix<T>> physGradsU;
            physGradsU.resize(nQuPoints);

            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformGradients(k, basisGradsU, physGradsU[k]);
                solUGrads[k].noalias() = solActUCoeffs * physGradsU[k].transpose();
            }

            // Evaluate turbulent viscosity at element nodes
            evalRANSDiffusionCoefficient(solUGrads, quNodes, basisRefs, geoEval);
        }
        //-------------------------------------------------

        T h_advection = 0.;
        T h_diffusion = 0.;
        if (m_dirElemLength)
        {
            gsMatrix<T> solGrads;
            gsMatrix<T> physGrad;
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEval.transformGradients(k, basisGradsU, physGrad);
                solGrads.noalias() = solActUCoeffs * physGrad.transpose();

                gsVector<T> advection = m_solUVals.col(k);

                T gradProduct = advection.transpose() * physGrad.col(k);
                if (advection.norm() == 0.)
                    h_advection += 0.;
                else
                    h_advection += math::abs(gradProduct / advection.norm());

                gsVector<T> adv = (m_solUVals.col(k)).transpose() * solGrads / (m_solUVals.col(k)).norm();
                if (adv.norm() == 0.)
                    h_diffusion += 0.;
                else
                {
                    h_diffusion += (
                                   adv.transpose() * physGrad.col(k) / adv.norm()
                                   ).norm();
                }
            }
            if (h_advection == 0.)
                h_advection = m_elementLength;
            else
                h_advection = 2. / h_advection;
            if (h_diffusion == 0.)
                h_diffusion = m_elementLength;
            else
                h_diffusion = 2. / h_diffusion;
        }

        m_pStabEvaluator->initAtElement(m_solUVals, m_diffusionCoeff);
        m_pStabEvaluator->setSUPGvars(m_tauStabType, m_degU, m_timeStep, m_elementLength, m_r,
                                      m_dirElemLength, h_advection, h_diffusion);
    }

    void evalRANSDiffusionCoefficient(std::vector<gsMatrix<T> >& solUGrads, gsMatrix<T>& quNodes, gsBasisRefs<T> const & basisRefs, gsGeometryEvaluator<T> & geoEval)
    {
        GISMO_ENSURE(m_pTMsolver != NULL, "uwbINSSUPGBlockVisitor: No turbulent model set!");
        std::string sTMEvaluatorType = m_pTMsolver->getTMEvaluator();
        GISMO_ASSERT(sTMEvaluatorType != "", "No evaluator type set.");

        typename uwbTMEvaluator<T>::uPtr evaluatorTM = uwbTMEvaluator<T>::make(sTMEvaluatorType);

        evaluatorTM->initialize(m_viscosity, quNodes.cols());
        evaluatorTM->setKOmegaVariant(sTMEvaluatorType);

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

        evaluatorTM->initAtElement(solUGrads, solKOVals, solKOGrads);
        evaluatorTM->evalQuantities_turbViscosity();
        gsVector<T> turbViscosityVals = evaluatorTM->getTurbViscosity();

        for (int i = 0; i < quNodes.cols(); i++)
            m_diffusionCoeff(0, i) = m_viscosity + turbViscosityVals(i);
    }

    void setTauStabType(const int tauStabType, bool unsteady, const T timeStep = 0., const int r = 2)
    {
        m_tauStabType = tauStabType;
        m_timeStep = timeStep;
        if (!unsteady)
            m_timeStep = 0.;
        m_r = r;
    }

    void setTurbulenceSolver(uwbTMSolverBase<T> * pTMsolver) { m_pTMsolver = pTMsolver; m_bRANSsupg = true; }

protected:
    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

protected:
    int m_degU;
    int m_tauStabType;
    T m_timeStep;
    int m_r;

    uwbTMSolverBase<T> * m_pTMsolver;
    bool m_bRANSsupg;
    bool m_bAverageTau;
    bool m_dirElemLength;

    // Basis values
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    gsMatrix<T> m_solUVals;
    gsMatrix<T> m_physGradU;

    uwbStabilizationEvaluator<T>* m_pStabEvaluator;

    using Base::m_viscosity;
    using Base::m_patchIndex;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_diffusionCoeff;

    using uwbVisitorBase<T>::m_elementLength;
};


// ============================================================= BLOCK A_SUPG ============================================================= //

template <class T>
class uwbINSSUPGBlockAVisitor : public uwbINSSUPGBlockVisitor<T>
{

public:
    typedef uwbINSSUPGBlockVisitor<T> Base;

public:

    uwbINSSUPGBlockAVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    virtual void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2_u;
        basisRefs.front().deriv2_into(quNodes, deriv2_u);

        m_basisDataU.push_back(deriv2_u);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        //const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];
        const gsMatrix<T> & bHessian_u = m_basisDataU[2]; //-> obsahuje i smisenou druhou derivaci

        const index_t nQuPoints = quWeights.rows();

        gsVector<T> tau_s;
        tau_s.setZero(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k)
            tau_s(k) = m_pStabEvaluator->getTauS(k);

        if (Base::m_bAverageTau)
        {
            T sumTau = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
                sumTau += tau_s(k);
            tau_s.setConstant(nQuPoints, sumTau / nQuPoints);
        }
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, m_physGradU);
            geoEval.transformLaplaceHgrad(k, basisGradsU, bHessian_u, m_physLaplacianU); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n

            // Local block
            //m_localMat.noalias() -= weight * m_viscosity * tau_s(k) * (m_physLaplacianU.transpose() * (m_solUVals.col(k).transpose() * m_physGradU));
            m_localMat.noalias() -= weight * m_viscosity * tau_s(k) * ((m_physGradU.transpose() * m_solUVals.col(k)) * m_physLaplacianU);
        }
    }

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
                        //sysBlock.coeffRef(ii, jj) += m_localMat(j, i);
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(s*usz + ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                    }
                }
            }
        }
    }

protected:

    gsMatrix<T> m_localMat; // Local matrix

    // Basis values
    gsMatrix<T> m_physLaplacianU;

    // members from uwbINSSUPGBlockVisitor
    using Base::m_activesU;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_physGradU;
    using Base::m_degU;
    using Base::m_pStabEvaluator;

protected:
    // members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_dim;
    using uwbINSBlockVisitor<T>::m_patchIndex;
};


// ============================================================= BLOCK B_SUPG ============================================================= //

template <class T>
class uwbINSSUPGBlockBVisitor : public uwbINSSUPGBlockVisitor<T>
{

public:
    typedef uwbINSSUPGBlockVisitor<T> Base;

public:

    uwbINSSUPGBlockBVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        // Evaluate basis functions on element nodes
        basisRefs.back().active_into(quNodes.col(0), m_activesP);
        basisRefs.back().deriv_into(quNodes, m_basisGradsP);

        /*const index_t numActU = m_activesU.rows();
        index_t nQuPoints = quNodes.cols();
        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        //std::vector<gsMatrix<T> > solUGrads(nQuPoints);
        m_solUGrads.resize(nQuPoints);
        std::vector< gsMatrix<T>> physGradsU;
        physGradsU.resize(nQuPoints);
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, basisGradsU, physGradsU[k]);
            m_solUGrads[k].noalias() = solActUCoeffs * physGradsU[k].transpose();
        }*/
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();
        const index_t numActP = m_activesP.rows();

        m_localMat.resize(m_dim);
        //m_localRhsLSIC.resize(m_dim);
        for (index_t i = 0; i != m_dim; ++i)
        {
            m_localMat[i].setZero(numActP, numActU);//local_B_SUPG_i
            //m_localRhsLSIC[i].setZero(numActU, 1);
        }
                
        //const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        const index_t nQuPoints = quWeights.rows();

        gsVector<T> tau_s;
        tau_s.setZero(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k)
            tau_s(k) = m_pStabEvaluator->getTauS(k);

        if (Base::m_bAverageTau)
        {
            T sumTau = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
                sumTau += tau_s(k);
            tau_s.setConstant(nQuPoints, sumTau / nQuPoints);
        }
        //T divergenceU;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, m_physGradU);
            geoEval.transformGradients(k, m_basisGradsP, m_physGradP);

            //T tau_s = m_pStabEvaluator->getTauS(k);

            for (index_t i = 0; i != m_dim; ++i) 
            {
                m_localMat[i].noalias() += weight * tau_s(k) *
                                           (m_physGradP.row(i).transpose() * (m_solUVals.col(k).transpose() * m_physGradU));
            }

            /*divergenceU = 0.;
            for (index_t i = 0; i != m_dim; ++i)
                divergenceU += m_solUGrads[k].coeff(i,i);

            //least-squares
            for (index_t i = 0; i != m_dim; ++i)
                m_localRhsLSIC[i].noalias() -= weight * tau_s(k) *
                                               (divergenceU * m_physGradU.row(i)).transpose();*/

        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        const index_t usz = m_Umap.freeSize();

        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        // Local Dofs to global dofs
        m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
        const index_t numActU = m_activesU.rows();

        for (index_t i = 0; i < numActU; ++i)
        {
            const int ii = m_activesU(i);
            if (m_Umap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActP; ++j) // Build B-part of the matrix
                {
                    const int jj = m_activesP(j);
                    if (m_Pmap.is_free_index(jj))
                    {
                        for (index_t s = 0; s != m_dim; ++s)
                            sysBlock.coeffRef(jj, s * usz + ii) += m_localMat[s](j, i);
                    }
                    else //m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        for (index_t s = 0; s < m_dim; ++s)
                            rhs(s * usz + ii, 0) -= m_localMat[s](i, j) * eliminatedDofs[1](bb, 0);
                    }
                }//end for j
            }
            //else {
            // not dependent on 'u' and assuming 'p = 0' (if selected) 
            //}
        }//end for i

        //LSIC
        /*for (index_t i = 0; i < numActU; ++i)
        {
            const int ii = m_activesU(i);
            if (m_Umap.is_free_index(ii))
            {
                for (index_t s = 0; s != m_dim; ++s)
                    rhs(ii + s*usz, 0) -= m_localRhsLSIC[s](i, 0);
            }
        }*/
    }

protected:

    std::vector<gsMatrix<T> > m_localMat;
    //std::vector<gsMatrix<T> > m_localRhsLSIC;
    //std::vector<gsMatrix<T> > m_solUGrads;

    // Basis values
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsP;
    gsMatrix<T> m_physGradP;

    // members from uwbINSSUPGBlockVisitor
    using Base::m_activesU;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_physGradU;
    using Base::m_degU;
    using Base::m_pStabEvaluator;

protected:
    // members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_Pmap;
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_dim;
    using uwbINSBlockVisitor<T>::m_patchIndex;

};//end B_SUPG


// ============================================================= BLOCK N_SUPG ============================================================= //

template <class T>
class uwbINSSUPGBlockNVisitor : public uwbINSSUPGBlockVisitor<T>
{
public:
    typedef uwbINSSUPGBlockVisitor<T> Base;

public:

    uwbINSSUPGBlockNVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        //const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        const index_t nQuPoints = quWeights.rows();

        gsVector<T> tau_s;
        tau_s.setZero(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k)
            tau_s(k) = m_pStabEvaluator->getTauS(k);

        if (Base::m_bAverageTau)
        {
            T sumTau = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
                sumTau += tau_s(k);
            tau_s.setConstant(nQuPoints, sumTau / nQuPoints);
        }

        gsMatrix<T> solActUCoeffs(m_dim, numActU);
        for (int j = 0; j < numActU; j++)
            solActUCoeffs.col(j) = m_solU.coefficientVector(m_patchIndex).row(m_activesU(j)).transpose();

        //gsMatrix<T> solUGrad;
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, m_physGradU);

            //solUGrad.noalias() = solActUCoeffs * m_physGradU.transpose();

            //T tau_s = m_pStabEvaluator->getTauS(k);

            m_localMat.noalias() += weight * tau_s(k) * ((m_solUVals.col(k).transpose() * m_physGradU).transpose()) * (m_solUVals.col(k).transpose() * m_physGradU);

            /*m_localMat.noalias() -= weight * tau_s(k) *
                                    basisValsU.col(k) *
                                    ((m_solUVals.col(k).transpose() * solUGrad.transpose()) * m_physGradU);*/
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

protected:

    gsMatrix<T> m_localMat; // Local matrix

    // members from uwbINSSUPGBlockVisitor
    using Base::m_activesU;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_physGradU;
    using Base::m_degU;
    using Base::m_pStabEvaluator;

protected:
    // members from uwbINSBlockVisitor
    using Base::m_viscosity;
    using Base::m_Umap;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_patchIndex;
};

// ====================================================== M_SUPG ========================================================== //

template <class T>
class uwbINSSUPGBlockMVisitor : public uwbINSSUPGBlockVisitor<T>
{
public:
    typedef uwbINSSUPGBlockVisitor<T> Base;

public:

    uwbINSSUPGBlockMVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        const index_t nQuPoints = quWeights.rows();

        gsVector<T> tau_s;
        tau_s.setZero(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k)
            tau_s(k) = m_pStabEvaluator->getTauS(k);

        if (Base::m_bAverageTau)
        {
            T sumTau = 0.;
            for (index_t k = 0; k < nQuPoints; ++k)
                sumTau += tau_s(k);
            tau_s.setConstant(nQuPoints, sumTau / nQuPoints);
        }

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, m_physGradU);

            //T tau_s = m_pStabEvaluator->getTauS(k);
 
            m_localMat.noalias() += weight * tau_s(k) * (m_physGradU.transpose() * m_solUVals.col(k)) * basisValsU.col(k).transpose();

            /*m_localMat.noalias() -= weight * tau_s(k) *
                                    basisValsU.col(k) *
                                    (m_solUVals.col(k).transpose() * m_physGradU);*/
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
                        //sysBlock.coeffRef(ii, jj) += m_localMat(j, i);
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                        {
                            rhs(s*usz + ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                        }
                                
                    }
                }//end for j
            }//end if free_index
        }//end for i

    }//end local2Global

protected:

    gsMatrix<T> m_localMat; // Local matrix

    // members from uwbINSSUPGBlockVisitor
    using Base::m_activesU;
    using Base::m_basisDataU;
    using Base::m_solUVals;
    using Base::m_physGradU;
    using Base::m_degU;
    using Base::m_pStabEvaluator;
        
protected:
    // members from uwbINSBlockVisitor
    using Base::m_viscosity;
    using Base::m_Umap;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_patchIndex;
};

} // namespace gismo
