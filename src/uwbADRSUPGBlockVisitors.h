/** @file uwbADRSUPGBlockVisitors.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbADRBlockVisitors.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbADRSUPGBlockVisitor : public uwbADRBlockVisitor<T>
{

public:
    typedef uwbADRBlockVisitor<T> Base;

public:
    uwbADRSUPGBlockVisitor(gsDofMapper& dofMapper, std::string evaluatorType) :
        Base(dofMapper, evaluatorType), m_deg(0), m_tauStabType(1), m_timeStep(0.)
    { }

    ~uwbADRSUPGBlockVisitor() { }

    virtual inline void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV2 | NEED_2ND_DER;
    }

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        Base::evaluate(basisRefs, geoEval, quNodes);

        gsMatrix<T> deriv2KO;
        basisRefs.front().deriv2_into(quNodes, deriv2KO);
        m_basisData.push_back(deriv2KO);

        m_deg = basisRefs.front().maxDegree();

        getADREvaluator()->setSUPGvars(m_tauStabType, m_deg, m_timeStep, m_elementLength);
    }

    void setSUPG(const int tauStabType = 1, const T timeStep = 0.) { m_tauStabType = tauStabType; m_timeStep = timeStep; }

protected:
    int m_deg;
    int m_tauStabType;
    T m_timeStep;

    using Base::m_diffusionCoeff;
    using Base::m_advectionCoeff;
    using Base::m_patchIndex;
    using Base::m_actives;
    using Base::m_basisData;
    using Base::m_elementLength;
    using Base::m_dim;
    using Base::getADREvaluator;
};

// ============================================================= BLOCK M ============================================================= //
template <class T>
class uwbMassMatrixSUPGVisitor : public uwbADRSUPGBlockVisitor<T>
{
public:
    typedef uwbADRSUPGBlockVisitor<T> Base;

public:

    uwbMassMatrixSUPGVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];
        const gsMatrix<T> & basisGrads = m_basisData[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);
            T tau_s = this->getADREvaluator()->getTauS(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGrads, m_physGrad);

            gsVector<T> advection = this->getADREvaluator()->getAdvectionCoefficient(k);
            m_localMat.noalias() += weight * tau_s * ((m_physGrad.transpose() * advection) * basisVals.col(k).transpose());
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        // Local Dofs to global dofs
        m_mapper.localToGlobal(m_actives, m_patchIndex, m_actives);
        const index_t numAct = m_actives.rows();

        for (index_t i = 0; i < numAct; ++i)
        {
            const int ii = m_actives(i);
            if (m_mapper.is_free_index(ii))
            {
                for (index_t j = 0; j < numAct; ++j)
                {
                    const int jj = m_actives(j);
                    if (m_mapper.is_free_index(jj))
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    else // m_mapper.is_boundary_index(jj)
                    {
                        const int bb = m_mapper.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, 0);
                    }
                }
            }
        }
    }

protected:
    gsMatrix<T> m_localMat; // Local matrix
    gsMatrix<T> m_physGrad;

    //members from uwbADRBlockVisitor<T>
    using uwbADRBlockVisitor<T>::m_mapper;
    using uwbADRBlockVisitor<T>::m_patchIndex;
    using uwbADRBlockVisitor<T>::m_actives;
    using uwbADRBlockVisitor<T>::m_basisData;
    using uwbADRBlockVisitor<T>::m_advectionCoeff;
};


// ============================ MATRICES N_ADR advection --- A_ADR diffusion --- R_ADR reaction ============================================================= //
template <class T>
class uwbADRSUPGVisitor : public uwbADRSUPGBlockVisitor<T>
{

public:
    typedef uwbADRSUPGBlockVisitor<T> Base;

public:

    uwbADRSUPGVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];
        const gsMatrix<T> & basisGrads = m_basisData[1];
        const gsMatrix<T> & basisHessian = m_basisData[2];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);
            T tau_s = this->getADREvaluator()->getTauS(k);
            gsVector<T> advection = this->getADREvaluator()->getAdvectionCoefficient(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGrads, m_physGrad);
            // compute Laplatians
            geoEval.transformLaplaceHgrad(k, basisGrads, basisHessian, m_physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n

            m_localMat.noalias() += weight * tau_s * ((m_physGrad.transpose() * advection) * (advection.transpose() * m_physGrad))   //advection
                                 -  weight * tau_s * this->getADREvaluator()->getDiffusionCoefficient(k) * ((m_physGrad.transpose() * advection) * m_physLaplacian)
                                 +  weight * tau_s * this->getADREvaluator()->getReactionCoefficient(k) * ((m_physGrad.transpose() * advection) * basisVals.col(k).transpose());//reaction

            //part with the derivative of the diffusion coefficient
            //now case with the constant coefficients
            //m_localMat.noalias() -= weight * tau_s * ((m_physGrad.transpose() * m_advectionCoeff) * (solDiffusionCoeffGrad * m_physGrad));
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        // Local Dofs to global dofs
        m_mapper.localToGlobal(m_actives, m_patchIndex, m_actives);
        const index_t numAct = m_actives.rows();

        for (index_t i = 0; i < numAct; ++i)
        {
            const int ii = m_actives(i);
            if (m_mapper.is_free_index(ii))
            {
                for (index_t j = 0; j < numAct; ++j)
                {
                    const int jj = m_actives(j);
                    if (m_mapper.is_free_index(jj))
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    else // m_mapper.is_boundary_index(jj)
                    {
                        const int bb = m_mapper.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, 0);
                    }
                }
            }
        }
    }

protected:
    gsMatrix<T> m_physGrad;
    gsMatrix<T> m_physLaplacian;
    gsMatrix<T> m_localMat;

    //members from uwbADRBlockVisitor<T>
    using uwbADRBlockVisitor<T>::m_basisData;
    using uwbADRBlockVisitor<T>::m_actives;
    using uwbADRBlockVisitor<T>::m_mapper;
    using uwbADRBlockVisitor<T>::m_patchIndex;
    using uwbADRBlockVisitor<T>::m_advectionCoeff;
};

} // namespace gismo

