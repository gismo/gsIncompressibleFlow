/** @file uwbADRCrosswindVisitors.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbADRBlockVisitors.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbADRCrosswindVisitor : public uwbADRBlockVisitor<T>
{

public:
    typedef uwbADRBlockVisitor<T> Base;

public:
    uwbADRCrosswindVisitor(gsDofMapper& dofMapper, std::string evaluatorType) :
        Base(dofMapper, evaluatorType), m_deg(0), m_crosswindType(0), m_tauType(1), m_timeStep(0.), m_bOldSolSet(false)
    { }

    ~uwbADRCrosswindVisitor() { }

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

        gsMatrix<T> deriv2;
        basisRefs.front().deriv2_into(quNodes, deriv2);
        m_basisData.push_back(deriv2);

        m_deg = basisRefs.front().maxDegree();

        const index_t numAct = m_actives.rows(); // number of active basis functions

        gsMatrix<T> solActCoeffs;
        solActCoeffs.setZero(1, numAct);
        for (int j = 0; j < numAct; j++)
            solActCoeffs.col(j) = m_sol.coefficientVector(m_patchIndex).row(m_actives(j)).transpose();

        gsMatrix<T> solVals = m_sol.value(quNodes, m_patchIndex);
        gsMatrix<T> solValsOld = solVals;
        if (m_bUnsteady)
        {
            GISMO_ASSERT(m_bOldSolSet, "Old time step solution not set in Crosswind visitor.");
            solValsOld = m_oldSol.value(quNodes, m_patchIndex);
        }

        const gsMatrix<T> & bGrads = m_basisData[1];
        gsMatrix<T> physGrad;
        gsMatrix<T> physLaplacian;

        index_t nQuPoints = quNodes.cols();
        std::vector<gsMatrix<T> > solGrads(nQuPoints);
        gsVector<T> residual(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, bGrads, physGrad);
            solGrads[k].noalias() = solActCoeffs * physGrad.transpose();

            gsVector<T> advection = this->getADREvaluator()->getAdvectionCoefficient(k);
            T diffusion = this->getADREvaluator()->getDiffusionCoefficient(k);
            T reaction = this->getADREvaluator()->getReactionCoefficient(k);

            geoEval.transformLaplaceHgrad(k, bGrads, deriv2, physLaplacian); //[B1_xx + B1_yy + B1_zz, ..., Bi_xx + Bi_yy + Bi_zz, ...] - 1 x n
            gsMatrix<T> solLaplacian = solActCoeffs * physLaplacian.transpose(); //scalar for 1 unknown ADR variable

            residual(k) = solGrads[k].row(0) * advection
                        - diffusion * solLaplacian(0, 0)
                        + reaction * solVals(0, k);//.col(k);

            if (m_bUnsteady)
                residual(k) += 1./m_timeStep * (solVals(0, k) - solValsOld(0, k));
                        //1./m_timeStep * (solVals.col(k) - solValsOld.col(k));
        }

        getADREvaluator()->setCrosswindVars(m_crosswindType, m_deg, m_timeStep, solGrads, residual, m_elementLength);
        getADREvaluator()->setTauType(m_tauType);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisGrads = m_basisData[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);
            T cwStabParam = this->getADREvaluator()->getCrosswindStabParam(k);
            gsMatrix<T> proj = this->getADREvaluator()->getCrosswindProjection(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGrads, m_physGrad);

            m_localMat.noalias() += weight * cwStabParam * ((proj * m_physGrad).transpose() * m_physGrad);
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

    void setCrosswind(const int crosswindType, bool unsteady = false, const T timeStep = 0.)
    {
        m_crosswindType = crosswindType;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setCrosswind(const int crosswindType, const int tauType = 1, bool unsteady = false, const T timeStep = 0.)
    {
        m_crosswindType = crosswindType;
        m_tauType = tauType;
        m_bUnsteady = unsteady;
        m_timeStep = timeStep;
        if (m_bUnsteady)
            GISMO_ASSERT(timeStep > 0, "Incomplete or wrong setting of the unsteady problem.");
    }

    void setOldSolutionField(gsField<T>& solution)
    {
        m_oldSol = solution;
        m_bOldSolSet = true;
    }

protected:
    int m_deg;
    int m_crosswindType;
    int m_tauType;
    T m_timeStep;
    bool m_bUnsteady;
    bool m_bOldSolSet;

    gsMatrix<T> m_localMat; // Local matrix
    gsMatrix<T> m_physGrad;
    gsMatrix<T> m_physLaplacian;

    using Base::m_sol;
    using Base::m_oldSol;
    using Base::m_diffusionCoeff;
    using Base::m_advectionCoeff;
    using Base::m_patchIndex;
    using Base::m_actives;
    using Base::m_basisData;
    using Base::m_elementLength;
    using Base::m_dim;
    using Base::getADREvaluator;
    using Base::m_mapper;
};

} // namespace gismo

