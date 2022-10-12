/** @file uwbADRBlockVisitors.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbADREvaluators.h"

namespace gismo
{
// ============================================================= PARENT ============================================================= //
template <class T>
class uwbADRBlockVisitor
{
protected:

    const gsDofMapper & m_mapper;
    gsField<T> m_sol, m_oldSol;
    index_t m_dim;
    int m_patchIndex;
    T m_diffusionCoeff;
    gsVector<T> m_advectionCoeff;
    T m_reactionCoeff;
    bool m_bLinCoeffsSet;
    bool m_bNonlinADcoeffsSet;
    bool m_bNonlinADRcoeffsSet;
    bool m_bSolutionSet;
    T m_elementLength;

    gsMatrix<index_t> m_actives;
    std::vector<gsMatrix<T> > m_basisData;

    std::string m_evaluatorType;
    uwbADREvaluator<T>* m_pEvaluator;

    gsField<T> m_advectionField;
    gsField<T> m_diffusionField;
    gsField<T> m_reactionField;
    gsMultiPatch<T> m_patchesCoeffs;
    gsMultiBasis<T> m_basesCoeffs;

public:
    uwbADRBlockVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : m_mapper(dofMapper), m_evaluatorType(evaluatorType)
    {
        m_bLinCoeffsSet = false;
        m_bNonlinADcoeffsSet = false;
        m_bNonlinADRcoeffsSet = false;
        m_bSolutionSet = false;
        m_pEvaluator = NULL;
        m_elementLength = 0.;
    }

    ~uwbADRBlockVisitor()
    {
        if (m_pEvaluator)
        {
            delete m_pEvaluator;
            m_pEvaluator = NULL;
        }
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions & options,
        gsQuadRule<T>    & rule,
        unsigned         & evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basisRefs.dim();
        m_patchIndex = patchIndex;

        rule = gsGaussRule<T>(basis, options.quA, options.quB);

        if (m_evaluatorType == "linConstCoeffs")
            GISMO_ASSERT(m_bLinCoeffsSet, "No ADR coefficients set in the visitor.");
        else if (m_evaluatorType == "nonlinCoeffsField")
            GISMO_ASSERT(m_bNonlinADcoeffsSet || m_bNonlinADRcoeffsSet, "No ADR coefficients set in the visitor.");
        else
            GISMO_ERROR("Wrong or no ADR evaluator set in the visitor. 'linNonConstCoeffs' and 'nonlinCoeffsBurgers' not finished yet.");

        m_pEvaluator = new uwbADREvaluator<T>(rule.numNodes(), m_dim);

        initializeSpecific(evFlags);

        if (getADREvaluator() != NULL)
            getADREvaluator()->initialize();
    }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>    & rule,
        unsigned         & evFlags)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basisRefs.dim();
        m_patchIndex = patchIndex;

        gsVector<index_t> numQuadNodes(m_dim);
        for (int i = 0; i < m_dim; ++i)
            numQuadNodes[i] = basis.maxDegree() + 1; 

        rule = gsGaussRule<T>(numQuadNodes);

        if (m_evaluatorType == "linConstCoeffs")
            GISMO_ASSERT(m_bLinCoeffsSet, "No ADR coefficients set in the visitor.");
        else if (m_evaluatorType == "nonlinCoeffsField")
            GISMO_ASSERT(m_bNonlinADcoeffsSet || m_bNonlinADRcoeffsSet, "No ADR coefficients set in the visitor.");
        else
            GISMO_ERROR("Wrong or no ADR evaluator set in the visitor. 'linNonConstCoeffs' and 'nonlinCoeffsBurgers' not finished yet.");

        m_pEvaluator = new uwbADREvaluator<T>(rule.numNodes(), m_dim);

        initializeSpecific(evFlags);

        if (getADREvaluator() != NULL)
            getADREvaluator()->initialize();
    }

    /*virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    { GISMO_NO_IMPLEMENTATION }*/

    virtual inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.front().active_into(quNodes.col(0), m_actives);
        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisData);
        geoEval.evaluateAt(quNodes);

        if (m_evaluatorType == "linConstCoeffs")
            getADREvaluator()->initAtElement(m_diffusionCoeff, m_advectionCoeff, m_reactionCoeff);
        else if (m_evaluatorType == "nonlinCoeffsField" && m_bNonlinADcoeffsSet)
        {
            gsMatrix<T> advectionVals = m_advectionField.value(quNodes, m_patchIndex);
            gsMatrix<T> diffusionVals = m_diffusionField.value(quNodes, m_patchIndex);
            getADREvaluator()->initAtElement(advectionVals, diffusionVals);
        }
        else if (m_evaluatorType == "nonlinCoeffsField" && m_bNonlinADRcoeffsSet)
        {
            gsMatrix<T> advectionVals = m_advectionField.value(quNodes, m_patchIndex);
            gsMatrix<T> diffusionVals = m_diffusionField.value(quNodes, m_patchIndex);
            gsMatrix<T> reactionVals = m_reactionField.value(quNodes, m_patchIndex);
            getADREvaluator()->initAtElement(advectionVals, diffusionVals, reactionVals);
        }
        else
            GISMO_ERROR("Wrong or no ADR evaluator set in the visitor.");

    }

    virtual inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    virtual inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    { GISMO_NO_IMPLEMENTATION }

    void setADRCoefficients(T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff)
    {
        m_diffusionCoeff = diffusionCoeff;
        m_advectionCoeff = advectionCoeff;
        m_reactionCoeff = reactionCoeff;
        m_bLinCoeffsSet = true;
    }

    void setAdvectionDiffusionCoeffFields(gsField<T>& advectionField, gsField<T>& diffusionField)
    {
        m_advectionField = advectionField;
        m_diffusionField = diffusionField;
        m_bNonlinADcoeffsSet = true;
    }
    void setADRcoeffFields(gsField<T>& advectionField, gsField<T>& diffusionField, gsField<T>& reactionField)
    {
        m_advectionField = advectionField;
        m_diffusionField = diffusionField;
        m_reactionField = reactionField;
        m_bNonlinADRcoeffsSet = true;
    }
    void setCoeffGeometry(const gsMultiPatch<T>& patchesCoeffs, const gsMultiBasis<T>& basesCoeffs)
    { m_patchesCoeffs = patchesCoeffs; m_basesCoeffs = basesCoeffs; }

    void setCurrentSolution(gsField<T>& solution)
    {
        m_sol = solution;
        this->m_bSolutionSet = true;
    }

    void setElementLength(T elementLength) { m_elementLength = elementLength; }

    virtual inline void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }

    uwbADREvaluator<T>* getADREvaluator()
    {
        GISMO_ASSERT(m_evaluatorType != "", "No evaluator type set in the visitor.");
        return m_pEvaluator;
    }

};

// ============================================================= mass matrix ============================================================= //
template <class T>
class uwbMassMatrixVisitor : public uwbADRBlockVisitor<T>
{
public:
    typedef uwbADRBlockVisitor<T> Base;

public:

    uwbMassMatrixVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    inline void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            m_localMat.noalias() += weight * (basisVals.col(k) * basisVals.col(k).transpose());
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

    //members from uwbADRBlockVisitor<T>
    using Base::m_mapper;
    using Base::m_patchIndex;
    using Base::m_actives;
    using Base::m_basisData;
};

// ================================ Matrices N_ADR advection -- A_ADR diffusion -- R_ADR reaction ============================================================= //
template <class T>
class uwbADRVisitor : public uwbADRBlockVisitor<T>
{

public:
    typedef uwbADRBlockVisitor<T> Base;

public:

    uwbADRVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];
        const gsMatrix<T> & bGrads = m_basisData[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, m_physGrad);

            m_localMat.noalias() += weight * (basisVals.col(k) * ((getADREvaluator()->getAdvectionCoefficient(k)).transpose() * m_physGrad))     //advection
                                 +  weight * getADREvaluator()->getDiffusionCoefficient(k) * (m_physGrad.transpose() * m_physGrad)            //diffusion
                                 +  weight * getADREvaluator()->getReactionCoefficient(k) * (basisVals.col(k) * basisVals.col(k).transpose()); //reaction
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
    gsMatrix<T> m_localMat; // Local matrix

    using Base::m_basisData;
    using Base::m_actives;
    using Base::m_mapper;
    using Base::m_patchIndex;
    using Base::m_advectionCoeff;

    using Base::getADREvaluator;
};

// ============================================================= mass matrix ============================================================= //
template <class T>
class uwbMassMatrixFCTLowOrderVisitor : public uwbADRBlockVisitor<T>
{
public:
    typedef uwbADRBlockVisitor<T> Base;

public:

    uwbMassMatrixFCTLowOrderVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    inline void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_localMat.setZero(numAct, numAct);
        gsMatrix<T> tempMat;
        tempMat.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            tempMat.noalias() += weight * (basisVals.col(k) * basisVals.col(k).transpose());
        }
        for (index_t i = 0; i < numAct; i++)
            for (index_t j = 0; j < numAct; j++)
                m_localMat.coeffRef(i, i) += tempMat(i, j);
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

    //members from uwbADRBlockVisitor<T>
    using Base::m_mapper;
    using Base::m_patchIndex;
    using Base::m_actives;
    using Base::m_basisData;
};


// ========== FCT stabilization -- low order ======= Matrices N_ADR advection -- A_ADR diffusion -- R_ADR reaction ============================================================= //
template <class T>
class uwbADRFCTLowOrderVisitor : public uwbADRBlockVisitor<T>
{

public:
    typedef uwbADRBlockVisitor<T> Base;

public:

    uwbADRFCTLowOrderVisitor(gsDofMapper& dofMapper, std::string evaluatorType) : Base(dofMapper, evaluatorType) { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numAct = m_actives.rows();

        m_mLocalAdvection.setZero(numAct, numAct);
        m_mLocalDiffusion.setZero(numAct, numAct);
        m_mLocalReaction.setZero(numAct, numAct);

        gsMatrix<T> mNumDiffusion;
        mNumDiffusion.setZero(numAct, numAct);

        const gsMatrix<T> & basisVals = m_basisData[0];
        const gsMatrix<T> & bGrads = m_basisData[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, m_physGrad);

            m_mLocalAdvection.noalias() -= weight * (basisVals.col(k) * ((getADREvaluator()->getAdvectionCoefficient(k)).transpose() * m_physGrad));     //advection
            m_mLocalDiffusion.noalias() += weight * getADREvaluator()->getDiffusionCoefficient(k) * (m_physGrad.transpose() * m_physGrad);             //diffusion
            m_mLocalReaction.noalias() += weight * getADREvaluator()->getReactionCoefficient(k) * (basisVals.col(k) * basisVals.col(k).transpose()); //reaction
        }

        for (index_t i = 0; i < numAct; i++)
        {
            T numDiffDiag = 0.;
            for (index_t j = 0; j < numAct; j++)
            {
                if (j != i)
                {
                    mNumDiffusion.coeffRef(i, j) = math::max(-m_mLocalAdvection(i, j), -m_mLocalAdvection(j, i));
                    mNumDiffusion.coeffRef(i, j) = math::max(mNumDiffusion(i, j), 0.);
                    numDiffDiag -= mNumDiffusion(i, j);
                }
            }
            mNumDiffusion.coeffRef(i, i) = numDiffDiag;
        }
        m_mLocalAdvection = -m_mLocalAdvection - mNumDiffusion;

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
                        sysBlock.coeffRef(ii, jj) += m_mLocalAdvection(i, j) + m_mLocalDiffusion(i, j) + m_mLocalReaction(i, j);
                    else // m_mapper.is_boundary_index(jj)
                    {
                        const int bb = m_mapper.global_to_bindex(jj);
                        rhs(ii, 0) -= m_mLocalAdvection(i, j) * eliminatedDofs[0](bb, 0)
                                   -  m_mLocalDiffusion(i, j) * eliminatedDofs[0](bb, 0)
                                   -  m_mLocalReaction(i, j) * eliminatedDofs[0](bb, 0);
                    }
                }
            }
        }
    }

protected:
    gsMatrix<T> m_physGrad;
    gsMatrix<T> m_mLocalAdvection; // Local matrix advection
    gsMatrix<T> m_mLocalDiffusion; // Local matrix diffusion
    gsMatrix<T> m_mLocalReaction; // Local matrix reaction

    using Base::m_basisData;
    using Base::m_actives;
    using Base::m_mapper;
    using Base::m_patchIndex;
    using Base::m_advectionCoeff;

    using Base::getADREvaluator;
};

} // namespace gismo

