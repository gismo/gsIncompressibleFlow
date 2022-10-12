/** @file uwbINSBlockVisitors.h

    Author(s): J. Sourek, H. Hornikova
*/

#pragma once

namespace gismo
{
// ========================================================== SUPER CLASS ========================================================== //
template <class T>
class uwbVisitorBase
{
public:
    uwbVisitorBase()
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
class uwbINSBlockVisitor : public uwbVisitorBase<T>
{
protected:

    const T m_viscosity;
    const gsDofMapper & m_Umap;
    const gsDofMapper & m_Pmap;
    gsField<T> m_solU;
    index_t m_dim; // Velocity vector dimension
    int m_patchIndex;
    gsMatrix<T> m_diffusionCoeff;

public:
    uwbINSBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        m_viscosity(viscosity),
        m_Umap(dofMappers.front()),
        m_Pmap(dofMappers.back())
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
        this->m_bSolutionSet = true;
    }

protected:
    virtual inline void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }
};

// ============================================================= BLOCK A ============================================================= //
template <class T>
class uwbINSBlockAsymVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockAsymVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        basisRefs.front().deriv_into(quNodes, m_basisGradsU);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);//local_A

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, m_basisGradsU, m_physGradU);

            // Local block A
            m_localMat.template triangularView<Eigen::Upper>() += weight * m_viscosity * (m_physGradU.transpose() * m_physGradU);
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
                            rhs(ii + s*usz, 0) -= // assuming single rhs
                            m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_basisGradsU;

    gsMatrix<T> m_physGradU;
    gsMatrix<T> m_localMat; // Local matrix

protected:
    using Base::m_viscosity;
    using Base::m_Umap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK B ============================================================= //
template <class T>
class uwbINSBlocksBVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlocksBVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }


    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), m_activesP);

        basisRefs.front().deriv_into(quNodes, m_basisGradsU);
        basisRefs.back().eval_into(quNodes, m_basisValsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();
        const index_t numActP = m_activesP.rows();

        m_localMat.resize(m_dim);
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].setZero(numActP, numActU);//local_B_i

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, m_basisGradsU, m_physGradU);

            // Local blocks B_i
            for (index_t i = 0; i != m_dim; ++i)
                m_localMat[i].noalias() += weight * (m_basisValsP.col(k) * m_physGradU.row(i));
        }

    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {

        const index_t usz = m_Umap.freeSize();
        const index_t ps = m_dim*usz;

        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

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
                            sysBlock.coeffRef(jj, ii + s * usz) += m_localMat[s](j, i);
                    }
                    else //m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        for (index_t s = 0; s < m_dim; ++s)
                            rhs(s * usz + ii, 0) += m_localMat[s](j, i) * eliminatedDofs[1](bb, 0);
                    }
                }

            }
            else // m_Umap.is_boundary_index(ii)
            {
                const int bb = m_Umap.global_to_bindex(ii);
                for (index_t k = 0; k < numActP; ++k)
                {
                    const int kk = m_activesP(k);
                    if (m_Pmap.is_free_index(kk)) 
                    {
                        T tmp = m_localMat[0](k, i)*eliminatedDofs[0](bb, 0);
                        for (index_t s = 1; s != m_dim; ++s)
                            tmp += m_localMat[s](k, i) * eliminatedDofs[0](bb, s);
                        rhs(ps + kk, 0) -= tmp;// assuming single rhs
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesU;
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsU;
    gsMatrix<T> m_basisValsP;

    gsMatrix<T> m_physGradU;
    std::vector<gsMatrix<T> > m_localMat;

protected:

    using Base::m_Umap;
    using Base::m_Pmap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK C ============================================================= //
template <class T>
class uwbINSBlocksCVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlocksCVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }


    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element  nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        basisRefs.back().active_into(quNodes.col(0), m_activesP);

        basisRefs.front().eval_into(quNodes, m_basisValsU);
        basisRefs.back().deriv_into(quNodes, m_basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();
        const index_t numActP = m_activesP.rows();

        m_localMat.resize(m_dim);
        for (index_t i = 0; i != m_dim; ++i)
            m_localMat[i].setZero(numActP, numActU);//local_B_i

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, m_basisGradsP, m_physGradP);

            // Local blocks B_i
            for (index_t i = 0; i != m_dim; ++i)
                m_localMat[i].noalias() += weight * (m_physGradP.row(i).transpose() * m_basisValsU.col(k).transpose());
        }

    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {

        const index_t usz = m_Umap.freeSize();
        const index_t ps = m_dim*usz;

        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

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
                            sysBlock.coeffRef(jj, ii + s * usz) += m_localMat[s](j, i);
                    }
                    else //m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        for (index_t s = 0; s < m_dim; ++s)
                            rhs(s * usz + ii, 0) += m_localMat[s](j, i) * eliminatedDofs[1](bb, 0);
                    }
                }

            }
            else // m_Umap.is_boundary_index(ii)
            {
                const int bb = m_Umap.global_to_bindex(ii);
                for (index_t k = 0; k < numActP; ++k)
                {
                    const int kk = m_activesP(k);
                    if (m_Pmap.is_free_index(kk)) {
                        T tmp = m_localMat[0](k, i)*eliminatedDofs[0](bb, 0);
                        for (index_t s = 1; s != m_dim; ++s)
                            tmp += m_localMat[s](k, i) * eliminatedDofs[0](bb, s);
                        rhs(ps + kk, 0) -= tmp;// assuming single rhs
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesU;
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsP;
    gsMatrix<T> m_basisValsU;

    gsMatrix<T> m_physGradP;
    std::vector<gsMatrix<T> > m_localMat;

protected:

    using Base::m_Umap;
    using Base::m_Pmap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK N ============================================================= //
template <class T>
class uwbINSBlockNVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockNVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        //const index_t numActU = m_activesU.rows();

        basisRefs.front().evalAllDers_into(quNodes, 1, m_basisDataU);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);
    }

    inline void assemble( gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);

        const gsMatrix<T> & basisValsU = m_basisDataU[0];
        const gsMatrix<T> & basisGradsU = m_basisDataU[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsU, m_physGradU);

            // Local block
            m_localMat.noalias() += weight * (basisValsU.col(k) * (m_solUVals.col(k).transpose() * m_physGradU));
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
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    }
                    else // m_Umap.is_boundary_index(jj)
                    {
                        const int bb = m_Umap.global_to_bindex(jj);
                        for (index_t s = 0; s != m_dim; ++s)
                            rhs(ii + s*usz, 0) -= m_localMat(i, j) * eliminatedDofs[0](bb, s); // assuming single rhs
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    gsMatrix<T> solActUCoeffs;
    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_physGradU;
    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_Umap;
    using Base::m_solU;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK M ============================================================= //
template <class T>
class uwbINSBlockMVisitor : public uwbINSBlockVisitor<T>
{
public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockMVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.front().active_into(quNodes.col(0), m_activesU);
        // Evaluate basis functions on element  nodes
        m_basisDataU.resize(1);

        basisRefs.front().eval_into(quNodes, m_basisDataU.front());

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);
    }

    inline void assemble( gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localMat.setZero(numActU, numActU);//local_A

        const gsMatrix<T> & basisValsU = m_basisDataU.front();

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            m_localMat.template triangularView<Eigen::Upper>() += weight * (basisValsU.col(k) * basisValsU.col(k).transpose());
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
                            rhs(ii + s*usz, 0) -= // assuming single rhs
                            m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;

    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_Umap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK Ap ============================================================= //
template <class T>
class uwbINSBlockApVisitor : public uwbINSBlockVisitor<T>
{
public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockApVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.back().active_into(quNodes.col(0), m_activesP);
        // Evaluate basis functions on element  nodes

        basisRefs.back().deriv_into(quNodes, m_basisGradsP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);
    }

    inline void assemble( gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActP = m_activesP.rows();

        m_localMat.setZero(numActP, numActP);//local_A

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            geoEval.transformGradients(k, m_basisGradsP, m_physGradP);

            m_localMat.template triangularView<Eigen::Upper>() += weight * (m_physGradP.transpose() * m_physGradP);
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {
        
        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        for (index_t i = 0; i < numActP; ++i)
        {
            const int ii = m_activesP(i);
            if (m_Pmap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActP; ++j) // Build A-part of the matrix
                {
                    const int jj = m_activesP(j);
                    if (m_Pmap.is_free_index(jj))
                    {
                        if (j >= i) {
                            sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                        }
                    }
                    else // m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[1](bb, 0);
                    }
                }
            }
        }
         
    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsP;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_Pmap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK Np ============================================================= //
template <class T>
class uwbINSBlockNpVisitor : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockNpVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.back().active_into(quNodes.col(0), m_activesP);

        //const index_t numActU = m_activesU.rows();

        basisRefs.back().evalAllDers_into(quNodes, 1, m_basisDataP);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActP = m_activesP.rows();

        m_localMat.setZero(numActP, numActP);

        const gsMatrix<T> & basisValsP = m_basisDataP[0];
        const gsMatrix<T> & basisGradsP = m_basisDataP[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, basisGradsP, m_physGradP);

            // Local block
            m_localMat.noalias() += weight * (basisValsP.col(k) * (m_solUVals.col(k).transpose() * m_physGradP));
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {

        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        for (index_t i = 0; i < numActP; ++i)
        {
            const int ii = m_activesP(i);
            if (m_Pmap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActP; ++j)
                {
                    const int jj = m_activesP(j);
                    if (m_Pmap.is_free_index(jj))
                    {
                        sysBlock.coeffRef(ii, jj) += m_localMat(i, j);
                    }
                    else // m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(i, j) * eliminatedDofs[1](bb, 0);
                    }
                }
            }
        }

    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesP;
    std::vector<gsMatrix<T> > m_basisDataP;
    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_Pmap;
    using Base::m_solU;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= BLOCK Mp ============================================================= //
template <class T>
class uwbINSBlockMpVisitor : public uwbINSBlockVisitor<T>
{
public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockMpVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE;
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        basisRefs.back().active_into(quNodes.col(0), m_activesP);
        // Evaluate basis functions on element  nodes

        basisRefs.back().eval_into(quNodes, bVals_p);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActP = m_activesP.rows();

        m_localMat.setZero(numActP, numActP);//local_A

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            m_localMat.template triangularView<Eigen::Upper>() += weight * (bVals_p.col(k) * bVals_p.col(k).transpose());
        }
    }

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    {

        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        for (index_t i = 0; i < numActP; ++i)
        {
            const int ii = m_activesP(i);
            if (m_Pmap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActP; ++j) // Build A-part of the matrix
                {
                    const int jj = m_activesP(j);
                    if (m_Pmap.is_free_index(jj))
                    {
                        if (j >= i) {
                            sysBlock.coeffRef(std::min(ii, jj), std::max(ii, jj)) += m_localMat(i, j);
                        }
                    }
                    else // m_Pmap.is_boundary_index(jj)
                    {
                        const int bb = m_Pmap.global_to_bindex(jj);
                        rhs(ii, 0) -= m_localMat(std::min(i, j), std::max(i, j)) * eliminatedDofs[1](bb, 0);
                    }
                }
            }
        }

    }

protected:

    // Basis values
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> bVals_p;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat; // Local matrix

protected:

    using Base::m_Pmap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

// ============================================================= RHS ============================================================= //
template <class T>
class uwbINSRhsVisitor : public uwbINSBlockVisitor<T>
{
public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSRhsVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_pRhsFcn(NULL)
    { }

    void setRhsFunction(const gsFunction<T>& rhsFcn)
    { m_pRhsFcn = &rhsFcn; }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE;

        GISMO_ASSERT(m_pRhsFcn != NULL, "No rhs function set in uwbINSRhsVisitor!");
    }

    inline void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes)
    {
        // Evaluate basis functions on element nodes
        basisRefs.front().active_into(quNodes.col(0), m_activesU);

        basisRefs.front().eval_into(quNodes, m_basisValsU);

        // Evaluate Geometry on element nodes
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand side at the geometry points
        m_pRhsFcn->eval_into( geoEval.values(), m_rhsVals ); //PERFORMANCE LOSS IF OPENMP ENABLED
    }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActU = m_activesU.rows();

        m_localRhsU.setZero(numActU, m_rhsVals.rows());

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian.
            const T weight = quWeights(k) * geoEval.measure(k);

            // Right-hand side
            m_localRhsU.noalias() += weight * (m_basisValsU.col(k) *  m_rhsVals.col(k).transpose());
        }
    }

    inline void localToGlobal( gsMatrix<T> & rhs )
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
                for (index_t s = 0; s != m_dim; ++s) 
                    rhs(ii + s*usz, 0) += m_localRhsU(i, s);
            }
        }
    }

protected:

    const gsFunction<T>* m_pRhsFcn;

    // Basis values
    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_basisValsU;
    gsMatrix<T> m_rhsVals;

    gsMatrix<T> m_localRhsU;

protected:

    using Base::m_Umap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

};

} // namespace gismo

