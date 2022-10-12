/** @file uwbINSdgBlockVisitorsBnd.h

Author(s) : H.Hornikova
*/

#pragma once
#include "uwbINSBlockVisitors.h"

namespace gismo {

// ============================================================= PARENT ============================================================= //

// boundary block visitor

template <class T>
class uwbINSBlockVisitorBnd : public uwbINSBlockVisitor<T>
{

public:
    typedef uwbINSBlockVisitor<T> Base;

public:

    uwbINSBlockVisitorBnd(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags,
        boxSide& side)
    {
        const gsBasis<T>& basis = basisRefs.front();

        m_dim = basis.dim();
        m_patchIndex = patchIndex;

        m_side = side;

        const int dir = m_side.direction();
        gsVector<int> numQuadNodes(basis.dim());
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = 2 * basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        this->initializeSpecific(evFlags);
    }

protected:

    gsVector<T> m_unormal;
    boxSide m_side;

    // members from uwbINSBlockVisitor
    using Base::m_dim;
    using Base::m_patchIndex;

}; // class uwbINSBlockVisitorBnd

   // ============================================================= Robin BCs - PCD ============================================================= //

template <class T>
class uwbINSBlockVisitorRobinPCD : public uwbINSBlockVisitorBnd<T>
{

public:
    typedef uwbINSBlockVisitorBnd<T> Base;

public:

    uwbINSBlockVisitorRobinPCD(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const & bases,
        gsGeometryEvaluator<T> & geoEval,
        const gsMatrix<T>      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        bases.back().active_into(quNodes.col(0), m_actives);
        const index_t numActive = m_actives.rows();

        bases.back().eval_into(quNodes, m_bVals);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate solution on element nodes
        m_solUVals = m_solU.value(quNodes, m_patchIndex);

        // Initialize local matrix/rhs
        m_localMat.setZero(numActive, numActive);
    }

    inline void assemble(gsDomainIterator<T>    & element,
        gsGeometryEvaluator<T> & geoEval,
        const gsVector<T>      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, m_side, m_unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * m_unormal.norm();

            // Compute the unit normal vector 
            m_unormal.normalize();

            T coef = - m_solUVals.col(k).transpose() * m_unormal;

            m_localMat.noalias() += weight * coef * (m_bVals.col(k) * m_bVals.col(k).transpose());
        }
    }

    inline void localToGlobal(gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs)
    {

        // Local Dofs to global dofs
        m_Pmap.localToGlobal(m_actives, m_patchIndex, m_actives);
        const index_t numActP = m_actives.rows();

        for (index_t i = 0; i < numActP; ++i)
        {
            const int ii = m_actives(i);
            if (m_Pmap.is_free_index(ii))
            {
                for (index_t j = 0; j < numActP; ++j)
                {
                    const int jj = m_actives(j);
                    if (m_Pmap.is_free_index(jj))
                        blockMatrix.coeffRef(ii, jj) += m_localMat(i, j);
                }
            }
        }
    }

protected:
    gsMatrix<index_t> m_actives;
    gsMatrix<T> m_bVals, m_solUVals, m_localMat;

    // members from uwbINSBlockVisitorBnd
    using Base::m_side;
    using Base::m_unormal;

    // members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_solU;
    using uwbINSBlockVisitor<T>::m_Pmap;    

}; // class uwbINSBlockVisitorRobinPCD

// ============================================================= Parent - Nitsche ============================================================= //

template <class T>
class uwbINSBlockVisitorNitsche : public uwbINSBlockVisitorBnd<T>
{

public:
    typedef uwbINSBlockVisitorBnd<T> Base;

public:

    uwbINSBlockVisitorNitsche(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags,
        const boundary_condition<T> & dirCond)
    {
        Base::initialize(basisRefs, patchIndex, rule, evFlags, dirCond.side());

        const gsBasis<T>& basis = basisRefs.front();

        const int deg = basis.maxDegree();
        m_penalty = (deg + basis.dim()) * (deg + 1) * T(2.5);
        m_pDirFcn = dirCond.function().get();
    }

protected:

    T m_penalty;
    const gsFunction<T>* m_pDirFcn;

}; // class uwbINSBlockVisitorNitsche


// ============================================================= Velocity block - Nitsche ============================================================= //

template <class T>
class uwbINSBlockVisitorNitscheA : public uwbINSBlockVisitorNitsche<T>
{

public:
    typedef uwbINSBlockVisitorNitsche<T> Base;

public:

    uwbINSBlockVisitorNitscheA(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const & bases,
        gsGeometryEvaluator<T> & geoEval,
        const gsMatrix<T>      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        bases.front().active_into(quNodes.col(0), m_actives);
        const index_t numActive = m_actives.rows();

        // Evaluate basis values and derivatives on element
        bases.front().evalAllDers_into(quNodes, 1, m_basisData);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Dirichlet data
        m_pDirFcn->eval_into(geoEval.values(), m_dirData);

        // Initialize local matrix/rhs
        m_localMat.setZero(numActive, numActive);
        m_localRhs.setZero(numActive, m_pDirFcn->targetDim());

    }

    inline void assemble(gsDomainIterator<T>    & element,
        gsGeometryEvaluator<T> & geoEval,
        const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & bGrads = m_basisData[1];
        

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            const gsMatrix<T> bVals = m_basisData[0].col(k);

            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, m_side, m_unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * m_unormal.norm();

            // Compute the unit normal vector 
            m_unormal.normalize();

            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, m_physGrad);

            // Get penalty parameter
            const T h = element.getCellSize();
            const T gamma = m_penalty / (0 != h ? h : 1);

            // Sum up quadrature point evaluations
            m_localRhs.noalias() -= weight * m_viscosity * ((m_physGrad.transpose() * m_unormal - gamma * bVals)
                * m_dirData.col(k).transpose());

            m_localMat.noalias() -= weight * m_viscosity * (bVals * m_unormal.transpose() * m_physGrad
                + (bVals * m_unormal.transpose() * m_physGrad).transpose()
                - gamma * bVals * bVals.transpose());
        }
    }

    void localToGlobal( gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs)
    {
        const index_t usz = m_Umap.freeSize();
        m_Umap.localToGlobal(m_actives, m_patchIndex, m_actives);
        const index_t numActive = m_actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t i = 0; i != numActive; ++i)
        {
            const unsigned ii = m_actives(i);

            for (index_t s = 0; s != m_dim; s++)
                blockRhs(ii + s*usz, 0) += m_localRhs(i,s);

            for (index_t j = 0; j != numActive; ++j)
            {
                const unsigned jj = m_actives(j);
                blockMatrix.coeffRef(ii, jj) += m_localMat(i, j);
            }
        }

    }

    protected:
        gsMatrix<index_t> m_actives;
        std::vector<gsMatrix<T> > m_basisData;
        gsMatrix<T> m_physGrad;
        gsMatrix<T> m_dirData;

        gsMatrix<T> m_localMat;
        gsMatrix<T> m_localRhs;

        // members from uwbINSBlockVisitorNitsche
        using Base::m_penalty;
        using Base::m_pDirFcn;

        // members from uwbINSBlockVisitorBnd
        using uwbINSBlockVisitorBnd<T>::m_side;
        using uwbINSBlockVisitorBnd<T>::m_unormal;

        // members from uwbINSBlockVisitor
        using uwbINSBlockVisitor<T>::m_patchIndex;
        using uwbINSBlockVisitor<T>::m_viscosity;
        using uwbINSBlockVisitor<T>::m_Umap;
        using uwbINSBlockVisitor<T>::m_dim;

}; // class uwbINSBlockVisitorNitscheA

// ============================================================= Pressure block - Nitsche ============================================================= //

template <class T>
class uwbINSBlockVisitorNitscheB : public uwbINSBlockVisitorNitsche<T>
{
public:
    typedef uwbINSBlockVisitorNitsche<T> Base;

public:

    uwbINSBlockVisitorNitscheB(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_JACOBIAN;
    }

    // Evaluate on element.
    inline void evaluate(gsBasisRefs<T> const & bases,
        gsGeometryEvaluator<T> & geoEval,
        const gsMatrix<T>      & quNodes)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        bases.front().active_into(quNodes.col(0), m_activesU);
        bases.back().active_into(quNodes.col(0), m_activesP);
        const index_t numActU = m_activesU.rows();
        const index_t numActP = m_activesP.rows();

        bases.front().eval_into(quNodes, m_bValsU);
        bases.back().eval_into(quNodes, m_bValsP);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate the Dirichlet data
        m_pDirFcn->eval_into(geoEval.values(), m_dirData);

        // Initialize local matrix/rhs
        m_localMatVec.resize(m_dim);
        m_localRhs.setZero(numActP, 1);

        for (index_t i = 0; i != m_dim; ++i)
            m_localMatVec[i].setZero(numActP, numActU); //local_B_i


    }

    inline void assemble(gsDomainIterator<T>    & element,
        gsGeometryEvaluator<T> & geoEval,
        const gsVector<T>      & quWeights)
    {

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            geoEval.outerNormal(k, m_side, m_unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * m_unormal.norm();

            // Compute the unit normal vector 
            m_unormal.normalize();

            m_localRhs.noalias() -= weight * m_bValsP.col(k) * (m_dirData.col(k).transpose() * m_unormal);

            // Sum up quadrature point evaluations
            for (index_t i = 0; i != m_dim; ++i)
            {
                const T w = weight * m_unormal(i);

                m_localMatVec[i].noalias() -= w * (m_bValsP.col(k) * m_bValsU.col(k).transpose());
            }
        }
    }

    void localToGlobal(gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs)
    {
        const index_t usz = m_Umap.freeSize();
        const index_t pshift = m_dim * usz;

        m_Umap.localToGlobal(m_activesU, m_patchIndex, m_activesU);
        const index_t numActU = m_activesU.rows();

        m_Pmap.localToGlobal(m_activesP, m_patchIndex, m_activesP);
        const index_t numActP = m_activesP.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t i = 0; i != numActP; ++i)
        {
            const unsigned ii = m_activesP(i);

            if (m_Pmap.is_free_index(ii))
            {

                blockRhs(pshift + ii, 0) += m_localRhs(i, 0);

                for (index_t j = 0; j != numActU; ++j)
                {
                    const unsigned jj = m_activesU(j);

                    for (index_t s = 0; s != m_dim; s++)
                        blockMatrix.coeffRef(ii, jj + s * usz) += m_localMatVec[s](i, j);
                }
            }
        }

    }

protected:
    gsMatrix<index_t> m_activesU, m_activesP;
    gsMatrix<T> m_bValsU, m_bValsP;
    gsMatrix<T> m_dirData;

    std::vector<gsMatrix<T> > m_localMatVec;
    gsMatrix<T> m_localRhs;

    // members from uwbINSBlockVisitorNitsche
    using Base::m_pDirFcn;

    // members from uwbINSBlockVisitorBnd
    using uwbINSBlockVisitorBnd<T>::m_side;
    using uwbINSBlockVisitorBnd<T>::m_unormal;

    // members from uwbINSBlockVisitor
    using uwbINSBlockVisitor<T>::m_patchIndex;
    using uwbINSBlockVisitor<T>::m_viscosity;
    using uwbINSBlockVisitor<T>::m_Umap;
    using uwbINSBlockVisitor<T>::m_Pmap;
    using uwbINSBlockVisitor<T>::m_dim;

}; // class uwbINSBlockVisitorNitscheB


} // end namespace gismo
