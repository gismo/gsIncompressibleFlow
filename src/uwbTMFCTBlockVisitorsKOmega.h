/** @file uwbTMFCTBlockVisitorsKOmega.h

    Author(s): E. Turnerova

*/

#pragma once
#include "uwbINSBlockVisitors.h"
#include "uwbTMBlockVisitorsKOmega.h"
#include "uwbTMEvaluators.h"

namespace gismo
{
// ============================================================= BLOCK FCT M ============================================================= //
template <class T>
class uwbTMFCTBlockMVisitorKOmega : public uwbTMBlockMVisitorKOmega<T>
{
public:
    typedef uwbTMBlockMVisitorKOmega<T> Base;

public:

    uwbTMFCTBlockMVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.setZero(numActKOmega, numActKOmega);

        gsMatrix<T> tempMat;
        tempMat.setZero(numActKOmega, numActKOmega);

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            tempMat.noalias() += weight * (m_basisValsKO.col(k) * m_basisValsKO.col(k).transpose());
        }

        for (index_t i = 0; i < numActKOmega; i++)
            for (index_t j = 0; j < numActKOmega; j++)
                m_localMat.coeffRef(i, i) += tempMat(i, j);
    }

protected:

    // Basis values
    using Base::m_activesKO;
    using Base::m_basisValsKO;
    using Base::m_localMat; // Local matrix
};


// ============================================================= BLOCK FCT N_TMkomega ============================================================= //
template <class T>
class uwbTMFCTBlockNVisitorKOmega : public uwbTMBlocksNVisitorKOmega<T>
{

public:
    typedef uwbTMBlocksNVisitorKOmega<T> Base;

public:

    uwbTMFCTBlockNVisitorKOmega(std::vector< gsDofMapper >& dofMappers, const T viscosity) :
        Base(dofMappers, viscosity) { }

    inline void assemble(gsGeometryEvaluator<T> & geoEval,
        gsVector<T> const      & quWeights)
    {
        const index_t numActKOmega = m_activesKO.rows();

        m_localMat.setZero(numActKOmega, numActKOmega);//local N_k(u^n) = local N_omega(u^n)

        gsMatrix<T> mNumDiffusion; //numerical diffusion
        mNumDiffusion.setZero(numActKOmega, numActKOmega);

        const gsMatrix<T> & basisValsKO = m_basisDataKO[0];
        const gsMatrix<T> & bGrads_komega = m_basisDataKO[1];

        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            const T weight = quWeights(k) * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads_komega, m_physGradKO);

            m_localMat.noalias() += weight * (basisValsKO.col(k) * (m_solUVals.col(k).transpose() * m_physGradKO));
        }

        for (index_t i = 0; i < numActKOmega; i++)
        {
            T numDiffDiag = 0.;
            for (index_t j = 0; j < numActKOmega; j++)
            {
                if (j != i)
                {
                    mNumDiffusion.coeffRef(i, j) = math::max(-m_localMat(i, j), -m_localMat(j, i));
                    mNumDiffusion.coeffRef(i, j) = math::max(mNumDiffusion(i, j), 0.);
                    numDiffDiag -= mNumDiffusion(i, j);
                }
            }
            mNumDiffusion.coeffRef(i, i) = numDiffDiag;
        }
        m_localMat = -m_localMat - mNumDiffusion;
    }

protected:

    // Basis values
    using Base::m_activesKO;
    using Base::m_basisDataKO;
    using Base::m_solUVals;
    using Base::m_physGradKO;
    using Base::m_localMat; // Local matrix
};

} // namespace gismo

