/** @file gsINSTermsALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{

/// @brief ALE convection term for incompressible Navier-Stokes: ((u - u_mesh) · ∇φ_trial) * φ_test
/// @tparam T coefficient type
template <class T>
class gsINSTerm_ALEConvection : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;
    
protected:
    using Base::m_solUVals;
    
    // Target dimension
    index_t m_tarDim;
    
    // Mesh velocity values at quadrature points
    gsMatrix<T> m_meshVelVals;
    
    // Mesh velocity field
    const gsField<T>* m_meshVelField;
    
public:
    /// @brief Constructor
    gsINSTerm_ALEConvection(const gsField<T>* meshVelField = nullptr, index_t tarDim = 2) 
    : Base(), m_meshVelField(meshVelField), m_tarDim(tarDim)
    { 
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_DERIV;
    }
    
    /// @brief Set mesh velocity field
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
    }
    
    /// @brief Set target dimension
    void setTargetDim(index_t tarDim) { m_tarDim = tarDim; }
    
    /// @brief Evaluate mesh velocity at quadrature points
    void computeMeshVelocity(const gsMapData<T>& mapData)
    {
        if (m_meshVelField != nullptr)
        {
            // Evaluate mesh velocity field at quadrature points
            m_meshVelVals = m_meshVelField->value(mapData.points, mapData.patchId);
        }
        else
        {
            // If no mesh velocity field is set, assume zero mesh velocity
            m_meshVelVals.setZero(m_tarDim, mapData.points.cols());
        }
    }
    
    /// @brief Compute coefficient times solution velocity
    virtual void computeCoeffSolU(const gsMapData<T>& mapData) override
    {
        Base::computeCoeffSolU(mapData);
    }
    
    /// @brief Assemble the ALE convection term
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        // If there is no evaluation point, return early
        if (mapData.points.cols()==0)
            return;

        // If derived data is empty, return directly
        if (trialFunData.size() < 2 || trialFunData[1].size()==0)
            return;

        const index_t numTestActive = testFunData[0].rows();
        const index_t numActive = trialFunData[0].rows();

        // Initialize local matrix
        localMat.setZero(numTestActive, numActive);

        // Calculate mesh velocity at quadrature points
        computeMeshVelocity(mapData);

        // Compute solution velocity values
        this->computeCoeffSolU(mapData);

        // Safety check: ensure dimensions match
        if (m_solUVals.rows() != m_meshVelVals.rows() ||
            m_solUVals.cols() != m_meshVelVals.cols())
        {
            gsWarn << "ALE convection: dimension mismatch between solution velocity ("
                   << m_solUVals.rows() << "x" << m_solUVals.cols()
                   << ") and mesh velocity ("
                   << m_meshVelVals.rows() << "x" << m_meshVelVals.cols()
                   << "). Skipping assembly.\n";
            return;
        }

        // Get relative velocity (u - u_mesh)
        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;
        
        // Assemble using transformed gradients
        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
        gsMatrix<T> trialFunPhysGrad;
        const index_t nQuPoints = quWeights.rows();
        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);
            transformGradients(mapData, k, trialFunData[1], trialFunPhysGrad);
            gsMatrix<T> convection = trialFunPhysGrad * relativeVel.col(k);
            localMat.noalias() += weight * (testFunData[0].col(k) * convection.transpose());
 
        }
    }
};

} // namespace gismo
