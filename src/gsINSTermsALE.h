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

        //Calculate mesh velocity at quadrature points
        computeMeshVelocity(mapData);

        // If there is no evaluation point, return early
        if (mapData.points.cols()==0)
            return;
        
        // Compute solution velocity values
        this->computeCoeffSolU(mapData);
        
        // Get relative velocity (u - u_mesh)
        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;
        
        const index_t numTestActive = testFunData[0].rows();
        const index_t numActive = trialFunData[0].rows();

        // 若派生数据为空，直接返回
        if (trialFunData[1].size()==0)
            return;
        
        // Initialize local matrix
        localMat.setZero(numTestActive, numActive);
        
        // Get physical gradients from mapData
        // mapData.values[0] contains the inverse Jacobian matrices
        const gsMatrix<T>& invJac = mapData.values[0];
        
        // Transform gradients and compute convection term
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            // Get Jacobian for this quadrature point
            gsMatrix<T> invJ = invJac.reshapeCol(k, mapData.dim.first, mapData.dim.second);
            
            // Transform trial function gradients to physical space
            gsMatrix<T> physGrad = trialFunData[1].reshapeCol(k, numActive, mapData.dim.first) * invJ.transpose();
            
            // Compute (u_rel · \nabla \phi_trial)
            gsMatrix<T> convection(numActive, 1);
            convection.setZero();
            
            for (index_t d = 0; d < m_tarDim; ++d)
            {
                convection += relativeVel(d, k) * physGrad.col(d);
            }
            
            // Add contribution: weight * (\phi_test * (u_rel · \nabla\phi_trial))
            localMat.noalias() += quWeights[k] * mapData.measures(k) * 
                                  testFunData[0].col(k) * convection.transpose();
        }
    }
};

} // namespace gismo