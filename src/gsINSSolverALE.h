/** @file gsINSSolverALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSAssemblerALE.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsIncompressibleFlow/src/gsBarrierPatchDynamic.h>
#include <gsModeling/gsBarrierPatch.h>

namespace gismo
{

/// @brief ALE solver for unsteady incompressible Navier-Stokes equations
/// @tparam T coefficient type
/// @tparam MatOrder sparse matrix storage order
template <class T = real_t, int MatOrder = RowMajor>
class gsINSSolverUnsteadyALE : public gsINSSolverUnsteady<T, MatOrder>
{
public:
    typedef gsINSSolverUnsteady<T, MatOrder> Base;
    
protected:
    using Base::m_paramsPtr;
    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_time;
    using Base::m_timeStepSize;
    
    // Flag to indicate if ALE is active
    bool m_isALEActive;
    
    // Mesh update callback function
    std::function<gsMatrix<T>(T)> m_meshUpdateFunc;
    
    // Flag to indicate if mesh optimization is enabled
    bool m_useMeshOptimization;
    
    // Flag to indicate if dynamic boundary mapping is enabled
    bool m_useDynamicBoundaryMapping;
    
    // gsBarrierPatch options
    gsOptionList m_meshOptOptions;
    
    // Store original mesh
    gsMultiPatch<T> m_originalMesh;
    
    // Store previous displacement for incremental update
    gsMatrix<T> m_previousDisp;
    
    // Rotation parameters (if known)
    T m_rotationPeriod;
    gsVector<T> m_rotationCenter;
    
public:
    /// @brief Constructor
    gsINSSolverUnsteadyALE(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    : Base(paramsPtr), m_isALEActive(false), m_useMeshOptimization(false),
      m_useDynamicBoundaryMapping(false), m_rotationPeriod(10.0) // Default 10s period
    {
        // Replace assembler with ALE version
        delete m_assemblerPtr;
        m_assemblerPtr = new gsINSAssemblerUnsteadyALE<T, MatOrder>(m_paramsPtr);
        
        Base::initMembers();
        
        // Initialize mesh optimization options
        m_meshOptOptions.addInt("Verbose", "Verbosity level for mesh optimization", 0);
        m_meshOptOptions.addInt("ParamMethod", "Parametrization method", 1);
        m_meshOptOptions.addInt("AAPreconditionType", "AA precondition type", 0);
        
        // Initialize rotation center to domain center
        m_rotationCenter.resize(2);
        m_rotationCenter << 0.5, 0.5;
    }
    
    /// @brief Destructor
    virtual ~gsINSSolverUnsteadyALE() { }
    
    /// @brief Activate/deactivate ALE formulation
    void setALEActive(bool active) 
    { 
        m_isALEActive = active;
        getALEAssembler()->setALEActive(active);
        
        if (active && m_originalMesh.nPatches() == 0)
        {
            // Store original mesh on first activation
            m_originalMesh = m_assemblerPtr->getPatches();
            m_previousDisp.setZero(getALEAssembler()->getUdofs() * m_originalMesh.targetDim(), 1);
        }
    }
    
    /// @brief Check if ALE is active
    bool isALEActive() const { return m_isALEActive; }
    
    /// @brief Set mesh update function
    /// @param[in] func function that returns mesh displacement given time t
    void setMeshUpdateFunction(std::function<gsMatrix<T>(T)> func)
    {
        m_meshUpdateFunc = func;
    }
    
    /// @brief Update mesh velocity field
    /// @param[in] meshVel mesh velocity coefficients
    void setMeshVelocity(const gsMatrix<T>& meshVel)
    {
        getALEAssembler()->setMeshVelocity(meshVel);
    }
    
    /// @brief Update mesh displacement and velocity
    /// @param[in] meshDisp new mesh displacement
    void updateMesh(const gsMatrix<T>& meshDisp)
    {
        getALEAssembler()->updateMesh(meshDisp);
    }
    
    /// @brief Get mesh velocity field
    gsField<T> getMeshVelocityField() const
    {
        return getALEAssembler()->getMeshVelocityField();
    }
    
    /// @brief Get mesh displacement field
    gsField<T> getMeshDisplacementField() const
    {
        return getALEAssembler()->getMeshDisplacementField();
    }
    
    /// @brief Enable/disable mesh optimization using gsBarrierPatch
    void setMeshOptimization(bool enable) { m_useMeshOptimization = enable; }
    
    /// @brief Check if mesh optimization is enabled
    bool isMeshOptimizationEnabled() const { return m_useMeshOptimization; }
    
    /// @brief Enable/disable dynamic boundary mapping for rotating domains
    void setDynamicBoundaryMapping(bool enable) { m_useDynamicBoundaryMapping = enable; }
    
    /// @brief Check if dynamic boundary mapping is enabled
    bool isDynamicBoundaryMappingEnabled() const { return m_useDynamicBoundaryMapping; }
    
    /// @brief Get mesh optimization options
    gsOptionList& getMeshOptOptions() { return m_meshOptOptions; }
    
    /// @brief Get mesh optimization options (const version)
    const gsOptionList& getMeshOptOptions() const { return m_meshOptOptions; }
    
    /// @brief Set rotation parameters for dynamic boundary mapping
    /// @param[in] period rotation period in seconds
    /// @param[in] center rotation center coordinates
    void setRotationParameters(T period, const gsVector<T>& center)
    {
        m_rotationPeriod = period;
        m_rotationCenter = center;
    }
    
    /// @brief Get current rotation angle based on time
    T getCurrentRotationAngle() const
    {
        return m_time * 2.0 * M_PI / m_rotationPeriod;
    }
    
    /// @brief Get ALE assembler
    gsINSAssemblerUnsteadyALE<T, MatOrder>* getALEAssembler() const
    {
        return dynamic_cast<gsINSAssemblerUnsteadyALE<T, MatOrder>*>(m_assemblerPtr);
    }
    
    /// @brief Get assembler (override to provide proper type)
    virtual gsINSAssemblerUnsteady<T, MatOrder>* getAssembler() const override
    {
        return dynamic_cast<gsINSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    }
    
    /// @brief Perform next iteration step with ALE
    virtual void nextIteration() override
    {
        if (m_isALEActive && m_meshUpdateFunc)
        {
            // Apply the displacement to the actual mesh geometry
            // This also updates the incremental displacement in the assembler
            applyMeshDisplacement();
            
            // Apply mesh optimization if enabled
            if (m_useMeshOptimization)
            {
                optimizeMesh();
            }
        }
        
        // Call base class iteration
        Base::nextIteration();
    }
    
protected:
    /// @brief Apply mesh displacement to the actual geometry
    void applyMeshDisplacement()
    {
        if (m_originalMesh.nPatches() == 0)
            return;
            
        // Get cumulative displacement from mesh update function
        gsMatrix<T> cumulativeDisp = m_meshUpdateFunc(m_time + m_timeStepSize);
        
        // Get current mesh from assembler
        gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());
        
        // Reset to original mesh and apply cumulative displacement
        for (size_t p = 0; p < patches.nPatches(); ++p)
        {
            // Reset to original mesh
            patches.patch(p).coefs() = m_originalMesh.patch(p).coefs();
            
            // Apply cumulative displacement
            const index_t nCoefs = patches.patch(p).coefsSize();
            const index_t dim = patches.patch(p).targetDim();
            
            // Extract displacement for this patch
            const gsDofMapper& mapper = getALEAssembler()->getMappers()[0];
            for (index_t i = 0; i < nCoefs; ++i)
            {
                if (mapper.is_free(i, p))
                {
                    index_t idx = mapper.index(i, p);
                    for (index_t d = 0; d < dim; ++d)
                    {
                        patches.patch(p).coef(i, d) += cumulativeDisp(idx + d * getALEAssembler()->getUdofs());
                    }
                }
            }
        }
        
        // Update incremental displacement for ALE assembler
        gsMatrix<T> incrementalDisp = cumulativeDisp - m_previousDisp;
        getALEAssembler()->updateMesh(incrementalDisp);
        m_previousDisp = cumulativeDisp;
    }
    
    /// @brief Optimize mesh using gsBarrierPatch
    void optimizeMesh()
    {
        try {
            // Get current mesh from assembler
            gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());
            
            // Use dynamic version if enabled and we have rotation
            if (m_useDynamicBoundaryMapping && m_meshUpdateFunc)
            {
                // Apply gsBarrierPatchDynamic for rotating domains
                gsBarrierPatchDynamic<2, T> opt;
                
                // Enable dynamic boundary mapping
                opt.setDynamicBoundaryMapping(true);
                
                // Set rotation parameters
                T rotationAngle = getCurrentRotationAngle();
                opt.setRotationAngle(rotationAngle);
                opt.setRotationCenter(m_rotationCenter);
                
                // Compute optimized mesh
                patches = opt.compute(patches, m_meshOptOptions);
                
                gsInfo << "Applied mesh optimization with dynamic boundary mapping\n";
            }
            else
            {
                // Use standard gsBarrierPatch
                gsBarrierPatch<2, T> opt(patches, false);
                opt.options() = m_meshOptOptions;
                opt.compute();
                patches = opt.result();
                
                gsInfo << "Applied standard mesh optimization\n";
            }
            
            // Update mesh in assembler
            getALEAssembler()->updateMeshAfterOptimization(patches);
        }
        catch (const std::exception& e) {
            gsWarn << "Mesh optimization failed: " << e.what() << "\n";
            gsWarn << "Continuing with unoptimized mesh.\n";
        }
    }
    
    /// @brief Returns the name of the class
    virtual std::string getName() override { return "gsINSSolverUnsteadyALE"; }
};

/// @brief Helper class for FSI simulations using ALE
/// This class provides utilities for coupling fluid and structure solvers
template <class T = real_t>
class gsFSIHelper
{
protected:
    // Fluid-structure interface patches
    std::vector<boxSide> m_fluidInterfaceSides;
    std::vector<boxSide> m_structInterfaceSides;
    
    // Interface mappers
    gsDofMapper m_fluidInterfaceMapper;
    gsDofMapper m_structInterfaceMapper;
    
public:
    /// @brief Constructor
    gsFSIHelper() { }
    
    /// @brief Set fluid-structure interface
    /// @param[in] fluidSides fluid domain interface sides
    /// @param[in] structSides structure domain interface sides
    void setInterface(const std::vector<boxSide>& fluidSides,
                     const std::vector<boxSide>& structSides)
    {
        m_fluidInterfaceSides = fluidSides;
        m_structInterfaceSides = structSides;
    }
    
    /// @brief Transfer displacement from structure to fluid mesh
    /// @param[in] structDisp structure displacement at interface
    /// @param[out] fluidMeshDisp fluid mesh displacement
    /// @param[in] fluidBasis fluid domain basis
    /// @param[in] structBasis structure domain basis
    void transferDisplacement(const gsMatrix<T>& structDisp,
                             gsMatrix<T>& fluidMeshDisp,
                             const gsMultiBasis<T>& fluidBasis,
                             const gsMultiBasis<T>& structBasis)
    {
        // TODO: Implement displacement transfer
        // This would involve:
        // 1. Evaluate structure displacement at interface quadrature points
        // 2. Project onto fluid mesh basis functions
        // 3. Extend displacement into fluid domain (e.g., using harmonic extension)
        
        GISMO_NO_IMPLEMENTATION
    }
    
    /// @brief Transfer traction from fluid to structure
    /// @param[in] fluidStress fluid stress tensor
    /// @param[in] fluidPressure fluid pressure
    /// @param[out] structTraction traction on structure interface
    /// @param[in] normal interface normal vector
    void transferTraction(const gsMatrix<T>& fluidStress,
                         const gsMatrix<T>& fluidPressure,
                         gsMatrix<T>& structTraction,
                         const gsMatrix<T>& normal)
    {
        // TODO: Implement traction transfer
        // Traction = -p*n + tau*n
        // where tau is the viscous stress tensor
        
        GISMO_NO_IMPLEMENTATION
    }
    
    /// @brief Compute mesh velocity from displacement history
    /// @param[in] dispNew new displacement
    /// @param[in] dispOld old displacement
    /// @param[in] dt time step
    /// @return mesh velocity
    static gsMatrix<T> computeMeshVelocity(const gsMatrix<T>& dispNew,
                                          const gsMatrix<T>& dispOld,
                                          T dt)
    {
        return (dispNew - dispOld) / dt;
    }
    
    /// @brief Extend interface displacement into domain using harmonic extension
    /// @param[in] interfaceDisp displacement at interface
    /// @param[out] domainDisp extended displacement in domain
    /// @param[in] basis domain basis
    /// @param[in] interfaceSides interface boundary sides
    void harmonicExtension(const gsMatrix<T>& interfaceDisp,
                          gsMatrix<T>& domainDisp,
                          const gsMultiBasis<T>& basis,
                          const std::vector<boxSide>& interfaceSides)
    {
        // TODO: Implement harmonic extension
        // Solve Laplace equation with interface displacement as BC
        
        GISMO_NO_IMPLEMENTATION
    }
};

} // namespace gismo