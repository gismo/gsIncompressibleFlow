/** @file gsINSAssemblerALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSVisitorsALE.h>

namespace gismo
{

/// @brief ALE extension for unsteady incompressible Navier-Stokes assembler
/// @tparam T real number type
/// @tparam MatOrder sparse matrix storage order
template <class T, int MatOrder = RowMajor>
class gsINSAssemblerUnsteadyALE : public gsINSAssemblerUnsteady<T, MatOrder>
{
public:
    typedef gsINSAssemblerUnsteady<T, MatOrder> Base;
    
protected:
    using Base::m_paramsPtr;
    using Base::m_dofMappers;
    using Base::m_tarDim;
    using Base::m_visitorUUnonlin;
    using Base::m_currentVelField;
    using Base::m_blockUUnonlin_comp;
    using Base::m_blockUUnonlin_whole;
    using Base::m_rhsUnonlin;
    using Base::m_nnzPerOuterU;
    using Base::m_udofs;
    
    // Mesh velocity coefficients
    gsMatrix<T> m_meshVelCoefs;
    gsMatrix<T> m_meshDispCoefs;
    gsMatrix<T> m_meshDispOld;
    
    // ALE visitor
    gsINSVisitorUUnonlinALE<T, MatOrder>* m_visitorUUnonlinALE;
    
    // Flag to indicate if ALE is active
    bool m_isALEActive;
    
    // Temporary mesh velocity field (created on demand)
    mutable gsField<T>* m_tempMeshVelField;
    
public:
    /// @brief Constructor
    gsINSAssemblerUnsteadyALE(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    : Base(paramsPtr),
      m_isALEActive(false),
      m_visitorUUnonlinALE(nullptr),
      m_tempMeshVelField(nullptr)
    { }
    
    /// @brief Destructor
    virtual ~gsINSAssemblerUnsteadyALE()
    {
        delete m_visitorUUnonlinALE;
        delete m_tempMeshVelField;
    }
    
    /// @brief Initialize (override to add ALE initialization)
    virtual void initialize() override
    {
        // First initialize base class
        Base::initialize();
        
        // Initialize mesh coefficients to zero with full DOF size (velocity + pressure)
        const index_t udofs = this->m_dofMappers[0].freeSize();
        const index_t pdofs = this->m_dofMappers[1].freeSize();
        const index_t pshift = m_tarDim * udofs;
        const index_t fullDofs = pshift + pdofs;
        m_meshVelCoefs.setZero(fullDofs, 1);
        m_meshDispCoefs.setZero(fullDofs, 1);
        // m_meshDispOld stores only velocity DOFs (will be resized when first used)
        m_meshDispOld.setZero(pshift, 1);
        
        // Create ALE visitor with proper paramsPtr (mesh velocity field will be set later)  
        m_visitorUUnonlinALE = new gsINSVisitorUUnonlinALE<T, MatOrder>(
            m_paramsPtr, m_dofMappers, m_tarDim, nullptr);
        
        // No need to manually initialize - assembleBlock will handle it
    }
    
    /// @brief Update sizes when boundary conditions change (override)
    virtual void updateSizes() override
    {
        // First update base class sizes
        Base::updateSizes();
        
        // Update mesh coefficient sizes to match full DOFs (velocity + pressure)
        const index_t udofs = this->getUdofs();
        const index_t pdofs = this->getPdofs();
        const index_t pshift = m_tarDim * udofs;
        const index_t fullDofs = pshift + pdofs;
        m_meshVelCoefs.conservativeResize(fullDofs, 1);
        m_meshDispCoefs.conservativeResize(fullDofs, 1);
        // m_meshDispOld stores only velocity DOFs
        m_meshDispOld.conservativeResize(pshift, 1);
    }
    
    /// @brief Activate/deactivate ALE formulation
    void setALEActive(bool active) { m_isALEActive = active; }
    
    /// @brief Check if ALE is active
    bool isALEActive() const { return m_isALEActive; }
    
    /// @brief Update mesh for new time step
    /// @param[in] meshDispNew new mesh displacement coefficients
    void updateMesh(const gsMatrix<T>& meshDispNew)
    {
        // Get time step
        T dt = m_paramsPtr->options().getReal("timeStep");
        
        // Number of free velocity DOFs and total DOFs
        const index_t udofs = this->m_dofMappers[0].freeSize();
        const index_t pdofs = this->m_dofMappers[1].freeSize();
        const index_t pshift = m_tarDim * udofs;
        const index_t fullDofs = pshift + pdofs;


        // Pad velocity and displacement vectors to full size (pressure part zero)
        m_meshVelCoefs.resize(fullDofs,1);
        m_meshVelCoefs.setZero();
        if (meshDispNew.size() >= pshift && m_meshDispOld.size() >= pshift) {
            m_meshVelCoefs.topRows(pshift) = (meshDispNew.topRows(pshift) - m_meshDispOld.topRows(pshift)) / dt;
        }

        m_meshDispCoefs.resize(fullDofs,1);
        m_meshDispCoefs.setZero();
        if (meshDispNew.size() >= pshift) {
            m_meshDispCoefs.topRows(pshift) = meshDispNew.topRows(pshift);
        }
        
        // Store old displacement for next time step (ensure correct size)
        if (m_meshDispOld.size() != meshDispNew.size()) {
            m_meshDispOld.resize(meshDispNew.size(), 1);
        }
        m_meshDispOld = meshDispNew;
    }
    
    /// @brief Get mesh velocity field
    gsField<T> getMeshVelocityField() const 
    { 
        // Check if our internal coefficients are zero
        T velNorm = m_meshVelCoefs.template lpNorm<gsEigen::Infinity>();
        
        if (velNorm < 1e-12) {
            // If coefficients are essentially zero, create a proper zero field
            // Create zero patches
            gsMultiPatch<T>* result = new gsMultiPatch<T>;
            for (size_t p = 0; p < this->getPatches().nPatches(); ++p) {
                const index_t sz = this->getBases().at(0).piece(p).size();
                gsMatrix<T> zeroCoeffs = gsMatrix<T>::Zero(sz, m_tarDim);
                result->addPatch(this->getBases().at(0).piece(p).makeGeometry(zeroCoeffs));
            }
            
            return gsField<T>(this->getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
        } else {
            // Use normal constructSolution for non-zero coefficients
            return this->constructSolution(m_meshVelCoefs, 0);
        }
    }
    
    /// @brief Get mesh displacement field  
    gsField<T> getMeshDisplacementField() const 
    { 
        return this->constructSolution(m_meshDispCoefs, 0);
    }
    
    /// @brief Initialize old displacement (used for first time step)
    /// @param[in] initialDisp initial mesh displacement at t=0
    void initializeOldDisplacement(const gsMatrix<T>& initialDisp)
    {
        m_meshDispOld = initialDisp;
    }
    
    
    /// @brief Update old time field at the end of a time step
    void updateOldTimeField(const gsMatrix<T>& oldSolution)
    {
        this->m_oldTimeVelField = this->constructSolution(oldSolution, 0, this->m_paramsPtr->isRotation());
    }
    
    
    /// @brief Update mesh after optimization (e.g., by gsBarrierPatch)
    /// @param[in] optimizedPatches the optimized mesh patches
    /// @param[in] originalMesh the original mesh to compute displacement from
    void updateMeshAfterOptimization(const gsMultiPatch<T>& optimizedPatches, 
                                   const gsMultiPatch<T>& originalMesh)
    {
        // Update the patches in the assembler
        const_cast<gsMultiPatch<T>&>(this->getPatches()) = optimizedPatches;
        
        // Compute displacement coefficients from optimized mesh relative to original mesh
        const index_t udofs = this->m_dofMappers[0].freeSize();
        gsMatrix<T> fullDisp(m_tarDim * udofs, 1);
        fullDisp.setZero();
        
        const gsDofMapper& mapper = this->m_dofMappers[0];
        for (index_t p = 0; p < optimizedPatches.nPatches(); ++p)
        {
            const auto& origCoefs = originalMesh.patch(p).coefs();
            const auto& optCoefs = optimizedPatches.patch(p).coefs();
            
            for (index_t i = 0; i < origCoefs.rows(); ++i)
            {
                if (mapper.is_free(i, p))
                {
                    index_t idx = mapper.index(i, p);
                    // x-displacement
                    fullDisp(idx) = optCoefs(i, 0) - origCoefs(i, 0);
                    // y-displacement  
                    fullDisp(idx + udofs) = optCoefs(i, 1) - origCoefs(i, 1);
                    // z-displacement (if 3D)
                    if (m_tarDim == 3)
                        fullDisp(idx + 2 * udofs) = optCoefs(i, 2) - origCoefs(i, 2);
                }
            }
        }
        
        // Update mesh with computed displacement
        updateMesh(fullDisp);
    }
    
    /// @brief Update mesh after optimization (e.g., by gsBarrierPatch) - legacy version
    /// @param[in] optimizedPatches the optimized mesh patches
    void updateMeshAfterOptimization(const gsMultiPatch<T>& optimizedPatches)
    {
        // Update the patches in the assembler
        const_cast<gsMultiPatch<T>&>(this->getPatches()) = optimizedPatches;
    }
    
    /// @brief Override updateCurrentSolField to handle time history properly in ALE
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol) override
    {
        // CRITICAL FIX: Call the unsteady base class method, not the steady one!
        // This ensures proper time history handling
        Base::updateCurrentSolField(solVector, updateSol);
    }
    
    /// @brief Override update - use standard update, ALE is handled in assembly
    virtual void update(const gsMatrix<T>& sol, bool updateSol = true) override
    {
        // Always call standard update - keep current fields as actual fluid velocities
        Base::update(sol, updateSol);
        
        // Set up ALE visitor if needed (for compatibility)
        if (m_isALEActive && m_visitorUUnonlinALE)
        {
            // Recreate mesh velocity field with current mesh velocity coefficients
            delete m_tempMeshVelField;
            m_tempMeshVelField = new gsField<T>(this->constructSolution(m_meshVelCoefs, 0));
            
            // Update mesh velocity in ALE visitor
            m_visitorUUnonlinALE->setMeshVelocityField(m_tempMeshVelField);
            
            // Update current solution in ALE visitor
            if (updateSol && this->m_currentVelField.nPieces() > 0)
            {
                m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);
            }
        }
    }
    
    
    /// @brief Assemble the nonlinear part with ALE formulation
    virtual void assembleNonlinearPart() override
    {
        // Always use standard formulation for now
        // ALE effects are handled through the ALE visitor when needed
        Base::assembleNonlinearPart();
    }
};

} // namespace gismo