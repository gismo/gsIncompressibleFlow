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
        // Initialize visitor so it's ready for assembly
        m_visitorUUnonlinALE->initialize();
        // Propagate current solution if available
        if (this->m_currentVelField.nPieces() > 0)
            m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);
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

    /// @brief Update the DOF mappers in all visitors (including ALE visitor)
    virtual void updateDofMappers() override
    {
        Base::updateDofMappers();
        if (m_visitorUUnonlinALE)
            m_visitorUUnonlinALE->updateDofMappers(this->m_dofMappers);
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

    /// @brief Set mesh velocity coefficients directly (top rows correspond to velocity DOFs)
    /// @param[in] meshVel mesh velocity coefficients for velocity DOFs
    void setMeshVelocity(const gsMatrix<T>& meshVel)
    {
        const index_t udofs = this->m_dofMappers[0].freeSize();
        const index_t pdofs = this->m_dofMappers[1].freeSize();
        const index_t pshift = m_tarDim * udofs;
        const index_t fullDofs = pshift + pdofs;

        m_meshVelCoefs.resize(fullDofs, 1);
        m_meshVelCoefs.setZero();
        if (meshVel.size() >= pshift)
            m_meshVelCoefs.topRows(pshift) = meshVel.topRows(pshift);
        else if (meshVel.size() > 0)
            m_meshVelCoefs.topRows(meshVel.rows()) = meshVel;
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
    
    /// @brief Override update: refresh Dirichlet data BEFORE assembly so u=w holds on moving walls
    virtual void update(const gsMatrix<T>& sol, bool updateSol = true) override
    {
        // If ALE is active, prepare mesh-velocity field for boundary updates and ALE terms
        if (m_isALEActive && m_visitorUUnonlinALE)
        {
            delete m_tempMeshVelField;
            m_tempMeshVelField = new gsField<T>(this->constructSolution(m_meshVelCoefs, 0));
            m_visitorUUnonlinALE->setMeshVelocityField(m_tempMeshVelField);

            // Refresh velocity Dirichlet values for Dirichlet sides
            // 1) Recompute ddof via standard path (interpolation/L2 proj/homogeneous)
            // 2) For homogeneous Dirichlet (function==nullptr), override with mesh velocity: u = w
            if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
            {
                // (1) recompute baseline Dirichlet values
                this->computeDirichletDofs(0, 0, this->m_ddof[0]);

                // (2) override homogeneous Dirichlet by mesh velocity (automatic no-slip on moving walls)
                const gsBoundaryConditions<T>& bcs = this->getBCs();
                const gsMultiBasis<T>& mbasis = this->getBases().at(0);
                const gsDofMapper& mapper = this->m_dofMappers[0];

                for (typename gsBoundaryConditions<T>::const_iterator it = bcs.dirichletBegin(); it != bcs.dirichletEnd(); ++it)
                {
                    if (it->unknown() != 0) continue; // only velocity
                    if (!it->isHomogeneous()) continue; // only homogeneous entries are auto-set to w

                    const index_t patchID = it->patch();
                    const gsBasis<T>& basis = mbasis.piece(patchID);

                    // Boundary DOF indices on this side (in patch coefficient numbering)
                    gsMatrix<index_t> boundary = basis.boundary(it->side());

                    // Evaluate mesh velocity directly at the boundary basis anchors
                    typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
                    gsMatrix<T> anch = h->anchors();
                    gsMatrix<T> wVals = m_tempMeshVelField->value(anch, patchID);

                    // Interpolate onto boundary basis anchors to get boundary coefficients
                    typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(wVals);
                    const gsMatrix<T>& dVals = geo->coefs();

                    // Write into eliminated Dirichlet vector with size safety
                    const index_t nrows = static_cast<index_t>(dVals.rows());
                    const index_t nbnd  = static_cast<index_t>(boundary.size());
                    if (nrows != nbnd)
                        gsWarn << "[ALE/Asm] Boundary coefficients count mismatch: dVals.rows()="
                               << nrows << ", boundary.size()=" << nbnd << ". Clamping.\n";
                    const index_t ncopy = std::min(nrows, nbnd);
                    for (index_t i = 0; i < ncopy; ++i)
                    {
                        const index_t ii = mapper.bindex(boundary.at(i), patchID);
                        if (ii >= 0 && ii < this->m_ddof[0].rows())
                            this->m_ddof[0].row(ii) = dVals.row(i);
                    }
                }
            }
        }

        // Update current solution field used e.g. for convection linearization
        Base::updateCurrentSolField(sol, updateSol);

        // Reassemble both linear and nonlinear parts to incorporate updated Dirichlet data
        this->assembleLinearPart();
        this->assembleNonlinearPart();
        this->fillSystem();

        // Update visitor with current solution field (optional, for compatibility)
        if (m_isALEActive && m_visitorUUnonlinALE && this->m_currentVelField.nPieces() > 0)
            m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);
    }
    
    
    /// @brief Assemble the nonlinear part; if ALE active and mesh is moving, use ALE visitor ((u-w)·∇u)
    virtual void assembleNonlinearPart() override
    {
        // matrix and rhs cleaning
        this->m_blockUUnonlin_comp.resize(this->m_udofs, this->m_udofs);
        this->m_blockUUnonlin_whole.resize(this->m_pshift, this->m_pshift);
        this->m_rhsUnonlin.setZero();

        // If ALE is active but mesh velocity is essentially zero, fall back to standard formulation
        bool useALE = m_isALEActive && m_visitorUUnonlinALE;
        if (useALE)
        {
            // Check infinity norm of mesh velocity coefficients (only velocity part)
            const index_t udofs = this->m_dofMappers[0].freeSize();
            const index_t pshift = m_tarDim * udofs;
            T velInf = (m_meshVelCoefs.size() >= pshift) ? m_meshVelCoefs.topRows(pshift).template lpNorm<gsEigen::Infinity>() : T(0);
            if (velInf < static_cast<T>(1e-14))
                useALE = false;
        }

        if (useALE)
        {
            // Assemble ALE convection into component block (udofs x udofs)
            // and let fillGlobalMat_UU replicate across components (robust path)
            this->m_blockUUnonlin_comp.reserve(gsVector<index_t>::Constant(this->m_blockUUnonlin_comp.outerSize(), this->m_nnzPerOuterU));
            this->assembleBlock(*m_visitorUUnonlinALE, 0, this->m_blockUUnonlin_comp, this->m_rhsUnonlin);
        }
        else
        {
            // Fallback to standard formulation
            Base::assembleNonlinearPart();
        }
    }
};

} // namespace gismo
