/** @file gsINSSolverALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSAssemblerALE.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsIncompressibleFlow/src/gsBarrierPatchDynamic.h>
#include <gsModeling/gsBarrierPatch.h>
#include <gsElasticity/src/gsALE.h>

namespace gismo
{

/// @brief Boundary mesh velocity function for FSI interface
/// Computes mesh velocity as w = (X_cur - X_prev) / dt
/// IMPORTANT: Input is PHYSICAL coordinates (not parametric)
/// This is required for proper Dirichlet BC evaluation on moving boundaries
template <class T = real_t>
class gsFSIMeshVelocityFunction : public gsFunction<T>
{
protected:
    const gsGeometry<T>* m_cur;   ///< Current geometry
    const gsGeometry<T>* m_prev;  ///< Previous geometry
    T m_dt;                        ///< Time step

public:
    /// @brief Constructor
    gsFSIMeshVelocityFunction(const gsGeometry<T>* cur, const gsGeometry<T>* prev, T dt)
    : m_cur(cur), m_prev(prev), m_dt(dt) {}

    /// @brief Set current geometry
    void setCurrent(const gsGeometry<T>* cur) { m_cur = cur; }

    /// @brief Set previous geometry
    void setPrevious(const gsGeometry<T>* prev) { m_prev = prev; }

    /// @brief Set time step
    void setDt(T dt) { m_dt = dt; }

    /// @brief Domain dimension (physical space)
    short_t domainDim() const { return 2; }

    /// @brief Target dimension (velocity vector)
    short_t targetDim() const { return 2; }

    /// @brief Evaluate mesh velocity at PHYSICAL points
    /// @param x Physical coordinates (2 x N)
    /// @param result Mesh velocity at these points (2 x N)
    void eval_into(const gsMatrix<T>& x, gsMatrix<T>& result) const
    {
        const index_t N = x.cols();
        result.resize(2, N);

        for (index_t i = 0; i < N; ++i)
        {
            // Input is physical coordinate, need to invert to parametric space
            gsMatrix<T> phys(2, 1);
            phys.col(0) = x.col(i);

            // Invert current geometry to find parametric coordinates
            gsMatrix<T> uv;
            m_cur->invertPoints(phys, uv);

            // Evaluate previous position at same parametric location
            gsMatrix<T> xPrev = m_prev->eval(uv);

            // Compute mesh velocity: w = (x_current - x_previous) / dt
            if (m_dt > 1e-14)
                result.col(i) = (phys.col(0) - xPrev.col(0)) / m_dt;
            else
                result.col(i).setZero();
        }
    }

    /// @brief Clone this function
    typename gsFunction<T>::uPtr clone() const
    {
        return memory::make_unique(new gsFSIMeshVelocityFunction<T>(*this));
    }
};

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

    // Flag to indicate if gsALE mesh deformation is enabled
    bool m_useGsALE;

    // gsALE method selection
    ale_method::method m_aleMethod;

    // gsBarrierPatch options
    gsOptionList m_meshOptOptions;

    // gsALE options
    gsOptionList m_gsALEOptions;

    // If true, freeze outer box boundary (xmin/xmax/ymin/ymax) during geometry update
    bool m_fixOuterBoundary;
    T m_fixOuterTol;
    
    // Store original mesh
    gsMultiPatch<T> m_originalMesh;
    
    // Rotation parameters (if known)
    T m_rotationPeriod;
    gsVector<T> m_rotationCenter;

    // FSI interface sides that should remain free during mesh optimization
    std::vector<patchSide> m_fsiInterfaceSides;

    // Boundary displacement field for gsALE (stores interface displacements)
    gsMultiPatch<T> m_boundaryDisplacement;

    // Boundary interface for gsALE (maps FSI interface sides)
    gsBoundaryInterface m_aleInterface;

    // Storage for boundary mesh velocity functions (for automatic no-slip on FSI interface)
    std::vector<std::shared_ptr<gsFunction<T>>> m_fsiBoundaryVelFuncs;

    // Pointer to boundary conditions (to allow automatic setup)
    gsBoundaryConditions<T>* m_bcInfoPtr;

    // Previous geometry for computing mesh velocity
    gsMultiPatch<T> m_prevGeometry;

public:
    /// @brief Constructor
    gsINSSolverUnsteadyALE(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    : Base(paramsPtr), m_isALEActive(false), m_useMeshOptimization(false),
      m_useDynamicBoundaryMapping(false), m_useGsALE(false),
      m_aleMethod(ale_method::TINE), m_rotationPeriod(10.0), m_bcInfoPtr(nullptr),
      m_fixOuterBoundary(false), m_fixOuterTol(static_cast<T>(1e-8))
    {
        delete m_assemblerPtr;
        m_assemblerPtr = new gsINSAssemblerUnsteadyALE<T, MatOrder>(m_paramsPtr);
        Base::initMembers();

        // Mesh optimization options
        m_meshOptOptions.addInt("Verbose", "Verbosity level", 0);
        m_meshOptOptions.addInt("ParamMethod", "Parametrization method", 1);
        m_meshOptOptions.addInt("AAPreconditionType", "AA precondition type", 0);

        // gsALE options
        m_gsALEOptions.addReal("PoissonsRatio", "Poisson's ratio", 0.3);
        m_gsALEOptions.addReal("LocalStiff", "Local stiffening degree", 2.3);
        m_gsALEOptions.addSwitch("Check", "Check mesh bijectivity", true);

        m_rotationCenter.resize(2);
        m_rotationCenter << 0.5, 0.5;
    }
    
    /// @brief Destructor
    virtual ~gsINSSolverUnsteadyALE() { }
    
    /// @brief Activate/deactivate ALE formulation
    void setALEActive(bool active) {
        m_isALEActive = active;
        getALEAssembler()->setALEActive(active);
        if (active && m_originalMesh.nPatches() == 0) {
            m_originalMesh = m_assemblerPtr->getPatches();
            if (m_meshUpdateFunc) {
                gsMatrix<T> initialDisp = m_meshUpdateFunc(m_time);
                getALEAssembler()->initializeOldDisplacement(initialDisp);
            }
        }
    }

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

        // If using gsALE, also update boundary displacement field
        if (m_useGsALE && m_boundaryDisplacement.nPatches() > 0)
        {
            updateBoundaryDisplacement(meshDisp);
        }
    }

    /// @brief Update boundary displacement field for gsALE from interface displacement
    /// @param[in] interfaceDisp interface displacement vector
    /// @param[in] dofLocations mapping from interface points to DOF locations (patch, coefIdx)
    void updateBoundaryDisplacementFromInterface(
        const gsMatrix<T>& interfaceDisp,
        const std::vector<std::pair<index_t,index_t>>& dofLocations)
    {
        if (!m_useGsALE || m_boundaryDisplacement.nPatches() == 0)
            return;

        const index_t N = interfaceDisp.cols();
        GISMO_ASSERT(N == static_cast<index_t>(dofLocations.size()),
                     "Mismatch between interface points and DOF locations");

        // Set boundary displacement coefficients directly
        for (index_t k = 0; k < N; ++k)
        {
            const index_t p = dofLocations[k].first;
            const index_t coefIdx = dofLocations[k].second;

            gsMatrix<T>& coefs = m_boundaryDisplacement.patch(p).coefs();
            coefs(coefIdx, 0) = interfaceDisp(0, k);  // x displacement
            coefs(coefIdx, 1) = interfaceDisp(1, k);  // y displacement
        }
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

    /// @brief Enable/disable gsALE mesh deformation
    /// @param[in] enable enable/disable gsALE
    /// @param[in] method ALE method to use (default: TINE)
    void setGsALE(bool enable, ale_method::method method = ale_method::TINE)
    {
        m_useGsALE = enable;
        m_aleMethod = method;
        if (enable)
        {
            gsInfo << "[gsINSSolverALE] gsALE mesh deformation enabled with method " << method << "\n";
            // When using gsALE, typically disable gsBarrierPatch optimization
            if (m_useMeshOptimization)
            {
                gsWarn << "[gsINSSolverALE] gsBarrierPatch optimization is also enabled. "
                       << "Consider disabling one method for consistency.\n";
            }
        }
    }

    /// @brief Check if gsALE mesh deformation is enabled
    bool isGsALEEnabled() const { return m_useGsALE; }

    /// @brief Enable/disable dynamic boundary mapping for rotating domains
    void setDynamicBoundaryMapping(bool enable) { m_useDynamicBoundaryMapping = enable; }

    /// @brief Check if dynamic boundary mapping is enabled
    bool isDynamicBoundaryMappingEnabled() const { return m_useDynamicBoundaryMapping; }

    /// @brief Get mesh optimization options
    gsOptionList& getMeshOptOptions() { return m_meshOptOptions; }

    /// @brief Get mesh optimization options (const version)
    const gsOptionList& getMeshOptOptions() const { return m_meshOptOptions; }

    /// @brief Get gsALE options
    gsOptionList& getGsALEOptions() { return m_gsALEOptions; }

    /// @brief Get gsALE options (const version)
    const gsOptionList& getGsALEOptions() const { return m_gsALEOptions; }

    /// @brief Freeze outer (box) boundary during geometry update
    void setFixOuterBoundary(bool enable, T tol = static_cast<T>(1e-8))
    {
        m_fixOuterBoundary = enable;
        m_fixOuterTol = tol;
    }
    
    /// @brief Set rotation parameters for dynamic boundary mapping
    /// @param[in] period rotation period in seconds
    /// @param[in] center rotation center coordinates
    void setRotationParameters(T period, const gsVector<T>& center)
    {
        m_rotationPeriod = period;
        m_rotationCenter = center;
    }

    /// @brief Set FSI interface sides that should remain free during mesh optimization
    /// @param[in] interfaceSides vector of patchSide objects representing FSI interface
    void setFSIInterfaceSides(const std::vector<patchSide>& interfaceSides)
    {
        m_fsiInterfaceSides = interfaceSides;
        gsInfo << "[gsINSSolverALE] Set " << m_fsiInterfaceSides.size() << " FSI interface sides\n";

        // Also setup gsALE interface if using gsALE
        if (m_useGsALE)
        {
            setupGsALEInterface();
        }
    }

    /// @brief Automatically setup no-slip boundary conditions on FSI interface
    /// This method creates boundary mesh velocity functions and adds Dirichlet BC u=w
    /// on the FSI interface sides. Must be called AFTER geometry is set up.
    /// @param[in] bcInfo reference to boundary conditions object (will be modified)
    /// @param[in] dt time step size (for mesh velocity computation)
    void setupFSINoSlipBC(gsBoundaryConditions<T>& bcInfo, T dt)
    {
        if (m_fsiInterfaceSides.empty())
        {
            gsWarn << "[gsINSSolverALE] No FSI interface sides specified. "
                   << "Call setFSIInterfaceSides() first.\n";
            return;
        }

        // Store pointer to boundary conditions for later updates
        m_bcInfoPtr = &bcInfo;

        // Get current geometry patches and store as previous (initially same)
        const gsMultiPatch<T>& patches = m_assemblerPtr->getPatches();
        m_prevGeometry = patches;  // Store initial geometry

        // Create boundary mesh velocity functions for each FSI interface side
        m_fsiBoundaryVelFuncs.clear();
        m_fsiBoundaryVelFuncs.reserve(m_fsiInterfaceSides.size());

        for (const auto& ps : m_fsiInterfaceSides)
        {
            // Create a custom function that computes mesh velocity from geometry change
            // This function will be updated each time step in nextIteration()
            auto meshVelFunc = std::make_shared<gsFSIMeshVelocityFunction<T>>(
                &patches.patch(ps.patch), &m_prevGeometry.patch(ps.patch), dt);

            m_fsiBoundaryVelFuncs.push_back(meshVelFunc);

            // Add Dirichlet BC: u = w (no-slip on moving boundary)
            bcInfo.addCondition(ps.patch, ps.side(), condition_type::dirichlet,
                               meshVelFunc.get(), 0);
        }

        gsInfo << "[gsINSSolverALE] Automatically applied Dirichlet BC u=w (no-slip) on "
               << m_fsiInterfaceSides.size() << " FSI interface sides\n";
    }

    /// @brief Update FSI boundary mesh velocity functions
    /// This should be called each time step before solving to update the mesh velocity
    /// @param[in] dt time step size
    void updateFSIBoundaryVelocity(T dt)
    {
        if (m_fsiBoundaryVelFuncs.empty() || m_fsiInterfaceSides.empty())
            return;

        const gsMultiPatch<T>& curGeo = m_assemblerPtr->getPatches();

        // Update each boundary velocity function with current and previous geometry
        for (size_t i = 0; i < m_fsiBoundaryVelFuncs.size(); ++i)
        {
            const index_t p = m_fsiInterfaceSides[i].patch;

            // Cast to the specific function type to access update methods
            auto* velFunc = dynamic_cast<gsFSIMeshVelocityFunction<T>*>(
                m_fsiBoundaryVelFuncs[i].get());

            if (velFunc)
            {
                velFunc->setCurrent(&curGeo.patch(p));
                velFunc->setPrevious(&m_prevGeometry.patch(p));
                velFunc->setDt(dt);
            }
        }
    }

    /// @brief Store current geometry as previous (for next time step)
    /// This should be called after each successful time step
    void updatePreviousGeometry()
    {
        if (m_prevGeometry.nPatches() > 0)
        {
            m_prevGeometry = m_assemblerPtr->getPatches();
        }
    }

    /// @brief Initialize gsALE boundary displacement field
    /// This should be called after the mesh is set up and before starting time integration
    void initializeGsALE()
    {
        if (!m_useGsALE) return;

        const gsMultiPatch<T>& patches = m_assemblerPtr->getPatches();

        // Create boundary displacement field (initially zero, same structure as geometry)
        m_boundaryDisplacement.clear();
        for (index_t p = 0; p < patches.nPatches(); ++p)
        {
            typename gsGeometry<T>::uPtr patch = patches.patch(p).clone();
            patch->coefs().setZero();
            m_boundaryDisplacement.addPatch(std::move(patch));
        }
        m_boundaryDisplacement.computeTopology();

        gsInfo << "[gsINSSolverALE] Initialized gsALE with " << m_boundaryDisplacement.nPatches()
               << " patches\n";
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
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
        
        // CRITICAL: Follow exact same sequence as non-ALE solver, with FSI BC refresh before assembly

        // Ensure FSI boundary velocity functions reflect current and previous geometry
        // This must happen BEFORE updating the assembler so that Dirichlet u=w is assembled
        updateFSIBoundaryVelocity(this->m_timeStepSize);

        // Step 1: Apply ALE mesh displacement for the CURRENT time first
        // This ensures the system is assembled on the updated geometry
        if (m_isALEActive && m_meshUpdateFunc)
        {
            applyMeshDisplacement();
            // If geometry changed here, refresh FSI boundary velocity again
            updateFSIBoundaryVelocity(this->m_timeStepSize);

            // Choose mesh deformation method
            if (m_useGsALE)
            {
                // Use gsALE mesh deformation
                applyGsALEDeformation();
            }
            else if (m_useMeshOptimization)
            {
                optimizeMesh();
            }
            gsWriteOutputLine(this->m_outFile, "[ALE] Mesh update done.", this->m_fileOutput, this->m_dispOutput);
        }

        gsWriteOutputLine(this->m_outFile, "[ALE] Updating assembler...", this->m_fileOutput, this->m_dispOutput);
        // Step 2: Update assembler on the updated geometry
        this->updateAssembler();
        
        // Step 3: Initialize if first iteration (must be AFTER updateAssembler)
        if (!this->m_iterationNumber)
            this->initIteration();
        
        // Step 4: Start solving with current solution as initial guess
        gsMatrix<T> tmpSolution = this->m_solution;
        
        this->applySolver(tmpSolution);
        this->writeSolChangeRelNorm(this->m_solution, tmpSolution);
        
        // Step 5: Picard iterations (exactly like non-ALE solver)
        index_t picardIter = 0;
        T relNorm = this->solutionChangeRelNorm(this->m_solution, tmpSolution);
        
        gsWriteOutputLine(this->m_outFile, "        [u, p] Picard's iterations...", this->m_fileOutput, this->m_dispOutput);
        
        while((relNorm > this->m_innerTol) && (picardIter < this->m_innerIter))
        {
            gsWriteOutput(this->m_outFile, "         ", this->m_fileOutput, this->m_dispOutput);
            
            gsMatrix<T> oldSol = tmpSolution;
            
            this->updateAssembler(tmpSolution, false);
            this->applySolver(tmpSolution);
            this->writeSolChangeRelNorm(oldSol, tmpSolution);
            
            relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
            picardIter++;
        }
        
        // Step 6: Update solution and time (exactly like non-ALE solver)
        this->m_solution = tmpSolution;

        // CRITICAL FIX: Explicitly update time history for next iteration
        // This must be done AFTER m_solution is updated
        getALEAssembler()->updateCurrentSolField(tmpSolution, true);

        this->m_time += this->m_timeStepSize;
        this->m_avgPicardIter += picardIter;
        this->m_iterationNumber++;

        // Note: Previous geometry is NOT updated here automatically
        // User must call updatePreviousGeometry() after checkpoint is accepted
    }
    
protected:
    /// @brief Update boundary displacement field from full mesh displacement vector
    /// @param[in] meshDisp full mesh displacement vector (dim * udofs, 1)
    void updateBoundaryDisplacement(const gsMatrix<T>& meshDisp)
    {
        if (m_boundaryDisplacement.nPatches() == 0)
            return;

        const index_t dim = m_assemblerPtr->getPatches().geoDim();
        const index_t udofs = getALEAssembler()->getUdofs();
        const gsDofMapper& uMap = getALEAssembler()->getMappers()[0];

        // Convert mesh displacement vector to boundary displacement field
        for (index_t p = 0; p < m_boundaryDisplacement.nPatches(); ++p)
        {
            gsMatrix<T>& coefs = m_boundaryDisplacement.patch(p).coefs();
            const index_t nCoefs = coefs.rows();

            for (index_t i = 0; i < nCoefs; ++i)
            {
                if (!uMap.is_free(i, p)) continue;
                const index_t gi = uMap.index(i, p);

                if (gi < udofs)
                {
                    coefs(i, 0) = meshDisp(gi, 0);         // x displacement
                    coefs(i, 1) = meshDisp(gi + udofs, 0); // y displacement
                }
            }
        }
    }

    /// @brief Setup gsALE boundary interface from FSI interface sides
    void setupGsALEInterface()
    {
        m_aleInterface = gsBoundaryInterface();
        for (const auto& ps : m_fsiInterfaceSides)
        {
            // Map each FSI interface side to itself for gsALE
            m_aleInterface.addInterfaceSide(ps.patch, ps.side(), ps.patch, ps.side());
        }
        gsInfo << "[gsINSSolverALE] Setup gsALE interface with " << m_fsiInterfaceSides.size()
               << " sides\n";
    }

    /// @brief Apply gsALE mesh deformation
    void applyGsALEDeformation()
    {
        if (m_originalMesh.nPatches() == 0 || m_boundaryDisplacement.nPatches() == 0)
        {
            gsWarn << "[gsINSSolverALE] gsALE not properly initialized. Call initializeGsALE() first.\n";
            return;
        }

        try
        {
            // Get current mesh from assembler (this is the deformed mesh from applyMeshDisplacement)
            gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());

            // Create gsALE object with current boundary displacement and interface
            gsALE<T> meshDeformer(patches, m_boundaryDisplacement, m_aleInterface, m_aleMethod);

            // Set gsALE options
            meshDeformer.options().setReal("PoissonsRatio",
                m_gsALEOptions.getReal("PoissonsRatio"));
            meshDeformer.options().setReal("LocalStiff",
                m_gsALEOptions.getReal("LocalStiff"));
            meshDeformer.options().setSwitch("Check",
                m_gsALEOptions.getSwitch("Check"));

            // Update mesh using gsALE
            index_t badPatch = meshDeformer.updateMesh();

            if (badPatch != -1)
            {
                gsWarn << "[gsINSSolverALE] gsALE mesh deformation failed! Bad patch: "
                       << badPatch << "\n";
                gsWarn << "Continuing with current mesh state.\n";
                return;
            }

            // Extract ALE displacement field from gsALE
            gsMultiPatch<T> aleDisplacement;
            meshDeformer.constructSolution(aleDisplacement);

            // Update the mesh in assembler's patches (already updated by gsALE::updateMesh)
            // The geometry control points in 'patches' are now modified

            gsInfo << "[gsINSSolverALE] Applied gsALE mesh deformation with method "
                   << m_aleMethod << "\n";
        }
        catch (const std::exception& e)
        {
            gsWarn << "[gsINSSolverALE] gsALE deformation failed: " << e.what() << "\n";
            gsWarn << "Continuing with current mesh.\n";
        }
    }

    /// @brief Apply mesh displacement to the actual geometry
    void applyMeshDisplacement()
    {
        if (m_originalMesh.nPatches() == 0)
            return;
            
        // IMPORTANT: Use m_time (current time) not m_time + m_timeStepSize
        // The mesh should be at the position corresponding to the current time step
        // being solved. The time will be incremented AFTER the solution is computed.
        gsMatrix<T> cumulativeDisp = m_meshUpdateFunc(m_time);
        
        // CRITICAL FIX: Check if displacement is essentially zero
        T dispNorm = cumulativeDisp.template lpNorm<gsEigen::Infinity>();
        if (dispNorm < 1e-12) {
            // For zero displacement, just update the ALE assembler without touching the mesh
            // This avoids disrupting the assembler's internal state
            getALEAssembler()->updateMesh(cumulativeDisp);
            return;
        }
        // Get current mesh from assembler
        gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());

        // Build a displacement field from the provided coefficients (use current cumulativeDisp)
        // Pad with zeros to full size (pressure part zero)
        const index_t udofs = getALEAssembler()->getUdofs();
        const index_t pdofs = getALEAssembler()->getPdofs();
        const index_t dim   = m_assemblerPtr->getPatches().geoDim();
        const index_t pshift = dim * udofs;
        gsMatrix<T> dispFull(pshift + pdofs, 1);
        dispFull.setZero();
        if (cumulativeDisp.size() >= pshift)
            dispFull.topRows(pshift) = cumulativeDisp.topRows(pshift);
        gsField<T> dispField = getALEAssembler()->constructSolution(dispFull, 0);

        // Reset to original mesh and apply displacement evaluated at geometry anchors
        // Optionally clamp displacement to zero on the outer box boundaries
        gsMatrix<T> bb;
        m_originalMesh.boundingBox(bb);
        const T xmin = bb(0,0), xmax = bb(0,1);
        const T ymin = bb(1,0), ymax = bb(1,1);

        for (size_t p = 0; p < patches.nPatches(); ++p)
        {
            // Reset to original mesh
            patches.patch(p).coefs() = m_originalMesh.patch(p).coefs();

            // Evaluate displacement at the (parametric) anchors of the geometry basis
            const index_t dim = patches.patch(p).targetDim();
            const gsBasis<T> & geoBasis = patches.patch(p).basis();
            gsMatrix<T> anchors = geoBasis.anchors(); // parametric anchors (paramDim x nCoefs)

            if (anchors.cols() == 0)
                continue;

            gsMatrix<T> U; // dim x nCoefs
            U = dispField.value(anchors, static_cast<index_t>(p));

            if (m_fixOuterBoundary)
            {
                // Evaluate physical coordinates of anchors on original mesh
                gsMatrix<T> XY;
                m_originalMesh.patch(static_cast<index_t>(p)).eval_into(anchors, XY); // dim x nCoefs

                for (index_t i = 0; i < XY.cols(); ++i)
                {
                    const T x = XY(0,i);
                    const T y = XY(1,i);
                    const bool onLeft   = std::abs(x - xmin) <= m_fixOuterTol;
                    const bool onRight  = std::abs(x - xmax) <= m_fixOuterTol;
                    const bool onBottom = std::abs(y - ymin) <= m_fixOuterTol;
                    const bool onTop    = std::abs(y - ymax) <= m_fixOuterTol;
                    if (onLeft || onRight || onBottom || onTop)
                    {
                        // Zero displacement on outer boundary anchors
                        for (index_t d = 0; d < dim; ++d)
                            U(d,i) = 0;
                    }
                }
            }

            // Add displacement to geometry control points
            for (index_t i = 0; i < patches.patch(p).coefsSize(); ++i)
            {
                for (index_t d = 0; d < dim; ++d)
                    patches.patch(p).coef(i, d) += U(d, i);
            }
        }

        // Update incremental displacement for ALE assembler (keeps m_meshVelCoefs consistent)
        getALEAssembler()->updateMesh(cumulativeDisp);
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
                // Use global mode: all patches optimized together with free interfaces
                // patchWise=true: each patch independently with HLBFGS, interfaces fixed
                // patchWise=false: global optimization with HLBFGS, all patches together, interfaces free
                gsBarrierPatch<2, T> opt(patches, true);  // false = global mode (interfaces free)
                opt.options() = m_meshOptOptions;
                opt.compute();
                patches = opt.result();

                gsInfo << "Applied global mesh optimization with HLBFGS (interfaces free)\n";
            }
            
            // Update mesh in assembler with displacement consistent to original mesh
            getALEAssembler()->updateMeshAfterOptimization(patches, m_originalMesh);
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
