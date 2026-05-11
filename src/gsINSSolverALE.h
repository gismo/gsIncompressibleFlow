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
#include <cmath>

namespace gismo
{

/// @brief Compute fluid traction at a parametric interface point.
///
/// The pressure field is assumed to be kinematic pressure. The returned
/// traction is physical stress times the reference normal. If
/// pullBackToReference is true, the Cauchy stress is pulled back with
/// J * sigma * F^{-T}.
template <class T = real_t>
gsVector<T> computeFluidTractionAtPoint(
    const gsField<T>& velField,
    const gsField<T>& presField,
    const gsMultiPatch<T>& refPatches,
    const gsMultiPatch<T>& curPatches,
    const gsVector<T>& uv,
    index_t patchId,
    boxSide side,
    T density,
    T viscosity,
    bool includeViscous,
    bool pullBackToReference,
    const gsVector<T>* normalTargetPoint = nullptr)
{
    const index_t dim = refPatches.geoDim();
    gsMatrix<T> uvM(uv.rows(), 1);
    uvM.col(0) = uv;
    const gsMatrix<T> I = gsMatrix<T>::Identity(dim, dim);

    gsMapData<T> mdRef(NEED_GRAD_TRANSFORM | NEED_MEASURE);
    mdRef.patchId = patchId;
    mdRef.points = uvM;
    mdRef.side = side;
    refPatches.patch(patchId).computeMap(mdRef);

    gsMapData<T> mdCur(NEED_GRAD_TRANSFORM | NEED_MEASURE);
    mdCur.patchId = patchId;
    mdCur.points = uvM;
    mdCur.side = side;
    curPatches.patch(patchId).computeMap(mdCur);

    gsVector<T> nref = mdRef.outNormal(0);
    const T nlen = nref.norm();
    if (nlen > T(0))
        nref /= nlen;

    if (normalTargetPoint && normalTargetPoint->rows() == dim)
    {
        gsMatrix<T> xRef;
        refPatches.patch(patchId).eval_into(uvM, xRef);
        if (nref.dot((*normalTargetPoint) - xRef.col(0)) < T(0))
            nref = -nref;
    }

    const gsMatrix<T> pk = presField.value(uvM, patchId);
    const T pPhys = density * pk(0, 0);

    const gsMatrix<T> F = mdCur.jacobian(0) * mdRef.jacobian(0).cramerInverse();
    const gsMatrix<T> invF = F.cramerInverse();
    const gsMatrix<T> invFT = invF.transpose();
    const T J = F.determinant();

    gsMapData<T> mdVel(NEED_DERIV);
    mdVel.patchId = patchId;
    mdVel.points = uvM;
    velField.igaFunction(patchId).computeMap(mdVel);
    const gsMatrix<T> gradXU = mdVel.jacobian(0) * mdRef.jacobian(0).cramerInverse();

    const gsMatrix<T> sigmaP = pPhys * I;
    gsMatrix<T> sigmaV(dim, dim);
    sigmaV.setZero();
    if (includeViscous && density * viscosity > T(0))
        sigmaV -= (density * viscosity) * (gradXU * invF + invFT * gradXU.transpose());

    const gsMatrix<T> sigma = sigmaP + sigmaV;
    gsVector<T> traction(dim);
    if (pullBackToReference)
        traction = J * sigma * invFT * nref;
    else
        traction = sigma * nref;

    for (index_t d = 0; d < dim; ++d)
        if (!std::isfinite(traction(d)))
            traction(d) = T(0);

    return traction;
}

/// @brief Boundary mesh velocity function for FSI interfaces.
/// Input points are physical coordinates.
///
/// Dimension (2D or 3D) is detected automatically from the current-geometry
/// pointer at construction time. A callable with `cur=nullptr` is still
/// supported as long as the dim is passed explicitly; otherwise 2 is assumed
/// until a geometry is set via setCurrent().
template <class T = real_t>
class gsFSIMeshVelocityFunction : public gsFunction<T>
{
protected:
    const gsGeometry<T>* m_cur;   ///< Current geometry
    const gsGeometry<T>* m_prev;  ///< Previous geometry
    T m_dt;                        ///< Time step
    short_t m_dim;                 ///< Spatial dimension (2 or 3)

public:
    /// @brief Constructor
    /// @param[in] cur   current geometry (used for dim auto-detection if dim<0)
    /// @param[in] prev  previous geometry
    /// @param[in] dt    time step
    /// @param[in] dim   spatial dimension; if <=0, taken from cur->geoDim()
    gsFSIMeshVelocityFunction(const gsGeometry<T>* cur, const gsGeometry<T>* prev,
                              T dt, short_t dim = -1)
    : m_cur(cur), m_prev(prev), m_dt(dt),
      m_dim( dim > 0 ? dim : (cur ? static_cast<short_t>(cur->geoDim()) : short_t(2)) )
    {}

    /// @brief Set current geometry (also refreshes auto-detected dim)
    void setCurrent(const gsGeometry<T>* cur)
    {
        m_cur = cur;
        if (cur) m_dim = static_cast<short_t>(cur->geoDim());
    }

    /// @brief Set previous geometry
    void setPrevious(const gsGeometry<T>* prev) { m_prev = prev; }

    /// @brief Set time step
    void setDt(T dt) { m_dt = dt; }

    /// @brief Domain dimension (physical space)
    short_t domainDim() const { return m_dim; }

    /// @brief Target dimension (velocity vector)
    short_t targetDim() const { return m_dim; }

    /// @brief Evaluate mesh velocity at PHYSICAL points
    /// @param x Physical coordinates (dim x N)
    /// @param result Mesh velocity at these points (dim x N)
    void eval_into(const gsMatrix<T>& x, gsMatrix<T>& result) const
    {
        const index_t N = x.cols();
        result.resize(m_dim, N);

        for (index_t i = 0; i < N; ++i)
        {
            gsMatrix<T> phys(m_dim, 1);
            phys.col(0) = x.col(i);

            gsMatrix<T> uv;
            m_cur->invertPoints(phys, uv);

            gsMatrix<T> xPrev = m_prev->eval(uv);

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
    using Base::m_iterationNumber;
    using Base::m_time;
    using Base::m_timeStepSize;
    
    bool m_isALEActive;

    std::function<gsMatrix<T>(T)> m_meshUpdateFunc;

    bool m_useMeshOptimization;

    bool m_useDynamicBoundaryMapping;

    // gsBarrierPatch options
    gsOptionList m_meshOptOptions;

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

    // Boundary mesh velocity functions used for u=w on FSI interfaces.
    std::vector<std::shared_ptr<gsFunction<T>>> m_fsiBoundaryVelFuncs;

    // Boundary condition storage owned by the caller.
    gsBoundaryConditions<T>* m_bcInfoPtr;

    // Previous geometry for computing mesh velocity
    gsMultiPatch<T> m_prevGeometry;

public:
    struct ALEState
    {
        T time = T(0);
        T timeStepSize = T(0);
        T avgPicardIter = T(0);
        unsigned iterationNumber = 0;
        gsMultiPatch<T> prevGeometry;
        typename gsINSAssemblerUnsteadyALE<T, MatOrder>::ALEState assemblerState;
    };

    /// @brief Constructor
    gsINSSolverUnsteadyALE(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    : Base(paramsPtr), m_isALEActive(false), m_useMeshOptimization(false),
      m_useDynamicBoundaryMapping(false), m_rotationPeriod(10.0), m_bcInfoPtr(nullptr),
      m_fixOuterBoundary(false), m_fixOuterTol(static_cast<T>(1e-8))
    {
        delete m_assemblerPtr;
        m_assemblerPtr = new gsINSAssemblerUnsteadyALE<T, MatOrder>(m_paramsPtr);
        Base::initMembers();

        // Mesh optimization options
        m_meshOptOptions.addInt("Verbose", "Verbosity level", 0);
        m_meshOptOptions.addInt("ParamMethod", "Parametrization method", 1);
        m_meshOptOptions.addInt("AAPreconditionType", "AA precondition type", 0);

        // Detect spatial dimension from the PDE's geometry so that both 2D and
        // 3D problems get a sensible rotation center. Default: domain bbox center.
        const short_t dim = static_cast<short_t>(
            m_assemblerPtr->getPatches().geoDim());
        m_rotationCenter.setZero(dim);
        if (dim >= 2)
        {
            m_rotationCenter(0) = static_cast<T>(0.5);
            m_rotationCenter(1) = static_cast<T>(0.5);
        }
        // 3D: z stays at 0 (user should override via setRotationParameters)
    }

    ALEState saveALEState() const
    {
        ALEState state;
        state.time = m_time;
        state.timeStepSize = m_timeStepSize;
        state.avgPicardIter = this->m_avgPicardIter;
        state.iterationNumber = m_iterationNumber;
        state.prevGeometry = m_prevGeometry;
        if (getALEAssembler())
            state.assemblerState = getALEAssembler()->saveALEState();
        return state;
    }

    void restoreALEState(const ALEState& state)
    {
        m_time = state.time;
        m_timeStepSize = state.timeStepSize;
        this->m_avgPicardIter = state.avgPicardIter;
        m_iterationNumber = state.iterationNumber;
        m_prevGeometry = state.prevGeometry;
        m_paramsPtr->options().setReal("timeStep", m_timeStepSize);
        if (getALEAssembler())
            getALEAssembler()->restoreALEState(state.assemblerState);
    }

    /// @brief Get spatial dimension of the ALE problem (2 or 3)
    short_t getDim() const
    {
        return static_cast<short_t>(m_assemblerPtr->getPatches().geoDim());
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
    }

    /// @brief Set no-slip boundary conditions on FSI interfaces.
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

        m_bcInfoPtr = &bcInfo;

        const gsMultiPatch<T>& patches = m_assemblerPtr->getPatches();
        m_prevGeometry = patches;

        m_fsiBoundaryVelFuncs.clear();
        m_fsiBoundaryVelFuncs.reserve(m_fsiInterfaceSides.size());

        for (const auto& ps : m_fsiInterfaceSides)
        {
            auto meshVelFunc = std::make_shared<gsFSIMeshVelocityFunction<T>>(
                &patches.patch(ps.patch), &m_prevGeometry.patch(ps.patch), dt);

            m_fsiBoundaryVelFuncs.push_back(meshVelFunc);

            bcInfo.addCondition(ps.patch, ps.side(), condition_type::dirichlet,
                               meshVelFunc.get(), 0);
        }

        gsInfo << "[gsINSSolverALE] Applied Dirichlet BC u=w on "
               << m_fsiInterfaceSides.size() << " FSI interface sides\n";
    }

    /// @brief Update FSI boundary mesh velocity functions
    /// @param[in] dt time step size
    void updateFSIBoundaryVelocity(T dt)
    {
        if (m_fsiBoundaryVelFuncs.empty() || m_fsiInterfaceSides.empty())
            return;

        const gsMultiPatch<T>& curGeo = m_assemblerPtr->getPatches();

        for (size_t i = 0; i < m_fsiBoundaryVelFuncs.size(); ++i)
        {
            const index_t p = m_fsiInterfaceSides[i].patch;

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
    void updatePreviousGeometry()
    {
        if (m_prevGeometry.nPatches() > 0)
        {
            m_prevGeometry = m_assemblerPtr->getPatches();
        }
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

        this->m_timeStepSize = this->m_paramsPtr->options().getReal("timeStep");

        updateFSIBoundaryVelocity(this->m_timeStepSize);

        if (m_isALEActive && m_meshUpdateFunc)
        {
            applyMeshDisplacement();
            updateFSIBoundaryVelocity(this->m_timeStepSize);

            if (m_useMeshOptimization)
                optimizeMesh();
            this->m_paramsPtr->logger() << "[ALE] Mesh update done.\n";
        }

        this->m_paramsPtr->logger() << "[ALE] Updating assembler...\n";
        this->updateAssembler();

        if (!this->m_iterationNumber)
            this->initIteration();

        gsMatrix<T> tmpSolution = this->m_solution;

        this->applySolver(tmpSolution);
        this->writeSolChangeRelNorm(this->m_solution, tmpSolution, std::string("[ALE]"));

        index_t picardIter = 0;
        T relNorm = this->solutionChangeRelNorm(this->m_solution, tmpSolution);

        this->m_paramsPtr->logger() << "        [u, p] Picard's iterations...\n";

        while((relNorm > this->m_innerTol) && (picardIter < this->m_innerIter))
        {
            this->m_paramsPtr->logger() << "         ";

            gsMatrix<T> oldSol = tmpSolution;

            this->updateAssembler(tmpSolution, false);
            this->applySolver(tmpSolution);
            this->writeSolChangeRelNorm(oldSol, tmpSolution, std::string("[ALE]"));

            relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
            picardIter++;
        }

        this->m_solution = tmpSolution;

        getALEAssembler()->updateCurrentSolField(tmpSolution, true);

        this->m_time += this->m_timeStepSize;
        this->m_avgPicardIter += picardIter;
        this->m_iterationNumber++;

        // Previous geometry is advanced after the coupling checkpoint is accepted.
    }
    
protected:
    /// @brief Apply mesh displacement to the actual geometry
    void applyMeshDisplacement()
    {
        if (m_originalMesh.nPatches() == 0)
            return;

        gsMatrix<T> cumulativeDisp = m_meshUpdateFunc(m_time);

        T dispNorm = cumulativeDisp.template lpNorm<gsEigen::Infinity>();
        if (dispNorm < 1e-12) {
            getALEAssembler()->updateMesh(cumulativeDisp);
            return;
        }

        gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());

        const index_t udofs = getALEAssembler()->getUdofs();
        const index_t pdofs = getALEAssembler()->getPdofs();
        const index_t dim   = m_assemblerPtr->getPatches().geoDim();
        const index_t pshift = dim * udofs;
        gsMatrix<T> dispFull(pshift + pdofs, 1);
        dispFull.setZero();
        if (cumulativeDisp.size() >= pshift)
            dispFull.topRows(pshift) = cumulativeDisp.topRows(pshift);
        gsField<T> dispField = getALEAssembler()->constructSolution(dispFull, 0);

        gsMatrix<T> bb;
        m_originalMesh.boundingBox(bb);
        const index_t meshDim = static_cast<index_t>(m_originalMesh.geoDim());

        for (size_t p = 0; p < patches.nPatches(); ++p)
        {
            patches.patch(p).coefs() = m_originalMesh.patch(p).coefs();

            const index_t dim = patches.patch(p).targetDim();
            const gsBasis<T> & geoBasis = patches.patch(p).basis();
            gsMatrix<T> anchors = geoBasis.anchors();

            if (anchors.cols() == 0)
                continue;

            gsMatrix<T> U; // dim x nCoefs
            U = dispField.value(anchors, static_cast<index_t>(p));

            if (m_fixOuterBoundary)
            {
                gsMatrix<T> XYZ;
                m_originalMesh.patch(static_cast<index_t>(p)).eval_into(anchors, XYZ);

                for (index_t i = 0; i < XYZ.cols(); ++i)
                {
                    bool onOuter = false;
                    for (index_t d = 0; d < meshDim && !onOuter; ++d)
                    {
                        const T c    = XYZ(d, i);
                        const T cMin = bb(d, 0);
                        const T cMax = bb(d, 1);
                        if (std::abs(c - cMin) <= m_fixOuterTol ||
                            std::abs(c - cMax) <= m_fixOuterTol)
                            onOuter = true;
                    }
                    if (onOuter)
                    {
                        // Zero displacement on outer boundary anchors (all components)
                        for (index_t d = 0; d < dim; ++d)
                            U(d, i) = 0;
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
    
    /// @brief Optimize mesh using gsBarrierPatch (dispatches 2D / 3D at runtime)
    void optimizeMesh()
    {
        try {
            // Get current mesh from assembler
            gsMultiPatch<T>& patches = const_cast<gsMultiPatch<T>&>(m_assemblerPtr->getPatches());

            const short_t dim = static_cast<short_t>(patches.geoDim());
            GISMO_ASSERT(dim == 2 || dim == 3,
                "gsINSSolverUnsteadyALE::optimizeMesh only supports 2D or 3D");

            // Use dynamic version if enabled and we have rotation
            if (m_useDynamicBoundaryMapping && m_meshUpdateFunc)
            {
                const T rotationAngle = getCurrentRotationAngle();

                if (dim == 2)
                {
                    gsBarrierPatchDynamic<2, T> opt;
                    opt.setDynamicBoundaryMapping(true);
                    opt.setRotationAngle(rotationAngle);
                    opt.setRotationCenter(m_rotationCenter);
                    patches = opt.compute(patches, m_meshOptOptions);
                }
                else // dim == 3
                {
                    gsBarrierPatchDynamic<3, T> opt;
                    opt.setDynamicBoundaryMapping(true);
                    opt.setRotationAngle(rotationAngle);
                    opt.setRotationCenter(m_rotationCenter);
                    patches = opt.compute(patches, m_meshOptOptions);
                }

                gsInfo << "Applied mesh optimization with dynamic boundary mapping (dim="
                       << dim << ")\n";
            }
            else
            {
                // Global mode: all patches optimized together with free interfaces.
                // patchWise=false: global optimization with HLBFGS, interfaces free
                if (dim == 2)
                {
                    gsBarrierPatch<2, T> opt(patches, false);
                    opt.options() = m_meshOptOptions;
                    opt.compute();
                    patches = opt.result();
                }
                else // dim == 3
                {
                    gsBarrierPatch<3, T> opt(patches, false);
                    opt.options() = m_meshOptOptions;
                    opt.compute();
                    patches = opt.result();
                }

                gsInfo << "Applied global mesh optimization with HLBFGS (interfaces free, dim="
                       << dim << ")\n";
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

} // namespace gismo
