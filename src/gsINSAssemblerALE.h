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
    using Base::m_tarDim;
    using Base::m_visitorUUnonlin;
    using Base::m_currentVelField;
    using Base::m_blockUUnonlin_comp;
    using Base::m_blockUUnonlin_whole;
    using Base::m_rhsUnonlin;
    using Base::m_nnzPerOuterU;
    using Base::m_udofs;
    using Base::m_pdofs;
    using Base::m_pshift;
    using Base::m_dofs;
    using Base::m_solution;

    // Mesh velocity coefficients
    gsMatrix<T> m_meshVelCoefs;
    gsMatrix<T> m_meshDispCoefs;
    gsMatrix<T> m_meshDispOld;

    gsSparseMatrix<T, MatOrder> m_blockTimeDiscrOld;
    bool m_hasOldMassMatrix;

    // ALE visitor
    gsINSVisitorUUnonlinALE<T, MatOrder>* m_visitorUUnonlinALE;

    // Flag to indicate if ALE is active
    bool m_isALEActive;

    // Temporary mesh velocity field (created on demand)
    mutable gsField<T>* m_tempMeshVelField;

    bool m_hasSUPG;
    gsINSVisitorPU_SUPG_ALE<T, MatOrder>* m_visitorSUPG_PU_ALE;
    gsINSVisitorUU_SUPGTimeDiscrALE<T, MatOrder>* m_visitorSUPG_UUtd_ALE;
    gsSparseMatrix<T, MatOrder> m_blockSUPG_UP;
    gsSparseMatrix<T, MatOrder> m_blockSUPG_UUtd_comp;
    gsMatrix<T> m_rhsSUPG_UP, m_rhsSUPG_UUtd;

    // Old-time velocity DOFs for SUPG time discretization RHS
    gsMatrix<T> m_oldTimeSolU;

    gsSparseMatrix<T, MatOrder> m_blockGradDiv;
    bool m_hasGradDiv;
    gsMapData<T> m_gradDivMapData;

public:
    /// @brief Constructor
    gsINSAssemblerUnsteadyALE(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    : Base(paramsPtr),
      m_isALEActive(false),
      m_visitorUUnonlinALE(nullptr),
      m_tempMeshVelField(nullptr),
      m_hasOldMassMatrix(false),
      m_hasSUPG(false),
      m_visitorSUPG_PU_ALE(nullptr),
      m_visitorSUPG_UUtd_ALE(nullptr),
      m_hasGradDiv(false)
    { }

    /// @brief Destructor
    virtual ~gsINSAssemblerUnsteadyALE()
    {
        delete m_visitorUUnonlinALE;
        delete m_tempMeshVelField;
        delete m_visitorSUPG_PU_ALE;
        delete m_visitorSUPG_UUtd_ALE;
    }

    struct ALEState
    {
        gsMatrix<T> meshVelCoefs;
        gsMatrix<T> meshDispCoefs;
        gsMatrix<T> meshDispOld;
        gsMatrix<T> oldTimeSolU;
    };

    ALEState saveALEState() const
    {
        ALEState state;
        state.meshVelCoefs = m_meshVelCoefs;
        state.meshDispCoefs = m_meshDispCoefs;
        state.meshDispOld = m_meshDispOld;
        state.oldTimeSolU = m_oldTimeSolU;
        return state;
    }

    void restoreALEState(const ALEState& state)
    {
        m_meshVelCoefs = state.meshVelCoefs;
        m_meshDispCoefs = state.meshDispCoefs;
        m_meshDispOld = state.meshDispOld;
        m_oldTimeSolU = state.oldTimeSolU;
    }

    /// @brief Initialize (override to add ALE initialization)
    virtual void initialize() override
    {
        Base::initialize();

        const index_t udofs = this->getMapper(0).freeSize();
        const index_t pdofs = this->getMapper(1).freeSize();
        const index_t pshift = m_tarDim * udofs;
        const index_t fullDofs = pshift + pdofs;
        m_meshVelCoefs.setZero(fullDofs, 1);
        m_meshDispCoefs.setZero(fullDofs, 1);
        m_meshDispOld.setZero(pshift, 1);

        m_visitorUUnonlinALE = new gsINSVisitorUUnonlinALE<T, MatOrder>(
            m_paramsPtr, m_tarDim, nullptr);
        m_visitorUUnonlinALE->initialize();
        if (this->m_currentVelField.nPieces() > 0)
            m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);

        m_hasSUPG = m_paramsPtr->options().getSwitch("SUPG_NS");

        if (m_hasSUPG)
        {
            m_visitorSUPG_PU_ALE = new gsINSVisitorPU_SUPG_ALE<T, MatOrder>(
                m_paramsPtr, m_tarDim, nullptr);
            m_visitorSUPG_PU_ALE->initialize();

            m_visitorSUPG_UUtd_ALE = new gsINSVisitorUU_SUPGTimeDiscrALE<T, MatOrder>(
                m_paramsPtr, m_tarDim, nullptr);
            m_visitorSUPG_UUtd_ALE->initialize();

            m_blockSUPG_UP.resize(m_pshift, m_pdofs);
            m_blockSUPG_UUtd_comp.resize(m_udofs, m_udofs);
            m_rhsSUPG_UP.setZero(m_dofs, 1);
            m_rhsSUPG_UUtd.setZero(m_pshift, 1);
        }

        if (m_hasSUPG)
            m_oldTimeSolU.setZero(m_pshift, 1);

        m_hasGradDiv = (m_paramsPtr->options().getReal("graddiv.gamma") > T(0));
    }

    /// @brief Assemble grad-div stabilization: gamma * (div u)(div v) as a whole UU block.
    void assembleGradDiv()
    {
        const T gamma = m_paramsPtr->options().getReal("graddiv.gamma");
        if (gamma <= T(0)) return;

        m_blockGradDiv.resize(m_pshift, m_pshift);
        m_blockGradDiv.reserve(gsVector<index_t>::Constant(m_pshift, this->m_nnzPerOuterU * m_tarDim));

        const gsMultiBasis<T>& velBasis = this->getBases().at(0);
        const gsMultiPatch<T>& geo      = this->getPatches();
        const gsDofMapper& mapper        = this->getMapper(0);

        gsMatrix<T> quNodes, bGrads, physGrads;
        gsVector<T> quWeights;
        gsMatrix<index_t> actives;
        index_t prevPatch = -1;
        gsGaussRule<T> rule;

        typename gsBasis<T>::domainIter domIt    = velBasis.domain()->beginAll();
        typename gsBasis<T>::domainIter domItEnd = velBasis.domain()->endAll();

        for (; domIt < domItEnd; ++domIt)
        {
            const index_t p = domIt.patch();

            if (p != prevPatch)
            {
                const gsBasis<T>& bp = velBasis.piece(p);
                gsVector<index_t> numNodes(m_tarDim);
                for (short_t d = 0; d < m_tarDim; ++d)
                    numNodes[d] = bp.degree(d) + 1;
                rule = gsGaussRule<T>(numNodes);
                prevPatch = p;
            }

            const gsBasis<T>& bp = velBasis.piece(p);
            rule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);

            m_gradDivMapData.flags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            m_gradDivMapData.points = quNodes;
            geo.patch(p).computeMap(m_gradDivMapData);

            bp.deriv_into(quNodes, bGrads);
            bp.active_into(quNodes.col(0), actives);
            const index_t nAct = actives.rows();

            for (index_t q = 0; q < quWeights.rows(); ++q)
            {
                const T w = gamma * quWeights(q) * m_gradDivMapData.measure(q);
                transformGradients(m_gradDivMapData, q, bGrads, physGrads);
                // physGrads: (nAct × tarDim)

                for (index_t a = 0; a < nAct; ++a)
                {
                    if (!mapper.is_free(actives(a, 0), p)) continue;
                    const index_t ag = mapper.index(actives(a, 0), p);

                    for (index_t b = 0; b < nAct; ++b)
                    {
                        if (!mapper.is_free(actives(b, 0), p)) continue;
                        const index_t bg = mapper.index(actives(b, 0), p);

                        for (short_t i = 0; i < m_tarDim; ++i)
                            for (short_t j = 0; j < m_tarDim; ++j)
                                m_blockGradDiv.coeffRef(ag + i * m_udofs,
                                                        bg + j * m_udofs)
                                    += w * physGrads(a, i) * physGrads(b, j);
                    }
                }
            }
        }
        m_blockGradDiv.makeCompressed();
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

        // SUPG blocks are managed by the ALE visitors.
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
        const index_t udofs = this->getMapper(0).freeSize();
        const index_t pdofs = this->getMapper(1).freeSize();
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
        const index_t udofs = this->getMapper(0).freeSize();
        const index_t pdofs = this->getMapper(1).freeSize();
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
        T velNorm = m_meshVelCoefs.template lpNorm<gsEigen::Infinity>();

        if (velNorm < 1e-12) {
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
        // Compute displacement coefficients from optimized mesh relative to original mesh
        const index_t udofs = this->getMapper(0).freeSize();
        gsMatrix<T> fullDisp(m_tarDim * udofs, 1);
        fullDisp.setZero();

        const gsDofMapper& mapper = this->getMapper(0);
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
        GISMO_UNUSED(optimizedPatches);
        // No direct patch replacement here; solver already updates patches.
    }

    /// @brief Update current solution field using the unsteady base history.
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol) override
    {
        Base::updateCurrentSolField(solVector, updateSol);
    }

    /// @brief Override update: refresh Dirichlet data BEFORE assembly so u=w holds on moving walls
    virtual void update(const gsMatrix<T>& sol, bool updateSol = true) override
    {
        // If ALE is active, prepare mesh-velocity field for boundary updates and ALE terms
        if (m_isALEActive && m_visitorUUnonlinALE)
        {
            if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
                this->computeDirichletDofs(0, 0, this->m_ddof[0]);

            delete m_tempMeshVelField;
            m_tempMeshVelField = new gsField<T>(this->constructSolution(m_meshVelCoefs, 0));
            m_visitorUUnonlinALE->setMeshVelocityField(m_tempMeshVelField);
        }

        Base::updateCurrentSolField(sol, updateSol);

        this->assembleLinearPart();

        // Assemble grad-div stabilization (linear term, but geometry changes per step)
        if (m_hasGradDiv)
            this->assembleGradDiv();

        this->fillBaseSystem();

        this->assembleNonlinearPart();

        this->fillSystem();

        if (m_isALEActive && m_visitorUUnonlinALE && this->m_currentVelField.nPieces() > 0)
            m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);

        // Store old-time velocity DOFs for next SUPG time discretization RHS
        if (m_hasSUPG)
            m_oldTimeSolU = m_solution.topRows(this->m_pshift);
    }


    /// @brief Advance ALE time-step data.
    void advanceTimeStep() { /* No-op for now */ }

    /// @brief Get stored old mass matrix.
    const gsSparseMatrix<T, MatOrder>& getBlockTimeDiscrOld() const { return m_blockTimeDiscrOld; }

    /// @brief Set stored old mass matrix.
    void setBlockTimeDiscrOld(const gsSparseMatrix<T, MatOrder>& mat)
    {
        m_blockTimeDiscrOld = mat;
        m_hasOldMassMatrix = (mat.nonZeros() > 0);
    }

    /// @brief Check if old mass matrix is available
    bool hasOldMassMatrix() const { return m_hasOldMassMatrix; }

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
            const index_t udofs = this->getMapper(0).freeSize();
            const index_t pshift = m_tarDim * udofs;
            T velInf = (m_meshVelCoefs.size() >= pshift) ? m_meshVelCoefs.topRows(pshift).template lpNorm<gsEigen::Infinity>() : T(0);
            if (velInf < static_cast<T>(1e-14))
                useALE = false;
        }

        if (useALE)
        {
            m_visitorUUnonlinALE->setParamsPtr(m_paramsPtr);

            if (this->m_currentVelField.nPieces() > 0)
                m_visitorUUnonlinALE->setCurrentSolution(this->m_currentVelField);
            else
            {
                gsWarn << "[ALE] No current velocity field set, falling back to standard formulation\n";
                Base::assembleNonlinearPart();
                return;
            }

            if (m_tempMeshVelField)
                m_visitorUUnonlinALE->setMeshVelocityField(m_tempMeshVelField);

            this->m_blockUUnonlin_comp.reserve(gsVector<index_t>::Constant(this->m_blockUUnonlin_comp.outerSize(), this->m_nnzPerOuterU));
            this->assembleBlock(*m_visitorUUnonlinALE, 0, this->m_blockUUnonlin_comp, this->m_rhsUnonlin);

            // SUPG convection + viscous terms are assembled in m_visitorUUnonlinALE above.
            // Now assemble SUPG PressureGrad and TimeDiscr blocks using ALE visitors.
            if (m_hasSUPG && m_visitorSUPG_PU_ALE && m_visitorSUPG_UUtd_ALE)
            {
                if (m_tempMeshVelField)
                {
                    m_visitorSUPG_PU_ALE->setMeshVelocityField(m_tempMeshVelField);
                    m_visitorSUPG_UUtd_ALE->setMeshVelocityField(m_tempMeshVelField);
                }
                if (this->m_currentVelField.nPieces() > 0)
                {
                    m_visitorSUPG_PU_ALE->setCurrentSolution(this->m_currentVelField);
                    m_visitorSUPG_UUtd_ALE->setCurrentSolution(this->m_currentVelField);
                }

                // SUPG PressureGrad -> m_blockSUPG_UP (ALE member)
                m_blockSUPG_UP.resize(this->m_pshift, this->m_pdofs);
                m_rhsSUPG_UP.setZero(this->m_dofs, 1);
                m_blockSUPG_UP.reserve(gsVector<index_t>::Constant(
                    m_blockSUPG_UP.outerSize(), this->m_nnzPerOuterU));
                this->assembleBlock(*m_visitorSUPG_PU_ALE, 0,
                                    m_blockSUPG_UP, m_rhsSUPG_UP);

                // SUPG TimeDiscr -> m_blockSUPG_UUtd_comp (ALE member)
                m_blockSUPG_UUtd_comp.resize(this->m_udofs, this->m_udofs);
                m_rhsSUPG_UUtd.setZero(this->m_pshift, 1);
                m_blockSUPG_UUtd_comp.reserve(gsVector<index_t>::Constant(
                    m_blockSUPG_UUtd_comp.outerSize(), this->m_nnzPerOuterU));
                this->assembleBlock(*m_visitorSUPG_UUtd_ALE, 0,
                                    m_blockSUPG_UUtd_comp, m_rhsSUPG_UUtd);
            }
        }
        else
        {
            Base::assembleNonlinearPart();
        }
    }

    /// @brief Override fillSystem to add ALE SUPG blocks and grad-div.
    virtual void fillSystem() override
    {
        Base::fillSystem();

        // Add ALE SUPG blocks
        if (m_hasSUPG)
        {
            if (m_blockSUPG_UP.nonZeros() > 0)
                this->fillGlobalMat_UP(this->m_matrix, m_blockSUPG_UP);

            if (m_blockSUPG_UUtd_comp.nonZeros() > 0)
                this->fillGlobalMat_UU(this->m_matrix, m_blockSUPG_UUtd_comp);

            if (!this->m_matrix.isCompressed())
                this->m_matrix.makeCompressed();

            // RHS: SUPG time discretization from eliminated DOFs
            this->m_rhs.topRows(this->m_pshift) += m_rhsSUPG_UUtd;

            // RHS: SUPG time discretization old-time contribution
            if (m_blockSUPG_UUtd_comp.nonZeros() > 0)
            {
                gsSparseMatrix<T, MatOrder> supgTdExpanded(this->m_pshift, this->m_pshift);
                supgTdExpanded.reserve(gsVector<index_t>::Constant(
                    supgTdExpanded.outerSize(), this->m_nnzPerOuterU));
                this->fillGlobalMat_UU(supgTdExpanded, m_blockSUPG_UUtd_comp);
                supgTdExpanded.makeCompressed();
                this->m_rhs.topRows(this->m_pshift) += supgTdExpanded * m_oldTimeSolU;
            }

            // RHS: SUPG pressure gradient from Dirichlet BCs
            this->m_rhs += m_rhsSUPG_UP;
        }

        // Add grad-div stabilization (whole UU block)
        if (m_hasGradDiv && m_blockGradDiv.nonZeros() > 0)
        {
            this->fillGlobalMat_UU(this->m_matrix, m_blockGradDiv);
            if (!this->m_matrix.isCompressed())
                this->m_matrix.makeCompressed();
        }
    }
};

} // namespace gismo
