/** @file gsINSMeshDeformer.h

    @brief Mesh deformation solver for incompressible flow using gsElasticity::gsALE

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li (2025)
*/

#pragma once

#include <gsElasticity/src/gsALE.h>
#include <gsElasticity/src/gsBaseUtils.h>
#include <gsCore/gsMultiPatch.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsMemory.h>

namespace gismo
{

/// @brief Enumeration of mesh deformation methods
/// These correspond to gsElasticity::ale_method
namespace mesh_deform_method
{
    enum method
    {
        HE = 0,          ///< Harmonic Extension (Laplace equation)
        IHE = 1,         ///< Incremental Harmonic Extension
        LE = 2,          ///< Linear Elasticity
        ILE = 3,         ///< Incremental Linear Elasticity
        TINE = 4,        ///< Tangential-Incremental-Nonlinear-Elasticity (Neo-Hookean)
        TINE_StVK = 5,   ///< TINE with St. Venant-Kirchhoff material
        BHE = 6,         ///< Biharmonic Extension
        IBHE = 7         ///< Incremental Biharmonic Extension
    };
}

/// @brief Mesh deformer for incompressible flow problems
///
/// This class wraps gsElasticity::gsALE to provide mesh deformation
/// for fluid domains in FSI or moving boundary problems. It solves
/// a PDE (Laplace, elasticity, or biharmonic) to propagate boundary
/// displacement into the domain interior.
///
/// @tparam T coefficient type
template <class T>
class gsINSMeshDeformer
{
public:
    /// @brief Constructor
    /// @param[in,out] fluidGeometry reference to fluid domain geometry (will be modified)
    /// @param[in] method mesh deformation method
    gsINSMeshDeformer(gsMultiPatch<T>& fluidGeometry,
                      mesh_deform_method::method method = mesh_deform_method::IHE)
    : m_geometry(fluidGeometry),
      m_method(method),
      m_isInitialized(false)
    {
        // Store original geometry
        m_originalGeometry = fluidGeometry;
    }

    /// @brief Destructor
    ~gsINSMeshDeformer() { }

    /// @brief Set up interface boundaries for mesh deformation
    /// @param[in] movingBoundaries list of patch sides that move with prescribed displacement
    ///
    /// Example: For FSI, these are the fluid-structure interface boundaries
    void setMovingBoundaries(const std::vector<patchSide>& movingBoundaries)
    {
        m_movingBoundaries = movingBoundaries;

        // Build a persistent interface mapping (displacement sides A -> geometry sides B)
        // Use one-to-one mapping on the same patch/side by default.
        m_interface = gsBoundaryInterface();
        for (const auto &ps : m_movingBoundaries)
            m_interface.addInterfaceSide(ps.patch, ps.side(), ps.patch, ps.side());
    }

    /// @brief Initialize the mesh deformer
    ///
    /// This creates the internal gsALE solver. Call this after setting
    /// moving boundaries and before the first mesh update.
    void initialize()
    {
        if (m_isInitialized)
            return;

        // Create a dummy displacement field (same structure as geometry)
        m_boundaryDisplacement.clear();
        for (size_t p = 0; p < m_geometry.nPatches(); ++p)
        {
            m_boundaryDisplacement.addPatch(m_geometry.patch(p).clone());
            m_boundaryDisplacement.patch(p).coefs().setZero();
        }

        // Ensure we have an interface; if not provided, default to mapping all boundaries
        if (m_interface.sidesA.empty())
        {
            for (index_t p = 0; p < static_cast<index_t>(m_geometry.nPatches()); ++p)
            {
                for (index_t s = 1; s <= 2 * m_geometry.parDim(); ++s)
                {
                    m_interface.addInterfaceSide(p, boundary::side(s), p, boundary::side(s));
                }
            }
        }

        // Map method enum to gsElasticity::ale_method
        ale_method::method aleMethod = static_cast<ale_method::method>(m_method);

        // Create gsALE solver with a persistent interface mapping
        m_aleSolver.reset(new gsALE<T>(m_geometry, m_boundaryDisplacement, m_interface, aleMethod));

        // Set default options
        m_aleSolver->options() = defaultOptions();

        m_isInitialized = true;

        gsInfo << "gsINSMeshDeformer initialized with method: " << methodName() << "\n";
    }

    /// @brief Update mesh by prescribing displacement on moving boundaries
    /// @param[in] boundaryDisplacementValues displacement values at moving boundaries
    ///
    /// This solves the mesh deformation PDE and updates the geometry.
    ///
    /// @note The input should be a gsMultiPatch where only the patches
    /// corresponding to moving boundaries have meaningful values.
    void updateMesh(const gsMultiPatch<T>& boundaryDisplacementValues)
    {
        if (!m_isInitialized)
            initialize();

        // Update boundary displacement
        m_boundaryDisplacement = boundaryDisplacementValues;

        // Solve mesh deformation
        index_t result = m_aleSolver->updateMesh();

        if (result >= 0)
        {
            gsInfo << "Mesh deformation check: " << result << " elements with negative Jacobian\n";
            if (result > 0)
                gsWarn << "Warning: Mesh has inverted elements!\n";
        }

        // Update geometry with deformed mesh
        gsMultiPatch<T> meshDisplacement;
        m_aleSolver->constructSolution(meshDisplacement);

        // Apply displacement to geometry
        for (size_t p = 0; p < m_geometry.nPatches(); ++p)
        {
            m_geometry.patch(p).coefs() = m_originalGeometry.patch(p).coefs() +
                                          meshDisplacement.patch(p).coefs();
        }
    }

    /// @brief Get current mesh displacement field
    gsMultiPatch<T> getMeshDisplacement() const
    {
        if (!m_isInitialized)
        {
            gsWarn << "gsINSMeshDeformer not initialized, returning zero displacement\n";
            gsMultiPatch<T> zero = m_geometry;
            for (size_t p = 0; p < zero.nPatches(); ++p)
                zero.patch(p).coefs().setZero();
            return zero;
        }

        gsMultiPatch<T> meshDisp;
        m_aleSolver->constructSolution(meshDisp);
        return meshDisp;
    }

    /// @brief Reset mesh to original configuration
    void resetToOriginal()
    {
        m_geometry = m_originalGeometry;

        // Reset ALE solver
        if (m_isInitialized)
        {
            for (size_t p = 0; p < m_boundaryDisplacement.nPatches(); ++p)
                m_boundaryDisplacement.patch(p).coefs().setZero();
        }
    }

    /// @brief Get options for mesh deformation solver
    gsOptionList& options()
    {
        if (!m_isInitialized)
            initialize();
        return m_aleSolver->options();
    }

    /// @brief Get default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addReal("PoissonsRatio", "Poisson's ratio for elasticity methods", 0.4);
        opt.addReal("LocalStiff", "Local stiffening parameter (0=uniform, >0=stiffer near boundaries)", 0.0);
        opt.addSwitch("Check", "Check mesh quality (Jacobian positivity)", true);
        opt.addInt("NumIter", "Number of nonlinear iterations for TINE methods", 5);
        return opt;
    }

    /// @brief Get the name of current deformation method
    std::string methodName() const
    {
        switch(m_method)
        {
            case mesh_deform_method::HE: return "Harmonic Extension (HE)";
            case mesh_deform_method::IHE: return "Incremental Harmonic Extension (IHE)";
            case mesh_deform_method::LE: return "Linear Elasticity (LE)";
            case mesh_deform_method::ILE: return "Incremental Linear Elasticity (ILE)";
            case mesh_deform_method::TINE: return "TINE (Neo-Hookean)";
            case mesh_deform_method::TINE_StVK: return "TINE (St. Venant-Kirchhoff)";
            case mesh_deform_method::BHE: return "Biharmonic Extension (BHE)";
            case mesh_deform_method::IBHE: return "Incremental Biharmonic Extension (IBHE)";
            default: return "Unknown";
        }
    }

    /// @brief Check if initialized
    bool isInitialized() const { return m_isInitialized; }

    /// @brief Save current state (for time step rollback)
    void saveState()
    {
        if (m_isInitialized)
            m_aleSolver->saveState();
    }

    /// @brief Recover saved state
    void recoverState()
    {
        if (m_isInitialized)
            m_aleSolver->recoverState();
    }

    /// @brief Get access to underlying gsALE solver (advanced use)
    gsALE<T>* getALESolver() { return m_aleSolver.get(); }

protected:
    /// Reference to fluid geometry (will be modified)
    gsMultiPatch<T>& m_geometry;

    /// Original undeformed geometry
    gsMultiPatch<T> m_originalGeometry;

    /// Current boundary displacement
    gsMultiPatch<T> m_boundaryDisplacement;

    /// List of moving boundaries
    std::vector<patchSide> m_movingBoundaries;

    /// Persistent interface mapping (disp sides A -> geom sides B)
    gsBoundaryInterface m_interface;

    /// Mesh deformation method
    mesh_deform_method::method m_method;

    /// Internal gsALE solver
    typename memory::unique_ptr<gsALE<T>> m_aleSolver;

    /// Initialization flag
    bool m_isInitialized;
};

} // namespace gismo
