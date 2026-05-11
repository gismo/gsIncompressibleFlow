/** @file gsINSVisitorsALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTermsALE.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSTermsSUPG.h>

namespace gismo
{

/// @brief ALE visitor for nonlinear velocity-velocity terms in incompressible Navier-Stokes
/// @tparam T coefficient type
/// @tparam MatOrder sparse matrix storage order
template <class T, int MatOrder = RowMajor>
class gsINSVisitorUUnonlinALE : public gsINSVisitorUU<T, MatOrder>
{
public:
    typedef gsINSVisitorUU<T, MatOrder> Base;

protected:
    using Base::m_terms;

    // Target dimension
    index_t m_tarDim;

    // Mesh velocity field (only for ALE convection term)
    const gsField<T>* m_meshVelField;

    // ALE convection term pointer
    gsINSTerm_ALEConvection<T>* m_aleConvectionTerm;

    // SUPG convection ALE term pointer (for mesh velocity propagation)
    gsINSTerm_SUPGConvectionALE<T>* m_supgConvectionTerm;

    // SUPG viscous ALE term pointer (for mesh velocity propagation)
    gsINSTerm_SUPGViscousALE<T>* m_supgViscousTerm;

public:
    /// @brief Constructor with paramsPtr
    gsINSVisitorUUnonlinALE(typename gsFlowSolverParams<T>::Ptr paramsPtr,
                           index_t targetDim = 2,
                           const gsField<T>* meshVelField = nullptr)
    : Base(paramsPtr),
      m_tarDim(targetDim),
      m_meshVelField(meshVelField),
      m_aleConvectionTerm(nullptr),
      m_supgConvectionTerm(nullptr),
      m_supgViscousTerm(nullptr)
    {
        // Terms will be created in initialize() via defineTerms()
    }

    /// @brief Destructor
    ~gsINSVisitorUUnonlinALE()
    {
        // Terms are deleted by base class
    }

protected:
    /// @brief Define ALE convection terms.
    virtual void defineTerms() override
    {
        m_aleConvectionTerm = new gsINSTerm_ALEConvection<T>(m_meshVelField, m_tarDim);
        this->m_terms.push_back(m_aleConvectionTerm);

        if (this->m_paramsPtr->options().getSwitch("SUPG_NS"))
        {
            T visc = this->m_paramsPtr->getPde().viscosity();
            T dt   = this->m_paramsPtr->options().getReal("timeStep");

            m_supgConvectionTerm = new gsINSTerm_SUPGConvectionALE<T>(
                m_meshVelField, visc, dt, m_tarDim);
            this->m_terms.push_back(m_supgConvectionTerm);
        }
    }

public:
    /// @brief Set mesh velocity field.
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
        if (m_aleConvectionTerm)
            m_aleConvectionTerm->setMeshVelocityField(meshVelField);
        if (m_supgConvectionTerm)
            m_supgConvectionTerm->setMeshVelocityField(meshVelField);
        if (m_supgViscousTerm)
            m_supgViscousTerm->setMeshVelocityField(meshVelField);
    }

    /// @brief Update params pointer (useful for thread-private visitor copies)
    void setParamsPtr(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    {
        this->m_paramsPtr = paramsPtr;
    }

    /// @brief Set current solution
    void setCurrentSolution(const gsField<T>& solution)
    {
        for (size_t i = 0; i < m_terms.size(); ++i)
        {
            gsFlowTermNonlin<T>* nonlinTerm = dynamic_cast<gsFlowTermNonlin<T>*>(m_terms[i]);
            if (nonlinTerm)
                nonlinTerm->setCurrentSolution(const_cast<gsField<T>&>(solution));
        }
    }
};

// ===================================================================================================================

/// @brief ALE visitor for SUPG pressure gradient in the PU block:
///        tau_SUPG * ((u-w) . grad phi_i) * (d psi_j / d x_d)
/// Inherits from gsINSVisitorPU for correct localToGlobal mapping.
/// @tparam T coefficient type
/// @tparam MatOrder sparse matrix storage order
template <class T, int MatOrder = RowMajor>
class gsINSVisitorPU_SUPG_ALE : public gsINSVisitorPU<T, MatOrder>
{
public:
    typedef gsINSVisitorPU<T, MatOrder> Base;

protected:
    using Base::m_terms;

    index_t m_tarDim;
    const gsField<T>* m_meshVelField;
    gsINSTerm_SUPGPressureGradALE<T>* m_supgPressGradTerm;

public:
    /// @brief Constructor
    gsINSVisitorPU_SUPG_ALE(typename gsFlowSolverParams<T>::Ptr paramsPtr,
                             index_t targetDim = 2,
                             const gsField<T>* meshVelField = nullptr)
    : Base(paramsPtr),
      m_tarDim(targetDim),
      m_meshVelField(meshVelField),
      m_supgPressGradTerm(nullptr)
    { }

    ~gsINSVisitorPU_SUPG_ALE()
    {
        // Terms are deleted by base class
    }

protected:
    virtual void defineTerms() override
    {
        T visc = this->m_paramsPtr->getPde().viscosity();
        T dt   = this->m_paramsPtr->options().getReal("timeStep");

        m_supgPressGradTerm = new gsINSTerm_SUPGPressureGradALE<T>(
            m_meshVelField, visc, dt, m_tarDim);
        this->m_terms.push_back(m_supgPressGradTerm);
    }

public:
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
        if (m_supgPressGradTerm)
            m_supgPressGradTerm->setMeshVelocityField(meshVelField);
    }

    void setCurrentSolution(const gsField<T>& solution)
    {
        for (size_t i = 0; i < m_terms.size(); ++i)
        {
            gsFlowTermNonlin<T>* nonlinTerm = dynamic_cast<gsFlowTermNonlin<T>*>(m_terms[i]);
            if (nonlinTerm)
                nonlinTerm->setCurrentSolution(const_cast<gsField<T>&>(solution));
        }
    }
};

// ===================================================================================================================

/// @brief ALE visitor for SUPG time discretization in the UU block:
///        (1/dt) * tau_SUPG * ((u-w) . grad phi_i) * phi_j
/// @tparam T coefficient type
/// @tparam MatOrder sparse matrix storage order
template <class T, int MatOrder = RowMajor>
class gsINSVisitorUU_SUPGTimeDiscrALE : public gsINSVisitorUU<T, MatOrder>
{
public:
    typedef gsINSVisitorUU<T, MatOrder> Base;

protected:
    using Base::m_terms;

    index_t m_tarDim;
    const gsField<T>* m_meshVelField;
    gsINSTerm_SUPGTimeDiscrALE<T>* m_supgTimeDiscrTerm;

public:
    /// @brief Constructor
    gsINSVisitorUU_SUPGTimeDiscrALE(typename gsFlowSolverParams<T>::Ptr paramsPtr,
                                     index_t targetDim = 2,
                                     const gsField<T>* meshVelField = nullptr)
    : Base(paramsPtr),
      m_tarDim(targetDim),
      m_meshVelField(meshVelField),
      m_supgTimeDiscrTerm(nullptr)
    { }

    ~gsINSVisitorUU_SUPGTimeDiscrALE()
    {
        // Terms are deleted by base class
    }

protected:
    virtual void defineTerms() override
    {
        m_supgTimeDiscrTerm = new gsINSTerm_SUPGTimeDiscrALE<T>(
            m_meshVelField,
            this->m_paramsPtr->getPde().viscosity(),
            this->m_paramsPtr->options().getReal("timeStep"),
            m_tarDim);
        this->m_terms.push_back(m_supgTimeDiscrTerm);
    }

public:
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
        if (m_supgTimeDiscrTerm)
            m_supgTimeDiscrTerm->setMeshVelocityField(meshVelField);
    }

    void setCurrentSolution(const gsField<T>& solution)
    {
        for (size_t i = 0; i < m_terms.size(); ++i)
        {
            gsFlowTermNonlin<T>* nonlinTerm = dynamic_cast<gsFlowTermNonlin<T>*>(m_terms[i]);
            if (nonlinTerm)
                nonlinTerm->setCurrentSolution(const_cast<gsField<T>&>(solution));
        }
    }
};

} // namespace gismo
