/** @file gsINSAssemblerSteady.h
    
    @brief Steady incompressible Navier-Stokes assembler. 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSAssemblerBase.h>

namespace gismo
{

/// @brief Assembler of linear systems arising from discretization of the steady incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template<class T>
class gsINSAssemblerSteady : public gsINSAssemblerBase<T>
{

public:
    typedef gsINSAssemblerBase<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_baseMatrix;
    gsMatrix<T> m_baseRhs;
    gsSparseMatrix<T> m_matrix;
    gsMatrix<T> m_rhs;
    bool m_bMatrixReady;
    bool m_bRhsReady;

protected: // *** Base class members ***

    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;

public: // *** Constructor/destructor ***

    gsINSAssemblerSteady(gsINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~gsINSAssemblerSteady()
    {
    }


public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize()
    {
        m_blockAssembler.assembleLinearStokesPart();
        m_blockAssembler.assembleBlockNpattern();

        fillBase();

        m_bInitialized = true;
    }

    /// @brief Assemble the nonlinear part of the problem.
    virtual void updateAssembly()
    {
        Base::updateAssembly();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    /// @brief Re-initialize all members.
    virtual void reinitMembers() { initMembers(); }

    /// @brief Fill the base matrix (the linear part).
    virtual void fillBase();

    /// @brief Fill the global matrix.
    virtual void fillMatrix();

    /// @brief Fill the right-hand side.
    virtual void fillRhs()
    {
        m_rhs.noalias() = m_baseRhs + m_blockAssembler.getRhsN();

        m_bRhsReady = true;
    }

public: // *** Getters/setters ***

    /// @brief Returns the assembled matrix.
    virtual const gsSparseMatrix<T> & matrix() const
    {
        GISMO_ASSERT(m_bMatrixReady, "Matrix not ready, update() must be called first");
        return m_matrix;
    }

    /// @brief Returns the assembled right-hand side.
    virtual const gsMatrix<T> & rhs() const
    {
        GISMO_ASSERT(m_bRhsReady, "Rhs not ready, update() must be called first");
        return m_rhs;
    }

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSAssemblerSteady.hpp)
#endif
