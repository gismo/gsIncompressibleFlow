/** @file gsINSAssemblerUnsteady.h

    @brief Unsteady incompressible Navier-Stokes assembler. 

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

/// @brief Assembler of linear systems arising from discretization of the unsteady incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template<class T>
class gsINSAssemblerUnsteady : public gsINSAssemblerBase<T>
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
    T m_timeStepSize;

protected: // *** Base class members ***

    using Base::m_blockAssembler;
    using Base::m_solution;
    using Base::m_bInitialized;

public: // *** Constructor/destructor ***

    gsINSAssemblerUnsteady(gsINSSolverParams<T>& params) : Base(params)
    {
        m_timeStepSize = params.options().getReal("timeStep");
        initMembers();
    }

    virtual ~gsINSAssemblerUnsteady()
    {
    }

public: // *** Member functions ***

    /// @brief Initialize the assembler.
   virtual void initialize()
    {
        m_blockAssembler.assembleLinearStokesPart();
        m_blockAssembler.assembleBlockNpattern();

        fillBase();

        if (!m_baseMatrix.isCompressed())
            m_baseMatrix.makeCompressed();

        m_bInitialized = true;
    }

    /// @brief Initialize the assembler with initial condition given by coefficient vector.
    /// @param[in] initialConditionCoeffVector coefficient vector of the initial condition
    virtual void initialize(gsMatrix<T> & initialConditionCoeffVector)
    {
        initialize();
        Base::update(initialConditionCoeffVector);
    }

    /// @brief Assemble the nonlinear part of the problem.
    virtual void updateAssembly()
    {
        Base::updateAssembly();
    }

    /// @brief Assemble the nonlinear part of the problem for inner Picard linearization process.
    virtual void updatePicardAssembly()
    {
        m_blockAssembler.assembleNonlinearPart();
    }

    /**
     * @brief Fill the pressure Poisson and convection blocks.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[out]  Ap          a reference to the pressure Poisson matrix to be filled
     * @param[out]  Fp          a reference to the pressure convection matrix to be filled
     * @param[in]   bcType      boundary conditions strategy for PCD (see the documentation of applyPCDboundaryConditions for description)
     * @param[in]   assembAp    assemble the pressure Poisson matrix
     * @param[in]   assembFp    assemble the part of the pressure convection matrix which is equal to the pressure Poisson matrix 
     * @param[in]   lumping     use mass lumping
     */ 
    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        Fp = (1. / m_timeStepSize) * m_blockAssembler.getBlockMp();

     
        Base::fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
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
    void fillRhs();


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

    /// @brief Returns the time step size.
    const T getTimeStepSize() const { return m_timeStepSize; }

    /// @brief Change the time step size.
    /// @param[in] timeStepSize the new time step size
    virtual void changeTimeStep(const T timeStepSize);

}; // class gsINSAssemblerUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSAssemblerUnsteady.hpp)
#endif