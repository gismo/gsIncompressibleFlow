/** @file gsINSAssemblerBase.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include <gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsINSBlockAssembler.h>
#include <gsIncompressibleFlow/src/gsINSSolverParams.h>

namespace gismo
{

/// @brief A base class for the assemblers of linear systems arising from discretization of the incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template<class T>
class gsINSAssemblerBase
{

protected: // *** Class members ***

    gsINSBlockAssembler<T> m_blockAssembler;
    gsMatrix<T> m_solution;
    bool m_bInitialized;

public: // *** Constructor/destructor ***

    /// @brief  Constructor of the object.
    /// @param params a set of all needed parameters
    gsINSAssemblerBase(gsINSSolverParams<T>& params) :
        m_blockAssembler(params)
    {
    }

    virtual ~gsINSAssemblerBase()
    {
    }

protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers() { GISMO_NO_IMPLEMENTATION }

    /// @brief Re-initialize all members.
    virtual void reinitMembers() { initMembers(); }

public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize() { GISMO_NO_IMPLEMENTATION }

    /**
     * @brief Update the assembler with a new computed solution (e.g. in a new linearization step or time step).
     * 
     * The new solution is saved into the member \c m_solution.
     * 
     * @param[in] solVector a reference to the new solution vector
    */
    virtual void update(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);

        updateAssembly();

        fillMatrix();
        fillRhs();
    }

    /**
     * @brief Update the assembler with a new computed solution (e.g. in a new linearization step or time step).
     * 
     * The new solution is \b not saved into the member \c m_solution.
     * 
     * @param[in] solVector a reference to the new solution vector
    */
    virtual void updatePicard(const gsMatrix<T> & solVector)
    {
        GISMO_ASSERT(m_bInitialized, "Assembler must be initialized first, call initialize()");

        m_blockAssembler.updateCurrentSolField(solVector, false);

        updatePicardAssembly();

        fillMatrix();
        fillRhs();
    }

    /// @brief Assemble the Stokes matrix.
    /// @param[out] stokesMatrix    a reference to a sparse matrix to be filled
    /// @param[out] stokesRhs       a reference to a RHS vector to be filled
    virtual void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        m_blockAssembler.fillStokesSystem_into(stokesMatrix, stokesRhs);
    }

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @return     the resulting solution field
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk) const
    {
        return m_blockAssembler.constructSolution(solVector, unk);
    }

    /// @brief Compute flow rate through a side of a given patch.
    /// @param[in] patch        the given patch ID
    /// @param[in] side         the given patch side
    /// @param[in] solution     solution vector to compute the flow rate from
    virtual T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const
    {
        return m_blockAssembler.computeFlowRate(patch, side, solution);
    }

    /// @brief Fix zero pressure on a given patch side.
    /// @param[in] patch  the given patch ID
    /// @param[in] side   the given patch side
    void addPressureOutletCondition(int patch, boxSide side);

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown (0 - velocity, 1 - pressure)
    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        m_blockAssembler.markDofsAsEliminatedZeros(boundaryDofs, unk);

        reinitMembers();
    }

    /**
     * @brief Identify pressure boundary DOF IDs and assemble the Robin boundary block if needed.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[in] bndIn     a list of pairs (patch, side) corresponding to inflow boundary
     * @param[in] bndOut    a list of pairs (patch, side) corresponding to outflow boundary
     * @param[in] bndWall   a list of pairs (patch, side) corresponding to solid wall
     * @param[in] bcType    boundary conditions strategy for PCD (see the documentation of applyPCDboundaryConditions for description)
     */ 
    void preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
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
    virtual void fillPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType, bool assembAp, bool assembFp, bool lumping);

    /**
     * @brief Get the pressure Poisson and convection blocks for the specified boundary parts.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[out]  Ap          a reference to the pressure Poisson matrix to be filled
     * @param[out]  Fp          a reference to the pressure convection matrix to be filled
     * @param[in]   bcType      boundary conditions strategy for PCD (see the documentation of applyPCDboundaryConditions for description)
     * @param[in]   bndIn       a list of pairs (patch, side) corresponding to inflow boundary
     * @param[in]   bndOut      a list of pairs (patch, side) corresponding to outflow boundary
     * @param[in]   bndWall     a list of pairs (patch, side) corresponding to solid wall
     * @param[in]   assembAp    assemble the pressure Poisson matrix
     * @param[in]   assembFp    assemble the part of the pressure convection matrix which is equal to the pressure Poisson matrix 
     * @param[in]   lumping     use mass lumping
     */ 
    void getPCDblocks(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType, bool assembAp, bool assembFp, bool lumping)
    {
        m_blockAssembler.preparePCDboundary(bndIn, bndOut, bndWall, bcType);
        this->fillPCDblocks(Ap, Fp, bcType, assembAp, assembFp, lumping);
    }


protected: // *** Member functions ***

    /// @brief Assemble the nonlinear part of the problem.
    virtual void updateAssembly()
    {
        m_blockAssembler.assembleNonlinearPart();
    }

    /// @brief Assemble the nonlinear part of the problem for inner Picard linearization process in the unsteady case.
    virtual void updatePicardAssembly()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Fill the base matrix (the linear part).
    virtual void fillBase() { GISMO_NO_IMPLEMENTATION }

    /// @brief Fill the global matrix.
    virtual void fillMatrix() { GISMO_NO_IMPLEMENTATION }

    /// @brief Fill the right-hand side.
    virtual void fillRhs() { GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    /// @brief Returns the assembled matrix.
    virtual const gsSparseMatrix<T> & matrix() const { GISMO_NO_IMPLEMENTATION }

    /// @brief Returns the assembled right-hand side.
    virtual const gsMatrix<T> & rhs() const { GISMO_NO_IMPLEMENTATION }

    /// @brief Returns the velocity mass matrix.
    virtual const gsSparseMatrix<T>& getVelocityMassMatrix()
    { 
        return m_blockAssembler.getBlockM();
    }

    /// @brief Returns the pressure mass matrix.
    virtual const gsSparseMatrix<T>& getPressureMassMatrix()
    {
        return m_blockAssembler.getBlockMp();
    }

    /// @brief Set a given solution vector as current solution.
    /// @param[in] solVector the given solution
    void setSolution(gsMatrix<T> solVector)
    {
        m_solution = solVector;
        m_blockAssembler.updateCurrentSolField(solVector, true);
    }

    /// @brief Returns the total number of DOFs (the matrix size).
    virtual int numDofs() const { return m_blockAssembler.numDofs(); }

    /// @brief  Returns the number of velocity DOFs (one velocity component).
    virtual int getUdofs() const { return m_blockAssembler.getUdofs(); }

    /// @brief Returns the number of pressure DOFs.
    virtual int getPdofs() const { return m_blockAssembler.getPdofs(); }

    /// @brief Returns the DOF shift of pressure (i.e. the total number of velocity DOFs).
    virtual int getPshift() const { return m_blockAssembler.getPshift(); }

    int getTarDim() const { return m_blockAssembler.getTarDim(); }
    real_t getViscosity() const { return m_blockAssembler.getViscosity(); }
    bool isInitialized() { return m_bInitialized; }
    const gsINSBlockAssembler<T>& getBlockAssembler() const { return m_blockAssembler; }
    gsINSBlockAssembler<T>& getBlockAssembler() { return m_blockAssembler; }
    virtual const gsMatrix<T>& getSolution() const { return m_solution; }
    //bool isRotation() const { return m_blockAssembler.isRotation(); }
    

};

} //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSAssemblerBase.hpp)
#endif