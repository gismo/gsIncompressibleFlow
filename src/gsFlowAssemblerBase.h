/** @file gsFlowAssemblerBase.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>

namespace gismo
{

/// @brief  A           base class for all assemblers in gsIncompressibleFlow.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsFlowAssemblerBase
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    index_t m_dofs;
    short_t m_tarDim;
    bool m_isInitialized;
    bool m_isBaseReady;
    bool m_isSystemReady;
    std::vector<gsDofMapper> m_dofMappers;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution;


public: // *** Constructor/destructor ***

    gsFlowAssemblerBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
    { }

    virtual ~gsFlowAssemblerBase()
    { }


protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateSizes()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Update the DOF mappers in all visitors (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateDofMappers()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  basisID     the index of the basis corresponding to \a unk (same as \a unk in this case)
    /// @param[in]  ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofs(const index_t unk, const index_t basisID, gsMatrix<T>& ddofVector);

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using interpolation.
    /// @param[in]  unk          the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper       reference to the DOF mapper for \a unk
    /// @param[in]  mbasis       reference to the basis corresponding to \a unk
    /// @param[out] ddofVector   reference to the vector where computed coefficients will be stored
    void computeDirichletDofsIntpl(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using L2-projection.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper      reference to the DOF mapper for \a unk
    /// @param[in]  mbasis      reference to the basis corresponding to \a unk
    /// @param[out] ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofsL2Proj(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);

    /// @brief Assemble a matrix block.
    /// @param[in]  visitor     visitor for the required block
    /// @param[in]  testBasisID ID of the test basis
    /// @param[out] block       the resulting matrix block
    /// @param[out] blockRhs    right-hand side for the matrix block (arising from eliminated Dirichlet DOFs)
    void assembleBlock(gsFlowVisitor<T, MatOrder>& visitor, index_t testBasisID, gsSparseMatrix<T, MatOrder>& block, gsMatrix<T>& blockRhs);

    /// @brief Assemble the right-hand side.
    /// @param[in]  visitor     visitor for the right-hand side
    /// @param[in]  testBasisID ID of the test basis
    /// @param[out] rhs         the resulting right-hand side vector
    void assembleRhs(gsFlowVisitor<T, MatOrder>& visitor, index_t testBasisID, gsMatrix<T>& rhs);

    /// @brief Assemble the linear part of the problem.
    virtual void assembleLinearPart()
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Assemble the nonlinear part of the problem.
    virtual void assembleNonlinearPart()
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Assemble all that needs to be updated in each nonlinear iteration.
    virtual void updateAssembly();


    /// @brief Update current solution field stored in the assembler.
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
    {GISMO_NO_IMPLEMENTATION}


public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize();

    /// @brief Update the assembler in new nonlinear iteration.
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk) const
    {GISMO_NO_IMPLEMENTATION}


public: // *** Getters/setters ***

    /// @brief Returns the number of degrees of freedom (DOFs).
    index_t numDofs() const 
    { 
        GISMO_ASSERT(m_dofs > 0, "Something went wrong, number of DOFs is zero!");
        return m_dofs; 
    }

    /// @brief Returns the target dimension.
    short_t getTarDim() const { return m_tarDim; }

    /// @brief Returns true if the assembler has been initialized.
    bool isInitialized() { return m_isInitialized; }

    /**
     * @brief Returns a const reference to the DOF mappers.
     *
     * In the case of velocity and pressure, the mapper for velocity is stored first, the mapper for pressure is second.
     */
    const std::vector<gsDofMapper>& getMappers() const { return m_dofMappers; }

    /**
     * @brief Returns a const reference to the vectors of coefficients at the Dirichlet boundaries.
     *
     * In the case of velocity and pressure, the vector of velocity coefficients is stored first, the vector of pressure coefficients is second.
     */
    const std::vector<gsMatrix<T> >& getDirichletDofs() const { return m_ddof; }

    /// @brief Returns a const reference to the current computed solution.
    const gsMatrix<T>& getSolution() const { return m_solution; }

    /// @brief Returns a const reference to the multipatch representing the computational domain.
    const gsMultiPatch<T>& getPatches() const { return m_paramsPtr->getPde().patches(); }

    /**
     * @brief Returns a reference to the discretization bases.
     *
     * In the case of velocity and pressure, the velocity basis is stored first, the  pressure basis is second.
     * 
     * There is also a const version returning a const reference.
     */
    std::vector< gsMultiBasis<T> >& getBases() { return m_paramsPtr->getBases(); }
    const std::vector< gsMultiBasis<T> >& getBases() const { return m_paramsPtr->getBases(); }

    /// @brief Returns a const reference to the boundary conditions.
    const gsBoundaryConditions<T>& getBCs() const { return m_paramsPtr->getBCs(); }

    /// @brief Returns a pointer to the right-hand-side function.
    const gsFunction<T>* getRhsFcn() const { return m_paramsPtr->getPde().rhs(); }

    /// @brief Returns the assembler options.
    gsAssemblerOptions getAssemblerOptions() const { return m_paramsPtr->assemblerOptions(); }
    
    /// @brief Returns the flow solver option list.
    gsOptionList options() const { return m_paramsPtr->options(); }

    /// @brief Returns the assembled matrix.
    virtual const gsSparseMatrix<T, MatOrder>& matrix() const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled matrix for unknown with index \a unk (e.g. from two-equation turbulence models).
    /// @param[in] unk index of the unknown
    virtual const gsSparseMatrix<T, MatOrder>& matrix(index_t unk) const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the mass matrix for unknown with index \a unk. There is also a const version.
    /// @param[in] unkID index of the unknown
    virtual gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID)
    {GISMO_NO_IMPLEMENTATION}
    
    virtual const gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID) const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled right-hand side.
    virtual const gsMatrix<T>& rhs() const
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Returns the assembled right-hand side for unknown with index \a unk (e.g. from two-equation turbulence models).
    /// @param[in] unk index of the unknown
    virtual const gsMatrix<T>& rhs(index_t unk) const
    {GISMO_NO_IMPLEMENTATION}

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowAssemblerBase.hpp)
#endif