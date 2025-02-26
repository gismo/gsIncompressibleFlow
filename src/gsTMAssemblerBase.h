/** @file gsINSAssembler.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
//#include <gsIncompressibleFlow/src/gsTMVisitors.h>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for incompressible Navier-Stokes assemblers.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMAssemblerBase: public gsFlowAssemblerBase<T, MatOrder>
{

public:
    typedef gsFlowAssemblerBase<T, MatOrder> Base;


protected: // *** Class members ***

    T m_viscosity;
    gsSparseMatrix<T, MatOrder> m_baseMatrix, m_matrix;
    gsMatrix<T> m_baseRhs, m_rhs;

    gsBoundaryConditions<T> m_bc;
    std::vector<gsMultiBasis<T> > m_bases;
    //std::vector<gsDofMapper> m_dofMappersTM;
    index_t m_nnzPerRowTM;
    index_t numTMvars;
    std::vector<index_t> m_kdofs;
    gsField<T> m_currentSolField;
    //index_t m_TMdofsAll;
    bool m_bInitialized;
    gsMatrix<T> m_RANSsolution;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_dofs;
    using Base::m_tarDim;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_isBaseReady;
    using Base::m_isSystemReady;

public: // *** Base class member functions ***

    //using Base::getBases;
    using Base::getPatches;
    using Base::getAssemblerOptions;
    //using Base::getBCs;

public: // *** Constructor/destructor ***

    gsTMAssemblerBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { 
        initMembers();
    }

    virtual ~gsTMAssemblerBase()
    { }

protected: // *** Member functions ***

    /// @brief Initialize the class members.
    virtual void initMembers();

    /// @brief Update sizes of members (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateSizes();

    /// @brief Update the current solution field stored in the assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner nonlinear iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol);

    /// @brief Fill the velocity-velocity block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    velocity-velocity block (either for one velocity component or the whole block for all components)
    void fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat, const index_t unk);

public: // *** Member functions ***

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk);

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk) const;

public: // *** Getters/setters ***

    /// @brief Returns the viscosity value.
    T getViscosity() const { return m_viscosity; }

    /**
     * @brief Returns a reference to the discretization bases.
     *
     * In the case of velocity and pressure, the velocity basis is stored first, the  pressure basis is second.
     * 
     * There is also a const version returning a const reference.
     */
    virtual std::vector< gsMultiBasis<T> >& getBases() { return m_bases; }
    virtual const std::vector< gsMultiBasis<T> >& getBases() const { return m_bases; }

    virtual const gsBoundaryConditions<T>& getBCs() const { return m_bc; }

    /// @brief Returns the assembled matrix.
    virtual const gsSparseMatrix<T, MatOrder>& matrix() const
    {
        GISMO_ASSERT(m_isSystemReady, "Matrix not ready, is fillGlobalSyst in the solver params true? If so, update() must be called first.");
        return m_matrix;
    }

    /// @brief Returns the assembled right-hand side.
    virtual const gsMatrix<T>& rhs() const
    {
        GISMO_ASSERT(m_isSystemReady, "Rhs not ready, is fillGlobalSyst in the solver params true? If so, update() must be called first.");
        return m_rhs;
    }

    /// @brief Returns a const reference to the current computed solution.
    const gsMatrix<T>& getSolution() const { return m_solution; }

    void setRANSsolution(gsMatrix<T>& sol) { m_RANSsolution = sol; }

    /// @brief Returns the number of DOFs for the i-th turbulent unknown.
    index_t numDofsUnk(size_t i);

    /// @brief Returns the total number of DOFs (the matrix size).
    index_t numDofs() const { return m_dofs; }

    bool isInitialized() { return m_bInitialized; }

}; 


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMAssemblerBase.hpp)
#endif