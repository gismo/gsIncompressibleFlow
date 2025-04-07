/** @file gsTMAssemblerBase.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for turbulence models assemblers.
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
    index_t m_nnzPerRowTM;
    index_t numTMvars;
    std::vector<index_t> m_kdofs;
    gsField<T> m_currentSolField;
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
    using Base::m_isInitialized;

public: // *** Base class member functions ***

    using Base::getBases;
    using Base::getPatches;
    using Base::getAssemblerOptions;
    using Base::getBCs;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
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

public: // *** Member functions ***

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk);

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  unk         the considered unknown
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    gsField<T> constructSolution(int unk) const
    { return constructSolution(m_solution, unk); }
    
    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk) const;

public: // *** Getters/setters ***

    /// @brief Returns the kinematic viscosity value.
    T getViscosity() const { return m_viscosity; }

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

}; 


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMAssemblerBase.hpp)
#endif