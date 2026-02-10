/** @file gsTMAssemblerSST.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsTMAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsTMVisitors.h>
#include <gsIncompressibleFlow/src/gsTMModels.h>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A class for the assembler of the k-omega SST turbulence model.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMAssemblerSST: public gsTMAssemblerBase<T, MatOrder>
{

public:
    typedef gsTMAssemblerBase<T, MatOrder> Base;


protected: // *** Class members ***

    gsTMVisitorLinearSST<T, MatOrder> m_visitorLinearSST_K, m_visitorLinearSST_O;
    gsTMVisitorTimeIterationSST<T, MatOrder> m_visitorTimeIterationSST_K, m_visitorTimeIterationSST_O;
    gsTMVisitorNonlinearSST<T, MatOrder> m_visitorNonlinearSST_K, m_visitorNonlinearSST_O;
    gsTMVisitorSSTTCSDStabilization_time<T, MatOrder> m_visitorSST_TCSD_time_K, m_visitorSST_TCSD_time_O;
    gsTMVisitorSSTTCSDStabilization_advection<T, MatOrder> m_visitorSST_TCSD_advection_K, m_visitorSST_TCSD_advection_O;

    gsSparseMatrix<T, MatOrder> m_blockLinearK, m_blockLinearO, m_blockTimeIterationK, m_blockTimeIterationO, m_blockNonlinearK, m_blockNonlinearO;
    gsSparseMatrix<T, MatOrder> m_matSST_TCSD_time_K, m_matSST_TCSD_time_O, m_matSST_TCSD_advection_K, m_matSST_TCSD_advection_O;
    gsMatrix<T> m_rhsLinearK, m_rhsLinearO, m_rhsTimeIterationK, m_rhsTimeIterationO, m_rhsNonlinearK, m_rhsNonlinearO;
    gsMatrix<T> m_rhsSST_TCSD_time_K, m_rhsSST_TCSD_time_O, m_rhsSST_TCSD_advection_K, m_rhsSST_TCSD_advection_O;

    gsField<T> m_currentFieldK, m_currentFieldO;
    gsField<T> m_oldTimeFieldK, m_oldTimeFieldO;

    gsField<T> m_distanceField;
    
    typename gsTMModelData<T>::Ptr m_TMModelPtr; 

    real_t m_kin, m_kwall, m_oin, m_owall;

    
protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_dofs;
    using Base::m_tarDim;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_isBaseReady;
    using Base::m_isSystemReady;

    using Base::m_viscosity;
    using Base::m_baseMatrix;
    using Base::m_matrix;
    using Base::m_baseRhs;
    using Base::m_rhs;
    using Base::m_bc;
    using Base::m_nnzPerRowTM;
    using Base::m_numVars;
    using Base::m_kdofs;
    using Base::m_isInitialized;

public: // *** Base class member functions ***

    using Base::getBasis;
    using Base::getPatches;
    using Base::getAssemblerOptions;
    using Base::getBCs;
    using Base::constructSolution;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMAssemblerSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr):
    Base(paramsPtr), m_TMModelPtr(TMModelPtr)
    { 
        m_numVars = 2;
        initMembers();
    }

    virtual ~gsTMAssemblerSST()
    { }

protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change after constructing the assembler).
    virtual void updateSizes();

    /// @brief Update the current solution field stored in the assembler.
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner nonlinear iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T>& solVector, bool updateSol);

    /// @brief Assemble the linear part of the problem.
    virtual void assembleLinearPart();

    /// @brief Assemble the non-linear part of the problem.
    virtual void assembleNonlinearPart();

    /// @brief Fill the block matrix into the global matrix.
    /// @param globalMat[out]   global matrix
    /// @param sourceMat[in]    block matrix (for one turbulence variable)
    void fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat, index_t unk);

    /// @brief Fill the linear part of the global matrix and right-hand side.
    virtual void fillBaseSystem();

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();


public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize();

    /// @brief Update the assembler in new nonlinear iteration.
    /// @param solVector 
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);

    
public: // *** Getters/setters ***

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

    /// @brief Returns the block of the linear system corresponding to the first turbuelnce variable.
    /// @param[in] linPartOnly if true, returns only the linear part of this block
    virtual gsSparseMatrix<T, MatOrder> getBlockKK(bool linPartOnly = false);

    /// @brief Returns the block of the linear system corresponding to the second turbuelnce variable.
    /// @param[in] linPartOnly if true, returns only the linear part of this block
    virtual gsSparseMatrix<T, MatOrder> getBlockOO(bool linPartOnly = false);

    /// @brief /// @brief Returns the part of the right-hand side corresponding to the first turbuelnce variable.
    virtual gsMatrix<T> getRhsK() const;

    /// @brief /// @brief Returns the part of the right-hand side corresponding to the second turbuelnce variable.
    virtual gsMatrix<T> getRhsO() const;

}; // gsTMAssemblerSST

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMAssemblerSST.hpp)
#endif