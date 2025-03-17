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

/// @brief              A base class for incompressible Navier-Stokes assemblers.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMAssemblerSST: public gsTMAssemblerBase<T, MatOrder>
{

public:
    typedef gsTMAssemblerBase<T, MatOrder> Base;


protected: // *** Class members ***

    //gsBoundaryConditions<T> m_bc;
    gsTMVisitorLinearSST<T, MatOrder> m_visitorLinearSST_K, m_visitorLinearSST_O;
    gsTMVisitorTimeIterationSST<T, MatOrder> m_visitorTimeIterationSST_K, m_visitorTimeIterationSST_O;
    gsTMVisitorNonlinearSST<T, MatOrder> m_visitorNonlinearSST_K, m_visitorNonlinearSST_O;

    gsSparseMatrix<T, MatOrder> m_blockLinearK, m_blockLinearO, m_blockTimeIterationK, m_blockTimeIterationO, m_blockNonlinearK, m_blockNonlinearO;
    gsMatrix<T> m_rhsLinearK, m_rhsLinearO, m_rhsTimeIterationK, m_rhsTimeIterationO, m_rhsNonlinearK, m_rhsNonlinearO;

    gsField<T> m_currentFieldK, m_currentFieldO;
    gsField<T> m_oldTimeFieldK, m_oldTimeFieldO;

    gsField<T> m_distanceField;
    //gsDistanceField<T> dField;

    //typename SSTModel<T>::tmePtr m_SSTPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr; 

    real_t m_kin, m_kwall, m_oin, m_owall;

    //bool m_isMassMatReady;
    //std::vector< gsSparseMatrix<T, MatOrder> > m_massMatBlocks;
    
protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_dofs;
    using Base::m_tarDim;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_isBaseReady;
    using Base::m_isSystemReady;

    using Base::m_viscosity;
    using Base::m_baseMatrix;
    using Base::m_matrix;
    using Base::m_baseRhs;
    using Base::m_rhs;
    using Base::m_bases;
    using Base::m_bc;
    using Base::m_nnzPerRowTM;
    using Base::numTMvars;
    using Base::m_kdofs;
    using Base::m_isInitialized;

public: // *** Base class member functions ***

    using Base::getBases;
    using Base::getPatches;
    using Base::getAssemblerOptions;
    using Base::getBCs;
    using Base::constructSolution;

public: // *** Constructor/destructor ***

    gsTMAssemblerSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr):
    Base(paramsPtr), m_TMModelPtr(TMModelPtr)
    { 
        initMembers();
    }

    virtual ~gsTMAssemblerSST()
    { }

protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateSizes();

    /// @brief Update the DOF mappers in all visitors (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateDofMappers();

    /// @brief Update the current solution field stored in the assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner nonlinear iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T>& solVector, bool updateSol);

    /// @brief Assemble the linear part of the problem.
    virtual void assembleLinearPart();

    /// @brief Assemble the linear part of the problem.
    virtual void assembleNonlinearPart();

    /// @brief Fill the velocity-velocity block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    velocity-velocity block (either for one velocity component or the whole block for all components)
    void fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat, index_t unk);

    /// @brief Fill the linear part of the global matrix and right-hand side.https://www.twitch.tv/mikeses
    virtual void fillBaseSystem();

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();


    // --- PCD member functions ---

    // void findPressureBoundaryIDs();

    // void findPressureBoundaryPartIDs(std::vector<std::pair<int, boxSide> > bndPart, std::vector<index_t>& idVector);


public: // *** Member functions ***

    /// @brief Initialize the assembler.
    virtual void initialize();

    /// @brief 
    /// @param solVector 
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown
    //virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk);

    /// @brief Fill the matrix and right-hand side for the Stokes problem.
    //virtual void fillStokesSystem(gsSparseMatrix<T, MatOrder>& stokesMat, gsMatrix<T>& stokesRhs);

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    //virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk) const;

    
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

    // void setSSTModelEvaluator(typename SSTModel<T>::Ptr SSTPtr)
    // {
    //     m_SSTPtr = SSTPtr;
    // }

    //typename SSTModel<T>::tmePtr getSSTModelEvaluator() { return m_SSTPtr; }

    /// @brief Returns the number of DOFs for the i-th unknown.
    /// @param[in] i    0 - velocity, 1 - pressure
    //index_t numDofsUnk(index_t i);

    /*
    /// @brief  Returns the number of velocity DOFs (one velocity component).
    index_t getUdofs() const { return m_udofs; }

    /// @brief Returns the number of pressure DOFs.
    index_t getPdofs() const { return m_pdofs; }

    /// @brief Returns the DOF shift of pressure (i.e. the total number of velocity DOFs).
    index_t getPshift() const { return m_pshift; }
    */

    /// @brief Returns the velocity-velocity block of the linear system.
    /// @param[in] linPartOnly if true, returns only the linear part of the velocity-velocity block
    virtual gsSparseMatrix<T, MatOrder> getBlockKK(bool linPartOnly = false);

    /// @brief Returns the velocity-velocity block of the linear system.
    /// @param[in] linPartOnly if true, returns only the linear part of the velocity-velocity block
    virtual gsSparseMatrix<T, MatOrder> getBlockOO(bool linPartOnly = false);

    /*
    /// @brief Returns the diagonal block of velocity-velocity block for i-th component.
    virtual gsSparseMatrix<T, MatOrder> getBlockUUcompDiag(index_t i = 0)
    { 
        GISMO_ASSERT(i >= 0 && i < m_tarDim, "Component index out of range.");
        return getBlockUU().block(i * m_udofs, i * m_udofs, m_udofs, m_udofs); 
    }


    /// @brief Returns the velocity-pressure block of the linear system.
    const gsSparseMatrix<T, MatOrder>& getBlockUP() const
    { return m_blockUP; }

    /// @brief Returns the pressure-velocity block of the linear system.
    gsSparseMatrix<T, MatOrder> getBlockPU() const
    { return (-1.0)*gsSparseMatrix<T, MatOrder>(m_blockUP.transpose()); }

    /// @brief Returns the part of velocity-pressure block for i-th velocity component.
    virtual gsSparseMatrix<T, MatOrder> getBlockUPcomp(index_t i) const
    { 
        GISMO_ASSERT(i >= 0 && i < m_tarDim, "Component index out of range.");
        return getBlockUP().middleRows(i * m_udofs, m_udofs);
    }

    /// @brief Returns part of pressure-velocity block for i-th velocity component.
    virtual gsSparseMatrix<T, MatOrder> getBlockPUcomp(index_t i) const
    { 
        GISMO_ASSERT(i >= 0 && i < m_tarDim, "Component index out of range.");
        return (-1.0)*gsSparseMatrix<T, MatOrder>(getBlockUPcomp(i).transpose());
    }

    /// @brief Returns the mass matrix for unknown with index \a unk.  There is also a const version.
    /// @param[in] unkID index of the unknown (0 - velocity, 1 - pressure)
    virtual gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID)
    { 
        GISMO_ASSERT(unkID == 0 || unkID == 1, "unkID must be 0 (velocity) or 1 (pressure).");
        GISMO_ASSERT(m_isMassMatReady, "Mass matrices not assembled in gsINSAssembler.");
        return m_massMatBlocks[unkID];
    }

    virtual const gsSparseMatrix<T, MatOrder>& getMassMatrix(index_t unkID) const
    { 
        GISMO_ASSERT(unkID == 0 || unkID == 1, "unkID must be 0 (velocity) or 1 (pressure).");
        GISMO_ASSERT(m_isMassMatReady, "Mass matrices not assembled in gsINSAssembler.");
        return m_massMatBlocks[unkID];
    }
    */

    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual gsMatrix<T> getRhsK() const;

    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual gsMatrix<T> getRhsO() const;

}; // gsTMAssemblerSST

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMAssemblerSST.hpp)
#endif