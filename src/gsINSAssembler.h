/** @file gsINSAssembler.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for incompressible Navier-Stokes assemblers.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsINSAssembler: public gsFlowAssemblerBase<T, MatOrder>
{

public:
    typedef gsFlowAssemblerBase<T, MatOrder> Base;


protected: // *** Class members ***

    T m_viscosity;
    gsSparseMatrix<T, MatOrder> m_baseMatrix, m_matrix;
    gsMatrix<T> m_baseRhs, m_rhs;

    index_t m_udofs;
    index_t m_pdofs;
    index_t m_pshift;
    index_t m_nnzPerRowU, m_nnzPerRowP;
    gsINSVisitorUUlin<T, MatOrder> m_visitorUUlin;
    gsINSVisitorUUnonlin<T, MatOrder> m_visitorUUnonlin;
    gsINSVisitorPU_withUPrhs<T, MatOrder> m_visitorUP;
    gsINSVisitorRhsU<T, MatOrder> m_visitorF;
    gsINSVisitorRhsP<T, MatOrder> m_visitorG;
    gsSparseMatrix<T, MatOrder> m_blockUUlin_comp, m_blockUUnonlin_comp, m_blockUP;
    gsSparseMatrix<T, MatOrder> m_blockUUlin_whole, m_blockUUnonlin_whole;
    gsMatrix<T> m_rhsUlin, m_rhsUnonlin, m_rhsBtB, m_rhsFG;
    gsField<T>  m_currentVelField, m_currentPresField;

    bool m_isMassMatReady;
    std::vector< gsSparseMatrix<T, MatOrder> > m_massMatBlocks;

    // PCD members
    // std::vector<index_t> m_presInIDs, m_presOutIDs, m_presWallIDs
    // std::vector< gsSparseMatrix<T, MatOrder> > m_pcdBlocks(3); // [laplaceP, robinBCblock, convectionP]


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_dofs;
    using Base::m_tarDim;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_isBaseReady;
    using Base::m_isSystemReady;

public: // *** Base class member functions ***

    using Base::getBases;
    using Base::getPatches;
    using Base::getAssemblerOptions;
    using Base::getBCs;

public: // *** Constructor/destructor ***

    gsINSAssembler(const gsFlowSolverParams<T>& params):
    Base(params)
    { }

    virtual ~gsINSAssembler()
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
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol);

    /// @brief Assemble the linear part of the problem.
    virtual void assembleLinearPart();

    /// @brief Assemble the linear part of the problem.
    virtual void assembleNonlinearPart();

    /// @brief Fill the velocity-velocity block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    velocity-velocity block (either for one velocity component or the whole block for all components)
    void fillGlobalMat_UU(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat);

    /// @brief Fill the velocity-pressure block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    velocity-pressure block
    void fillGlobalMat_UP(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat);

    /// @brief Fill the pressure-velocity block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    pressure-velocitye block
    void fillGlobalMat_PU(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat);

    /// @brief Fill the pressure-pressure block into the global saddle-point matrix.
    /// @param globalMat[out]   global saddle-point matrix
    /// @param sourceMat[in]    pressure-pressure block
    void fillGlobalMat_PP(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat);

    /// @brief Fill the linear part of the global matrix and right-hand side.
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
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk);

    /// @brief Fill the matrix and right-hand side for the Stokes problem.
    virtual void fillStokesSystem(gsSparseMatrix<T, MatOrder>& stokesMat, gsMatrix<T>& stokesRhs);

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, index_t unk) const;

    /// @brief Compute flow rate through a side of a given patch.
    /// @param[in] patch        the given patch ID
    /// @param[in] side         the given patch side
    /// @param[in] solution     solution vector to compute the flow rate from
    T computeFlowRate(index_t patch, boxSide side, gsMatrix<T> solution) const;


public: // *** Getters/setters ***

    /// @brief Returns the viscosity value.
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

    /// @brief Returns the number of DOFs for the i-th unknown.
    /// @param[in] i    0 - velocity, 1 - pressure
    index_t numDofsUnk(index_t i);

    /// @brief  Returns the number of velocity DOFs (one velocity component).
    index_t getUdofs() const { return m_udofs; }

    /// @brief Returns the number of pressure DOFs.
    index_t getPdofs() const { return m_pdofs; }

    /// @brief Returns the DOF shift of pressure (i.e. the total number of velocity DOFs).
    index_t getPshift() const { return m_pshift; }

    /// @brief Returns the velocity-velocity block of the linear system.
    /// @param[in] linPartOnly if true, returns only the linear part of the velocity-velocity block
    virtual gsSparseMatrix<T, MatOrder> getBlockUU(bool linPartOnly = false);

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

    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual gsMatrix<T> getRhsU() const;

    /// @brief /// @brief Returns part of the right-hand side for i-th velocity component.
    virtual gsMatrix<T> getRhsUcomp(index_t i) const
    { 
        GISMO_ASSERT(i >= 0 && i < m_tarDim, "Component index out of range.");
        return getRhsU().middleRows(i * m_udofs, m_udofs);
    }

    /// @brief Returns the pressure part of the right-hand side.
    gsMatrix<T> getRhsP() const;

}; // gsINSAssembler


// ===================================================================================================================


/// @brief  The steady incompressible Navier--Stokes assembler.
/// @tparam T real number type
template<class T, int MatOrder>
class gsINSAssemblerSteady: public gsINSAssembler<T, MatOrder>
{

public:
    typedef gsINSAssembler<T, MatOrder> Base;


public: // *** Constructor/destructor ***

    gsINSAssemblerSteady(const gsFlowSolverParams<T>& params): 
    Base(params)
    {
        Base::initMembers();
    }

    virtual ~gsINSAssemblerSteady()
    { }

}; // gsINSAssemblerSteady


// ===================================================================================================================


/// @brief  The unsteady incompressible Navier--Stokes assembler.
/// @tparam T real number type
template<class T, int MatOrder>
class gsINSAssemblerUnsteady: public gsINSAssembler<T, MatOrder>
{

public:
    typedef gsINSAssembler<T, MatOrder> Base;


protected: // *** Class members ***

    gsINSVisitorUUtimeDiscr<T, MatOrder> m_visitorTimeDiscr;
    gsSparseMatrix<T, MatOrder> m_blockTimeDiscr;
    gsMatrix<T> m_rhsTimeDiscr;
    gsField<T> m_oldTimeVelField;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_pshift;
    using Base::m_nnzPerRowU;
    using Base::m_dofMappers;
    using Base::m_solution;
    using Base::m_baseMatrix;
    using Base::m_rhs;
    using Base::m_currentVelField;


public: // *** Constructor/destructor ***

    gsINSAssemblerUnsteady(const gsFlowSolverParams<T>& params): 
    Base(params)
    {
        initMembers();
    }

    virtual ~gsINSAssemblerUnsteady()
    { }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change, e.g. after markDofsAsEliminatedZeros()).
    virtual void updateSizes();

    /// @brief Update the DOF mappers in all visitors (e.g. after markDofsAsEliminatedZeros() ).
    virtual void updateDofMappers();

    /// @brief Update the current solution field stored in the assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner nonlinear iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol);

    /// @brief Assemble the linear part of the matrix.
    virtual void assembleLinearPart();

    /// @brief Fill the linear part of the global matrix and right-hand side.
    virtual void fillBaseSystem();

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();


public: // *** Member functions ***

    /// @brief 
    /// @param solVector 
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);


public: // *** Getters/setters ***

    /// @brief Returns the velocity-velocity block of the linear system.
    virtual gsSparseMatrix<T, MatOrder> getBlockUU()
    { return Base::getBlockUU() + m_blockTimeDiscr; }


    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual gsMatrix<T> getRhsU() const
    { return Base::getRhsU() + m_rhsTimeDiscr; }


}; //gsINSAssemblerUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSAssembler.hpp)
#endif