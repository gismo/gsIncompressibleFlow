/** @file gsINSAssembler.h
    
    @brief 
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova (Hornikova)
 */

#pragma once

#include<gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>

namespace gismo
{

template<class T>
class gsINSAssemblerBase
{

protected: // *** Class members ***

    gsINSSolverParams<T> m_params;
    int m_dofs;
    int m_tarDim;
    T m_viscosity;
    bool m_isInitialized;
    bool m_isSystemReady;
    std::vector<gsDofMapper> m_dofMappers;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution;
    gsSparseMatrix<T, RowMajor> m_baseMatrix, m_matrix;
    gsMatrix<T> m_baseRhs, m_rhs;

public: // *** Constructor/destructor ***

    gsINSAssemblerBase(const gsINSSolverParams<T>& params):
    m_params(params)
    {
        initMembers();
    }

    virtual ~gsINSAssemblerBase()
    { }

protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers();

    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  basisID     the index of the basis corresponding to \a unk (same as \a unk in this case)
    /// @param[in]  ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofs(const int unk, const int basisID, gsMatrix<T>& ddofVector);


    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using interpolation.
    /// @param[in]  unk          the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper       reference to the DOF mapper for \a unk
    /// @param[in]  mbasis       reference to the basis corresponding to \a unk
    /// @param[out] ddofVector   reference to the vector where computed coefficients will be stored
    void computeDirichletDofsIntpl(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);


    /// @brief Compute the coefficients of the basis functions at the Dirichlet boundaries using L2-projection.
    /// @param[in]  unk         the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in]  mapper      reference to the DOF mapper for \a unk
    /// @param[in]  mbasis      reference to the basis corresponding to \a unk
    /// @param[out] ddofVector  reference to the vector where computed coefficients will be stored
    void computeDirichletDofsL2Proj(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector);


    /// @brief 
    /// @param visitor 
    /// @param basisID 
    /// @param block 
    /// @param blockRhs 
    void assembleBlock(gsINSVisitor<T>& visitor, index_t testBasisID, gsSparseMatrix<T, RowMajor>& block, gsMatrix<T>& blockRhs);


    /// @brief 
    /// @param visitor 
    /// @param testBasisID 
    /// @param blockRhs 
    void assembleRhs(gsINSVisitor<T>& visitor, index_t testBasisID, gsMatrix<T>& rhs);

    /// @brief Assemble the linear part of the matrix.
    virtual void assembleLinearPart() {GISMO_NO_IMPLEMENTATION}


    /// @brief Assemble the nonlinear part of the matrix.
    virtual void assembleNonlinearPart() {GISMO_NO_IMPLEMENTATION}


    /// @brief Assemble the nonlinear part of the problem.
    virtual void updateAssembly();


    /// @brief Update the current solution field stored in the block assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector obtained from the linear system
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
    {GISMO_NO_IMPLEMENTATION}


    /// @brief Fill the linear part of the global matrix and right-hand side.
    virtual void fillBaseSystem() {GISMO_NO_IMPLEMENTATION}


    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem() {GISMO_NO_IMPLEMENTATION}


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
    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk);


    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk) const
    {GISMO_NO_IMPLEMENTATION}
    

public: // *** Getters/setters ***

    /// @brief Returns the number of degrees of freedom (DOFs).
    int numDofs() const 
    { 
        GISMO_ASSERT(m_dofs > 0, "Something went wrong, number of DOFs is zero!");
        return m_dofs; 
    }

    /// @brief Returns the target dimension.
    int getTarDim() const { return m_tarDim; }

    /// @brief Returns the viscosity value.
    T getViscosity() const { return m_viscosity; }

    bool isInitialized() { return m_isInitialized; }

    /// @brief Returns the assembled matrix.
    const gsSparseMatrix<T, RowMajor>& matrix() const
    {
        GISMO_ASSERT(m_isSystemReady, "Matrix not ready, update() must be called first.");
        return m_matrix;
    }

    /// @brief Returns the assembled right-hand side.
    const gsMatrix<T>& rhs() const
    {
        GISMO_ASSERT(m_isSystemReady, "Rhs not ready, update() must be called first.");
        return m_rhs;
    }

    /**
     * @brief Returns a const reference to the DOF mappers.
     *
     * The mapper for velocity is stored first, the mapper for pressure is second.
     */
    const std::vector<gsDofMapper>& getMappers() const { return m_dofMappers; }

    /**
     * @brief Returns a const reference to the vectors of coefficients at the Dirichlet boundaries.
     *
     * The vector of velocity coefficients is stored first, the vector of pressure coefficients is second.
     */
    const std::vector<gsMatrix<T> >& getDirichletDofs() const { return m_ddof; }

    /// @brief Returns a const reference to the current computed solution.
    const gsMatrix<T>& getSolution() const { return m_solution; }

    /// @brief Returns a const reference to the multipatch representing the computational domain.
    const gsMultiPatch<T>& getPatches() const { return m_params.getPde().patches(); }
    
    /**
     * @brief Returns a reference to the discretization bases.
     *
     * The velocity basis is stored first, the  pressure basis is second.
     * 
     * There is also a const version returning a const reference.
     */
    std::vector< gsMultiBasis<T> >& getBases() { return m_params.getBases(); }
    const std::vector< gsMultiBasis<T> >& getBases() const { return m_params.getBases(); }
    
    /// @brief Returns a const reference to the boundary conditions.
    const gsBoundaryConditions<T>& getBCs() const { return m_params.getBCs(); }
    
    /// @brief Returns a pointer to the right-hand-side function.
    const gsFunction<T>* getRhsFcn() const { return m_params.getPde().rhs(); }

    /// @brief Returns the assembler options.
    gsAssemblerOptions  getAssemblerOptions() const { return m_params.assemblerOptions(); }
    
    /// @brief Returns the INS solver option list.
    gsOptionList    options() const { return m_params.options(); }

    /// @brief Check if the solved problem is unsteady.
    bool isUnsteady() const { return m_params.options().getSwitch("unsteady"); }

};

// ===================================================================================================================

template<class T>
class gsINSAssembler: public gsINSAssemblerBase<T>
{

public:
    typedef gsINSAssemblerBase<T> Base;


protected: // *** Class members ***

    int m_udofs;
    int m_pdofs;
    int m_pshift;
    int m_nnzPerRowU, m_nnzPerRowP;
    gsINSVisitorUUlin<T> m_visitorUUlin;
    gsINSVisitorUUnonlin<T> m_visitorUUnonlin;
    gsINSVisitorPU<T> m_visitorUP;
    gsINSVisitorRhsU<T> m_visitorF;
    gsINSVisitorRhsP<T> m_visitorG;
    gsSparseMatrix<T, RowMajor> m_blockUUlin, m_blockUUnonlin, m_blockUP, m_blockPU;
    gsMatrix<T> m_rhsUlin, m_rhsUnonlin, m_rhsBtB, m_rhsFG;
    gsField<T>  m_currentVelField, m_currentPresField, m_oldTimeVelField;

protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_dofs;
    using Base::m_tarDim;
    using Base::m_viscosity;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_baseMatrix;
    using Base::m_matrix;
    using Base::m_baseRhs;
    using Base::m_rhs;
    using Base::m_isSystemReady;


public: // *** Constructor/destructor ***

    gsINSAssembler(const gsINSSolverParams<T>& params): 
    Base(params)
    {
        initMembers();
    }

    virtual ~gsINSAssembler()
    { }

protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();


    /// @brief Assemble the linear part of the matrix.
    virtual void assembleLinearPart();


    /// @brief Assemble the linear part of the matrix.
    virtual void assembleNonlinearPart();


    /// @brief Update the current solution field stored in the block assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector obtained from the linear system
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol);


    /// @brief Fill the linear part of the global matrix and right-hand side.
    virtual void fillBaseSystem();


    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();


public: // *** Member functions ***


    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    virtual gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk) const;


    /// @brief Compute flow rate through a side of a given patch.
    /// @param[in] patch        the given patch ID
    /// @param[in] side         the given patch side
    /// @param[in] solution     solution vector to compute the flow rate from
    T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const;

public: // *** Getters/setters ***

    /// @brief  Returns the number of velocity DOFs (one velocity component).
    int getUdofs() const { return m_udofs; }

    /// @brief Returns the number of pressure DOFs.
    int getPdofs() const { return m_pdofs; }

    /// @brief Returns the DOF shift of pressure (i.e. the total number of velocity DOFs).
    int getPshift() const { return m_pshift; }

    /// @brief Returns the velocity-velocity block of the linear system.
    virtual const gsSparseMatrix<T, RowMajor> getBlockUU() const
    { return m_blockUUlin + m_blockUUnonlin; }

    /// @brief Returns the velocity-pressure block of the linear system.
    const gsSparseMatrix<T, RowMajor>& getBlockUP() const
    { return m_blockUP; }

    /// @brief Returns the pressure-velocity block of the linear system.
    const gsSparseMatrix<T, RowMajor> getBlockPU() const
    { return gsSparseMatrix<T, RowMajor>(m_blockUP.transpose()); }

    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual const gsMatrix<T> getRhsU() const
    { 
        gsMatrix<T> rhsUpart = (m_rhsFG + m_rhsBtB).topRows(m_pshift);

        return (rhsUpart + m_rhsUlin + m_rhsUnonlin);
    }

    /// @brief Returns the pressure part of the right-hand side.
    const gsMatrix<T> getRhsP() const
    { return (m_rhsFG + m_rhsBtB).bottomRows(m_pdofs); }

};

// ===================================================================================================================

template<class T>
class gsINSAssemblerUnsteady: public gsINSAssembler<T>
{

public:
    typedef gsINSAssembler<T> Base;


protected: // *** Class members ***

    gsINSVisitorUUtimeDiscr<T> m_visitorTimeDiscr;
    gsSparseMatrix<T, RowMajor> m_blockTimeDiscr;
    gsMatrix<T> m_rhsTimeDiscr;

    // int m_udofs;
    // int m_pdofs;
    // int m_pshift;
    // int m_nnzPerRowU, m_nnzPerRowP;
    // gsINSVisitorUUlin<T> m_visitorUUlin;
    // gsINSVisitorUUnonlin<T> m_visitorUUnonlin;
    // gsINSVisitorPU<T> m_visitorUP;
    // gsINSVisitorRhsU<T> m_visitorF;
    // gsINSVisitorRhsP<T> m_visitorG;
    // gsSparseMatrix<T, RowMajor> m_blockUUlin, m_blockUUnonlin, m_blockUP, m_blockPU;
    // gsMatrix<T> m_rhsUlin, m_rhsUnonlin, m_rhsBtB, m_rhsFG;
    // gsField<T>  m_currentVelField, m_currentPresField, m_oldTimeVelField;

protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_pshift;
    using Base::m_nnzPerRowU;
    using gsINSAssemblerBase<T>::m_solution;
    using gsINSAssemblerBase<T>::m_baseMatrix;
    using gsINSAssemblerBase<T>::m_rhs;
    using gsINSAssemblerBase<T>::m_isInitialized;


public: // *** Constructor/destructor ***

    gsINSAssemblerUnsteady(const gsINSSolverParams<T>& params): 
    Base(params)
    {
        initMembers();
    }

    virtual ~gsINSAssemblerUnsteady()
    { }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();


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
    virtual const gsSparseMatrix<T, RowMajor> getBlockUU() const
    { return Base::getBlockUU() + m_blockTimeDiscr; }


    /// @brief /// @brief Returns the velocity part of the right-hand side.
    virtual const gsMatrix<T> getRhsU() const
    { return Base::getRhsU() + m_rhsTimeDiscr; }

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSAssembler.hpp)
#endif