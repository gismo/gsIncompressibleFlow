/** @file gsINSBlockAssembler.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include<gsCore/gsField.h>
#include <gsIncompressibleFlow/src/gsINSBlockAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSBlockVisitorsBnd.h>

#ifdef _OPENMP 
#include <omp.h>
#endif

namespace gismo
{

/// @brief Assembler of matrix blocks needed for the incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template<class T>
class gsINSBlockAssembler : public gsINSBlockAssemblerBase<T>
{

public:
    typedef gsINSBlockAssemblerBase<T> Base;

protected: // *** Class members ***

    int m_udofs;
    int m_pdofs;
    int m_pshift;
    int m_nonZerosPerColU, m_nonZerosPerColP;
    gsSparseMatrix<T> m_blockA;
    std::vector< gsSparseMatrix<T> > m_blockB;
    std::vector< gsSparseMatrix<T> > m_blockMinusBT;
    std::vector< gsSparseMatrix<T> > m_blockCT;
    gsSparseMatrix<T> m_blockN;
    gsSparseMatrix<T> m_blockNpattern;
    gsSparseMatrix<T> m_blockM;
    gsSparseMatrix<T> m_blockMinv;
    gsSparseMatrix<T> m_blockAp;
    gsSparseMatrix<T> m_blockNp;
    gsSparseMatrix<T> m_blockMp;
    //gsSparseMatrix<T> m_blockPdg;
    gsSparseMatrix<T> m_blockRobin;
    gsMatrix<T> m_rhsA;
    gsMatrix<T> m_rhsB;
    gsMatrix<T> m_rhsC;
    gsMatrix<T> m_rhsN;
    gsMatrix<T> m_rhsM;
    gsMatrix<T> m_rhsAp;
    gsMatrix<T> m_rhsMp;
    //gsMatrix<T> m_rhsPdg;
    gsMatrix<T> m_rhsF;
    gsField<T>  m_currentVelField, m_currentPresField, m_oldTimeVelField;

    // PCD assembly
    std::vector<index_t> m_presInIDs, m_presOutIDs, m_presWallIDs;
    bool m_bPCDbndPrepared;

protected: // *** Base class members ***

    using Base::m_dofs;
    using Base::m_numThreads;
    using Base::m_tarDim;
    using Base::m_bUnsteady;
    using Base::m_params;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_viscosity;

public: // *** Constructor/destructor ***

    /// @brief  Constructor of the object.
    /// @param params a set of all needed parameters
    gsINSBlockAssembler(const gsINSSolverParams<T>& params) :
        Base(params)
    {
        initMembers();
    }

    virtual ~gsINSBlockAssembler()
    {}

protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

public: // *** Member functions ***

    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[out] result      the resulting solution as a gsMultiPatch object
    /// @param[in]  unk         the considered unknown
    void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result, int unk) const;
    
    /// @brief Construct solution from computed solution vector for unknown \a unk.
    /// @param[in]  solVector   the solution vector obtained from the linear system
    /// @param[in]  unk         the considered unknown
    /// @return     the resulting solution field
    gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk) const;

    /// @brief Update the current solution field stored in the block assembler (used as convection coefficient).
    /// @param[in] solVector    new solution vector obtained from the linear system
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol);

    /// @brief Compute flow rate through a side of a given patch.
    /// @param[in] patch        the given patch ID
    /// @param[in] side         the given patch side
    /// @param[in] solution     solution vector to compute the flow rate from
    T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const;

    /// @brief Assemble the linear part of the Stokes linear system.
    void assembleLinearStokesPart();

    /// @brief Assemble the nonlinear part of the Navier-Stokes linear system.
    void assembleNonlinearPart();

    /// @brief Assemble velocity mass matrix.
    void assembleMassMatrix();
    
    /// @brief Assemble pressure mass matrix.
    void assemblePressureMassMatrix();
    
    //void assembleBlocksC();

    /** @brief Assemble block N (multiplied by 0) to get its sparsity pattern.
     * 
     *  Since block N arises from the nonlinear convective term, it is reassembled
     *  in each step of the linearization method. The block assembled here is used
     *  for exact memory allocation in the global matrix.
     * 
     */
    void assembleBlockNpattern();

    /// @brief Assemble all blocks that do not change during the linearization steps.
    void assembleAllLinearBlocks()
    {
        assembleLinearStokesPart();
        assemblePressurePoisson();
    }

    /// @brief Assemble the Stokes matrix from individual blocks assembled before.
    /// @param[out] stokesMatrix    a reference to a sparse matrix to be filled
    /// @param[out] stokesRhs       a reference to a RHS vector to be filled
    void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const;

    /**
     * @brief Assemble the Poisson matrix for the pressure basis.
     * 
     * This is needed for the PCD block preconditioner.
     */ 
    void assemblePressurePoisson();

    /**
     * @brief Assemble the convection matrix for the pressure basis.
     * 
     * This is needed for the PCD block preconditioner.
     */ 
    void assemblePressureConvection();

    /**
     * @brief Assemble block arising from Robin boundary conditions for a Poisson problem on the pressure space.
     * 
     * This is needed for the PCD block preconditioner.
     */ 
    void assembleApRobinBlock(std::vector<std::pair<int, boxSide> > bndPart);

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
    void preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType = 0);

    /**
     * @brief Returns the Poisson matrix for the pressure basis.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * The matrix \f$ A_p \f$ is either assembled or defined as \f$ B \hat{M}_u^{-1} B^T \ \f$, where \f$ \hat{M}_u_\f$ is an
     * approximation of the velocity mass matrix. 
     * 
     * The function assemblePressurePoisson() is called here, if was not called before.
     * 
     *  @param[in] assemb   assemble the matrix
     *  @param[in] lumping  use mass lumping (in the case assemb = false)
     */ 
    gsSparseMatrix<T> getPressurePoissonMatrix(bool assemb = true, bool lumping = false);

    /**
     * @brief Assembles and returns the convection matrix for the pressure basis.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     */ 
    gsSparseMatrix<T> getPressureConvectionMatrix()
    {
        assemblePressureConvection();

        return m_blockNp;
    }

    /**
     * @brief Returns the pressure mass matrix.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * The function assemblePressureMassMatrix() is called here, if was not called before.
     * 
     */ 
    gsSparseMatrix<T> getPressureMassMatrix()
    {
        if (!m_blockMp.nonZeros())
            assemblePressureMassMatrix();

        return m_blockMp;
    }


    // bcType references:
    // 0 - Kay, Loghin (original PCD)
    // 1 - Elman, Tuminaro variant 1
    // 2 - Elman, Tuminaro variant 2
    // 3 - Blechta Y-variant

    /**
     * @brief Apply boundary modifications to the pressure Poisson and convection matrices.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * Admissible values of \a bcType and their meaning:
     * 0 - Dirichlet inflow + Neumann outflow and walls for both \f$ Ap, Fp \f$
     * 1 - Robin inflow + Neumann outflow and walls for \f$ Fp \f$, nothing for \f$ Ap\f$ (\f$ Ap \f$ should not be assembled, i.e., \f$ B \hat{M}_u^{-1} B^T \ \f$)
     * 2 - Robin inflow + Dirichlet outflow + Neumann walls for \f$ Fp \f$, nothing for \f$ Ap \f$ (\f$ Ap\f$ should not be assembled)
     * 3 - Robin inflow + Dirichlet outflow + Neumann walls for \f$ Fp \f$, Dirichlet outflow + Neumann inflow and walls for \f$ Ap \f$ (\f$ Ap \f$ should be assembled)
     * 4 - no boundary modification
     * 
     * @param[in,out]   Ap      a reference to the pressure Poisson matrix to be modified
     * @param[in,out]   Fp      a reference to the pressure convection matrix to be modified
     * @param[in]       bcType  boundary conditions strategy
     */ 
    void applyPCDboundaryConditions(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType);

    //T computeDimensionlessWallDistance(gsMatrix<T> solution, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const;

protected: // *** Member functions ***

    /// @brief Set up the block visitor according to its type.
    /// @param[in,out] blockVisitor a pointer to the visitor to be set up
    virtual void blockSpecificSettings(gsINSVisitorBase<T>* blockVisitor)
    {
        std::vector<gsField<T> > sols;
        sols.push_back(m_currentVelField);
        blockVisitor->setCurrentSolution(sols);
    }

    /// @brief Set up the right-hand-side visitor according to its type.
    /// @param[in,out] rhsVisitor a pointer to the visitor to be set up
    virtual void rhsSpecificSettings(gsINSVisitorBase<T>* rhsVisitor)
    {
        gsINSRhsVisitor<T>* pVisitor = dynamic_cast<gsINSRhsVisitor<T>*>(rhsVisitor);

        if (pVisitor != NULL)
            pVisitor->setRhsFunction(getRhsFcn());
    }

    /**
     * @brief Identify pressure boundary DOF IDs.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[in] bndIn     a list of pairs (patch, side) corresponding to inflow boundary
     * @param[in] bndOut    a list of pairs (patch, side) corresponding to outflow boundary
     * @param[in] bndWall   a list of pairs (patch, side) corresponding to solid wall
     */ 
    void findPressureBoundaryIDs(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall)
    {
        findPressureBoundaryPartIDs(bndIn, m_presInIDs);
        findPressureBoundaryPartIDs(bndOut, m_presOutIDs);
        findPressureBoundaryPartIDs(bndWall, m_presWallIDs);
    }

    /**
     * @brief Identify pressure boundary DOF IDs for a given boundary part.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[in]   bndPart     a list of pairs (patch, side) of the required boundary part (inflow, outflow, wall)
     * @param[out]  idVector    a vector in which the found IDs will be stored
     */ 
    void findPressureBoundaryPartIDs(std::vector<std::pair<int, boxSide> > bndPart, std::vector<index_t>& idVector);

    /**
     * @brief Apply the Dirichlet boundary modification to the given block.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[in,out]   block           a reference to the block to be modified
     * @param[in]       bndPartIDs      the DOF IDs of the Dirichlet boundary part
     * @param[in]       scaleDiag       use scaling for the diagonal matrix entries corresponding to the Dirichlet boundary part (otherwise equal to 1)
     * @param[in]       exclInflowPart  exclude inflow boundary IDs from calculation of the scaling parameter
     */ 
    void applyPCDboundaryDirichletPart(gsSparseMatrix<T>& block, const std::vector<index_t>& bndPartIDs, bool scaleDiag = true, bool exclInflowPart = false);

    /**
     * @brief Apply the Robin boundary modification to the given block.
     * 
     * This is needed for the PCD block preconditioner.
     * 
     * @param[in,out]   block           a reference to the block to be modified
     */ 
    void applyPCDboundaryRobinPart(gsSparseMatrix<T>& block)
    {
        GISMO_ASSERT(m_bPCDbndPrepared, "Pressure boundary DOF IDs not prepared!");

        block += m_blockRobin;        
    }

    /// @brief Check if the given block is filled.
    void blockCheck(const gsSparseMatrix<T> & block) const
    {
        GISMO_ASSERT(block.nonZeros(), "The required block is empty!");
    }


public: // *** Getters/setters ***

    /// @brief  Returns the number of velocity DOFs (one velocity component).
    int getUdofs() const { return m_udofs; }

    /// @brief Returns the number of pressure DOFs.
    int getPdofs() const { return m_pdofs; }

    /// @brief Returns the DOF shift of pressure (i.e. the total number of velocity DOFs).
    int getPshift() const { return m_pshift; }

    T   getTimeStep() const { return m_params.options().getReal("timeStep"); }

    const gsSparseMatrix<T> & getBlockA() const { return m_blockA; }
    const std::vector< gsSparseMatrix<T> > & getBlockB() const { return m_blockB; }
    const std::vector< gsSparseMatrix<T> > & getBlockMinusBT() const { return m_blockMinusBT; }
    const gsSparseMatrix<T> & getBlockB(const int i) const { return m_blockB[i]; }
    const gsSparseMatrix<T> & getBlockMinusBT(const int i) const { return m_blockMinusBT[i]; }
    const gsSparseMatrix<T> & getBlockCT(const int i) const { return m_blockCT[i]; }
    const gsSparseMatrix<T> & getBlockN() const { return m_blockN; }
    const gsSparseMatrix<T> & getBlockNpattern() const { return m_blockNpattern; }
    const gsSparseMatrix<T> & getBlockM() const { return m_blockM; }
    const gsSparseMatrix<T> & getBlockAp() const { blockCheck(m_blockAp); return m_blockAp; }
    const gsSparseMatrix<T> & getBlockNp() const { blockCheck(m_blockNp); return m_blockNp; }
    const gsSparseMatrix<T> & getBlockMp() const { blockCheck(m_blockMp); return m_blockMp; }

    const gsMatrix<T> & getRhsA() const { return m_rhsA; }
    const gsMatrix<T> & getRhsB() const { return m_rhsB; }
    const gsMatrix<T> & getRhsC() const { return m_rhsC; }
    const gsMatrix<T> & getRhsN() const { return m_rhsN; }
    const gsMatrix<T> & getRhsM() const { return m_rhsM; }
    const gsMatrix<T> & getRhsAp() const { return m_rhsAp; }
    const gsMatrix<T> & getRhsMp() const { return m_rhsMp; }
    const gsMatrix<T> & getRhsF() const { return m_rhsF; }


public: // *** Base member functions ***

    using Base::isUnsteady;
    using Base::getPatches;
    using Base::getBases;
    using Base::getRhsFcn;
    using Base::getAssemblerOptions;

protected: // *** Base member functions ***

    using Base::computeDirichletDofs;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSBlockAssembler.hpp)
#endif