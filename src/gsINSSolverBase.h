/** @file gsINSSolverBase.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSSolverParams.h>

#include <gsCore/gsDebug.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief A base class for incompressible Navier-Stokes solvers.
/// @tparam T coefficient type
template<class T>
class gsINSSolverBase
{

protected: // *** Class members ***

    gsINSAssemblerBase<T>* m_pAssembler;
    gsMatrix<T> m_solution;
    unsigned m_iterationNumber;
    bool m_bPatternAnalyzed;
    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT, m_relNorm;

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver;
#else
    typename gsSparseSolver<T>::LU m_solver;
#endif


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverBase()
    {
        m_pAssembler = NULL;
    }

    virtual ~gsINSSolverBase()
    {
        if (m_pAssembler)
        {
            delete m_pAssembler;
            m_pAssembler = NULL;
        }
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    /// @brief Re-initialize all members.
    virtual void reinitMembers() { initMembers(); }

public: // *** Static functions ***

#ifdef GISMO_WITH_PARDISO
    /// @brief Setup the pardiso solver.
    /// @param[in,out] solver a reference to the pardiso solver
    static void pardisoSetup(typename gsSparseSolver<T>::PardisoLU& solver)
    {
        solver.setParam(7, 15);
        solver.setParam(9, 13);
        solver.setParam(12, 0);
    }
#endif

public: // *** Member functions ***

    /// @brief Initialize the solver.
    virtual void initialize() 
    { 
        if (!getAssembler()->isInitialized())
        {
            m_clock.restart();
            getAssembler()->initialize();
            m_assembT += m_clock.stop();
        }
    }

    /// @brief Prepare for the solution process.
    virtual void initIteration()
    {
        m_clock.restart();
        m_solver.analyzePattern(getAssembler()->matrix());
        m_solsetupT += m_clock.stop();
    }

    /// @brief Solve the linear system.
    /// @param[out] solution a reference to the vector, where the computed solution will be stored
    virtual void applySolver(gsMatrix<T>& solution);

    /// @brief Solve the linear system with underrelaxation.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    /// @param[in]  alpha_u     velocity relaxation parameter
    /// @param[in]  alpha_p     pressure relaxation parameter
    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p);

    /// @brief Perform next iteration step.
    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION }

    /// @brief Perform several iteration steps.
    /// @param[in] numberOfIterations the number of iterations to be performed
    virtual void nextIteration(const unsigned numberOfIterations)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        for (unsigned iter = 0; iter < numberOfIterations; iter++)
            nextIteration();
    }

    /// @brief Solve the incompressible Navier-Stokes problem.
    /// @param[in] maxIterations    the maximum number of linearization method iterations (in the steady case) or time steps (in the unsteady case)
    /// @param[in] epsilon          the stopping tolerance
    /// @param[in] minIterations    the minimum number of iterations/time steps
    void solve(const int maxIterations, const T epsilon = 1e-3, const int minIterations = 1);

    /// @brief Compute the Stokes problem.
    virtual void solveStokes();

    /// @brief Solve the generalized Stokes problem.
    virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Update the assembler with current solution.
    virtual void updateAssembler()
    {
        m_clock.restart();
        getAssembler()->update(m_solution);
        m_assembT += m_clock.stop();
    }

    /// @brief Compute and return the relative norm of the solution change.
    T solutionChangeRelNorm() const;

    /// @brief Compute and return the relative norm of the solution change given the two successive solutions.
    /// @param[in] solOld the old solution
    /// @param[in] solNew the new solution
    T solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const;

    /// @brief Compute and display the relative norm of the solution change given the two successive solutions.
    /// @param[in] solOld the old solution
    /// @param[in] solNew the new solution
    virtual void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const;

    /// @brief Compute and return the relative residual norm for the current solution.
    virtual T residualRelNorm() const 
    {
        return residualRelNorm(m_solution);
    }

    /// @brief Compute and return the relative residual norm for the given solution.
    virtual T residualRelNorm(const gsMatrix<T>& solution) const
    {
        gsMatrix<T> residual = getAssembler()->rhs() - getAssembler()->matrix() * solution;
        T resNorm = residual.norm() / getAssembler()->rhs().norm();

        return resNorm;
    }

    /// @brief Fix zero pressure on a given patch side.
    /// @param[in] patch  the given patch ID
    /// @param[in] side   the given patch side
    void addPressureOutletCondition(int patch, boxSide side)
    {
        getAssembler()->addPressureOutletCondition(patch, side);
        reinitMembers();
    }

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown (0 - velocity, 1 - pressure)
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        getAssembler()->markDofsAsEliminatedZeros(boundaryDofs, unk);
        reinitMembers();
    }

    /// @brief Construct solution field for the unknown \a unk for the current solution vector.
    /// @param unk the considered unknown (0 - velocity, 1 - pressure)
    /// @return 
    gsField<T> constructSolution(int unk) const
    {
        return getAssembler()->constructSolution(m_solution, unk);
    }


public: // *** Getters/setters ***

    /// @brief Set a given solution vector as current solution.
    /// @param[in] solVector the given solution
    virtual void setSolution(const gsMatrix<T> & solVector) { m_solution = solVector; }

    /// @brief Returns a pointer to the assembler.
    virtual gsINSAssemblerBase<T>* getAssembler() const { return m_pAssembler; }
    
    /// @brief Returns the current solution vector.
    const gsMatrix<T> & getSolution() const { return m_solution; }

    /// @brief Returns the current iteration number.
    unsigned getIterationNumber() const { return m_iterationNumber; }

    /// @brief Returns the total number of DOFs (the matrix size).
    int numDofs() const { return getAssembler()->numDofs(); }

    /// @brief Returns the total time spent on matrix assembly.
    virtual const T getAssemblyTime() const { return m_assembT; }

    /// @brief Returns the total time spent on linear solver setup.
    virtual const T getSolverSetupTime() const { return m_solsetupT; }

    /// @brief Returns the total time spent on solving of the linear systems.
    virtual const T getSolveTime() const { return m_solveT; }

}; //gsINSSolverBase


// ===========================================================================


/// @brief A base class for incompressible Navier-Stokes solvers with iterative solvers for the linear system.
/// @tparam T           coefficient type
/// @tparam SolverType  type of the linear solver
template<class T, class SolverType>
class gsINSSolverBaseIter
{

protected: // *** Class members ***

    std::map<std::string, gsSparseMatrix<> > m_matrices;
    typename gsINSPreconditioner<T>::Ptr m_pPrec;
    std::string m_precType;
    gsOptionList m_precOpt;
    int m_maxLinIter;
    T m_linTol;
    gsINSAssemblerBase<T>& m_assemblerRef;
    std::vector<index_t> m_linIterVector;
    gsSparseMatrix<T> m_matGammaPart;
    gsMatrix<T> m_rhsGammaPart;

    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT;

public: // *** Constructor/destructor ***

    // /// @brief Constructor.
    // gsINSSolverBaseIter()
    // {}

    /// @brief Constructor of the object.
    /// @param params[in]       list of optional parameters 
    /// @param assembler[in]    a reference to the assembler
    gsINSSolverBaseIter(gsINSSolverParams<T>& params, gsINSAssemblerBase<T>& assembler) :
        m_assemblerRef(assembler)
    {
        m_maxLinIter = params.options().getInt("maxIt_lin");
        m_linTol = params.options().getReal("tol_lin");
        m_precType = params.options().getString("precType");
        m_precOpt = params.precOptions();
        
        initMembers();
    }

    virtual ~gsINSSolverBaseIter()
    {
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    /// @brief Re-initialize all members.
    virtual void reinitMembers() { initMembers(); }

    /// @brief Initialize matrices needed for the preconditioner.
    void initPrecMat();

    /// @brief Update matrices needed for the preconditioner.
    void updatePrecMat();

public: // *** Member functions ***

    /// @brief Initialize the solver.
    /// @param[in] solution a const reference to the current solution vector
    virtual void initialize(const gsMatrix<T>& solution);

    /// @brief Prepare for the solution process.
    virtual void initIteration()
    {}

    /// @brief Solve the linear system.
    /// @param[out] solution a reference to the vector, where the computed solution will be stored
    virtual void applySolver(gsMatrix<T>& solution);

    /// @brief Solve the linear system with underrelaxation.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    /// @param[in]  alpha_u     velocity relaxation parameter
    /// @param[in]  alpha_p     pressure relaxation parameter
    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p);

    /// @brief Compute and return the solution of the Stokes problem.
    virtual gsMatrix<T> getStokesSolution();


public: // *** Getters/setters ***

    /// @brief Returns a vector of iteration counts of the linear solver during the solution process. 
    std::vector<index_t> getLinIterVector() const { return m_linIterVector; }

    /// @brief Returns the average number of iterations of the linear solver per Picard iteration.
    T getAvgLinIterations() const;

    /// @brief Returns the total time spent on matrix assembly.
    const T getAssemblyTime() const { return m_assembT; }

    /// @brief Returns the total time spent on linear solver setup.
    const T getSolverSetupTime() const { return m_solsetupT; }

    /// @brief Returns the total time spent on solving of the linear systems.
    const T getSolveTime() const { return m_solveT; }

}; //gsINSSolverBaseIter

} //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSSolverBase.hpp)
#endif
