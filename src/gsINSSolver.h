/** @file gsINSSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSSolverParams.h>

#include <gsCore/gsDebug.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief A base class for incompressible Navier-Stokes solvers.
/// @tparam T coefficient type
template<class T>
class gsINSSolver
{

protected: // *** Class members ***

    gsINSSolverParams<T> m_params;
    gsINSAssembler<T>* m_pAssembler;
    gsMatrix<T> m_solution;
    unsigned m_iterationNumber;
    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT, m_relNorm;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolver(gsINSSolverParams<T>& params): m_params(params)
    { }

    virtual ~gsINSSolver()
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

    /// @brief 
    void nextIteration_steady();


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
    virtual void initIteration() {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the linear system.
    /// @param[out] solution a reference to the vector, where the computed solution will be stored
    virtual void applySolver(gsMatrix<T>& solution) {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the linear system with underrelaxation.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    /// @param[in]  alpha_u     velocity relaxation parameter
    /// @param[in]  alpha_p     pressure relaxation parameter
    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
    {GISMO_NO_IMPLEMENTATION}

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
    virtual void solveStokes() { GISMO_NO_IMPLEMENTATION }

    /// @brief Solve the generalized Stokes problem.
    virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Update the assembler with current solution.
    /// @param updateSol 
    virtual void updateAssembler(bool updateSol = true)
    { getAssembler()->update(m_solution, updateSol); }

    /// @brief Update the assembler with a given solution.
    /// @param updateSol 
    virtual void updateAssembler(const gsMatrix<T>& sol, bool updateSol = true)
    { getAssembler()->update(sol, updateSol); }

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
    { GISMO_NO_IMPLEMENTATION }


    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk)
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
    gsINSAssembler<T>* getAssembler() const { return m_pAssembler; }
    
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

}; //gsINSSolver

// ===================================================================================================================

/// @brief 
/// @tparam T coefficient type
template<class T>
class gsINSSolverDirect : public gsINSSolver<T>
{

public:
    typedef gsINSSolver<T> Base;


protected: // *** Class members ***

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver;
#else
    typename gsSparseSolver<T>::LU m_solver;
#endif


protected: // *** Base class members / functions ***

    using Base::m_clock;
    using Base::m_assembT;
    using Base::m_solsetupT;
    using Base::m_solveT;
    using Base::getAssembler;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverDirect(gsINSSolverParams<T>& params): Base(params)
    { }


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

}; //gsINSSolverDirect

// ===================================================================================================================

/// @brief 
/// @tparam T coefficient type
template<class T>
class gsINSSolverDirectSteady : public gsINSSolverDirect<T>
{

public:
    typedef gsINSSolverDirect<T> Base;


protected: // *** Base class members ***

    using gsINSSolver<T>::m_pAssembler;
    using gsINSSolver<T>::m_solution;
    using gsINSSolver<T>::m_iterationNumber;
    using gsINSSolver<T>::m_params;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverDirectSteady(gsINSSolverParams<T>& params): Base(params)
    { 
        m_pAssembler = new gsINSAssembler<T>(params);

        Base::initMembers();
        m_params.options().setSwitch("unsteady", false);
    }


public: // *** Member functions ***

     /// @brief Perform next iteration step.
    virtual void nextIteration()
    { gsINSSolver<T>::nextIteration_steady(); }


    // /// @brief Compute the Stokes problem.
    // virtual void solveStokes();

    
    // /// @brief Solve the generalized Stokes problem.
    // virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    // { GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    // gsINSAssembler<T>* getAssembler() const
    // {
    //     return dynamic_cast<gsINSAssembler<T>*>(m_pAssembler);
    // }


}; //gsINSSolverDirectSteady

// ===================================================================================================================

/// @brief 
/// @tparam T coefficient type
template<class T>
class gsINSSolverDirectUnsteady : public gsINSSolverDirect<T>
{

public:
    typedef gsINSSolverDirect<T> Base;


protected: // *** Class members ***

    T m_time, m_timeStepSize;
    T m_innerIter, m_avgPicardIter;
    T m_innerTol;

protected: // *** Base class members ***

    using Base::m_solution;
    using Base::m_iterationNumber;
    using gsINSSolver<T>::m_pAssembler;
    using gsINSSolver<T>::m_params;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverDirectUnsteady(gsINSSolverParams<T>& params): Base(params)
    { 
        m_pAssembler = new gsINSAssemblerUnsteady1<T>(params);

        initMembers();
        m_params.options().setSwitch("unsteady", true);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();


public: // *** Member functions ***

     /// @brief Perform next iteration step.
    virtual void nextIteration();
    

    // /// @brief Compute the Stokes problem.
    // virtual void solveStokes();

    
    // /// @brief Solve the generalized Stokes problem.
    // virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    // { GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    gsINSAssemblerUnsteady1<T>* getAssembler() const
    {
        return dynamic_cast<gsINSAssemblerUnsteady1<T>*>(m_pAssembler);
    }

    // @brief Returns the elapsed simulation time.
    T getSimulationTime() const { return m_time; }

    /// @brief Returns the average number of Picard iterations per time step.
    T getAvgPicardIterations() const { return m_avgPicardIter / m_iterationNumber; }


}; //gsINSSolverDirectUnsteady

} //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSSolver.hpp)
#endif