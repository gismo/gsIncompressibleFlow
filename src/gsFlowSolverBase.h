/** @file gsFlowSolverBase.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for all flow solvers in gsIncompressibleFlow.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsFlowSolverBase
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr;
    gsFlowLinSystSolver<T, MatOrder>* m_linSolverPtr;

    gsMatrix<T> m_solution;
    unsigned m_iterationNumber;
    gsStopwatch m_clock;
    T m_initAssembT, m_assembT, m_relNorm;


public: // *** Constructor/destructor ***

    gsFlowSolverBase(const gsFlowSolverParams<T>& params):
    gsFlowSolverBase(memory::make_shared_not_owned(&params))
    { }

    gsFlowSolverBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
    {
        m_assemblerPtr = NULL;

        #ifdef _OPENMP
        omp_set_num_threads(m_paramsPtr->options().getInt("numThreads"));
        #endif
    }

    virtual ~gsFlowSolverBase()
    {
        if (m_assemblerPtr)
        {
            delete m_assemblerPtr;
            m_assemblerPtr = NULL;
        }

        if (m_linSolverPtr)
        {
            delete m_linSolverPtr;
            m_linSolverPtr = NULL;
        }
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    /// @brief Update sizes of members (when DOF numbers change after constructing the solver).
    virtual void updateSizes();

    /// @brief Print the options used by the solver to the logger.
    void printOptions()
    {
        std::stringstream sstr;
        sstr << "\n-----------------------------------\n";
        sstr << "Flow solver options:\n";
        sstr << m_paramsPtr->options() << "-----------------------------------\n";
        m_paramsPtr->logger().log(sstr.str(), true); // true = log to file only
    }


public: // *** Member functions ***

    /// @brief Initialize the solver.
    virtual void initialize();

    /// @brief Prepare for the solution process.
    virtual void initIteration()
    { getLinSolver()->setupSolver(getAssembler()->matrix()); }

    /// @brief Prepare for the solution process.
    /// @param[in] mat the linear system matrix
    virtual void initIteration(const gsSparseMatrix<T, MatOrder>& mat)
    { getLinSolver()->setupSolver(mat); }

    /// @brief Solve the linear system.
    /// @param[out] solution a reference to the vector, where the computed solution will be stored
    virtual void applySolver(gsMatrix<T>& solution)
    { getLinSolver()->applySolver(getAssembler()->matrix(), getAssembler()->rhs(), solution); }

    /// @brief Perform next iteration step.
    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION }

    /// @brief Perform several iteration steps.
    /// @param[in] numberOfIterations the number of iterations to be performed
    void nextIteration(const unsigned numberOfIterations);
    
    /// @brief Solve the incompressible Navier-Stokes problem.
    /// @param[in] maxIterations    the maximum number of linearization method iterations (in the steady case) or time steps (in the unsteady case)
    /// @param[in] epsilon          the stopping tolerance
    /// @param[in] minIterations    the minimum number of iterations/time steps
    void solve(const int maxIterations, const T epsilon = 1e-3, const int minIterations = 1);

    /// @brief Update the assembler with a given solution.
    /// @param[in] updateSol save the given solution into m_solution in the assembler
    virtual void updateAssembler(const gsMatrix<T>& sol, bool updateSol = true);

    /// @brief Update the assembler with current solution.
    /// @param[in] updateSol save current solution into m_solution in the assembler
    virtual void updateAssembler(bool updateSol = true)
    { updateAssembler(m_solution, updateSol); }

    /// @brief Compute and return the relative norm of the solution change.
    T solutionChangeRelNorm() const;

    /// @brief Compute and return the relative norm of the solution change given the two successive solutions.
    /// @param[in] solOld the old solution
    /// @param[in] solNew the new solution
    T solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const;

    /// @brief Compute and display the relative norm of the solution change given the two successive solutions.
    /// @param[in] solOld the old solution
    /// @param[in] solNew the new solution
    virtual void writeSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew, std::string solstr);

    /// @brief Compute and return the relative residual norm for the current solution.
    virtual T residualRelNorm() const 
    { return residualRelNorm(m_solution); }

    /// @brief Compute and return the relative residual norm for the given solution.
    virtual T residualRelNorm(const gsMatrix<T>& solution) const
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Construct solution field for the unknown \a unk for the current solution vector.
    /// @param[in] unk          the considered unknown (0 - velocity, 1 - pressure)
    /// @param[in] customSwitch a switch to be used for any purpose by derived classes
    virtual gsField<T> constructSolution(int unk, bool customSwitch = false) const
    { return getAssembler()->constructSolution(m_solution, unk, customSwitch); }

    /// @brief Check values of jacobian near boundaries of all patches. 
    /// @param npts[in] number of evaluation points along a patch side in each direction
    /// @param dist[in] distance from patch boundaries (in the parametric space)
    /// @param tol[in]  critical (positive) value of jacobian to throw warning (too close to zero)
    /// @return error code (0 = all jacobians OK, 1 = jacobian close to zero in some points, -1 = jacobian negative in some points)
    int checkGeoJacobian(int npts = -1, T dist = -1, T tol = -1);
    
public: // *** Getters/setters ***

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsFlowSolverBase"; }

    /// @brief Retrurns the solver parameters.
    virtual typename gsFlowSolverParams<T>::Ptr getParams() { return m_paramsPtr; }

    /// @brief Returns a pointer to the assembler.
    virtual gsFlowAssemblerBase<T, MatOrder>* getAssembler() const { return m_assemblerPtr; }

    /// @brief Returns a pointer to the linear system solver.
    gsFlowLinSystSolver<T, MatOrder>* getLinSolver() const { return m_linSolverPtr; }

    /// @brief Returns the current solution vector.
    const gsMatrix<T>& getSolution() const { return m_solution; }

    /// @brief Set a given solution vector as current solution.
    /// @param[in] solVector the given solution
    virtual void setSolution(const gsMatrix<T> & solVector) { m_solution = solVector; }

    /// @brief Returns the current iteration number.
    unsigned getIterationNumber() const { return m_iterationNumber; }

    /// @brief Returns the total number of DOFs (the matrix size).
    int numDofs() const { return getAssembler()->numDofs(); }

    /// @brief Returns the time of the initial matrix assembly.
    virtual const T getInitAssemblyTime() const { return m_initAssembT; }

    /// @brief Returns the total time spent on matrix assembly.
    virtual const T getAssemblyTime() const { return m_initAssembT + m_assembT; }

    /// @brief Returns the total time spent on linear solver setup.
    virtual const T getSolverSetupTime() const { return getLinSolver()->getSolverSetupTime(); }

    /// @brief Returns the total time spent on solving of the linear systems.
    virtual const T getSolveTime() const { return getLinSolver()->getSolveTime(); }


}; // gsFlowSolverBase

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowSolverBase.hpp)
#endif