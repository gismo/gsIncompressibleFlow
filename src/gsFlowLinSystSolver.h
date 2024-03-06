/** @file gsFlowLinSystSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

namespace gismo
{

/// @brief  Interface for classes solving linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// @tparam T real number type type
template<class T>
class gsFlowLinSystSolver
{

protected: // *** Class members ***

    const gsFlowSolverParams<T>& m_paramsRef;
    gsStopwatch m_clock;
    T m_setupT, m_solveT;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver(const gsFlowSolverParams<T>& params):
    m_paramsRef(params)
    {
        m_setupT = 0;
        m_solveT = 0;
    }

    virtual ~gsFlowLinSystSolver()
    { }


protected: // *** Member functions ***

    /// @brief Start measuring time (decides whether to use gsStopwatch or MPI_Wtime)
    T stopwatchStart();

    /// @brief Stop measuring time (decides whether to use gsStopwatch or MPI_Wtime)
    T stopwatchStop();


public: // *** Member functions ***

    /// @brief Setup the linear solver for a given matrix.
    virtual void setupSolver(const gsSparseMatrix<T>& mat)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the linear system.
    /// @param[in]  mat         a const reference to the system matrix
    /// @param[in]  rhs         a const reference to the system right-hand side
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the Navier--Stokes linear system with underrelaxation.
    /// @param[in]  mat         a const reference to the system matrix
    /// @param[in]  rhs         a const reference to the system right-hand side
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    /// @param[in]  alpha_u     velocity relaxation parameter
    /// @param[in]  alpha_p     pressure relaxation parameter
    /// @param[in]  usize       size of the velocity part of the system
    /// @param[in]  pdofs       number of pressure DOFs
    virtual void applySolver(const gsSparseMatrix<T>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p, index_t usize, index_t pdofs);


public: // *** Getters ***

    /// @brief Returns the total time spent on linear solver setup.
    virtual const T getSolverSetupTime() const { return m_setupT; }

    /// @brief Returns the total time spent on solving of the linear systems.
    virtual const T getSolveTime() const { return m_solveT; }


}; // gsFlowLinSystSolver

// ===================================================================================================================

/// @brief Direct solver for linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// @tparam T   coefficient type
template<class T>
class gsFlowLinSystSolver_direct: public gsFlowLinSystSolver<T>
{

public:
    typedef gsFlowLinSystSolver<T> Base;


protected: // *** Class members ***

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver;
#else
    typename gsSparseSolver<T>::LU m_solver;
#endif


protected: // *** Base class members ***

    using Base::m_paramsRef;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::stopwatchStart;
    using Base::stopwatchStop;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_direct(const gsFlowSolverParams<T>& params):
    Base(params)
    {
        #ifdef GISMO_WITH_PARDISO
        pardisoSetup(m_solver);
        #endif
    }


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

    /// @brief Setup the linear solver for a given matrix.
    virtual void setupSolver(const gsSparseMatrix<T>& mat);

    /// @brief Solve the linear system.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution);


}; // gsFlowLinSystSolver_direct

// ===================================================================================================================

// ===================================================================================================================

template<class T>
gsFlowLinSystSolver<T>* createLinSolver(const gsFlowSolverParams<T>& params)
{
    std::string type = params.options().getString("linSolver");

    if (type == "direct")
        return new gsFlowLinSystSolver_direct<T>(params);
    // else if (type == _"iter")
    //     return new gsFlowLinSystSolver_iter<T>(params);
    // else if (type == _"petsc")
    //     return new gsFlowLinSystSolver_PETSc<T>(params);
    else
    {
        gsInfo << "Invalid linear system solver type, using direct.\n";
        return new gsFlowLinSystSolver_direct<T>(params);
    }

}

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowLinSystSolver.hpp)
#endif