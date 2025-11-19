/** @file gsFlowLinSystSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsSolver/gsGMRes.h>

namespace gismo
{

/// @brief              Interface for classes solving linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// @tparam T           real number type type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsFlowLinSystSolver
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    gsStopwatch m_clock;
    T m_setupT, m_solveT;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
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
    virtual void setupSolver(const gsSparseMatrix<T, MatOrder>& mat)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Setup the linear solver for a matrix given by its blocks.
    /// @param[in] matBlocks vector of saddle-point matrix blocks
    virtual void setupSolver(const std::vector< gsSparseMatrix<T, MatOrder> >& matBlocks)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the linear system.
    /// @param[in]  mat         a const reference to the system matrix
    /// @param[in]  rhs         a const reference to the system right-hand side
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the linear system.
    /// @param[in]  matBlocks vector of saddle-point matrix blocks
    /// @param[in]  rhsBlocks vector of right-hand side blocks
    /// @param[out] solution  a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const std::vector< gsSparseMatrix<T, MatOrder> >& matBlocks, const std::vector< gsMatrix<T> >& rhsBlocks, gsMatrix<T>& solution)
    {GISMO_NO_IMPLEMENTATION}

    /// @brief Solve the Navier--Stokes linear system with underrelaxation.
    /// @param[in]  mat         a const reference to the system matrix
    /// @param[in]  rhs         a const reference to the system right-hand side
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    /// @param[in]  alpha_u     velocity relaxation parameter
    /// @param[in]  alpha_p     pressure relaxation parameter
    /// @param[in]  usize       size of the velocity part of the system
    /// @param[in]  pdofs       number of pressure DOFs
    virtual void applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p, index_t usize, index_t pdofs);

    /// @brief Prints the linear iteration counts per call of \a applySolver() (only for iterative solvers, nothing happens here).
    virtual void reportLinIterations() { }


public: // *** Getters ***

    /// @brief Returns the total time spent on linear solver setup.
    virtual const T getSolverSetupTime() const { return m_setupT; }

    /// @brief Returns the total time spent on solving of the linear systems.
    virtual const T getSolveTime() const { return m_solveT; }


}; // gsFlowLinSystSolver

// ===================================================================================================================
// ===================================================================================================================

/// @brief Direct solver for linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsFlowLinSystSolver_direct: public gsFlowLinSystSolver<T, MatOrder>
{

public:
    typedef gsFlowLinSystSolver<T, MatOrder> Base;


protected: // *** Class members ***

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver;
#else
    typename gsSparseSolver<T>::LU m_solver;
#endif


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::stopwatchStart;
    using Base::stopwatchStop;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_direct(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    {
        #ifdef GISMO_WITH_PARDISO
        pardisoSetup<T>(m_solver);
        #endif
    }


public: // *** Member functions ***

    /// @brief Setup the linear solver for a given matrix.
    virtual void setupSolver(const gsSparseMatrix<T, MatOrder>& mat);

    /// @brief Solve the linear system.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution);


}; // gsFlowLinSystSolver_direct

// ===================================================================================================================

/// @brief G+Smo/Eigen iterative solver for linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @tparam SolverType  the G+Smp/Eigen iterative solver type
/// @ingroup IncompressibleFlow
template <class T, int MatOrder, class SolverType>
class gsFlowLinSystSolver_iter: public gsFlowLinSystSolver<T, MatOrder>
{

public:
    typedef gsFlowLinSystSolver<T, MatOrder> Base;


protected: // *** Class members ***

    typename gsLinearOperator<T>::Ptr m_precPtr;
    std::vector<index_t> m_linIterVector;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::stopwatchStart;
    using Base::stopwatchStop;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_iter(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


public: // *** Member functions ***

    /// @brief Setup the linear solver for a given matrix (nothing happens here).
    virtual void setupSolver(const gsSparseMatrix<T, MatOrder>& mat)
    { }

    /// @brief Setup the preconditioner for a given matrix.
    virtual void setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat);

    /// @brief Solve the linear system.
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution);

    /// @brief Prints the linear iteration counts per call of \a applySolver().
    virtual void reportLinIterations()
    {
        m_paramsPtr->logger() << "Iterations of linear solver for each call of applySolver():\n";
        for (size_t i = 0; i < m_linIterVector.size(); i++)
            m_paramsPtr->logger() << m_linIterVector[i] << ", ";

        m_paramsPtr->logger() << "\nAverage number of linear solver iterations per call of applySolver(): " << getAvgLinIterations() << "\n";
    }


public: // *** Getters/setters ***

    /// @brief Returns vector of iteration counts of the linear solver for each call of applySolver().
    std::vector<index_t> getLinIterVector() const { return m_linIterVector; }

    /// @brief Returns the average iteration count of the linear solver per applySolver() call.
    T getAvgLinIterations() const
    {
        index_t linIterSum = 0;

        for (size_t i = 0; i < m_linIterVector.size(); i++)
            linIterSum += m_linIterVector[i];

        return (T)linIterSum / m_linIterVector.size();
    }

}; // gsFlowLinSystSolver_iter

// ===================================================================================================================

/// @brief G+Smo/Eigen iterative solver for saddle-point linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase) with block preconditioners.
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @tparam SolverType  the G+Smp/Eigen iterative solver type
/// @ingroup IncompressibleFlow
template <class T, int MatOrder, class SolverType>
class gsFlowLinSystSolver_iterSP: public gsFlowLinSystSolver_iter<T, MatOrder, SolverType>
{

public:
    typedef gsFlowLinSystSolver_iter<T, MatOrder, SolverType> Base;


protected: // *** Class members ***

    std::string m_precType;
    gsOptionList m_precOpt;
    const gsINSAssembler<T, MatOrder>* m_assemblerPtr;
    std::map<std::string, gsSparseMatrix<T, MatOrder> > m_matrices;


protected: // *** Base class members ***

    using Base::m_precPtr;
    using Base::m_linIterVector;
    using Base::m_paramsPtr;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::stopwatchStart;
    using Base::stopwatchStop;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_iterSP(typename gsFlowSolverParams<T>::Ptr paramsPtr, const gsINSAssembler<T, MatOrder>* assemblerPtr):
    Base(paramsPtr), m_assemblerPtr(assemblerPtr)
    {
        m_precType = m_paramsPtr->options().getString("lin.precType");
        m_precOpt = m_paramsPtr->precOptions();
        m_precOpt.addInt("dim", "Problem dimension", m_assemblerPtr->getTarDim());
        m_precOpt.addReal("visc", "Viscosity", m_paramsPtr->getPde().viscosity());
        m_precOpt.addInt("udofs", "Number of velocity dofs", m_assemblerPtr->getUdofs());
        m_precOpt.addInt("pdofs", "Number of pressure dofs", m_assemblerPtr->getPdofs());
    }


public: // *** Member functions ***

    /// @brief Setup the preconditioner for a given matrix.
    virtual void setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat);

}; // gsFlowLinSystSolver_iterSP

// ===================================================================================================================

#ifdef gsPetsc_ENABLED

/// @brief PETSc solver for general linear systems with no specified block structure inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
/// Note: only MatOrder = RowMajor is supported
/// @tparam T           coefficient type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowLinSystSolver_PETSc: public gsFlowLinSystSolver<T, RowMajor>
{

public:
    typedef gsFlowLinSystSolver<T, RowMajor> Base;


protected: // *** Class members ***

    gsOptionList m_petscOpt;
    std::vector<index_t> m_linIterVector;
    std::pair<index_t, index_t> m_rowLocInfo, m_colLocInfo;

    Mat m_petscMat;
    Vec m_petscRhs, m_petscSol;
    KSP m_ksp;
    PC m_pc;
    real_t m_dataCopyT;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::stopwatchStart;
    using Base::stopwatchStop;
    using Base::setupSolver;
    using Base::applySolver;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_PETSc(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr), m_dataCopyT(0.0)
    { }

    ~gsFlowLinSystSolver_PETSc()
    {
        PetscCallVoid( MatDestroy(&m_petscMat) );
        PetscCallVoid( VecDestroy(&m_petscRhs) );
        PetscCallVoid( VecDestroy(&m_petscSol) );
        PetscCallVoid( KSPDestroy(&m_ksp) );
    }


public: // *** Member functions ***

    /// @brief Returns default PETSc options for the solver.
    virtual gsOptionList getDefaultOptions();

    /// @brief Setup the linear solver for a given matrix.
    virtual void setupSolver(const gsSparseMatrix<T, RowMajor>& mat);

    /// @brief Solve the linear system.
    /// @param[in]  mat         a const reference to the system matrix
    /// @param[in]  rhs         a const reference to the system right-hand side
    /// @param[out] solution    a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const gsSparseMatrix<T, RowMajor>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution);


    /// @brief Prints the linear iteration counts per call of \a applySolver().
    virtual void reportLinIterations()
    {
        m_paramsPtr->logger() << "Iterations of PETSc solver for each call of applySolver():\n";
        for (size_t i = 0; i < m_linIterVector.size(); i++)
            m_paramsPtr->logger() << m_linIterVector[i] << ", ";

        m_paramsPtr->logger() << "\nAverage number of linear solver iterations per call of applySolver(): " << getAvgLinIterations() << "\n";
    }


protected: // *** Member functions ***

    /// @brief  Apply options for PETSc.
    void applyOptions();

    /// @brief  Apply given options for PETSc.
    /// @param petscOpt option list
    void applyOptions(gsOptionList petscOpt);


public: // *** Getters/setters ***

    /// @brief Returns vector of iteration counts of the linear solver for each call of applySolver().
    std::vector<index_t> getLinIterVector() const { return m_linIterVector; }

    /// @brief Returns the average iteration count of the linear solver per applySolver() call.
    T getAvgLinIterations() const
    {
        index_t linIterSum = 0;

        for (size_t i = 0; i < m_linIterVector.size(); i++)
            linIterSum += m_linIterVector[i];

        return (T)linIterSum / m_linIterVector.size();
    }

    /// @brief Returns the total time spent on copying data to and from PETSc.
    virtual const T getDataCopyTime() const { return m_dataCopyT; }

}; // gsFlowLinSystSolver_PETSc


/**
 * @brief PETSc solver for saddle-point linear systems inside the incompressible flow solvers (classes derived from gsFlowSolverBase).
 * 
 * Assuming saddle-point linear systems from the incompressible flow problems, i.e., assuming the velocity blocks to have \a d blocks 
 * corresponding to velocity components, where \a d is the spatial dimension.
 * 
 * Note: only MatOrder = RowMajor is supported
 * 
 * @tparam T           coefficient type
 * @ingroup IncompressibleFlow
 */
template <class T>
class gsFlowLinSystSolver_PETSc_SP: public gsFlowLinSystSolver_PETSc<T>
{

public:
    typedef gsFlowLinSystSolver_PETSc<T> Base;


protected: // *** Class members ***

    std::pair<index_t, index_t> m_uLocInfo, m_pLocInfo;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_setupT;
    using Base::m_solveT;
    using Base::m_dataCopyT;
    using Base::m_linIterVector;
    using Base::m_petscMat;
    using Base::m_petscRhs;
    using Base::m_petscSol;
    using Base::m_ksp;
    using Base::m_pc;
    using Base::stopwatchStart;
    using Base::stopwatchStop;
    // using Base::setupSolver;
    // using Base::applySolver;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsFlowLinSystSolver_PETSc_SP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    {
        m_paramsPtr->options().setSwitch("fillGlobalSyst", false);
    }


protected: // *** Member functions ***

    /// @brief  Apply options for PETSc.
    void applyOptions();
    

public: // *** Member functions ***

    /// @brief Returns default PETSc options for the solver.
    virtual gsOptionList getDefaultOptions();

    /// @brief Setup the linear solver for a matrix given by its blocks.
    /// @param[in] matBlocks vector of saddle-point matrix blocks
    virtual void setupSolver(const std::vector< gsSparseMatrix<T, RowMajor> >& matBlocks);

    /// @brief Solve the linear system.
    /// @param[in]  matBlocks vector of saddle-point matrix blocks
    /// @param[in]  rhsBlocks vector of right-hand side blocks
    /// @param[out] solution  a reference to the vector, where the computed solution will be stored
    virtual void applySolver(const std::vector< gsSparseMatrix<T, RowMajor> >& matBlocks, const std::vector< gsMatrix<T> >& rhsBlocks, gsMatrix<T>& solution);
    
}; // gsFlowLinSystSolver_PETSc_SP

#endif

// ===================================================================================================================
// ===================================================================================================================

template <class T, int MatOrder, class SolverType = gsGMRes<T> >
gsFlowLinSystSolver<T, MatOrder>* createLinSolver(typename gsFlowSolverParams<T>::Ptr paramsPtr, const gsFlowAssemblerBase<T, MatOrder>* assemblerPtr = NULL)
{
    const gsINSAssembler<T, MatOrder>* INSassembPtr = dynamic_cast< const gsINSAssembler<T, MatOrder>* >(assemblerPtr);

    std::string type = paramsPtr->options().getString("lin.solver");

    if (type == "direct")
        return new gsFlowLinSystSolver_direct<T, MatOrder>(paramsPtr);
    else if (type == "iter" && INSassembPtr == NULL) // assembler is not an INS assembler
        return new gsFlowLinSystSolver_iter<T, MatOrder, SolverType>(paramsPtr);
    else if (type == "iter" && INSassembPtr != NULL) // assembler is an INS assembler
        return new gsFlowLinSystSolver_iterSP<T, MatOrder, SolverType>(paramsPtr, INSassembPtr);
    else if(type == "petsc")
    {
        if constexpr (MatOrder == RowMajor)
        {
            if (INSassembPtr == NULL) // assembler is not an INS assembler
                return new gsFlowLinSystSolver_PETSc<T>(paramsPtr);
            else // assembler is an INS assembler
                return new gsFlowLinSystSolver_PETSc_SP<T>(paramsPtr);
        }
        else
            GISMO_ERROR("Trying to use PETSc with ColMajor matrix ordering - currently not supported.");
    } 
    else
    {
        paramsPtr->logger() << "Invalid linear system solver type, using direct solver.\n";
        return new gsFlowLinSystSolver_direct<T, MatOrder>(paramsPtr);
    }
}

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowLinSystSolver.hpp)
#endif