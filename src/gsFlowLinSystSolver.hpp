/** @file gsFlowLinSystSolver.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>

namespace gismo
{

template<class T, int MatOrder>
T gsFlowLinSystSolver<T, MatOrder>::stopwatchStart()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsPtr->options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        m_clock.restart();

    return 0.0;
}


template<class T, int MatOrder>
T gsFlowLinSystSolver<T, MatOrder>::stopwatchStop()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsPtr->options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        return m_clock.stop();

}


template<class T, int MatOrder>
void gsFlowLinSystSolver<T, MatOrder>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p, index_t usize, index_t pdofs)
{
    GISMO_ASSERT(solution.rows() > 0, "The solution in applySolver() is empty!");

    gsMatrix<T> newSol;
    this->applySolver(mat, rhs, newSol);

    real_t time0 = stopwatchStart();
    solution.topRows(usize) = alpha_u * newSol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
    solution.bottomRows(pdofs) = alpha_p * newSol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
    real_t time1 = stopwatchStop();

    m_solveT += time1 - time0;
}


// ===================================================================================================================
// ===================================================================================================================


template<class T, int MatOrder>
void gsFlowLinSystSolver_direct<T, MatOrder>::setupSolver(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    m_solver.analyzePattern(mat);
    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}


template<class T, int MatOrder>
void gsFlowLinSystSolver_direct<T, MatOrder>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    real_t time0 = stopwatchStart();
    m_solver.factorize(mat);

    real_t time1 = stopwatchStop();
    
    solution = m_solver.solve(rhs);
    real_t time2 = stopwatchStop();

    m_setupT += time1 - time0;
    m_solveT += time2 - time1;
}

// ===================================================================================================================

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iter<T, MatOrder, SolverType>::setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    m_precPtr = gsIncompleteLUOp< gsSparseMatrix<T, MatOrder> >::make(mat);
    //m_precPtr = gsGaussSeidelOp< gsSparseMatrix<T, MatOrder>, gsGaussSeidel::symmetric >::make(mat);
    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iter<T, MatOrder, SolverType>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    real_t time0 = stopwatchStart();

    this->setupPreconditioner(mat);
    
    SolverType solver(mat, m_precPtr);
    solver.setMaxIterations(m_paramsPtr->options().getInt("lin.maxIt"));
    solver.setTolerance(m_paramsPtr->options().getReal("lin.tol"));

    real_t time1 = stopwatchStop();
    
    solver.solve(rhs, solution);

    real_t time2 = stopwatchStop();

    m_setupT += time1 - time0;
    m_solveT += time2 - time1;

    m_linIterVector.push_back(solver.iterations());
}

// ===================================================================================================================

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iterSP<T, MatOrder, SolverType>::setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    
    m_matrices.clear();
    m_matrices.insert(std::make_pair("matNS", mat));
    m_matrices.insert(std::make_pair("matMu", m_assemblerPtr->getMassMatrix(0)));
    m_matrices.insert(std::make_pair("matMp", m_assemblerPtr->getMassMatrix(1)));

    if (m_precType.substr(0, 3) == "PCD")
    {
        gsWarn << "PCD preconditioner not fully implemented yet, using LSC instead.";
        m_precType.replace(0, 3, "LSC");

        // gsSparseMatrix<T> Ap, Fp;
        // m_assemblerPtr->fillPCDblocks(Ap, Fp, m_precOpt.getInt("pcd_bcType"), m_precOpt.getSwitch("pcd_assembAp"), m_precOpt.getSwitch("pcd_assembFp"), m_precOpt.getSwitch("lumpingM"));
        // m_matrices.insert(std::make_pair("matFp", Fp));
        // m_matrices.insert(std::make_pair("matAp", Ap));
    }

    m_precPtr = gsINSPreconditioner<T, MatOrder>::make(m_precType, m_matrices, m_precOpt);

    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}


}