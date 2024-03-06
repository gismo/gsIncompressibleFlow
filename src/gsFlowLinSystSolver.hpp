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

template<class T>
T gsFlowLinSystSolver<T>::stopwatchStart()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsRef.options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        m_clock.restart();

    return 0.0;
}


template<class T>
T gsFlowLinSystSolver<T>::stopwatchStop()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsRef.options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        return m_clock.stop();

}


template<class T>
void gsFlowLinSystSolver<T>::applySolver(const gsSparseMatrix<T>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p, index_t usize, index_t pdofs)
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


template<class T>
void gsFlowLinSystSolver_direct<T>::setupSolver(const gsSparseMatrix<T>& mat)
{
    real_t time0 = stopwatchStart();
    m_solver.analyzePattern(mat);
    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}


template<class T>
void gsFlowLinSystSolver_direct<T>::applySolver(const gsSparseMatrix<T>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    real_t time0 = stopwatchStart();
    m_solver.factorize(mat);

    real_t time1 = stopwatchStop();
    
    solution = m_solver.solve(rhs);
    real_t time2 = stopwatchStop();

    m_setupT += time1 - time0;
    m_solveT += time2 - time1;
}


}