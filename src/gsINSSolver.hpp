/** @file gsINSSolver.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSSolver.h>

namespace gismo
{

template<class T>
void gsINSSolver<T>::initMembers()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
    m_iterationNumber = 0;
    m_relNorm = std::numeric_limits<T>::infinity();
}


template<class T>
void gsINSSolver<T>::nextIteration_steady()
{
    GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    updateAssembler();

    if (!m_iterationNumber)
        initIteration();

    applySolver(m_solution);

    m_iterationNumber++;
}


template<class T>
void gsINSSolver<T>::solve(const int maxIterations, const T epsilon, const int minIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
    int iter = 0;
    m_relNorm = solutionChangeRelNorm();

    while ((iter < minIterations) || ((m_relNorm > epsilon) && (iter < maxIterations)))
    {
        gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

        nextIteration();
        m_relNorm = solutionChangeRelNorm();
        gsInfo << " Solution change relative norm: " << m_relNorm << "\n";

        iter++;
    }
}


template<class T>
T gsINSSolver<T>::solutionChangeRelNorm() const
{
    T relNorm;

    if (m_iterationNumber)
    {
        gsMatrix<T> solChangeVector = getAssembler()->getSolution() - m_solution;
        relNorm = solChangeVector.norm() / m_solution.norm();

    }
    else
    {
        relNorm = std::numeric_limits<T>::infinity();
    }

    return relNorm;
}


template<class T>
T gsINSSolver<T>::solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsMatrix<T> solChangeVector = solOld - solNew;
    T relNorm = solChangeVector.norm() / solNew.norm();

    return relNorm;
}


template<class T>
void gsINSSolver<T>::dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsInfo << "     [u, p] solution change relative norm: ";

    for (int i = 0; i < solOld.cols(); i++)
        gsInfo << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

    gsInfo << "\n";
}

// ===================================================================================================================

template<class T>
void gsINSSolverDirect<T>::initMembers()
{
    Base::initMembers();

    #ifdef GISMO_WITH_PARDISO
        pardisoSetup(m_solver);
    #endif
}

template<class T>
void gsINSSolverDirect<T>::applySolver(gsMatrix<T>& solution)
{
    m_clock.restart();
    m_solver.factorize(getAssembler()->matrix());
    m_solsetupT += m_clock.stop();

    m_clock.restart();
    solution = m_solver.solve(getAssembler()->rhs());
    m_solveT += m_clock.stop();
}


template<class T>
void gsINSSolverDirect<T>::applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
{
    GISMO_ASSERT(solution.rows() > 0,"The solution in applySolver() is empty!");

    int usize = getAssembler()->getPshift();
    int pdofs = getAssembler()->getPdofs();

    m_clock.restart();
    m_solver.factorize(getAssembler()->matrix());
    m_solsetupT += m_clock.stop();

    m_clock.restart();
    gsMatrix<T> newsol = m_solver.solve(getAssembler()->rhs());
    solution.topRows(usize) = alpha_u * newsol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
    solution.bottomRows(pdofs) = alpha_p * newsol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
    m_solveT += m_clock.stop();
}

// ===================================================================================================================

template<class T>
void gsINSSolverDirectUnsteady<T>::initMembers()
{
    Base::initMembers();
    m_time = 0;
    m_timeStepSize = m_params.options().getReal("timeStep");
    m_innerIter = m_params.options().getInt("maxIt_picard");
    m_innerTol = m_params.options().getReal("tol_picard");
    m_avgPicardIter = 0;
}


template<class T>
void gsINSSolverDirectUnsteady<T>::nextIteration()
{
    GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    this->updateAssembler();

    if (!m_iterationNumber)
        this->initIteration();

    gsMatrix<T> tmpSolution = m_solution;

    this->applySolver(tmpSolution);

    this->dispSolChangeRelNorm(m_solution, tmpSolution);

    index_t picardIter = 0;
    T relNorm = this->solutionChangeRelNorm(m_solution, tmpSolution);

    gsInfo << "        [u, p] Picard's iterations...\n";
    while((relNorm > m_innerTol) && (picardIter < m_innerIter))
    {
        gsInfo << "         ";
        gsMatrix<T> oldSol = tmpSolution;

        this->updateAssembler(tmpSolution, false);

        this->applySolver(tmpSolution);

        this->dispSolChangeRelNorm(oldSol, tmpSolution);

        relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
        picardIter++;
    }
    
    m_solution = tmpSolution;

    m_time += m_timeStepSize;
    m_avgPicardIter += picardIter;
    m_iterationNumber++;
}

} //namespace gismo
