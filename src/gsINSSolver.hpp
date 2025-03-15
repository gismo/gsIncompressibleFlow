/** @file gsINSSolver.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSSolver.h>

namespace gismo
{

template<class T, int MatOrder>
void gsINSSolver<T, MatOrder>::solveStokes()
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
    gsWriteOutputLine(m_outFile, "Computing the steady Stokes problem...", m_fileOutput, m_dispOutput);

    gsSparseMatrix<T, MatOrder> stokesMat;
    gsMatrix<T> stokesRhs;

    getAssembler()->fillStokesSystem(stokesMat, stokesRhs);

    if (!m_iterationNumber)
        this->initIteration(stokesMat);

    this->getLinSolver()->applySolver(stokesMat, stokesRhs, m_solution);
}

// ===================================================================================================================

template<class T, int MatOrder>
void gsINSSolverSteady<T, MatOrder>::nextIteration()
{
    GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    this->updateAssembler();

    if (!m_iterationNumber)
        this->initIteration();

    this->applySolver(m_solution);

    m_iterationNumber++;
}

// ===================================================================================================================

template<class T, int MatOrder>
void gsINSSolverUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();
    m_time = 0;
    m_timeStepSize = m_paramsPtr->options().getReal("timeStep");
    m_innerIter = m_paramsPtr->options().getInt("nonlin.maxIt");
    m_innerTol = m_paramsPtr->options().getReal("nonlin.tol");
    m_avgPicardIter = 0;
}


template<class T, int MatOrder>
void gsINSSolverUnsteady<T, MatOrder>::plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::string fileNameSuffix, unsigned plotPts)
{
    int numPatches = m_paramsPtr->getPde().patches().nPatches();

    gsField<T> uSol = this->constructSolution(0);
    std::stringstream filenameU;
    filenameU << "velocity" + fileNameSuffix + "_" << m_iterationNumber << "it";
    gsWriteParaview<T>(uSol, filenameU.str(), plotPts);

    gsField<T> pSol = this->constructSolution(1);
    std::stringstream filenameP;
    filenameP << "pressure" + fileNameSuffix + "_" << m_iterationNumber << "it";
    gsWriteParaview<T>(pSol, filenameP.str(), plotPts);

    for (int p = 0; p < numPatches; p++)
    {
        std::stringstream fnU;
        fnU << filenameU.str() << p << ".vts";
        fileU << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnU.str() << "\"/>\n";

        std::stringstream fnP;
        fnP << filenameP.str() << p << ".vts";
        fileP << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnP.str() << "\"/>\n";
    }
}


template<class T, int MatOrder>
void gsINSSolverUnsteady<T, MatOrder>::nextIteration()
{
    GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    this->updateAssembler();

    if (!m_iterationNumber)
        this->initIteration();

    gsMatrix<T> tmpSolution = m_solution;

    this->applySolver(tmpSolution);

    this->writeSolChangeRelNorm(m_solution, tmpSolution, "[u, p]");

    index_t picardIter = 0;
    T relNorm = this->solutionChangeRelNorm(m_solution, tmpSolution);

    gsWriteOutputLine(m_outFile, "        [u, p] Picard's iterations...", m_fileOutput, m_dispOutput);

    while((relNorm > m_innerTol) && (picardIter < m_innerIter))
    {
        gsWriteOutput(m_outFile, "         ", m_fileOutput, m_dispOutput);

        gsMatrix<T> oldSol = tmpSolution;

        this->updateAssembler(tmpSolution, false);
        this->applySolver(tmpSolution);
        this->writeSolChangeRelNorm(oldSol, tmpSolution, "[u, p]");

        relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
        picardIter++;
    }
    
    m_solution = tmpSolution;

    m_time += m_timeStepSize;
    m_avgPicardIter += picardIter;
    m_iterationNumber++;
}


template<class T, int MatOrder>
void gsINSSolverUnsteady<T, MatOrder>::solveWithAnimation(const int totalIter, const int iterStep, std::string fileNameSuffix, const T epsilon, unsigned plotPts, const int minIterations)
{
    // prepare plotting
    std::string fileNameU = "velocity" + fileNameSuffix + "_animation.pvd";
    std::ofstream fileU(fileNameU.c_str());
    GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

    std::string fileNameP = "pressure" + fileNameSuffix + "_animation.pvd";
    std::ofstream fileP(fileNameP.c_str());
    GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

    startAnimationFile(fileU);
    startAnimationFile(fileP);

    plotCurrentTimeStep(fileU, fileP, fileNameSuffix, plotPts);

    for (int i = 0; i < totalIter; i += iterStep)
    {
        this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

        plotCurrentTimeStep(fileU, fileP, fileNameSuffix, plotPts);
    }

    endAnimationFile(fileU);
    endAnimationFile(fileP);
}

// from old version of the solver (TODO):

// template<class T>
// void gsINSSolverUnsteady<T>::initGeneralizedStokesSolution(gsSparseMatrix<T>& stokesMatrix, gsMatrix<T>& stokesRhs)
// {
//     GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

//     getAssembler()->fillStokesSystem(stokesMatrix, stokesRhs);

//     const int uDofs = getAssembler()->getUdofs();
//     const T invTimeStep = 1. / m_timeStepSize;

//     #pragma omp parallel for num_threads(getAssembler()->getBlockAssembler().getNumThreads())
//     for (index_t col = 0; col < uDofs; ++col)
//         for (typename gsSparseMatrix<T>::InnerIterator it(getAssembler()->getVelocityMassMatrix(), col); it; ++it)
//             for (index_t s = 0; s < getAssembler()->getTarDim(); ++s)
//                 stokesMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += invTimeStep * it.value();

//     m_clock.restart();
//     m_solver.analyzePattern(stokesMatrix);
//     m_solver.factorize(stokesMatrix);
//     m_solsetupT += m_clock.stop();
// }


// template<class T>
// void gsINSSolverUnsteady<T>::solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations)
// {
//     gsSparseMatrix<T> stokesMatrix;
//     gsMatrix<T> stokesRhs;

//     initGeneralizedStokesSolution(stokesMatrix, stokesRhs);

//     const int uDofs = getAssembler()->getUdofs();
//     const T invTimeStep = 1. / m_timeStepSize;
//     int iter = 0;
//     T relNorm = std::numeric_limits<T>::infinity();

//     while ((iter < minIterations) || ((relNorm > epsilon) && (iter < maxIterations)))
//     {
//         gsInfo << "Iteration number " << iter + 1 << "...";

//         gsMatrix<T> rhs = stokesRhs;
//         gsMatrix<T> newSol = m_solution;

//         for (index_t s = 0; s < getAssembler()->getTarDim(); ++s)
//             rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * getAssembler()->getVelocityMassMatrix() * m_solution.middleRows(s * uDofs, uDofs);

//         m_clock.restart();
//         newSol = m_solver.solve(rhs);
//         m_solveT += m_clock.stop();

//         relNorm = this->solutionChangeRelNorm(m_solution, newSol);
//         gsInfo << " Solution change relative norm: " << relNorm << "\n";

//         m_solution = newSol;
//         iter++;
//     }
// }

} //namespace gismo
