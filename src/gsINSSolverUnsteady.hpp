/** @file gsINSSolverUnsteady.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
#include <gsIO/gsWriteParaview.h>

namespace gismo
{

template<class T>
void gsINSSolverUnsteady<T>::plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, unsigned plotPts)
{
    int numPatches = getAssembler()->getBlockAssembler().getPatches().nPatches();

    gsField<T> uSol = this->constructSolution(0);
    std::stringstream filenameU;
    filenameU << "velocity_" << m_iterationNumber << "it";
    gsWriteParaview<T>(uSol, filenameU.str(), plotPts);

    gsField<T> pSol = this->constructSolution(1);
    std::stringstream filenameP;
    filenameP << "pressure_" << m_iterationNumber << "it";
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


template<class T>
void gsINSSolverUnsteady<T>::nextIteration()
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    this->updateAssembler();

    // in the first time step the pattern will get analyzed
    if (!m_iterationNumber)
        this->initIteration();

    gsMatrix<T> tmpSolution = m_solution;

    this->applySolver(tmpSolution);

    Base::dispSolChangeRelNorm(m_solution, tmpSolution);

    int iter = 0;
    T relNorm = this->solutionChangeRelNorm(m_solution, tmpSolution);

    gsInfo << "        [u, p] Picard's iterations...\n";
    while((relNorm > m_innerTol) && (iter < m_innerIter))
    {
        gsInfo << "         ";
        gsMatrix<T> oldSol = tmpSolution;

        this->updateAssemblerPicard(tmpSolution);

        this->applySolver(tmpSolution);

        dispSolChangeRelNorm(oldSol, tmpSolution);

        relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
        iter++;
    }

    m_avgPicardIter += iter;
    
    m_solution = tmpSolution;

    m_iterationNumber++;
    m_time += getAssembler()->getTimeStepSize();
}


template<class T>
void gsINSSolverUnsteady<T>::solveWithAnimation(const int totalIter, const int iterStep, const T epsilon, unsigned plotPts, const int minIterations)
{
    // prepare plotting
    std::string fileNameU = "velocity_animation.pvd";
    std::ofstream fileU(fileNameU.c_str());
    GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

    std::string fileNameP = "pressure_animation.pvd";
    std::ofstream fileP(fileNameP.c_str());
    GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

    startAnimationFile(fileU);
    startAnimationFile(fileP);

    plotCurrentTimeStep(fileU, fileP, plotPts);

    for (int i = 0; i < totalIter; i += iterStep)
    {
        this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

        plotCurrentTimeStep(fileU, fileP, plotPts);
    }

    endAnimationFile(fileU);
    endAnimationFile(fileP);
}


template<class T>
void gsINSSolverUnsteady<T>::initGeneralizedStokesSolution(gsSparseMatrix<T>& stokesMatrix, gsMatrix<T>& stokesRhs)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    getAssembler()->fillStokesSystem_into(stokesMatrix, stokesRhs);

    const int uDofs = getAssembler()->getUdofs();
    const T invTimeStep = 1. / m_timeStepSize;

    #pragma omp parallel for num_threads(getAssembler()->getBlockAssembler().getNumThreads())
    for (index_t col = 0; col < uDofs; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(getAssembler()->getVelocityMassMatrix(), col); it; ++it)
            for (index_t s = 0; s < getAssembler()->getTarDim(); ++s)
                stokesMatrix.coeffRef(it.row() + s * uDofs, it.col() + s * uDofs) += invTimeStep * it.value();

    m_clock.restart();
    m_solver.analyzePattern(stokesMatrix);
    m_solver.factorize(stokesMatrix);
    m_solsetupT += m_clock.stop();
}


template<class T>
void gsINSSolverUnsteady<T>::solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations)
{
    gsSparseMatrix<T> stokesMatrix;
    gsMatrix<T> stokesRhs;

    initGeneralizedStokesSolution(stokesMatrix, stokesRhs);

    const int uDofs = getAssembler()->getUdofs();
    const T invTimeStep = 1. / m_timeStepSize;
    int iter = 0;
    T relNorm = std::numeric_limits<T>::infinity();

    while ((iter < minIterations) || ((relNorm > epsilon) && (iter < maxIterations)))
    {
        gsInfo << "Iteration number " << iter + 1 << "...";

        gsMatrix<T> rhs = stokesRhs;
        gsMatrix<T> newSol = m_solution;

        for (index_t s = 0; s < getAssembler()->getTarDim(); ++s)
            rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * getAssembler()->getVelocityMassMatrix() * m_solution.middleRows(s * uDofs, uDofs);

        m_clock.restart();
        newSol = m_solver.solve(rhs);
        m_solveT += m_clock.stop();

        relNorm = this->solutionChangeRelNorm(m_solution, newSol);
        gsInfo << " Solution change relative norm: " << relNorm << "\n";

        m_solution = newSol;
        iter++;
    }
}


template<class T, class SolverType>
void gsINSSolverUnsteadyIter<T, SolverType>::solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations)
{
    gsSparseMatrix<T> stokesMatrix;
    gsMatrix<T> stokesRhs;

    Base::initGeneralizedStokesSolution(stokesMatrix, stokesRhs);

    const int uDofs = this->getAssembler()->getUdofs();
    const T invTimeStep = 1. / m_timeStepSize;
    int iter = 0;
    T relNorm = std::numeric_limits<T>::infinity();

    while ((iter < minIterations) || ((relNorm > epsilon) && (iter < maxIterations)))
    {
        gsInfo << "Iteration number " << iter + 1 << "...";

        gsMatrix<T> rhs = stokesRhs;
        gsMatrix<T> newSol = m_solution;

        for (index_t s = 0; s < this->getAssembler()->getTarDim(); ++s)
            rhs.middleRows(s * uDofs, uDofs).noalias() += invTimeStep * this->getAssembler()->getVelocityMassMatrix() * m_solution.middleRows(s * uDofs, uDofs);

        m_matrices.at("matNS") = stokesMatrix;

        if (m_precType.substr(0, 2) == "AL")
        {
            m_matrices.at("matNS") += m_matGammaPart;
            rhs += m_rhsGammaPart;
        }

        m_pPrec->update(m_matrices);
        
        SolverType solver(m_matrices.at("matNS"), m_pPrec);
        solver.setMaxIterations(m_maxLinIter);
        solver.setTolerance(m_linTol);
        solver.solve(rhs, newSol);

        relNorm = this->solutionChangeRelNorm(m_solution, newSol);
        gsInfo << " Solution change relative norm: " << relNorm << "\n";

        m_solution = newSol;
        iter++;
    }
}

} // namespace gismo
