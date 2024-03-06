/** @file gsFlowSolverBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>

namespace gismo
{

template<class T>
void gsFlowSolverBase<T>::initMembers()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
    m_iterationNumber = 0;
    m_relNorm = std::numeric_limits<T>::infinity();

    m_fileOutput = m_params.options().getSwitch("fileOutput");
    m_dispOutput = !m_params.options().getSwitch("quiet");

    if(m_fileOutput)
        createOutputFile();

    m_initAssembT = 0;
    m_assembT = 0;
}


template<class T>
void gsFlowSolverBase<T>::updateSizes()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
}


template<class T>
void gsFlowSolverBase<T>::createOutputFile()
{
    std::string fileName = m_params.options().getString("outFile");

    if (fileName == "")
        fileName = this->getName() + "_output.txt";

    m_outFile.open(fileName);

    std::stringstream output;
    output << "\n" << m_params.options() << "\n";
    gsWriteOutput(m_outFile, output.str(), m_fileOutput, false);
}


template<class T>
real_t gsFlowSolverBase<T>::stopwatchStart()
{

#ifdef GISMO_WITH_PETSC
    if (m_params.options().getSwitch("parallel"))
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
real_t gsFlowSolverBase<T>::stopwatchStop()
{

#ifdef GISMO_WITH_PETSC
    if (m_params.options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        return m_clock.stop();

}


template<class T>
void gsFlowSolverBase<T>::initialize()
{ 
    if (!getAssembler()->isInitialized())
    {
        real_t time0 = stopwatchStart();
        getAssembler()->initialize();
        real_t time1 = stopwatchStop();

        m_initAssembT += time1 - time0;
    }
}


template<class T>
void gsFlowSolverBase<T>::nextIteration(const unsigned numberOfIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    for (unsigned iter = 0; iter < numberOfIterations; iter++)
        nextIteration();
}


template<class T>
void gsFlowSolverBase<T>::solve(const int maxIterations, const T epsilon, const int minIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
    int iter = 0;
    m_relNorm = solutionChangeRelNorm();

    while ((iter < minIterations) || ((m_relNorm > epsilon) && (iter < maxIterations)))
    {
        m_outStream.str("");
        m_outStream << "Iteration number " << m_iterationNumber + 1 << "...";
        gsWriteOutput(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);

        nextIteration();
        m_relNorm = solutionChangeRelNorm();

        m_outStream.str("");
        m_outStream << " Solution change relative norm: " << m_relNorm;
        gsWriteOutputLine(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);

        iter++;
    }
}


template<class T>
void gsFlowSolverBase<T>::updateAssembler(const gsMatrix<T>& sol, bool updateSol)
{ 
    real_t time0 = stopwatchStart();
    getAssembler()->update(sol, updateSol);
    real_t time1 = stopwatchStop();

    m_assembT += time1 - time0;
}


template<class T>
T gsFlowSolverBase<T>::solutionChangeRelNorm() const
{
    T relNorm;

    if (m_iterationNumber)
        relNorm = solutionChangeRelNorm(getAssembler()->getSolution(), m_solution);
    else
        relNorm = std::numeric_limits<T>::infinity();

    return relNorm;
}


template<class T>
T gsFlowSolverBase<T>::solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsMatrix<T> solChangeVector = solOld - solNew;
    T relNorm = solChangeVector.norm() / solNew.norm();

    return relNorm;
}


template<class T>
void gsFlowSolverBase<T>::writeSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew)
{
    m_outStream.str("");
    m_outStream << "     [u, p] solution change relative norm: ";

    for (int i = 0; i < solOld.cols(); i++)
        m_outStream << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

    gsWriteOutputLine(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);
}


template<class T>
void gsFlowSolverBase<T>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk)
{
    getAssembler()->markDofsAsEliminatedZeros(boundaryDofs, unk);
    updateSizes();
}


} // namespace gismo