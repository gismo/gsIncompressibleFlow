/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

namespace gismo
{

template <class T, int MatOrder>
typename gsTMSolverBase<T, MatOrder>::tmPtr gsTMSolverBase<T, MatOrder>::make(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr)
{
    std::string turbModel = paramsPtr->options().getString("TM");
    if (turbModel == "SST") 
    {
        return gsTMSolverSST<T, MatOrder>::make(paramsPtr, TMModelPtr);
    }
    //elseif (m_paramsPtr->options().getSwitch("TM.eval") == "SA") 
    //{ }
    else 
    {
        gsInfo << "Unknown identifier of a turbulent model entered! Using k-omega SST model." << std::endl;
        return gsTMSolverSST<T, MatOrder>::make(paramsPtr, TMModelPtr);
    }
}


template<class T, int MatOrder>
void gsTMSolverBase<T, MatOrder>::initMembers()
{
    Base::initMembers();

    m_solution = getAssembler()->getSolution();
    
    m_TMtime = 0;
    m_TMtimeStepSize = m_paramsPtr->options().getReal("timeStep");
    m_TMinnerIter = m_paramsPtr->options().getInt("TM.maxIt");
    m_TMinnerTol = m_paramsPtr->options().getReal("TM.tol");
    m_TMavgPicardIter = 0;
}


template<class T, int MatOrder>
void gsTMSolverBase<T, MatOrder>::nextIteration()
{
    GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    this->updateAssembler();

    if (!m_iterationNumber)
        this->initIteration();

    gsMatrix<T> tmpSolution = m_solution;

    //this->applySolver(tmpSolution);

    //this->writeSolChangeRelNorm(m_solution, tmpSolution, "TM");

    index_t picardIter = 0;
    //T relNorm = this->solutionChangeRelNorm(m_solution, tmpSolution);
    T relNorm = 1.0;

    gsWriteOutputLine(m_outFile, "        TM solver Picard's iterations...", m_fileOutput, m_dispOutput);

    while((relNorm > m_TMinnerTol) && (picardIter < m_TMinnerIter))
    {
        gsWriteOutput(m_outFile, "         ", m_fileOutput, m_dispOutput);

        gsMatrix<T> oldSol = tmpSolution;
        //gsInfo << "Current solution sum: " << tmpSolution.sum() << std::endl;

        this->updateAssembler(tmpSolution, false);
        this->applySolver(tmpSolution);

        // bool cropped = false;
        // for (index_t i = 0; i < tmpSolution.rows(); i++)
        //     if (tmpSolution(i) < 0.0)
        //         {
        //             tmpSolution(i) = 0.0;
        //             cropped = true;
        //         }
        // //tmpSolution(i) = math::max(tmpSolution(i), 0.0);
        // if (cropped)
        //     gsInfo << "k, omega oriznuto nulou ..." << std::endl;

        this->writeSolChangeRelNorm(oldSol, tmpSolution, "TM");

        relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
        picardIter++;
    }
    
    m_solution = tmpSolution;

    m_TMtime += m_TMtimeStepSize;
    m_TMavgPicardIter += picardIter;
    m_iterationNumber++;
}


} // namespace gismo