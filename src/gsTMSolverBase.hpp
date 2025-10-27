/** @file gsTMSolverBase.hpp

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
        gsWarn << "Unknown identifier of a turbulent model entered! Using k-omega SST model." << std::endl;
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

    index_t picardIter = 0;
    T relNorm = 1.0;

    m_paramsPtr->logger() << "        TM solver Picard's iterations...\n";

    while((relNorm > m_TMinnerTol) && (picardIter < m_TMinnerIter))
    {
        m_paramsPtr->logger() << "         ";

        gsMatrix<T> oldSol = tmpSolution;
        
        this->updateAssembler(tmpSolution, false);
        this->applySolver(tmpSolution);

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