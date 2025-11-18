/** @file gsRANSSolverUnsteady.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsRANSSolverUnsteady.h>

namespace gismo
{

template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();

    getAssembler()->setTurbulenceSolver(m_TMsolverPtr);

    m_bComputeTMfirst = false;
    m_turbT = 0;

}


template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::ofstream& fileK, std::ofstream& fileO, std::ofstream& fileTV, std::string fileNamePrefix, unsigned plotPts)
{
    int numPatches = m_paramsPtr->getPde().patches().nPatches();

    gsField<T> uSol = this->constructSolution(0);
    std::stringstream filenameU;
    filenameU << fileNamePrefix + "_velocity_" << m_iterationNumber << "it";
    gsWriteParaview<T>(uSol, filenameU.str(), plotPts);

    gsField<T> pSol = this->constructSolution(1);
    std::stringstream filenameP;
    filenameP << fileNamePrefix + "_pressure_" << m_iterationNumber << "it";
    gsWriteParaview<T>(pSol, filenameP.str(), plotPts);

    gsField<T> kSol = constructSolution(2);
    std::stringstream filenameK;
    filenameK << fileNamePrefix + "_Ksol_" << m_iterationNumber << "it";
    gsWriteParaview<T>(kSol, filenameK.str(), plotPts);

    gsField<T> oSol = constructSolution(3);
    std::stringstream filenameO;
    filenameO << fileNamePrefix + "_Osol__" << m_iterationNumber << "it";
    gsWriteParaview<T>(oSol, filenameO.str(), plotPts);
    
    std::stringstream filenameTV;
    filenameTV << fileNamePrefix + "_turbViscosity_" << m_iterationNumber << "it";
    m_TMModelPtr->plotTurbulentViscosity(m_paramsPtr, filenameTV.str());
        
    for (int p = 0; p < numPatches; p++)
    {
        std::stringstream fnU;
        fnU << filenameU.str() << p << ".vts";
        fileU << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnU.str() << "\"/>\n";

        std::stringstream fnP;
        fnP << filenameP.str() << p << ".vts";
        fileP << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnP.str() << "\"/>\n";

        std::stringstream fnK;
        fnK << filenameK.str() << p << ".vts";
        fileK << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnK.str() << "\"/>\n";

        std::stringstream fnO;
        fnO << filenameO.str() << p << ".vts";
        fileO << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnO.str() << "\"/>\n";

        std::stringstream fnTV;
        fnTV << filenameTV.str() << p << ".vts";
        fileTV << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnTV.str() << "\"/>\n";
    }
}

template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::initialize()
{ 
    Base::initialize();
    m_TMsolverPtr->initialize();

    T viscSteady = 0.1;
    gsNavStokesPde<T> NSSteadyPde(m_paramsPtr->getPde().patches(), m_paramsPtr->getPde().bc(), m_paramsPtr->getPde().force(), viscSteady);
    gsFlowSolverParams<T> paramsSteady(NSSteadyPde, m_paramsPtr->getBases(), NULL, m_paramsPtr->getMpiComm());
    paramsSteady.copyAllOptionsFrom(*m_paramsPtr);
    gsINSSolverSteady<T, MatOrder> NSSteadySolver(paramsSteady);

    m_paramsPtr->logger() << "\n--------------------------\n";
    m_paramsPtr->logger() << "Computing initial solution for RANS by solving steady N-S problem for viscosity " << viscSteady << "\n";

    // gsStopwatch clock;

    m_paramsPtr->logger() << "Initialization...\n";
    NSSteadySolver.initialize();
    NSSteadySolver.solveStokes(); // start from Stokes solution
    NSSteadySolver.solve(m_paramsPtr->options().getInt("nonlin.maxIt"), m_paramsPtr->options().getReal("nonlin.tol"), 0);    

    // real_t totalT = clock.stop();

    // m_paramsPtr->logger() << "\nAssembly time (steady):" << NSSteadySolver.getAssemblyTime() << "\n";
    // m_paramsPtr->logger() << "Solve time (steady):" << NSSteadySolver.getSolveTime() << "\n";
    // m_paramsPtr->logger() << "Solver setup time (steady):" << NSSteadySolver.getSolverSetupTime() << "\n";
    // m_paramsPtr->logger() << "Total solveProblem time (steady):" << totalT << "\n\n";

    // gsField<> velocity = NSSteadySolver.constructSolution(0);
    // gsField<> pressure = NSSteadySolver.constructSolution(1);
    // int plotPts = 20000;

    // m_paramsPtr->logger() << "Plotting solution of steady N-S problem in Paraview...";
    // gsWriteParaview<T>(velocity, "NS_steady_velocity", plotPts, false);
    // gsWriteParaview<T>(pressure, "NS_steady_pressure", plotPts);
    // m_paramsPtr->logger() << " Done.\n";

    m_paramsPtr->logger() << "--------------------------\n";

    m_solution = NSSteadySolver.getSolution();
}

// upravit
template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::nextIteration()
{
    //GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    if (!m_bComputeTMfirst)
        Base::nextIteration();

    gsField<T> uSolField = this->constructSolution(0);
    m_paramsPtr->setVelocitySolution(uSolField);

    m_clock.restart();
    m_TMsolverPtr->nextIteration(); // update turbulence model
    m_turbT += m_clock.stop(); // NOTE: this is not exact, because there is some assembly and solver setup in the TM nextIteration() method

    if (m_bComputeTMfirst)
        Base::nextIteration();

        /*
        if ((this->m_iterationNumber % 250) == 0)
        {
            gsFileData<> fd;
            fd << m_solution;
            fd.save("RANS_solution_iter" + std::to_string(this->m_iterationNumber) + ".xml");

            gsFileData<> fd_TM;
            fd_TM << m_pTurbulenceSolver->getSolution();
            fd_TM.save("TM_solution_iter" + std::to_string(this->m_iterationNumber) + ".xml");
        }
        */
}


template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::solveWithAnimation(const int totalIter, const int iterStep, std::string fileNamePrefix, const T epsilon, unsigned plotPts, bool plotTurb, const int minIterations)
{
    if (!plotTurb)
            Base::solveWithAnimation(totalIter, iterStep, "", epsilon, plotPts);
    else
    {
        // prepare plotting
        std::string fileNameU = fileNamePrefix + "_velocity_animation.pvd";
        std::ofstream fileU(fileNameU.c_str());
        GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

        std::string fileNameP = fileNamePrefix + "_pressure_animation.pvd";
        std::ofstream fileP(fileNameP.c_str());
        GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

        std::string fileNameK = fileNamePrefix + "_Ksol_animation.pvd";
        std::ofstream fileK(fileNameK.c_str());
        GISMO_ASSERT(fileK.is_open(), "Error creating " << fileNameK);

        std::string fileNameO = fileNamePrefix + "_Osol_animation.pvd";
        std::ofstream fileO(fileNameO.c_str());
        GISMO_ASSERT(fileO.is_open(), "Error creating " << fileNameO);

        std::string fileNameTV = fileNamePrefix + "_turbViscosity_animation.pvd";
        std::ofstream fileTV(fileNameTV.c_str());
        GISMO_ASSERT(fileTV.is_open(), "Error creating " << fileNameTV);

        startAnimationFile(fileU);
        startAnimationFile(fileP);
        startAnimationFile(fileK);
        startAnimationFile(fileO);
        startAnimationFile(fileTV);

        plotCurrentTimeStep(fileU, fileP, fileK, fileO, fileTV, fileNamePrefix, plotPts);

        for (int i = 0; i < totalIter; i += iterStep)
        {
            this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

            plotCurrentTimeStep(fileU, fileP, fileK, fileO, fileTV, fileNamePrefix, plotPts);
        }

        endAnimationFile(fileU);
        endAnimationFile(fileP);
        endAnimationFile(fileK);
        endAnimationFile(fileO);
        endAnimationFile(fileTV);
    }
}

} //namespace gismo