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

    // inicializace turbulentniho solveru
    //m_TMsolverPtr = new gsTMSolverSST<T, MatOrder>(m_paramsPtr);

    // initialize turbulence solver
    /*if (!(m_pTurbulenceSolver->isInitialized()))
        {
            gsField<T> uSol = this->constructSolution(0);

            m_clock.restart();
            m_pTurbulenceSolver->initialize(uSol);
            m_turbT += m_clock.stop();
        }
    */
    
    //m_timeStepSize = m_params.options().getReal("timeStep");
}

// doplnit vykreslovani turbulentni viskozity a promennych turbulentniho modelu
template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::ofstream& fileN, std::ofstream& fileTM, std::string fileNameSuffix, unsigned plotPts)
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

    /*
    gsField<T> turbSol = m_TMsolverPtr->constructSolution();
    std::stringstream filenameTM;
    filenameTM << "TMsol" + fileNameSuffix + "_" << m_iterationNumber << "it";
    gsWriteParaview<T>(turbSol, filenameTM.str(), plotPts);

    std::stringstream filenameN;
    filenameN << "nuT_" + fileNameSuffix + "_" << m_iterationNumber << "it";
    plotTurbulentViscosity(filenameN.str(), plotPts);
    */

    for (int p = 0; p < numPatches; p++)
    {
        std::stringstream fnU;
        fnU << filenameU.str() << p << ".vts";
        fileU << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnU.str() << "\"/>\n";

        std::stringstream fnP;
        fnP << filenameP.str() << p << ".vts";
        fileP << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnP.str() << "\"/>\n";

        /*
        std::stringstream fnTM;
        fnTM << filenameTM.str() << p << ".vts";
        fileTM << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnTM.str() << "\"/>\n";

        std::stringstream fnN;
        fnN << filenameN.str() << p << ".vts";
        fileN << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fnN.str() << "\"/>\n";
        */
    }
}

// upravit
template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::nextIteration()
{
    //GISMO_ASSERT(this->getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    if (!m_bComputeTMfirst)
        Base::nextIteration();

    gsField<T> uSolField = this->constructSolution(0);
    m_paramsPtr->getSSTModel().setVelocitySolution(uSolField);

    m_clock.restart();
    m_TMsolverPtr->nextIteration(); // update turbulence model
    m_turbT += m_clock.stop(); // NOTE: this is not exact, because there is some assembly and solver setup in the TM nextIteration() method

    ////------------------------- PLOT -------------------------
    //std::string path = "D:/hhornik/gismo/motor/uwb-pilsen/outFiles/k-omega/";
    //std::stringstream filename;
    //filename << path << "k_omega_iter" << this->m_iterationNumber;
    //gsField<T> tmSol = m_pTurbulenceSolver->constructSolution();
    //gsWriteParaview<T>(tmSol, filename.str(), 50000);
    ////------------------------- PLOT -------------------------
        
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

// upravit, melo by stacit uz jen upravit plotCurrentTimeStep()
template<class T, int MatOrder>
void gsRANSSolverUnsteady<T, MatOrder>::solveWithAnimation(const int totalIter, const int iterStep, std::string fileNameSuffix, const T epsilon, unsigned plotPts, bool plotTurb, const int minIterations)
{
    if (!plotTurb)
            Base::solveWithAnimation(totalIter, iterStep, "", epsilon, plotPts);
    else
    {
        // prepare plotting
        std::string fileNameU = "velocity" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileU(fileNameU.c_str());
        GISMO_ASSERT(fileU.is_open(), "Error creating " << fileNameU);

        std::string fileNameP = "pressure" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileP(fileNameP.c_str());
        GISMO_ASSERT(fileP.is_open(), "Error creating " << fileNameP);

        std::string fileNameN = "nuT" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileN(fileNameN.c_str());
        GISMO_ASSERT(fileN.is_open(), "Error creating " << fileNameN);

        std::string fileNameTM = "TMsol" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileTM(fileNameTM.c_str());
        GISMO_ASSERT(fileTM.is_open(), "Error creating " << fileNameTM);

        startAnimationFile(fileU);
        startAnimationFile(fileP);
        startAnimationFile(fileN);
        startAnimationFile(fileTM);

        plotCurrentTimeStep(fileU, fileP, fileN, fileTM, fileNameSuffix, plotPts);

        for (int i = 0; i < totalIter; i += iterStep)
        {
            this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

            plotCurrentTimeStep(fileU, fileP, fileN, fileTM, fileNameSuffix, plotPts);
        }

        endAnimationFile(fileU);
        endAnimationFile(fileP);
        endAnimationFile(fileN);
        endAnimationFile(fileTM);
    }
}

} //namespace gismo