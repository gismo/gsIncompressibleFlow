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
void gsRANSSolverUnsteady<T, MatOrder>::plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::ofstream& fileK, std::ofstream& fileO, std::ofstream& fileTV, std::string fileNameSuffix, unsigned plotPts)
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

    gsField<T> kSol = constructSolution(2);
    std::stringstream filenameK;
    filenameK << "Ksol" + fileNameSuffix + "_" << m_iterationNumber << "it";
    gsWriteParaview<T>(kSol, filenameK.str(), plotPts);

    gsField<T> oSol = constructSolution(3);
    std::stringstream filenameO;
    filenameO << "Osol_" + fileNameSuffix + "_" << m_iterationNumber << "it";
    gsWriteParaview<T>(oSol, filenameO.str(), plotPts);
    
    std::stringstream filenameTV;
    filenameTV << "turbViscosity" + fileNameSuffix + "_" << m_iterationNumber << "it";
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

    // compute initial solution for RANS by solving steady N-S problem
    short_t dim = m_paramsPtr->getPde().patches().domainDim();
    gsFunctionExpr<T> f; // external force
    switch(dim)
    {
        case 2:
        default:
            f = gsFunctionExpr<>("0", "0", 2);
            break;

        case 3:
            f = gsFunctionExpr<>("0", "0", "0", 3);
            break;
    }

    // std::vector<gsMultiBasis<T> > discbases;
    // discbases.push_back(m_paramsPtr->getBases().at(0));
    // discbases.push_back(m_paramsPtr->getBases().at(1));
    gsNavStokesPde<T> NSSteadyPde(m_paramsPtr->getPde().patches(), m_paramsPtr->getPde().bc(), &f, 0.1);
    gsFlowSolverParams<T> paramsSteady(NSSteadyPde, m_paramsPtr->getBases());
    gsINSSolverSteady<T, MatOrder> NSSteadySolver(paramsSteady);

    gsInfo << "\n-------------------------------------------------------------------------------\n";
    gsInfo << "Computing initial solution for RANS by solving steady N-S problem for viscosity \n";

    gsStopwatch clock;

    gsInfo << "Initialization...\n";
    NSSteadySolver.initialize();

    gsInfo << "numDofs: " << NSSteadySolver.numDofs() << "\n";

    //if (m_paramsPtr->options().getSwitch("stokesInit"))
    NSSteadySolver.solveStokes();
    
    NSSteadySolver.solve(m_paramsPtr->options().getInt("nonlin.maxIt"), m_paramsPtr->options().getReal("nonlin.tol"), 0);    

    real_t totalT = clock.stop();

    gsInfo << "\nAssembly time (steady):" << NSSteadySolver.getAssemblyTime() << "\n";
    gsInfo << "Solve time (steady):" << NSSteadySolver.getSolveTime() << "\n";
    gsInfo << "Solver setup time (steady):" << NSSteadySolver.getSolverSetupTime() << "\n";
    gsInfo << "Total solveProblem time (steady):" << totalT << "\n\n";

    gsField<> velocity = NSSteadySolver.constructSolution(0);
    gsField<> pressure = NSSteadySolver.constructSolution(1);

    int plotPts = 20000;
 
    gsInfo << "Plotting solution of steady N-S problem in Paraview...";
    gsWriteParaview<T>(velocity, "NS_steady_velocity", plotPts, false);
    gsWriteParaview<T>(pressure, "NS_steady_pressure", plotPts);
    gsInfo << " Done.\n";
    gsInfo << "-------------------------------------------------------------------------------\n";

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

        std::string fileNameK = "Ksol" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileK(fileNameK.c_str());
        GISMO_ASSERT(fileK.is_open(), "Error creating " << fileNameK);

        std::string fileNameO = "Osol" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileO(fileNameO.c_str());
        GISMO_ASSERT(fileO.is_open(), "Error creating " << fileNameO);

        std::string fileNameTV = "turbViscosity" + fileNameSuffix + "_animation.pvd";
        std::ofstream fileTV(fileNameTV.c_str());
        GISMO_ASSERT(fileTV.is_open(), "Error creating " << fileNameTV);

        startAnimationFile(fileU);
        startAnimationFile(fileP);
        startAnimationFile(fileK);
        startAnimationFile(fileO);
        startAnimationFile(fileTV);

        plotCurrentTimeStep(fileU, fileP, fileK, fileO, fileTV, fileNameSuffix, plotPts);

        for (int i = 0; i < totalIter; i += iterStep)
        {
            this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

            plotCurrentTimeStep(fileU, fileP, fileK, fileO, fileTV, fileNameSuffix, plotPts);
        }

        endAnimationFile(fileU);
        endAnimationFile(fileP);
        endAnimationFile(fileK);
        endAnimationFile(fileO);
        endAnimationFile(fileTV);
    }
}

} //namespace gismo