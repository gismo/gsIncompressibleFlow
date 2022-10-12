/** @file uwbINSSolverUnsteady.h

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSSolverBase.h"
#include "uwbINSAssemblerUnsteady.h"
#include "uwbINSAssemblerUnsteadyPeriodic.h"

namespace gismo
{

template<class T>
class uwbINSSolverUnsteady : public uwbINSSolverBase<T>
{

public:
    typedef uwbINSSolverBase<T> Base;

public:

    uwbINSSolverUnsteady(uwbINSSolverParams<T>& params, bool createAssembler = true)
    {
        //create assembler
        if (createAssembler)
        {
            if (params.getBCs().numPeriodic())
                m_pAssembler = new uwbINSAssemblerUnsteadyPeriodic<T>(params);
            else
                m_pAssembler = new uwbINSAssemblerUnsteady<T>(params);


            initMembers();
        }

        m_innerIter = params.settings().get(constantsINS::unst_innerIt);
        m_innerTol = params.settings().get(constantsINS::unst_innerTol);
        m_avgPicardIter = 0;

        m_alpha_u = params.settings().get(constantsINS::alpha_u);
        m_alpha_p = params.settings().get(constantsINS::alpha_p);
    }

    virtual ~uwbINSSolverUnsteady()
    {
    }

protected:

    virtual void initMembers()
    {
        Base::initMembers();
        m_time = 0;
    }

public:
    void setInitialCondition(const gsMatrix<T> & solVector)
    {
        this->setSolution(solVector);

        m_iterationNumber = 0;
        m_time = 0;
    }

    void setStokesInitialCondition()
    {
        this->setStokesSolution();

        m_time = 0;
    }

    void changeTimeStepSize(const T timeStepSize)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->changeTimeStep(timeStepSize);
    }

    virtual void updateAssemblerPicard(const gsMatrix<T>& sol)
    {
        m_clock.restart();
        getAssembler()->updatePicard(sol);
        m_assembT += m_clock.stop();
    }

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        this->updateAssembler();

        // in the first time step the pattern will get analyzed
        if (!m_iterationNumber)
            this->initIteration();

        gsMatrix<T> tmpSolution = m_solution;

        this->applySolver(tmpSolution);
        //this->applySolver(tmpSolution, m_alpha_u, m_alpha_p);

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
            //this->applySolver(tmpSolution, m_alpha_u, m_alpha_p);

            dispSolChangeRelNorm(oldSol, tmpSolution);

            relNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
            iter++;
        }

        m_avgPicardIter += iter;
        
        m_solution = tmpSolution;

        m_iterationNumber++;
        m_time += getAssembler()->getTimeStepSize();
    }

    void nextTimeStep() { nextIteration(); }

    void nextTimeStep(const unsigned numberOfSteps)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        gsInfo << "Simulation ... \n";
        for (unsigned step = 0; step < numberOfSteps; step++)
        {
            nextIteration();
            gsInfo << "Step number " << m_iterationNumber << ", time " << m_time << " s, done. \n";
        }
    }

    void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     [u, p] Picard's solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << this->solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

        gsInfo << "\n";
    }

    void solveWithAnimation(const int totalIter, const int iterStep, const T epsilon = 1e-3, unsigned plotPts = 10000, const int minIterations = 1)
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

protected:
    void startAnimationFile(std::ofstream& file)
    {
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"Collection\" version=\"0.1\">";
        file << "<Collection>\n";
    }

    void endAnimationFile(std::ofstream& file)
    {
        file << "</Collection>\n";
        file << "</VTKFile>\n";
        file.close();
    }

    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, unsigned plotPts)
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

public:

    unsigned getTimeStepNumber() const { return m_iterationNumber; }
    T getSimulationTime() const { return m_time; }
    T getAvgPicardIterations() const { return m_avgPicardIter / m_iterationNumber; }

    virtual uwbINSAssemblerUnsteady<T>* getAssembler() const
    {
        uwbINSAssemblerUnsteadyPeriodic<T>* pAssembler = dynamic_cast<uwbINSAssemblerUnsteadyPeriodic<T>*>(m_pAssembler);

        if (pAssembler != NULL)
            return pAssembler;
        else
            return dynamic_cast<uwbINSAssemblerUnsteady<T>*>(m_pAssembler);
    }

protected:

    T m_time;
    T m_innerIter, m_avgPicardIter;
    T m_innerTol;

    real_t m_alpha_u, m_alpha_p;

    // members from uwbINSSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_clock;
    using Base::m_assembT;

}; // class uwbINSSolverUnsteady

} // namespace gismo
