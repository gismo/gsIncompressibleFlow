/** @file uwbADRSolverUnsteady.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRSolverBase.h"
#include "uwbADRAssemblerBase.h"

namespace gismo
{

template<class T>
class uwbADRSolverUnsteady : public uwbADRSolverBase<T>
{

public:
    typedef uwbADRSolverBase<T> Base;

public:

    uwbADRSolverUnsteady(uwbADRSolverParams<T>& params, bool createAssembler = true)
    {
        //create assembler
        if (createAssembler)
        {
                m_pAssembler = new uwbADRAssemblerBaseUnsteady<T>(params, 1);

            initMembers();
        }

        m_innerIter = params.settings().get(constantsADR::unst_innerIt);
        m_innerTol = params.settings().get(constantsADR::unst_innerTol);
    }

    virtual ~uwbADRSolverUnsteady()
    { }

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
        getAssembler()->setInitialCondition(solVector);

        m_iterationNumber = 0;
        m_time = 0;
    }

    void changeTimeStepSize(const T timeStepSize)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
        getAssembler()->changeTimeStep(timeStepSize);
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

    void solveWithAnimation(const int totalIter, const int iterStep, const T epsilon = 1e-3, unsigned plotPts = 10000)
    {
        // prepare plotting
        std::string fileName = "ADR_solAnimation.pvd";
        std::ofstream file(fileName.c_str());
        GISMO_ASSERT(file.is_open(), "Error creating " << fileName);

        startAnimationFile(file);

        plotCurrentTimeStep(file, plotPts);

        for (int i = 0; i < totalIter; i += iterStep)
        {
            this->solve(math::min(iterStep, totalIter), epsilon);

            plotCurrentTimeStep(file, plotPts);
        }

        endAnimationFile(file);
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

    void plotCurrentTimeStep(std::ofstream& file, unsigned plotPts)
    {
        int numPatches = getAssembler()->getBlockAssembler().getPatches().nPatches();

        gsField<T> sol = this->constructSolution();
        std::stringstream filename;
        filename << "ADR_sol" << m_iterationNumber << "it";
        gsWriteParaview<T>(sol, filename.str(), plotPts);

        for (int p = 0; p < numPatches; p++)
        {
            std::stringstream fn;
            fn << filename.str() << p << ".vts";
            file << "<DataSet timestep = \"" << m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fn.str() << "\"/>\n";
        }
    }

public:

    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION }
    unsigned getTimeStepNumber() const { return m_iterationNumber; }
    T getSimulationTime() const { return m_time; }

    virtual uwbADRAssemblerBaseUnsteady<T>* getAssembler() const { return dynamic_cast<uwbADRAssemblerBaseUnsteady<T>*>(m_pAssembler); }

protected:

    T m_time;
    T m_innerIter;
    T m_innerTol;

    // members from uwbADRSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;

}; // class uwbADRSolverUnsteady

} // namespace gismo
