/** @file uwbRANSSolver.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbINSSolverUnsteady.h"
#include "uwbTMSolverBase.h"
#include "uwbRANSAssembler.h"
#include "uwbRANSAssemblerPeriodic.h"

namespace gismo
{

template<class T>
class uwbRANSSolver : public uwbINSSolverUnsteady<T>
{

public:
    typedef uwbINSSolverUnsteady<T> Base;

    uwbRANSSolver(uwbINSSolverParams<T>& params) : Base(params, false)
    {
        // save pointer to turbulence solver
        m_pTurbulenceSolver = params.settings().getTurbulenceSolver();

        // create assembler
        if (!(params.getBCs().numPeriodic()))
            m_pAssembler = new uwbRANSAssembler<T>(params, m_pTurbulenceSolver);
        else
            m_pAssembler = new uwbRANSAssemblerPeriodic<T>(params, m_pTurbulenceSolver);

        Base::initMembers();

        m_bComputeTMfirst = false;
        m_turbT = 0;
    }

    virtual ~uwbRANSSolver()
    {  }

public:

    virtual void initialize()
    {
        Base::initialize();

        if (!(m_pTurbulenceSolver->isInitialized()))
        {
            gsField<T> uSol = this->constructSolution(0);

            m_clock.restart();
            m_pTurbulenceSolver->initialize(uSol);
            m_turbT += m_clock.stop();
        }
    }

    virtual void nextIteration()
    {
        if (!m_bComputeTMfirst)
            Base::nextIteration();

        gsField<T> uSol = this->constructSolution(0);
        m_pTurbulenceSolver->updateVelocitySolution(uSol);

        //gsMatrix<T> oldSol = m_pTurbulenceSolver->getSolution();
        //m_pTurbulenceSolver->getAssembler()->setOldSolution(oldSol);

        m_clock.restart();
        m_pTurbulenceSolver->nextIteration(); // update turbulence model
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

        /*if ((this->m_iterationNumber % 250) == 0)
        {
            gsFileData<> fd;
            fd << m_solution;
            fd.save("RANS_solution_iter" + std::to_string(this->m_iterationNumber) + ".xml");

            gsFileData<> fd_TM;
            fd_TM << m_pTurbulenceSolver->getSolution();
            fd_TM.save("TM_solution_iter" + std::to_string(this->m_iterationNumber) + ".xml");
        }*/
    }

    virtual uwbRANSAssembler<T>* getAssembler() const
    {
        uwbRANSAssemblerPeriodic<T>* pAssembler = dynamic_cast<uwbRANSAssemblerPeriodic<T>*>(m_pAssembler);

        if (pAssembler != NULL)
            return pAssembler;
        else
            return dynamic_cast<uwbRANSAssembler<T>*>(m_pAssembler);
    }

    void solveWithAnimation(const int totalIter, const int iterStep, const T epsilon = 1e-3, unsigned plotPts = 10000, bool plotTurb = false, const int minIterations = 1)
    {
        if (!plotTurb)
            Base::solveWithAnimation(totalIter, iterStep, epsilon, plotPts);
        else
        {
            // prepare plotting
            std::string fileName1 = "velocity_animation.pvd";
            std::string fileName2 = "TMSol_animation.pvd";
            std::string fileName3 = "nuT_animation.pvd";
            std::string fileName4 = "pressure_animation.pvd";
            std::ofstream file1(fileName1.c_str());
            std::ofstream file2(fileName2.c_str());
            std::ofstream file3(fileName3.c_str());
            std::ofstream file4(fileName4.c_str());
            GISMO_ASSERT(file1.is_open(), "Error creating " << fileName1);
            GISMO_ASSERT(file2.is_open(), "Error creating " << fileName2);
            GISMO_ASSERT(file3.is_open(), "Error creating " << fileName3);
            GISMO_ASSERT(file4.is_open(), "Error creating " << fileName4);

            startAnimationFile(file1);
            startAnimationFile(file2);
            startAnimationFile(file3);
            startAnimationFile(file4);

            plotCurrentTimeStep(file1, file2, file3, file4, plotPts);


            for (int i = 0; i < totalIter; i += iterStep)
            {
                this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

                plotCurrentTimeStep(file1, file2, file3, file4, plotPts);
            }

            /*for (int i = 0; i < 200; i += 10)
            {
                this->solve(math::min(10, totalIter), epsilon, minIterations);

                plotCurrentTimeStep(file1, file2, file3, file4, plotPts);
            }

            for (int i = 200; i < totalIter; i += iterStep)
            {
                this->solve(math::min(iterStep, totalIter), epsilon, minIterations);

                plotCurrentTimeStep(file1, file2, file3, file4, plotPts);
            }*/

            endAnimationFile(file1);
            endAnimationFile(file2);
            endAnimationFile(file3);
            endAnimationFile(file4);
        }
    }

    gsField<T> constructTurbulentViscositySol()
    {
        gsField<T> uSol = this->constructSolution(0);
        m_pTurbulenceSolver->updateVelocitySolution(uSol);
        return m_pTurbulenceSolver->constructTurbulentViscositySol();
    }

    void plotTurbulentViscosity(std::string const & fn, unsigned npts = 10000)
    {
        gsField<T> uSol = this->constructSolution(0);
        m_pTurbulenceSolver->updateVelocitySolution(uSol);
        m_pTurbulenceSolver->plotTurbulentViscosity(fn, npts);
    }

    gsField<T> constructTMCoefficientSol(std::string coeffType)
    {
        gsField<T> uSol = this->constructSolution(0);
        m_pTurbulenceSolver->updateVelocitySolution(uSol);
        return m_pTurbulenceSolver->constructCoefficientSol(coeffType);
    }

    const T getTurbModelTime() const { return m_turbT; }
    const T getTurbAssemblyTime() const { return m_pTurbulenceSolver->getAssemblyTime(); }
    const T getTurbSolverSetupTime() const { return m_pTurbulenceSolver->getSolverSetupTime(); }
    const T getTurbSolveTime() const { return m_pTurbulenceSolver->getSolveTime(); }

protected:

    void plotCurrentTimeStep(std::ofstream& file1, std::ofstream& file2, std::ofstream& file3, std::ofstream& file4, unsigned plotPts)
    {
        int numPatches = getAssembler()->getBlockAssembler().getPatches().nPatches();

        gsField<T> uSol = this->constructSolution(0);
        std::stringstream filenameU;
        filenameU << "velocity_" << this->m_iterationNumber << "it";
        gsWriteParaview<T>(uSol, filenameU.str(), plotPts);

        gsField<T> pSol = this->constructSolution(1);
        std::stringstream filenameP;
        filenameP << "pressure_" << this->m_iterationNumber << "it";
        gsWriteParaview<T>(pSol, filenameP.str(), plotPts);

        gsField<T> turbSol = m_pTurbulenceSolver->constructSolution();
        std::stringstream filenameTurb;
        filenameTurb << "TMSol_" << this->m_iterationNumber << "it";
        gsWriteParaview<T>(turbSol, filenameTurb.str(), plotPts);

        std::stringstream filenameTVisc;
        filenameTVisc << "nuT_" << this->m_iterationNumber << "it";
        plotTurbulentViscosity(filenameTVisc.str(), plotPts);

        for (int p = 0; p < numPatches; p++)
        {
            std::stringstream fn1;
            fn1 << filenameU.str() << p << ".vts";
            file1 << "<DataSet timestep = \"" << this->m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fn1.str() << "\"/>\n";

            std::stringstream fn2;
            fn2 << filenameTurb.str() << p << ".vts";
            file2 << "<DataSet timestep = \"" << this->m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fn2.str() << "\"/>\n";

            std::stringstream fn3;
            fn3 << filenameTVisc.str() << p << ".vts";
            file3 << "<DataSet timestep = \"" << this->m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fn3.str() << "\"/>\n";

            std::stringstream fn4;
            fn4 << filenameP.str() << p << ".vts";
            file4 << "<DataSet timestep = \"" << this->m_iterationNumber << "\" part = \"" << p << "\" file = \"" << fn4.str() << "\"/>\n";
        }
    }

public:
    void setComputationSequence(bool computeTMfirst) { m_bComputeTMfirst = computeTMfirst; }

protected:
    uwbTMSolverBaseUnsteady<T>* m_pTurbulenceSolver;
    T m_turbT;

    bool m_bComputeTMfirst;

    // members from uwbINSSolverUnsteady
    using Base::m_innerIter;
    using Base::m_innerTol;
    using Base::m_avgPicardIter;
    using Base::m_time;

    // member functions from uwbINSSolverUnsteady
    using Base::startAnimationFile;
    using Base::endAnimationFile;

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_pAssembler;
    using uwbINSSolverBase<T>::m_solution;
    using uwbINSSolverBase<T>::m_relNorm;
    using uwbINSSolverBase<T>::m_clock;

}; // class uwbRANSSolver

} // namespace gismo
