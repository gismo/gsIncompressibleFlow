/** @file uwbADRSolverNonlinear.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRSolverUnsteady.h"
#include "uwbADRAssemblerNonlinear.h"

namespace gismo
{

template<class T>
class uwbADRSolverNonlinear : public uwbADRSolverUnsteady<T>
{

public:
    typedef uwbADRSolverUnsteady<T> Base;

public:
    uwbADRSolverNonlinear(uwbADRSolverParams<T>& params) : Base(params, false)
    {
        gsInfo << "Solving Burger's equation...\n";
        params.settings().setADREvaluator("nonlinCoeffsBurgers");
        m_pAssembler = new uwbADRAssemblerNonlinear<T>(params);
        m_bUnsteady = params.settings().get(constantsADR::unsteady);

        Base::initMembers();
    }

    uwbADRSolverNonlinear(uwbADRSolverParams<T>& params, gsField<T>& advectionField,
                          gsField<T>& diffusionField) : Base(params, false)
    {
        gsInfo << "Solving nonlinear advection-diffusion equation...\n";
        params.settings().setADREvaluator("nonlinCoeffsField");
        m_pAssembler = new uwbADRAssemblerNonlinear<T>(params, advectionField, diffusionField);
        m_bUnsteady = params.settings().get(constantsADR::unsteady);

        Base::initMembers();
    }

    uwbADRSolverNonlinear(uwbADRSolverParams<T>& params, gsField<T>& advectionField,
                          gsField<T>& diffusionField, gsField<T>& reactionField) : Base(params, false)
    {
        gsInfo << "Solving nonlinear advection-diffusion-reaction equation...\n";
        params.settings().setADREvaluator("nonlinCoeffsField");
        m_pAssembler = new uwbADRAssemblerNonlinear<T>(params, advectionField, diffusionField, reactionField);
        m_bUnsteady = params.settings().get(constantsADR::unsteady);

        Base::initMembers();
    }

    uwbADRSolverNonlinear(uwbADRSolverParams<T>& params, T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff) : Base(params, false)
    {
        gsInfo << "Solving linear advection-diffusion-reaction equation...\n";
        params.settings().setADREvaluator("linConstCoeffs");
        m_pAssembler = new uwbADRAssemblerNonlinear<T>(params, diffusionCoeff, advectionCoeff, reactionCoeff);
        m_bUnsteady = params.settings().get(constantsADR::unsteady);

        Base::initMembers();
    }

    virtual ~uwbADRSolverNonlinear() { }

public:
    void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
    {
        if (m_bUnsteady)
        {
            gsInfo << "...solving unsteady case.\n";
            Base::solve(maxIterations, epsilon, minIterations);
        }
        else
        {
            gsInfo << "...solving steady case.\n";

            GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

            int iter = 0;
            T solRelNorm = std::numeric_limits<T>::infinity();
            gsInfo << "        Picard's iterations...\n";
            while ((iter < minIterations) || ((math::abs(solRelNorm) > m_innerTol) && (iter < m_innerIter)))
            {
                gsInfo << "         ";
                gsMatrix<T> oldSol = m_solution;
                getAssembler()->updatePicardSteady(m_solution);

                if(!iter)
                    m_solver.analyzePattern(getAssembler()->matrix());

                m_solver.factorize(getAssembler()->matrix());
                m_solution = m_solver.solve(getAssembler()->rhs());

                //------- positive variables -----
                //tmpSolution = tmpSolution.cwiseMax(1e-10);
                //------------------------------------------
                this->dispPicardSolChangeRelNorm(oldSol, m_solution);

                solRelNorm = this->solutionChangeRelNorm(oldSol, m_solution);
                iter++;
            }
        }
    }

    void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution);

        // in the first time step the pattern will get analyzed
        if (!m_iterationNumber)
            m_solver.analyzePattern(getAssembler()->matrix());

        gsMatrix<T> tmpSolution(getAssembler()->numVarDofs(), getAssembler()->getNumVar());

        m_solver.factorize(getAssembler()->matrix());
        tmpSolution = m_solver.solve(getAssembler()->rhs());

        //------- positive variables -----
        //tmpSolution = tmpSolution.cwiseMax(1e-10);
        //------------------------------------------

        this->dispSolChangeRelNorm(m_solution, tmpSolution);

        int iter = 0;
        T solRelNorm = this->solutionChangeRelNorm(m_solution, tmpSolution);

        int maxPicardIter = m_innerIter;
        gsInfo << "        Picard's iterations...\n";
        while ((math::abs(solRelNorm) > m_innerTol) && (iter < maxPicardIter))
        {
            gsInfo << "         ";
            getAssembler()->updatePicard(tmpSolution);

            gsMatrix<T> oldSol = tmpSolution;

            m_solver.factorize(getAssembler()->matrix());
            tmpSolution = m_solver.solve(getAssembler()->rhs());
       
            //------- positive variables -----
            //tmpSolution = tmpSolution.cwiseMax(1e-10);
            //------------------------------------------
            this->dispPicardSolChangeRelNorm(oldSol, tmpSolution);

            solRelNorm = this->solutionChangeRelNorm(oldSol, tmpSolution);
            iter++;
        }

        m_solution = tmpSolution;

        m_iterationNumber++;
        m_time += getAssembler()->getTimeStepSize();
    }

    void dispPicardSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     Picard's solution change relative norm: ";
        gsInfo << this->solutionChangeRelNorm(solOld, solNew) << ", ";
        gsInfo << "\n";
    }

public:
    virtual uwbADRAssemblerNonlinear<T>*  getAssembler() const { return dynamic_cast<uwbADRAssemblerNonlinear<T>*>(m_pAssembler); }

protected:
    bool m_bUnsteady;

    using Base::m_time;
    using Base::m_innerIter;
    using Base::m_innerTol;

    using uwbADRSolverBase<T>::m_solution;
    using uwbADRSolverBase<T>::m_pAssembler;
    using uwbADRSolverBase<T>::m_iterationNumber;
    using uwbADRSolverBase<T>::m_solver;

}; // class uwbADRSolverNonlinear

} // namespace gismo
