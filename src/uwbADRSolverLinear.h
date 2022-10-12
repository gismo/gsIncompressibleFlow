/** @file uwbADRSolverLinear.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRSolverUnsteady.h"
#include "uwbADRAssemblerLinear.h"
#include "uwbADRSolverParams.h"

namespace gismo
{

template<class T>
class uwbADRSolverLinear : public uwbADRSolverUnsteady<T>
{

public:
    typedef uwbADRSolverUnsteady<T> Base;

public:
    //constant coefficients
    uwbADRSolverLinear(uwbADRSolverParams<T>& params,
                       T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff,
                       bool createAssembler = true)  : Base(params, false)
    {
        gsInfo << "Solving linear advection-diffusion-reaction equation with constant coefficients...\n";
        GISMO_ASSERT(!(params.getBCs().numPeriodic()), "Use interface for periodic boundary conditions!");
        GISMO_ASSERT(!(params.settings().get(constantsADR::CROSSWIND)), "Nonlinear solver is necessary to use for crosswind stabilization!");
        GISMO_ASSERT(!(params.settings().get(constantsADR::isoArtificialDiffusion)),
                     "Nonlinear solver is necessary to use for isotropic artificial diffusion stabilization!");

        params.settings().setADREvaluator("linConstCoeffs");

        if (createAssembler)
            m_pAssembler = new uwbADRAssemblerLinear<T>(params, diffusionCoeff, advectionCoeff, reactionCoeff);
        Base::initMembers();

        m_bUnsteady = params.settings().get(constantsADR::unsteady);
    }

    //non-constant coefficients
    uwbADRSolverLinear(uwbADRSolverParams<T>& params,
                       bool createAssembler = true)  : Base(params, false)
    {
        gsInfo << "Solving linear advection-diffusion-reaction equation with non-constant coefficients...\n";
        GISMO_ASSERT(!(params.getBCs().numPeriodic()), "Use interface for periodic boundary conditions!");
        GISMO_ASSERT(!(params.settings().get(constantsADR::CROSSWIND)), "Nonlinear solver is necessary to use for crosswind stabilization!");
        GISMO_ASSERT(!(params.settings().get(constantsADR::isoArtificialDiffusion)),
                     "Nonlinear solver is necessary to use for isotropic artificial diffusion stabilization!");

        params.settings().setADREvaluator("linNonConstCoeffs");

        if (createAssembler)
            m_pAssembler = new uwbADRAssemblerLinear<T>(params);
        Base::initMembers();

        m_bUnsteady = params.settings().get(constantsADR::unsteady);
    }

    virtual ~uwbADRSolverLinear()
    { }

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

            getAssembler()->updateLinearSteady();

            m_solver.analyzePattern(getAssembler()->matrix());
            m_solver.factorize(getAssembler()->matrix());
            m_solution = m_solver.solve(getAssembler()->rhs());
        }
    }

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->updateLinear(m_solution);

        // in the first time step the pattern will get analyzed
        if (!m_iterationNumber)
            m_solver.analyzePattern(getAssembler()->matrix());

        m_solver.factorize(getAssembler()->matrix());
        m_solution = m_solver.solve(getAssembler()->rhs());

        //------- positive variables -----
        //tmpSolution = tmpSolution.cwiseMax(1e-10);
        //------------------------------------------
        m_iterationNumber++;
        m_time += getAssembler()->getTimeStepSize();
    }

public:
    virtual uwbADRAssemblerLinear<T>*  getAssembler() const { return dynamic_cast<uwbADRAssemblerLinear<T>*>(m_pAssembler); }
    bool isInitialized() { return getAssembler()->isInitialized(); }

protected:

    bool m_bUnsteady;

    using Base::m_time;

    using uwbADRSolverBase<T>::m_solution;
    using uwbADRSolverBase<T>::m_pAssembler;
    using uwbADRSolverBase<T>::m_iterationNumber;
    using uwbADRSolverBase<T>::m_solver;

}; // class uwbADRSolverLinear

} //namespace gismo
