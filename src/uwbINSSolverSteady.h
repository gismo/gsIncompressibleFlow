/** @file uwbINSSolverSteady.h

Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include "uwbINSSolverBase.h"
#include "uwbINSAssemblerSteady.h"
//#include "uwbINSSUPGAssemblerSteady.h"
#include "uwbINSAssemblerSteadyPeriodic.h"
#include <gsTrilinos/gsTrilinos.h>

namespace gismo
{

template<class T>
class uwbINSSolverSteady : public uwbINSSolverBase<T>
{

public:
    typedef uwbINSSolverBase<T> Base;

public:
    uwbINSSolverSteady(uwbINSSolverParams<T>& params)
    {
        //create assembler
        if (params.getBCs().numPeriodic())
            m_pAssembler = new uwbINSAssemblerSteadyPeriodic<T>(params);
        else
            m_pAssembler = new uwbINSAssemblerSteady<T>(params);

        Base::initMembers();

        m_alpha_u = params.settings().get(constantsINS::alpha_u);
        m_alpha_p = params.settings().get(constantsINS::alpha_p);
    }

    virtual ~uwbINSSolverSteady()
    {
    }

public:
    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        this->updateAssembler();

        // in the first iteration the pattern will get analyzed
        if (!m_iterationNumber)
            this->initIteration();

        this->applySolver(m_solution);
        //this->applySolver(m_solution, m_alpha_u, m_alpha_p);

        m_iterationNumber++;
    }

    virtual T residualRelNorm() const
    {
        gsMatrix<T> residual = getAssembler()->rhs() - getAssembler()->matrix() * m_solution;
        T resNorm = residual.norm() / getAssembler()->rhs().norm();

        return resNorm;
    }

    virtual uwbINSAssemblerSteady<T>* getAssembler() const
    {
        uwbINSAssemblerSteadyPeriodic<T>* pAssembler = dynamic_cast<uwbINSAssemblerSteadyPeriodic<T>*>(m_pAssembler);

        if (pAssembler != NULL)
            return pAssembler;
        else
            return dynamic_cast<uwbINSAssemblerSteady<T>*>(m_pAssembler);
    }

    void evalResiduum(const gsMatrix<T> & solVector, std::vector<T> & residuum)
    {
        gsInfo << "Evaluating residuum of steady N-S problem...\n";
        getAssembler()->evalResiduum(solVector, residuum);
    }

    void evalResiduum(std::vector<T> & residuum)
    {
        gsInfo << "Evaluating residuum of steady N-S problem...\n";
        getAssembler()->evalResiduum(m_solution, residuum);
    }

protected:
    real_t m_alpha_u, m_alpha_p;

    // members from uwbINSSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;

}; //uwbINSSolverSteady

} //namespace gismo
