/** @file uwbINSSolverSteadyIterative.h

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSSolverIterativeBase.h"
#include "uwbINSSolverSteady.h"

namespace gismo
{

template<class T, class SolverType>
class uwbINSSolverSteadyIterative : public uwbINSSolverSteady<T>, public uwbINSSolverIterativeBase<T, SolverType>
{

public:
    typedef uwbINSSolverSteady<T> Base;
    typedef uwbINSSolverIterativeBase<T, SolverType> BaseIter;

public:
    uwbINSSolverSteadyIterative(uwbINSSolverParams<T>& params) : 
        Base(params), BaseIter(params, *this->getAssembler())
    { }

    virtual ~uwbINSSolverSteadyIterative()
    {
    }

public:
    virtual void setStokesSolution()
    {
        m_solution = BaseIter::getStokesSolution();
        m_iterationNumber = 0;
    }

    virtual void initIteration()
    {
        BaseIter::initIteration();
    }

    virtual void applySolver(gsMatrix<T>& solution)
    {
        BaseIter::applySolver(solution);
    }

    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
    {
        BaseIter::applySolver(solution, alpha_u, alpha_p);
    }

    virtual void initialize()
    {
        Base::initialize();
        BaseIter::initialize(m_solution);
    }

    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        Base::markDofsAsEliminatedZeros(boundaryDofs, unk);
        BaseIter::reinitMembers();
    }

    const T getAssemblyTime() const 
    { 
        return (Base::getAssemblyTime() + BaseIter::getAssemblyTime());
    }

    const T getSolverSetupTime() const 
    { 
        return (Base::getSolverSetupTime() + BaseIter::getSolverSetupTime());
    }

    const T getSolveTime() const 
    { 
        return (Base::getSolveTime() + BaseIter::getSolveTime());
    }

protected:

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_solution;
    using uwbINSSolverBase<T>::m_iterationNumber;

}; //uwbINSSolverSteadyIterative

} //namespace gismo
