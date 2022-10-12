/** @file uwbINSSolverUnsteadyIterative.h

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSSolverIterativeBase.h"
#include "uwbINSSolverUnsteady.h"

namespace gismo
{

template<class T, class SolverType>
class uwbINSSolverUnsteadyIterative : public uwbINSSolverUnsteady<T>, public uwbINSSolverIterativeBase<T, SolverType>
{

public:
    typedef uwbINSSolverUnsteady<T> Base;
    typedef uwbINSSolverIterativeBase<T, SolverType> BaseIter;

public:

    uwbINSSolverUnsteadyIterative(uwbINSSolverParams<T>& params, bool createAssembler = true) :
        Base(params, createAssembler), BaseIter(params, *this->getAssembler())
    { }

    virtual ~uwbINSSolverUnsteadyIterative()
    {
    }

public:
    virtual void setStokesSolution()
    {
        m_solution = BaseIter::getStokesSolution();
        m_iterationNumber = 0;
    }

    void setStokesInitialCondition()
    {
        setStokesSolution();

        m_time = 0;
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

    // members from uwbINSSolverUnsteady
    using Base::m_avgPicardIter;
    using Base::m_time;

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_solution;
    using uwbINSSolverBase<T>::m_iterationNumber;

}; // class uwbINSSolverUnsteadyIterative

} // namespace gismo
