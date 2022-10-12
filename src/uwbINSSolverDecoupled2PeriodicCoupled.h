/** @file uwbINSSolverDecoupled2PeriodicCoupled.h

Author(s): H. Hornikova,  J. Sourek
*/

#pragma once

#include "uwbINSSolverDecoupled2PeriodicBase.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupled2PeriodicCoupled : public uwbINSSolverDecoupled2PeriodicBase<T>
{

public:
    typedef uwbINSSolverDecoupled2PeriodicBase<T> Base;

public:
    uwbINSSolverDecoupled2PeriodicCoupled(uwbINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSSolverDecoupled2PeriodicCoupled()
    { }

protected:
    void initMembers()
    {
        for (int s = 0; s < getAssembler()->getTarDim() - 1; s++)
        {
            #ifdef GISMO_WITH_PARDISO
            m_solver1.push_back(new typename gsSparseSolver<T>::PardisoLU);
            uwbINSSolverBase<T>::pardisoSetup(*m_solver1[s]);
            #else
            m_solver1.push_back(new typename gsSparseSolver<T>::LU);
            #endif 
        }
    }

    virtual void reinitMembers()
    {
        Base::initMembers();
        initMembers();
    }

public:

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution, m_solution2);

        int udofs = getAssembler()->getUdofs();
        short_t tarDim = getAssembler()->getTarDim();

        if (!m_iterationNumber)
        {
            switch (tarDim)
            {
            case 2:
                m_solver1[0]->analyzePattern(getAssembler()->matrix1());
                break;
            case 3:
#pragma omp parallel sections 
            {
#pragma omp section
            {
                m_solver1[0]->analyzePattern(getAssembler()->matrix1().block(0, 0, udofs, udofs));
            }
#pragma omp section
            {
                m_solver1[1]->analyzePattern(getAssembler()->matrix1().block(udofs, udofs, 2 * udofs, 2 * udofs));
            }
            }
            break;
            default:
                break;
            }
        }

        switch (tarDim)
        {
        case 2:
            m_solver1[0]->factorize(getAssembler()->matrix1());
            m_solution1 = m_solver1[0]->solve(getAssembler()->rhs1());
            break;
        case 3:
#pragma omp parallel sections 
        {
#pragma omp section
        {
            m_solver1[0]->factorize(getAssembler()->matrix1().block(0, 0, udofs, udofs));
            m_solution1.middleRows(0, udofs) = m_solver1[0]->solve(getAssembler()->rhs1().middleRows(0, udofs));
        }
#pragma omp section
        {
            m_solver1[1]->factorize(getAssembler()->matrix1().block(udofs, udofs, 2 * udofs, 2 * udofs));
            m_solution1.middleRows(udofs, 2 * udofs) = m_solver1[1]->solve(getAssembler()->rhs1().middleRows(udofs, 2 * udofs));
        }
        }
        break;
        default:
            break;
        }

        getAssembler()->updateRhs2(m_solution1);
        m_solution2 = m_solver2.solve(getAssembler()->rhs2());

        getAssembler()->updateRhs3(m_solution1, m_solution2);
        m_solution3 = m_solver3.front()->solve(getAssembler()->rhs3());

        getAssembler()->createSolVector_into(m_solution1, m_solution2, m_solution3, m_solution);

        m_iterationNumber++;
    }

protected:

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_solution;
    using uwbINSSolverBase<T>::m_iterationNumber;

    // members from uwbINSSolverDecoupled
    using uwbINSSolverDecoupledBase<T>::m_solution1;
    using uwbINSSolverDecoupledBase<T>::m_solution2;
    using uwbINSSolverDecoupledBase<T>::m_solution3;

    // members from uwbINSSolverDecoupled2PeriodicBase
    using Base::m_solver1;
    using Base::m_solver2;
    using Base::m_solver3;

public:
    // member function from uwbINSSolverDecoupled2PeriodicBase
    using Base::getAssembler;

}; //uwbINSSolverDecoupled2PeriodicCoupled

} //namespace gismo
