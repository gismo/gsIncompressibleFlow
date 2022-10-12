/** @file uwbINSSolverDecoupled1PeriodicCoupled.h

Author(s): H. Hornikova,  J. Sourek
*/

#pragma once

#include "uwbINSSolverDecoupled1PeriodicBase.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupled1PeriodicCoupled : public uwbINSSolverDecoupled1PeriodicBase<T>
{

public:
    typedef uwbINSSolverDecoupled1PeriodicBase<T> Base;

public:
    uwbINSSolverDecoupled1PeriodicCoupled(uwbINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSSolverDecoupled1PeriodicCoupled()
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

        for (int s = 0; s < getAssembler()->getTarDim() - 1; s++)
        {
            #ifdef GISMO_WITH_PARDISO
            m_solver3.push_back(new typename gsSparseSolver<T>::PardisoLU);
            uwbINSSolverBase<T>::pardisoSetup(*m_solver3[s]);
            #else
            m_solver3.push_back(new typename gsSparseSolver<T>::LU);
            #endif 
        }
    }

    virtual void reinitMembers()
    {
        Base::initMembers();
        initMembers();
    }

public:
    virtual void initialize()
    {
        Base::initialize();

        int udofs = getAssembler()->getUdofs();

        switch (getAssembler()->getTarDim())
        {
        case 2:
            m_solver3[0]->analyzePattern(getAssembler()->matrix3());
            m_solver3[0]->factorize(getAssembler()->matrix3());
            break;
        case 3:
#pragma omp parallel sections 
        {
#pragma omp section
        {
            m_solver3[0]->analyzePattern(getAssembler()->matrix3().block(0, 0, udofs, udofs));
            m_solver3[0]->factorize(getAssembler()->matrix3().block(0, 0, udofs, udofs));
        }
#pragma omp section
        {
            m_solver3[1]->analyzePattern(getAssembler()->matrix3().block(udofs, udofs, 2 * udofs, 2 * udofs));
            m_solver3[1]->factorize(getAssembler()->matrix3().block(udofs, udofs, 2 * udofs, 2 * udofs));
        }
        }
        break;
        default:
            GISMO_ERROR("Periodic conditions implemented only for 2D and 3D");
            break;
        }
    }

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution);

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

        getAssembler()->updateRhs3(m_solution2);

        switch (tarDim)
        {
        case 2:
            m_solution3 = m_solver3[0]->solve(getAssembler()->rhs3());
            break;
        case 3:
#pragma omp parallel sections 
        {
#pragma omp section
        {
            m_solution3.middleRows(0, udofs) = m_solver3[0]->solve(getAssembler()->rhs3().middleRows(0, udofs));
        }
#pragma omp section
        {
            m_solution3.middleRows(udofs, 2 * udofs) = m_solver3[1]->solve(getAssembler()->rhs3().middleRows(udofs, 2 * udofs));
        }
        }
        break;
        default:
            break;
        }

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

    // members from uwbINSSolverDecoupled1PeriodicBase
    using Base::m_solver1;
    using Base::m_solver2;
    using Base::m_solver3;

public:
    // member function from uwbINSSolverDecoupled1PeriodicBase
    using Base::getAssembler;

}; //uwbINSSolverDecoupled1PeriodicCoupled

} //namespace gismo
