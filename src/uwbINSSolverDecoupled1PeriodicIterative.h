/** @file uwbINSSolverDecoupled1PeriodicIterative.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSSolverDecoupled1PeriodicBase.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupled1PeriodicIterative : public uwbINSSolverDecoupled1PeriodicBase<T>
{

public:
    typedef uwbINSSolverDecoupled1PeriodicBase<T> Base;

public:
    uwbINSSolverDecoupled1PeriodicIterative(uwbINSSolverParams<T>& params) : Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSSolverDecoupled1PeriodicIterative()
    {
    }

protected:
    void initMembers()
    {
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
        {
            #ifdef GISMO_WITH_PARDISO
            m_solver1.push_back(new typename gsSparseSolver<T>::PardisoLU);
            uwbINSSolverBase<T>::pardisoSetup(*m_solver1[s]);
            #else
            m_solver1.push_back(new typename gsSparseSolver<T>::LU);
            #endif 
        }

        for (int s = 0; s < getAssembler()->getTarDim(); s++)
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

#pragma omp parallel for
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
        {
            m_solver3[s]->analyzePattern(getAssembler()->matrix3().block(s * udofs, s * udofs, udofs, udofs));
            m_solver3[s]->factorize(getAssembler()->matrix3().block(s * udofs, s * udofs, udofs, udofs));
        }
    }

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution);

        int udofs = getAssembler()->getUdofs();

        if (!m_iterationNumber)
        {
#pragma omp parallel for 
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
                m_solver1[s]->analyzePattern(getAssembler()->matrix1().block(s * udofs, s * udofs, udofs, udofs));
        }

#pragma omp parallel for
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
            m_solver1[s]->factorize(getAssembler()->matrix1().block(s * udofs, s * udofs, udofs, udofs));

        gsMatrix<T> rhs1 = getAssembler()->rhs1();
#pragma omp parallel for
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
        {
            for (int t = 0; t < getAssembler()->getTarDim(); t++)
            {
                if (s != t)
                    rhs1.middleRows(s * udofs, udofs) -= getAssembler()->matrix1().block(s * udofs, t * udofs, udofs, udofs) * m_solution.middleRows(t * udofs, udofs);
            }
        }

#pragma omp parallel for          
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
            m_solution1.middleRows(s * udofs, udofs) = m_solver1[s]->solve(rhs1.middleRows(s * udofs, udofs));

        //inner iterations
        int iter = 0;
        T relNorm = 1;

        for (iter = 0; iter < m_innerIterations; iter++)
        {

            gsMatrix<T> m_solution1_old = m_solution1;

            rhs1 = getAssembler()->rhs1();

#pragma omp parallel for  
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
            {
                for (int t = 0; t < getAssembler()->getTarDim(); t++)
                {
                    if (s != t)
                        rhs1.middleRows(s * udofs, udofs) -= getAssembler()->matrix1().block(s * udofs, t * udofs, udofs, udofs) * m_solution1.middleRows(t * udofs, udofs);
                }
            }

#pragma omp parallel for          
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
                m_solution1.middleRows(s * udofs, udofs) = m_solver1[s]->solve(rhs1.middleRows(s * udofs, udofs));

            relNorm = (m_solution1 - m_solution1_old).norm() / m_solution1.norm();

            if ((relNorm < m_maxNorm) || (!m_iterationNumber))
                break;
        }

        gsInfo << "iter sol1: " << iter + 1 << ", ";

        getAssembler()->updateRhs2(m_solution1);
        m_solution2 = m_solver2.solve(getAssembler()->rhs2());

        getAssembler()->updateRhs3(m_solution2);

#pragma omp parallel for          
        for (int s = 0; s < getAssembler()->getTarDim(); s++)
            m_solution3.middleRows(s * udofs, udofs) = m_solver3[s]->solve(getAssembler()->rhs3().middleRows(s * udofs, udofs));

        iter = 0;
        relNorm = 1;

        for (iter = 0; iter < m_innerIterations; ++iter)
        {
            gsMatrix<T> m_solution3_old = m_solution3;

            gsMatrix<T> rhs3 = getAssembler()->rhs3();

#pragma omp parallel for  
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
            {
                for (int t = 0; t < getAssembler()->getTarDim(); t++)
                {
                    if (s != t)
                        rhs3.middleRows(s * udofs, udofs) -= getAssembler()->matrix3().block(s * udofs, t * udofs, udofs, udofs) * m_solution3.middleRows(t * udofs, udofs);
                }
            }

#pragma omp parallel for          
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
                m_solution3.middleRows(s * udofs, udofs) = m_solver3[s]->solve(rhs3.middleRows(s * udofs, udofs));

            relNorm = (m_solution3 - m_solution3_old).norm() / m_solution3_old.norm();

            if (relNorm < m_maxNorm)
                break;
        }

        gsInfo << "iter sol3: " << iter + 1 << ", ";

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
    using uwbINSSolverDecoupledBase<T>::m_innerIterations;
    using uwbINSSolverDecoupledBase<T>::m_maxNorm;
    
    // members from uwbINSSolverDecoupled1PeriodicBase
    using Base::m_solver1;
    using Base::m_solver2;
    using Base::m_solver3;

public:
    // member function from uwbINSSolverDecoupled1PeriodicBase
    using Base::getAssembler;

}; //uwbINSSolverDecoupled1PeriodicIterative

} //namespace gismo
