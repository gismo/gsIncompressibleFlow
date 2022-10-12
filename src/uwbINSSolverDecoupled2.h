/** @file uwbINSSolverDecoupled2.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSSolverDecoupledBase.h"
#include "uwbINSAssemblerDecoupled2.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupled2 : public uwbINSSolverDecoupledBase<T>
{

public:
    typedef uwbINSSolverDecoupledBase<T> Base;

public:
    uwbINSSolverDecoupled2(uwbINSSolverParams<T>& params) : Base(params)
    {
        //create assembler
        m_pAssembler = new uwbINSAssemblerDecoupled2<T>(params);
        
        initMembers();
    }

    virtual ~uwbINSSolverDecoupled2()
    {
    }

protected:
    virtual void initMembers()
    {
        uwbINSSolverBase<T>::initMembers();

        m_solution1.setZero(getAssembler()->getPshift(), 1);
        m_solution2.setZero(getAssembler()->getPdofs(), 1);
        m_solution3.setZero(getAssembler()->getPdofs(), 1);

        #ifdef GISMO_WITH_PARDISO
        uwbINSSolverBase<T>::pardisoSetup(m_solver1);
        uwbINSSolverBase<T>::pardisoSetup(m_solver2);
        uwbINSSolverBase<T>::pardisoSetup(m_solver3);
        #endif
    }

public:
    virtual void initialize()
    {
        getAssembler()->initialize();

        #pragma omp parallel sections 
        {
        #pragma omp section
        {
            m_solver2.analyzePattern(getAssembler()->matrix2());
            m_solver2.factorize(getAssembler()->matrix2());
        }
        #pragma omp section
        {
            m_solver3.analyzePattern(getAssembler()->matrix3());
            m_solver3.factorize(getAssembler()->matrix3());
        }
        }
    }

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution, m_solution2);

        int udofs = getAssembler()->getUdofs();

        if (!m_iterationNumber)
            m_solver1.analyzePattern(getAssembler()->matrix1().block(0, 0, udofs, udofs));
        m_solver1.factorize(getAssembler()->matrix1().block(0, 0, udofs, udofs));

        gsMatrix<T> rhs1 = getAssembler()->rhs1();

        if (getAssembler()->isRotation())
        {
#pragma omp parallel for
            for (int s = 0; s < getAssembler()->getTarDim(); s++)
            {
                for (int t = 0; t < getAssembler()->getTarDim(); t++)
                {
                    if (s != t)
                        rhs1.middleRows(s * udofs, udofs) -= getAssembler()->matrix1().block(s * udofs, t * udofs, udofs, udofs) * m_solution.middleRows(t * udofs, udofs);
                }
            }
        }

#pragma omp parallel for          
        for (int s = 0; s < getAssembler()->getTarDim(); s++) {
            m_solution1.middleRows(s * udofs, udofs) = m_solver1.solve(rhs1.middleRows(s * udofs, udofs));
        }

        if (getAssembler()->isRotation())
        {
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
                {
                    m_solution1.middleRows(s * udofs, udofs) = m_solver1.solve(rhs1.middleRows(s * udofs, udofs));
                }

                relNorm = (m_solution1 - m_solution1_old).norm() / m_solution1.norm();

                if ((relNorm < m_maxNorm) || (!m_iterationNumber))
                {
                    break;
                }
            }

            gsInfo << "iter sol1: " << iter + 1 << ", ";
        }

        getAssembler()->updateRhs2(m_solution1);
        m_solution2 = m_solver2.solve(getAssembler()->rhs2());

        getAssembler()->updateRhs3(m_solution1, m_solution2);
        m_solution3 = m_solver3.solve(getAssembler()->rhs3());

        getAssembler()->createSolVector_into(m_solution1, m_solution2, m_solution3, m_solution);

        m_iterationNumber++;
    }

    virtual uwbINSAssemblerDecoupled2<T>* getAssembler()
    {
        return dynamic_cast<uwbINSAssemblerDecoupled2<T>* >(m_pAssembler);
    }

protected:

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver1, m_solver2, m_solver3;
    #else
    typename gsSparseSolver<T>::LU m_solver1, m_solver2, m_solver3;
    #endif

    // members from uwbINSSolverDecoupled
    using Base::m_solution1;
    using Base::m_solution2;
    using Base::m_solution3;
    using Base::m_innerIterations;
    using Base::m_maxNorm;

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_pAssembler;
    using uwbINSSolverBase<T>::m_solution;
    using uwbINSSolverBase<T>::m_iterationNumber;

}; //uwbINSSolverDecoupled2

} //namespace gismo
