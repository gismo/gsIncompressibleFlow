/** @file uwbINSSolverDecoupled2PeriodicBase.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSSolverDecoupledBase.h"
#include "uwbINSAssemblerDecoupled2Periodic.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupled2PeriodicBase : public uwbINSSolverDecoupledBase<T>
{

public:
    typedef uwbINSSolverDecoupledBase<T> Base;

public:
    uwbINSSolverDecoupled2PeriodicBase(uwbINSSolverParams<T>& params) : Base(params)
    {
        //create assembler
        m_pAssembler = new uwbINSAssemblerDecoupled2Periodic<T>(params);
        
        initMembers();
    }

    virtual ~uwbINSSolverDecoupled2PeriodicBase()
    {
        for (size_t s = 0; s < m_solver1.size(); s++)
        {
            delete m_solver1[s];
            m_solver1[s] = NULL;
        }
    }

protected:
    virtual void initMembers()
    {   
        uwbINSSolverBase<T>::initMembers();

        m_solution1.setZero(getAssembler()->getPshift(), 1);
        m_solution2.setZero(getAssembler()->getPdofs(), 1);
        m_solution3.setZero(getAssembler()->getPdofs(), 1);
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
            m_solver3.front()->analyzePattern(getAssembler()->matrix3());
            m_solver3.front()->factorize(getAssembler()->matrix3());
        }
        }
    }

    virtual uwbINSAssemblerDecoupled2Periodic<T>* getAssembler()
    {
        return dynamic_cast<uwbINSAssemblerDecoupled2Periodic<T>*>(m_pAssembler);
    }


protected:
    #ifdef GISMO_WITH_PARDISO
    std::vector< typename gsSparseSolver<T>::PardisoLU * > m_solver1, m_solver3;
    typename gsSparseSolver<T>::PardisoLU m_solver2;
    #else
    std::vector< typename gsSparseSolver<T>::LU * > m_solver1, m_solver3;
    typename gsSparseSolver<T>::LU m_solver2;
    #endif 

    // members from uwbINSSolverDecoupled
    using Base::m_solution1;
    using Base::m_solution2;
    using Base::m_solution3;

    // members from uwbINSSolverBase
    using uwbINSSolverBase<T>::m_pAssembler;

}; //uwbINSSolverDecoupled2PeriodicBase

} //namespace gismo
