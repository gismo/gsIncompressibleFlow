/** @file uwbINSSolverDecoupledBase.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#include "uwbINSSolverBase.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupledBase : public uwbINSSolverBase<T>
{

public:
    typedef uwbINSSolverBase<T> Base;

public:
    uwbINSSolverDecoupledBase(uwbINSSolverParams<T>& params)
    {
        m_innerIterations = params.settings().get(constantsINS::dec_innerIt);
        m_maxNorm = params.settings().get(constantsINS::dec_innerTol);
    }

    virtual ~uwbINSSolverDecoupledBase()
    {
    }

    void changeRelaxParametres(const T alpha_u, const T alpha_p)
    { this->getAssembler()->changeRelaxParametres(alpha_u, alpha_p); }

    void setInnerIterations(int iter) { m_innerIterations = iter; }
    void setMaxNorm(int maxNorm) { m_maxNorm = maxNorm; }

protected:

    int m_innerIterations;
    T m_maxNorm;
    gsMatrix<T> m_solution1, m_solution2, m_solution3;

}; //uwbINSSolverDecoupledBase

} //namespace gismo