
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{

template<class T>
class uwbINSPde : public gsPde<T>
{

public:
    uwbINSPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc,
        const gsFunction<T>&           rhs,
        const T                        viscosity
        ) : gsPde<T>(domain, bc), m_rhs(rhs), m_viscosity(viscosity)
    {
    }

    ~uwbINSPde( )
    { 
    }

    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Incompressible Navier-Stokes equation with:\n";
        os << "Source function f= " << m_rhs << ".\n";
        return os;
    }

    const gsFunction<T>& getRhs() const { return m_rhs; }

    T getViscosity() const { return m_viscosity; }

protected:
    const gsFunction<T>&  m_rhs;
    T m_viscosity;
}; // class uwbINSPde

} // namespace gismo
