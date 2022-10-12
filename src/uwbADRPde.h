
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{

template<class T>
class uwbADRPde : public gsPde<T>
{

public:
    uwbADRPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc
        //const gsFunction<T>&           rhs,
        ) : gsPde<T>(domain, bc)//, m_rhs(rhs)
    { }

    ~uwbADRPde( ) { }

    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Solving advection-diffusion-reaction equation:\n";
        return os;
    }

    //const gsFunction<T>& getRhs() const { return m_rhs; }

//protected:
    //const gsFunction<T>&  m_rhs;
}; // class uwbADRPde

} // namespace gismo
