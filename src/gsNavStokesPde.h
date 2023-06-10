/** @file gsNavStokesPde.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/

#pragma once

#include <gsPde/gsStokesPde.h>

namespace gismo
{

/** @brief
    An incompressible Navier-Stokes PDE.

    This class describes a Navier-Stokes PDE, with an arbitrary right-hand side
    function.

    @tparam T coefficient type

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsNavStokesPde : public gsStokesPde<T>
{

protected:

    typedef gsStokesPde<T> Base;
    gsNavStokesPde() { }

protected: // *** Base class members ***

    using Base::m_viscosity;
    using Base::m_force;
    using Base::m_source;


public: // *** Constructor/destructor ***

    gsNavStokesPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc,
        gsFunction<T>*                 force,
        const T                        viscosity)
        : gsStokesPde<T>(domain, bc, force, NULL, viscosity)
    { }

    gsNavStokesPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc,
        gsFunction<T>*                 force,
        gsFunction<T>*                 source,
        const T                        viscosity)
        : gsStokesPde<T>(domain, bc, force, source, viscosity)
    { }

    ~gsNavStokesPde( )
    { 
    }

public: // *** Member functions ***

    /// @brief Print a short description of the PDE.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << "Incompressible Navier-Stokes equation:\n"
           <<"u\u00B7\u2207u-\u03BD\u0394u-\u2207p = f,\n"
           <<" \u2207\u00B7u=0\n"
           <<"with:\n";
        os << "viscosity = " << m_viscosity << ".\n";
        if ( m_force )
        os<<"Force  function f = "<< *m_force <<".\n";
        if ( m_source )
        os<<"Source function g = "<< *m_source <<".\n";
        return os;
    }


}; // class gsNavStokesPde

} // namespace gismo
