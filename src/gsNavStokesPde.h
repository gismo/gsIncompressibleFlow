/** @file gsNavStokesPde.h

    @brief Describes the incompressible Navier-Stokes equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

#include <gsPde/gsStokesPde.h>

namespace gismo
{

/** @brief
    The incompressible Navier-Stokes PDE.

    This class describes the Navier-Stokes PDE, with an arbitrary right-hand side
    function.

    @tparam T coefficient type

    @ingroup Pde
    @ingroup pdeclass
    @ingroup IncompressibleFlow
 */

template<class T>
class gsNavStokesPde : public gsStokesPde<T>
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsNavStokesPde> Ptr;
    typedef memory::unique_ptr<gsNavStokesPde> uPtr;

protected:

    typedef gsStokesPde<T> Base;
    
protected: // *** Base class members ***

    using Base::m_viscosity;
    using Base::m_force;
    using Base::m_source;


public: // *** Constructor/destructor ***

    gsNavStokesPde() { }

    /// Copy constructor.
    gsNavStokesPde(const gsNavStokesPde & other) :
    gsNavStokesPde(other.m_domain, other.m_boundary_conditions, other.m_force, other.m_source, other.m_viscosity)
    { }

    /// @brief Constructor.
    /// @param[in] domain       multipatch computational domain
    /// @param[in] bc           boundary conditions
    /// @param[in] force        right-hand side function for momentum equations
    /// @param[in] viscosity    viscosity
    gsNavStokesPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc,
        const gsFunction<T>*           force,
        const T                        viscosity)
        : gsStokesPde<T>(domain, bc, force, NULL, viscosity)
    { }

    /// @brief Constructor.
    /// @param[in] domain       multipatch computational domain
    /// @param[in] bc           boundary conditions
    /// @param[in] force        right-hand side function for momentum equations
    /// @param[in] source       right-hand side function for continuity equation
    /// @param[in] viscosity    viscosity
    gsNavStokesPde(
        const gsMultiPatch<T>&         domain,
        const gsBoundaryConditions<T>& bc,
        const gsFunction<T>*           force,
        const gsFunction<T>*           source,
        const T                        viscosity)
        : gsStokesPde<T>(domain, bc, force, source, viscosity)
    { }

    ~gsNavStokesPde( )
    { 
    }

public: // *** Member functions ***

    /// Prints the object as a string.
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
