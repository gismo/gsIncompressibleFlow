/** @file gsRANSTerms.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once
#include <gismo.h>

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
    
namespace gismo
{

/// @brief      A class for integrals arising from symmetric gradient in RANS equations.
/// @tparam T   real number type
template <class T>
class gsRANSTerm_SymmetricGradient : public gsFlowTerm<T>
{

protected: // *** Class members ***

    real_t m_viscosity;
    gsVector<T> m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsRANSTerm_SymmetricGradient()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

public: // *** Getter/setters

    void setTurbulentViscosityVals(gsVector<T> turbViscosityVals) { m_turbViscosityVals = turbViscosityVals; }

    void setViscosity(real_t viscosity) { m_viscosity = viscosity; }

};    

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSTerms.hpp)
#endif