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

/**
 * @brief   A class for integrals arising from symmetric gradient in RANS equations.
 * 
 * Only upper off-diagonal blocks are assembled, the lower off-diagonal blocks are taken as transpositions of the corresponding upper off-diagonal ones.
 * This is suitable for element-by-element assembly.
 * 
 * @tparam T   real number type
 */
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
        this->m_trialFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

public: // *** Getter/setters

    void setTurbulentViscosityVals(gsVector<T> turbViscosityVals) { m_turbViscosityVals = turbViscosityVals; }

    void setViscosity(real_t viscosity) { m_viscosity = viscosity; }

}; // gsRANSTerm_SymmetricGradient


/**
 * @brief   A class for integrals arising from symmetric gradient in RANS equations.
 * 
 * All off-diagonal blocks are assembled, This is suitable for row-by-row assembly.
 * 
 * @tparam T   real number type
 */

template <class T>
class gsRANSTerm_SymmetricGradient_full : public gsRANSTerm_SymmetricGradient<T>
{

public:
    typedef gsRANSTerm_SymmetricGradient<T> Base;

protected: // *** Base class members ***

    using Base::m_viscosity;
    using Base::m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsRANSTerm_SymmetricGradient_full():
    Base()
    { }


public: // *** Member functions ***

    void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

}; // gsRANSTerm_SymmetricGradient_full

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSTerms.hpp)
#endif