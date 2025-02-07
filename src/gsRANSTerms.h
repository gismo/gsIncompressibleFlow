/** @file gsFlowTerms.h
    
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
    

/// @brief      A class for integrals of the form: test function gradient * shape function gradient.
/// @tparam T   real number type
template <class T>
class gsRANSTerm_SymmetricGradientDiag : public gsFlowTerm<T>
{

protected: // *** Class members ***

    real_t m_viscosity;
    gsVector<T> m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsRANSTerm_SymmetricGradientDiag(real_t viscosity, gsVector<T> turbViscVals)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_DERIV;

        m_viscosity = viscosity;
        m_turbViscosityVals = turbViscVals;
    }


public: // *** Member functions ***

    //virtual void gsRANSTerm_Diffusion<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMatA, std::vector<gsMatrix<T> > localMatDiag);
    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===========================================================================================

/// @brief      A class for integrals of the form: test function gradient * shape function gradient.
/// @tparam T   real number type
template <class T>
class gsRANSTerm_SymmetricGradientOffdiag : public gsFlowTerm<T>
{

protected: // *** Class members ***

    real_t m_viscosity;
    gsVector<T> m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsRANSTerm_SymmetricGradientOffdiag(real_t viscosity, gsVector<T> turbViscVals)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_DERIV;

        m_viscosity = viscosity;
        m_turbViscosityVals = turbViscVals;
    }


public: // *** Member functions ***

    //virtual void gsRANSTerm_Diffusion<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMatA, std::vector<gsMatrix<T> > localMatDiag);
    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};



} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSTerms.hpp)
#endif