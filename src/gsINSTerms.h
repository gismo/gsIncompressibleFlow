/** @file gsINSTerms.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once
#include <gismo.h>
#include <gsIncompressibleFlow/src/gsFlowTerms.h>

namespace gismo
{

/// @brief      A class for integrals of the form: pressure shape function value * velocity test function divergence.
/// @tparam T   real number type
template <class T>
class gsINSTerm_PvalUdiv : public gsFlowTerm<T> // order: shape, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_PvalUdiv()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_VALUE;
    }


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData0)
    { this->setConstCoeff(-1.0); } // -1 to get block -Bt 


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: velocity shape function divergence * pressure test function value.
/// @tparam T   real number type
template <class T>
class gsINSTerm_UdivPval : public gsFlowTerm<T> // order: shape, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_UdivPval()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: velocity solution * shape function gradient * test function value.
/// @tparam T   real number type
template <class T>
class gsINSTerm_UsolGradVal : public gsFlowTermNonlin<T> // order: shape, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_UsolGradVal()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSTerms.hpp)
#endif