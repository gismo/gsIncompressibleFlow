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
class gsINSTermPvalUdiv : public gsFlowTerm<T> // order: shape, test
{

public:
    typedef gsFlowTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermPvalUdiv()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_DERIV;
        Base::m_shapeFunFlags = NEED_VALUE;
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
class gsINSTermUdivPval : public gsFlowTerm<T> // order: shape, test
{

public:
    typedef gsFlowTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermUdivPval()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_VALUE;
        Base::m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: velocity solution * shape function gradient * test function value.
/// @tparam T   real number type
template <class T>
class gsINSTermUsolGradVal : public gsFlowTermNonlin<T> // order: shape, test
{

public:
    typedef gsFlowTermNonlin<T> Base;

protected: // *** Base class members ***

    using Base::m_solUVals;
    using gsFlowTerm<T>::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermUsolGradVal()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_VALUE;
        Base::m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSTerms.hpp)
#endif