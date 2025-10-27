/** @file gsTMTerms.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once
#include <gismo.h>

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsTMModels.h>

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
    
namespace gismo
{

/// @brief      A class for integrals of the form: velocity * test function gradient * shape function gradient
/// @tparam T   real number type
template <class T>
class gsTMTerm_VecCoeffGradVal : public gsFlowTermNonlin<T>
{

public: // *** Type definitions ***

    typedef gsFlowTermNonlin<T> Base;


protected: // *** Base class members ***

    using Base::m_currentSolU;
    using Base::m_isCurrentSolSet;
    using Base::m_solUVals;

public: // *** Constructor/destructor ***

    gsTMTerm_VecCoeffGradVal()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===========================================================================================================================

/// @brief      A class for integrals of the form: coefficient * test function gradient * shape function gradient,
/// @brief      where coefficient = k1 * turbViscosity + k2
/// @tparam T   real number type
template <class T>
class gsTMTerm_CoeffGradGrad : public gsFlowTerm_GradGrad<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr; 
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsTMTerm_CoeffGradGrad(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

};

// ===========================================================================================================================

/// @brief      A class for integrals of the form: vector coefficient * test function value * shape function value
/// @tparam T   real number type
template <class T>
class gsTMTerm_CoeffValVal : public gsFlowTerm_ValVal<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_ValVal<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;
    
protected: // *** Base class members ***

    using Base::m_coeff;
    
public: // *** Constructor/destructor ***

    gsTMTerm_CoeffValVal(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: coefficient * test function value * shape function gradient.
/// @brief      where coefficient = - 2 * (1 - F1) / omega^2 grad(k) * grad(omega)
/// @tparam T   real number type
template <class T>
class gsTMTerm_BlendCoeff : public gsFlowTerm<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;
    
protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsTMTerm_BlendCoeff(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ========================================================================================================================

/// @brief      A class for integrals of the form: test function value * rhs function value.
/// @tparam T   real number type
template <class T>
class gsTMTerm_BlendCoeffRhs : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;
    
protected: // *** Base class members ***

    using Base::m_rhsVals;

public: // *** Constructor/destructor ***

    gsTMTerm_BlendCoeffRhs(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_VALUE | NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ========================================================================================================================

/// @brief      A class for integrals of the form: test function value * rhs function value.
/// @tparam T   real number type
template <class T>
class gsTMTerm_ProductionRhs : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;
    
protected: // *** Base class members ***

    using Base::m_rhsVals;

public: // *** Constructor/destructor ***

    gsTMTerm_ProductionRhs(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_VALUE | NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMTerms.hpp)
#endif