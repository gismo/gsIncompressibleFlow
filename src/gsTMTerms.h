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

/// @brief      A class for integrals of the form: vector coefficient * shape function gradient * test function value
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
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

public: // *** Getter/setters

    //void setTurbulentViscosityVals(gsVector<T> turbViscosityVals) { m_turbViscosityVals = turbViscosityVals; }

};

// ===========================================================================================================================

/// @brief      A class for integrals of the form: coefficient * test function gradient * shape function gradient,
/// @brief      where coefficient = k1 * turbViscosity + k2
/// @tparam T   real number type
template <class T>
class gsTMTerm_CoeffGradGrad : public gsFlowTerm<T>
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;    
    real_t m_k1;
    real_t m_k2;
    real_t m_k3;
    gsVector<T> m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsTMTerm_CoeffGradGrad(typename gsFlowSolverParams<T>::Ptr paramsPtr, real_t k1, real_t k2, real_t k3) :
    m_paramsPtr(paramsPtr), m_k1(k1), m_k2(k2), m_k3(k3)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData)

    //virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

public: // *** Getter/setters

    //void setTurbulentViscosityVals(gsVector<T> turbViscosityVals) { m_turbViscosityVals = turbViscosityVals; }

};

// ===========================================================================================================================

/// @brief      A class for integrals of the form: vector coefficient * shape function gradient * test function value
/// @tparam T   real number type
template <class T>
class gsTMTerm_CoeffValVal : public gsFlowTerm_ValVal<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_ValVal<T> Base;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    index_t m_unknown;

protected: // *** Base class members ***

    //using Base::m_isCurrentSolSet;
    using Base::m_coeff;
    
public: // *** Constructor/destructor ***

    gsTMTerm_CoeffValVal(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_unknown(unk)
    {
        this->m_geoFlags = NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
        this->m_shapeFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    //virtual void computeCoeffSolU(const gsMapData<T>& mapData)
    //{ 
    //    GISMO_ASSERT(m_isCurrentSolSet, "No solution set in gsTMTerm_CoeffValVal.");
    //
    //    m_solVals.resize(mapData.dim.first, mapData.points.cols());
    //    m_solVals = m_currentSol.value(mapData.points, mapData.patchId);
    //}

    //virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

public: // *** Getter/setters

    //void setCurrentSolution(gsField<T>& solution)
    //{ 
    //    m_currentSol = solution;
    //    m_isCurrentSolSet = true;
    //}

};

// ===================================================================================================================
// ===================================================================================================================

/// @brief      A class for integrals of the form: test function value * rhs function value.
/// @tparam T   real number type
template <class T>
class gsTMTerm_BlendCoeffRhs : public gsFlowTerm_rhs<T>
{

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;

protected: // *** Base class members ***

    using Base::m_rhsVals;

public: // *** Constructor/destructor ***

    gsTMTerm_BlendCoeffRhs(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    m_paramsPtr(paramsPtr)
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

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    index_t m_unknown;

protected: // *** Base class members ***

    using Base::m_rhsVals;

public: // *** Constructor/destructor ***

    gsTMTerm_ProductionRhs(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    m_paramsPtr(paramsPtr), m_unknown(unk)
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