/** @file gsFlowTerms.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once
#include <gismo.h>

namespace gismo
{

/// @brief      A class computing individual terms of the weak formulation appearing in incompressible flow problems.
/// @tparam T   real number type
template <class T>
class gsFlowTerm
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowTerm> Ptr; 
    typedef memory::unique_ptr<gsFlowTerm> uPtr;


protected: // *** Class members ***

    unsigned m_geoFlags, m_testFunFlags, m_shapeFunFlags; // evaluation flags
    gsVector<T> m_coeff;


public: // *** Constructor/destructor ***

    gsFlowTerm()
    {
        m_geoFlags = 0;
        m_testFunFlags = 0;
        m_shapeFunFlags = 0;
    }

    virtual ~gsFlowTerm()
    {}


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
    { GISMO_NO_IMPLEMENTATION }

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
    { GISMO_NO_IMPLEMENTATION }

    void updateEvalFlags(unsigned& geoFlags, unsigned& testFunFlags, unsigned& shapeFunFlags)
    { 
        geoFlags |= m_geoFlags;
        testFunFlags |= m_testFunFlags;
        shapeFunFlags |= m_shapeFunFlags;
    }


protected: // *** Member functions ***

    void setConstCoeff(T value)
    {
        m_coeff.resize(1);
        m_coeff << value;
    }

    /**
     * @brief Evaluate the term coefficient.
     * 
     * The result is saved into m_coeff.
     * If the coefficient is constant, m_coeff is a vector of size 1.
     * If it is space-dependent, the size of m_coeff is equal to the number of evaluation points.
     * 
     * @param[in] mapData       geometry map information (including the evaluation points)
     */
    virtual void evalCoeff(const gsMapData<T>& mapData)
    { setConstCoeff(1.0); }

    virtual gsVector<T> getCoeffWeightsProduct(const gsVector<T>& quWeights);

};

// ===================================================================================================================

/// @brief      A class computing nonlinear terms of the weak formulation appearing in incompressible flow problems.
/// @tparam T   real number type
template <class T>
class gsFlowTermNonlin : public gsFlowTerm<T>
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowTermNonlin> Ptr; 
    typedef memory::unique_ptr<gsFlowTermNonlin> uPtr;

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;


protected: // *** Class members ***

    gsField<T> m_currentSolU;
    bool m_isCurrentSolSet;
    gsMatrix<T> m_solUVals;


public: // *** Constructor/destructor ***

    gsFlowTermNonlin()
    { }


public: // *** Member functions ***

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_currentSolU = solutions.front();
        m_isCurrentSolSet = true;
    }

    void setCurrentSolution(gsField<T>& solution)
    { 
        m_currentSolU = solution;
        m_isCurrentSolSet = true;
    }


protected: // *** Member functions ***

    virtual void computeCoeffSolU(const gsMapData<T>& mapData)
    { 
        GISMO_ASSERT(m_isCurrentSolSet, "No velocity solution set in gsFlowTermNonlin.");

        m_solUVals.resize(mapData.dim.first, mapData.points.cols());
        m_solUVals = m_currentSolU.value(mapData.points, mapData.patchId);
    }

};

// ===================================================================================================================
// ===================================================================================================================

/// @brief      A class for integrals of the form: test function value * shape function value.
/// @tparam T   real number type
template <class T>
class gsFlowTermValVal : public gsFlowTerm<T>
{

public:
    typedef gsFlowTerm<T> Base;


protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsFlowTermValVal()
    {
        Base::m_geoFlags = NEED_MEASURE;
        Base::m_testFunFlags = NEED_VALUE;
        Base::m_shapeFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: (1 / time step) * test function value * shape function value.
/// @tparam T   real number type
template <class T>
class gsFlowTermTimeDiscr : public gsFlowTermValVal<T>
{

public:
    typedef gsFlowTermValVal<T> Base;


protected: // *** Class members ***

    real_t m_timeStep;


public: // *** Constructor/destructor ***

    gsFlowTermTimeDiscr(real_t timeStep) :
    m_timeStep(timeStep)
    { }


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData0)
    { this->setConstCoeff(1./m_timeStep); }

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: test function gradient * shape function gradient.
/// @tparam T   real number type
template <class T>
class gsFlowTermGradGrad : public gsFlowTerm<T>
{

public:
    typedef gsFlowTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsFlowTermGradGrad()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_DERIV;
        Base::m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: viscosity * test function gradient * shape function gradient.
/// @tparam T   real number type
template <class T>
class gsFlowTermDiffusion : public gsFlowTermGradGrad<T>
{

public:
    typedef gsFlowTermGradGrad<T> Base;

protected: // *** Class members ***

    real_t m_viscosity;

public: // *** Constructor/destructor ***

    gsFlowTermDiffusion(real_t viscosity) :
    m_viscosity(viscosity)
    { }


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData)
    { this->setConstCoeff(m_viscosity); }

};

// ===================================================================================================================
// ===================================================================================================================

/// @brief      A class for integrals of the form: test function value * rhs function value.
/// @tparam T   real number type
template <class T>
class gsFlowTermRhs : public gsFlowTerm<T>
{

public:
    typedef gsFlowTerm<T> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;
    gsMatrix<T> m_rhsVals;


protected: // *** Base class members ***

    using Base::m_geoFlags;
    using Base::m_testFunFlags;


public: // *** Constructor/destructor ***

    gsFlowTermRhs(const gsFunction<T>* pRhsFun)
    {
        m_geoFlags = NEED_VALUE | NEED_MEASURE;
        m_testFunFlags = NEED_VALUE;
        m_pRhsFun = pRhsFun;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowTerms.hpp)
#endif