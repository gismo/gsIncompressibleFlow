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

    GISMO_CLONE_FUNCTION(gsRANSTerm_SymmetricGradient)


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

    GISMO_CLONE_FUNCTION(gsRANSTerm_SymmetricGradient_full)


public: // *** Member functions ***

    void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

}; // gsRANSTerm_SymmetricGradient_full

// ===================================================================================
// For SUPG stabilization

template <class T>
class gsRANSTerm_SG_SUPGstabilization_diffusion : public gsRANSTerm_SymmetricGradient<T>
{

public:
    typedef gsRANSTerm_SymmetricGradient<T> Base;

protected: // *** Class members ***

    std::vector< gsMatrix<T> > m_turbViscosityGrads;
    gsMatrix<T> m_tauS;
    gsMatrix<T> m_solUVals;

protected: // *** Base class members ***

    using Base::m_viscosity;
    using Base::m_turbViscosityVals;

public: // *** Constructor/destructor ***

    gsRANSTerm_SG_SUPGstabilization_diffusion() :
    Base()
    {  }

    GISMO_CLONE_FUNCTION(gsRANSTerm_SG_SUPGstabilization_diffusion)


public: // *** Member functions ***

    void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

public: // *** Getter/setters

    void setTauS(gsMatrix<T> tauS) { m_tauS = tauS;}

    gsMatrix<T> getTauS() { return m_tauS; }

    void setUSolVals(gsMatrix<T> mat) { m_solUVals = mat; }

    void setTurbulentViscosityGrads(std::vector< gsMatrix<T> > turbViscosityGrads) { m_turbViscosityGrads = turbViscosityGrads; }

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: viscosity * test function gradient * trial function gradient.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsRANSTerm_SUPGstabilization_pressure : public gsFlowTermNonlin<T>
{

public:
    typedef gsFlowTermNonlin<T> Base;    

protected: // *** Base class members ***
    
    using Base::m_currentSolU;
    using Base::m_isCurrentSolSet;
    using Base::m_solUVals;

protected: // *** Class members ***

    gsMatrix<T> m_tauS;

public: // *** Constructor/destructor ***

    gsRANSTerm_SUPGstabilization_pressure()
    { 
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_VALUE;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    GISMO_CLONE_FUNCTION(gsRANSTerm_SUPGstabilization_pressure)


protected: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat);

public: // *** Getters/setters ***

    void setTauS(gsMatrix<T> tauS) { m_tauS = tauS;}

    gsMatrix<T> getTauS() { return m_tauS; }

    void setUSolVals(gsMatrix<T> mat) { m_solUVals = mat; }
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSTerms.hpp)
#endif