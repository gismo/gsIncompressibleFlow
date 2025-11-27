/** @file gsFlowTerms.h

    @brief Base classes for individual terms in incompressible flow problems.

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
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowTerm> Ptr; 
    typedef memory::unique_ptr<gsFlowTerm> uPtr;


protected: // *** Class members ***

    unsigned m_geoFlags, m_testFunFlags, m_trialFunFlags; // evaluation flags
    gsVector<T> m_coeff;


public: // *** Constructor/destructor ***

    /// Constructor.
    gsFlowTerm()
    {
        m_geoFlags = 0;
        m_testFunFlags = 0;
        m_trialFunFlags = 0;
    }

    virtual ~gsFlowTerm()
    {}

    GISMO_UPTR_FUNCTION_PURE(gsFlowTerm, clone)


public: // *** Member functions ***

    /// @brief Assemble the current local matrix.
    /// @param[in]  mapData         geometry mapping information
    /// @param[in]  quWeights       quadrature weights
    /// @param[in]  testFunData     test basis data (0 - values, 1 - derivatives, 2 - 2nd derivatives)
    /// @param[in]  trialFunData    trial basis data (0 - values, 1 - derivatives, 2 - 2nd derivatives)
    /// @param[out] localMat        resulting local matrix
    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Assemble the current local matrices.
    /// @param[in]  mapData         geometry mapping information
    /// @param[in]  quWeights       quadrature weights
    /// @param[in]  testFunData     test basis data (0 - values, 1 - derivatives, 2 - 2nd derivatives)
    /// @param[in]  trialFunData    trial basis data (0 - values, 1 - derivatives, 2 - 2nd derivatives)
    /// @param[out] localMat        vector of resulting local matrices (e.g. one for each velocity component)
    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Add evaluation flags of this term to the given evaluation flags.
    /// @param[out] geoFlags        flags for geometry evaluation
    /// @param[out] testFunFlags    flags for test basis evaluation
    /// @param[out] trialFunFlags   flags for trial basis evaluation
    void updateEvalFlags(unsigned& geoFlags, unsigned& testFunFlags, unsigned& trialFunFlags)
    { 
        geoFlags |= m_geoFlags;
        testFunFlags |= m_testFunFlags;
        trialFunFlags |= m_trialFunFlags;
    }


protected: // *** Member functions ***

    /// @brief Set the term coefficient to a constant value.
    /// @param[in] value the constant coefficient value
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

    virtual gsVector<T> getCoeffGeoMapProduct(const gsMapData<T>& mapData);

};

// ===================================================================================================================

/// @brief      A class computing nonlinear terms of the weak formulation appearing in incompressible flow problems.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
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

    /// Constructor.
    gsFlowTermNonlin()
    { }

    GISMO_CLONE_FUNCTION(gsFlowTermNonlin)

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

/// @brief      A class for integrals of the form: test function value * trial function value.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_ValVal : public gsFlowTerm<T>
{

public: // *** Constructor/destructor ***

    gsFlowTerm_ValVal()
    {
        this->m_geoFlags = NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsFlowTerm_ValVal)


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: (1 / time step) * test function value * trial function value.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_TimeDiscr : public gsFlowTerm_ValVal<T>
{

protected: // *** Class members ***

    real_t m_timeStep;


public: // *** Constructor/destructor ***

    gsFlowTerm_TimeDiscr(real_t timeStep) :
    m_timeStep(timeStep)
    { }

    GISMO_CLONE_FUNCTION(gsFlowTerm_TimeDiscr)


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData0)
    { this->setConstCoeff(1./m_timeStep); }

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: test function gradient * trial function gradient.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_GradGrad : public gsFlowTerm<T>
{

public: // *** Constructor/destructor ***

    gsFlowTerm_GradGrad()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    GISMO_CLONE_FUNCTION(gsFlowTerm_GradGrad)


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: viscosity * test function gradient * trial function gradient.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_Diffusion : public gsFlowTerm_GradGrad<T>
{

protected: // *** Class members ***

    real_t m_viscosity;

public: // *** Constructor/destructor ***

    gsFlowTerm_Diffusion(real_t viscosity) :
    m_viscosity(viscosity)
    { }

    GISMO_CLONE_FUNCTION(gsFlowTerm_Diffusion)


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData)
    { this->setConstCoeff(m_viscosity); }

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: viscosity * test function gradient * trial function gradient.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_TCSDStabilization_time : public gsFlowTermNonlin<T>
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

    gsFlowTerm_TCSDStabilization_time()
    { 
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_VALUE;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsFlowTerm_TCSDStabilization_time)


protected: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

public: // *** Getters/setters ***

    void setTauS(gsMatrix<T> tauS) { m_tauS = tauS;}

    gsMatrix<T> getTauS() { return m_tauS; }

    void setUSolVals(gsMatrix<T> mat) { m_solUVals = mat; }
};

// ===================================================================================================================

/// @brief      A class for integrals of the form: viscosity * test function gradient * trial function gradient.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_TCSDStabilization_advection : public gsFlowTermNonlin<T>
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

    gsFlowTerm_TCSDStabilization_advection()
    { 
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_VALUE;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    GISMO_CLONE_FUNCTION(gsFlowTerm_TCSDStabilization_advection)


protected: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

public: // *** Getters/setters ***

    void setTauS(gsMatrix<T> tauS) { m_tauS = tauS;}

    gsMatrix<T> getTauS() { return m_tauS; }

    void setUSolVals(gsMatrix<T> mat) { m_solUVals = mat; }
};

// ===================================================================================================================
// ===================================================================================================================

/// @brief      A class for integrals of the form: test function value * rhs function value.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsFlowTerm_rhs : public gsFlowTerm<T>
{

protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;
    gsMatrix<T> m_rhsVals;

public: // *** Constructor/destructor ***

    gsFlowTerm_rhs()
    {
        this->m_geoFlags = NEED_VALUE | NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
    }

    gsFlowTerm_rhs(const gsFunction<T>* pRhsFun)
    {
        this->m_geoFlags = NEED_VALUE | NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
        m_pRhsFun = pRhsFun;
    }

    GISMO_CLONE_FUNCTION(gsFlowTerm_rhs)


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowTerms.hpp)
#endif