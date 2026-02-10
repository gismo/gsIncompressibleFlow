/** @file gsINSTerms.h
    
    @brief Classes for individual terms in incompressible Navier-Stokes equations.

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

/// @brief      A class for integrals of the form: pressure trial function value * velocity test function divergence.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_PvalUdiv : public gsFlowTerm<T> // PU: trial, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_PvalUdiv()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_PvalUdiv)


protected: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData0)
    { this->setConstCoeff(-1.0); } // -1 to get block -Bt 


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: velocity trial function divergence * pressure test function value.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_UdivPval : public gsFlowTerm<T> // UP: trial, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_UdivPval()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_DERIV;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_UdivPval)

public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: velocity solution * trial function gradient * test function value.
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_UsolGradVal : public gsFlowTermNonlin<T> // order: trial, test
{

public: // *** Constructor/destructor ***

    gsINSTerm_UsolGradVal()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_DERIV;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_UsolGradVal)


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

};


// ===================================================================================================================
// For weak imposition of Dirichlet boundary conditions
// ===================================================================================================================

/// @brief      A class for integrals arising from convective term: (u.n)*u*v for u being a test function value, v shape function value and up being a velocity solution from previous Picard's iteration
/// @tparam T   real number type
template <class T>
class gsINSTerm_CoeffUvalUval_WeakDirichlet : public gsFlowTerm_ValVal<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_ValVal<T> Base;

protected: // *** Class members ***

    gsField<T> m_currentSolU;
    bool m_isCurrentSolSet;
    gsMatrix<T> m_solUVals;
    gsMatrix<T> m_vals;
        
protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTerm_CoeffUvalUval_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_CoeffUvalUval_WeakDirichlet)

protected: // *** Member functions ***

    virtual void computeCoeffSolU(const gsMapData<T>& mapData)
    { 
        GISMO_ASSERT(m_isCurrentSolSet, "No velocity solution set in gsTMTerm_RhsVal_WeakDirichlet.");

        m_solUVals.resize(mapData.dim.first, mapData.points.cols());
        m_solUVals = m_currentSolU.value(mapData.points, mapData.patchId);
    }        

public: // *** Member functions ***

    void setCurrentSolution(gsField<T>& solution)
    { 
        m_currentSolU = solution;
        m_isCurrentSolSet = true;
    }

    void setBndVals(gsMatrix<T> vals)
    {
        m_vals = vals;
    }

    virtual void evalCoeff(const gsMapData<T>& mapData);
};

/// @brief      A class for integrals arising from convective term: (u.n)*u_d*v for u_d being a prescribed Dirichlet condition value, v shape function value and up being a velocity solution from previous Picard's iteration
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsUVal_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Class members ***

    gsField<T> m_currentSolU;
    bool m_isCurrentSolSet;
    gsMatrix<T> m_solUVals;
    
protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsUVal_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsUVal_WeakDirichlet)

protected: // *** Member functions ***

    virtual void computeCoeffSolU(const gsMapData<T>& mapData)
    { 
        GISMO_ASSERT(m_isCurrentSolSet, "No velocity solution set in gsINSTerm_RhsUVal_WeakDirichlet.");

        m_solUVals.resize(mapData.dim.first, mapData.points.cols());
        m_solUVals = m_currentSolU.value(mapData.points, mapData.patchId);
    }

public: // *** Member functions ***

    void setCurrentSolution(gsField<T>& solution)
    { 
        m_currentSolU = solution;
        m_isCurrentSolSet = true;
    }

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals representing penalty for velocity: u * v for u being a test function value, v shape function value
/// @tparam T   real number type
template <class T>
class gsINSTerm_CoeffUvalUvalPenalty_WeakDirichlet : public gsFlowTerm<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

protected: // *** Class members ***

    std::vector<real_t> m_penalties;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTerm_CoeffUvalUvalPenalty_WeakDirichlet(std::vector<real_t> penalties) :
    m_penalties(penalties)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_CoeffUvalUvalPenalty_WeakDirichlet)

protected: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat);

};

/// @brief      A class for integrals representing penalty for velocity: u_d * v for u_d being a prescribed Dirichlet condition value, v shape function value
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsUValPenalty_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Class members ***

    std::vector<real_t> m_penalties;
        
protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsUValPenalty_WeakDirichlet(std::vector<real_t> penalties) :
    m_penalties(penalties)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsUValPenalty_WeakDirichlet)

public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals arising from diffusion term: u.grad(v).n for u being a shape function value, grad(v) being test function gradient and n an outward normal vector
/// @tparam T   real number type
template <class T>
class gsINSTerm_CoeffUvalUdiv_WeakDirichlet : public gsFlowTerm<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

protected: // *** Class members ***

    real_t m_viscosity;

public: // *** Constructor/destructor ***

    gsINSTerm_CoeffUvalUdiv_WeakDirichlet(real_t viscosity) :
    m_viscosity(viscosity)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE | NEED_DERIV;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_CoeffUvalUdiv_WeakDirichlet)

public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData)
    { this->setConstCoeff(m_viscosity); } 

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

/// @brief      A class for integrals arising from diffusion term: u_d.grad(v).n for u_d being a prescribed Dirichlet condition value, grad(v) being test function gradient and n an outward normal vector
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsUdiv_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

protected: // *** Class members ***

    real_t m_viscosity;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsUdiv_WeakDirichlet(real_t viscosity) :
    m_viscosity(viscosity)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsUdiv_WeakDirichlet)

public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals of the form: q * (u.n) for q being a pressure test function value, u velocity trial function value and n outward normal vector
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_PvalUval_WeakDirichlet : public gsFlowTerm<T> // PU: test, trial
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;    

public: // *** Constructor/destructor ***

    gsINSTerm_PvalUval_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_PvalUval_WeakDirichlet)

public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat);

};

/// @brief      A class for integrals of the form: q * (u_d.n) for q being a pressure test function value, u_d prescribed Dirichlet condition value and n outward normal vector
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsPvalU_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsPvalU_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsPvalU_WeakDirichlet)

public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/*
/// @brief      A class for integrals of the form: p * (v.n) for p being a pressure trial function value, v velocity test function value and n outward normal vector
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_UvalPval_WeakDirichlet : public gsFlowTerm<T> // UP: trial, test
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

public: // *** Constructor/destructor ***

    gsINSTerm_UvalPval_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat);

};
*/

/// @brief      A class for integrals of the form: p_d * (v.n) for p_d being a prescribed Dirichlet condition value, v velocity test function value and n outward normal vector
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsUvalP_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsUvalP_WeakDirichlet()
    {
        this->m_geoFlags = NEED_MEASURE | NEED_OUTER_NORMAL;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsUvalP_WeakDirichlet)

public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

// ===================================================================================================================

/// @brief      A class for integrals representing penalty for pressure: p * q for q being a test function value, p shape function value
/// @tparam T   real number type
/// @ingroup IncompressibleFlow
template <class T>
class gsINSTerm_PvalPval_WeakDirichlet : public gsFlowTerm<T> // UP: trial, test
{

public: // *** Type definitions ***

    typedef gsFlowTerm<T> Base;

protected: // *** Class members ***

    std::vector<real_t> m_penalties;

public: // *** Constructor/destructor ***

    gsINSTerm_PvalPval_WeakDirichlet(std::vector<real_t> penalties) :
    m_penalties(penalties)
    {
        this->m_geoFlags = NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
        this->m_trialFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_PvalPval_WeakDirichlet)

public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat);

};

/// @brief      A class for integrals representing penalty for pressure: p_d * q for q being a test function value, p prescribed Dirichlet condition value
/// @tparam T   real number type
template <class T>
class gsINSTerm_RhsPvalP_WeakDirichlet : public gsFlowTerm_rhs<T>
{

public: // *** Type definitions ***

    typedef gsFlowTerm_rhs<T> Base;

protected: // *** Class members ***

    std::vector<real_t> m_penalties;
    
protected: // *** Base class members ***

    using Base::m_rhsVals;
    using Base::m_vals;

public: // *** Constructor/destructor ***

    gsINSTerm_RhsPvalP_WeakDirichlet(std::vector<real_t> penalties) :
    m_penalties(penalties)
    {
        this->m_geoFlags = NEED_MEASURE;
        this->m_testFunFlags = NEED_VALUE;
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_RhsPvalP_WeakDirichlet)

public: // *** Member functions ***

    virtual void evalCoeff(const gsMapData<T>& mapData);

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSTerms.hpp)
#endif