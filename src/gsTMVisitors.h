/** @file gsTMVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsTMTerms.h>
#include <gsIncompressibleFlow/src/gsTMModels.h>

namespace gismo
{

template <class T, int MatOrder>
class gsTMVisitorLinearSST : public gsFlowVisitor<T, MatOrder>
{
    
public:
    typedef gsFlowVisitor<T, MatOrder> Base;
    
protected: // *** Class members ***
    
    index_t m_unknown;
    
protected: // *** Base class members ***
    
    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_geoFlags;

    
public: // *** Constructor/destructor ***
    
    gsTMVisitorLinearSST() {}
    
    gsTMVisitorLinearSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    Base(paramsPtr), m_unknown(unk)
    { }
    
    
protected: // *** Member functions ***

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown;
        m_trialUnkID = m_unknown;
    }
    
    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );
    }

public: // *** Member functions *** 

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
    
};

// ==============================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorTimeIterationSST : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    typename gsTMModelData<T>::Ptr m_TMModelPtr;
    index_t m_unknown;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_geoFlags;

public: // *** Constructor/destructor ***

    gsTMVisitorTimeIterationSST() {}

    gsTMVisitorTimeIterationSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr, index_t unk) :
    Base(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    { }

    /// @brief Copy constructor.
    gsTMVisitorTimeIterationSST(gsTMVisitorTimeIterationSST<T, MatOrder> const & other) :
    gsTMVisitorTimeIterationSST()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();
    }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        // production term going to rhs
        m_terms.push_back( new gsTMTerm_ProductionRhs<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
    
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown;
        m_trialUnkID = m_unknown;
    }

public: // *** Member functions *** 

    void evaluate(index_t testFunID);
    
    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void assemble();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ====================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorNonlinearSST : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:

    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    
    typename gsTMModelData<T>::Ptr m_TMModelPtr;
    index_t m_unknown;
    real_t m_konst1, m_konst2, m_konst3;
    real_t m_konst;
    index_t m_numLhsTerms;
    index_t m_numRhsTerms;
    gsMatrix<T> m_solution;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_geoFlags;

public: // *** Constructor/destructor ***

    gsTMVisitorNonlinearSST() {}

    gsTMVisitorNonlinearSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr, index_t unk) :
    Base(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    { }

    /// @brief Copy constructor.
    gsTMVisitorNonlinearSST(gsTMVisitorNonlinearSST<T, MatOrder> const & other) :
    gsTMVisitorNonlinearSST()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();
    }

protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_numLhsTerms = 4;
        m_numRhsTerms = 0;
        
        // further, define all lhs terms first and all rhs terms then
        // LHS
        // diffusion term with coefficient (m_konst1 * F1 + m_konst2 * (1 - F1)) * turbulent viscosity + m_konst3
        m_terms.push_back( new gsTMTerm_CoeffGradGrad<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        // advection terrm
        m_terms.push_back( new gsTMTerm_VecCoeffGradVal<T>() );
        // nonlinear reaction term
        m_terms.push_back( new gsTMTerm_CoeffValVal<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        // blended term 2 * (1 - F1) * sigma0mega2 / omega * grad(k) * grad(omega) for omega equation
        m_terms.push_back( new gsTMTerm_BlendCoeff<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        
        // RHS
        //m_terms.push_back( new gsTMTerm_BlendCoeffRhs<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown; 
        m_trialUnkID = m_unknown;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void assemble();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setCurrentSolution(const gsMatrix<T>& solution) { m_solution = solution; }

};


// ================================================================================================================
// For T-CSD stabilization
//

template <class T, int MatOrder>
class gsTMVisitorSSTTCSDStabilization_time : public gsTMVisitorLinearSST<T, MatOrder>
{

public:
    typedef gsTMVisitorLinearSST<T, MatOrder> Base;

public:
    gsField<T> m_solution;
    real_t m_viscosity;
    gsVector<T> m_TurbulentViscosityVals;
    gsMatrix<T> m_tauS;

    typename gsTMModelData<T>::Ptr m_TMModelPtr;

protected: // *** Base class members ***

    using Base::m_unknown;
    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsTMVisitorSSTTCSDStabilization_time() {}

    gsTMVisitorSSTTCSDStabilization_time(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr, index_t unk) :
    Base(paramsPtr, unk), m_TMModelPtr(TMModelPtr)
    { 
        initMembers();
    }

    /// @brief Copy constructor.
    gsTMVisitorSSTTCSDStabilization_time(gsTMVisitorSSTTCSDStabilization_time<T, MatOrder> const & other) :
    gsTMVisitorSSTTCSDStabilization_time()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TCSDStabilization_time<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown;    // velocity
        m_trialUnkID = m_unknown;   // velocity
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setRANSsolution(gsField<T>& solution) { m_solution = solution;}

};

template <class T, int MatOrder>
class gsTMVisitorSSTTCSDStabilization_advection : public gsTMVisitorLinearSST<T, MatOrder>
{

public:
    typedef gsTMVisitorLinearSST<T, MatOrder> Base;

public:
    gsField<T> m_solution;
    real_t m_viscosity;
    gsVector<T> m_TurbulentViscosityVals;
    gsMatrix<T> m_tauS;

    typename gsTMModelData<T>::Ptr m_TMModelPtr;

protected: // *** Base class members ***

    using Base::m_unknown;
    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsTMVisitorSSTTCSDStabilization_advection() {}

    gsTMVisitorSSTTCSDStabilization_advection(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr, index_t unk) :
    Base(paramsPtr, unk), m_TMModelPtr(TMModelPtr)
    { 
        initMembers();
    }

    /// @brief Copy constructor.
    gsTMVisitorSSTTCSDStabilization_advection(gsTMVisitorSSTTCSDStabilization_advection<T, MatOrder> const & other) :
    gsTMVisitorSSTTCSDStabilization_advection()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of TM model
        if (other.m_TMModelPtr)
            m_TMModelPtr = typename gsTMModelData<T>::Ptr(other.m_TMModelPtr->clone().release());

        // create new terms
        // terms cannot be cloned here, because they store m_TMModelPtr
        m_terms.clear();
        defineTerms();

    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TCSDStabilization_advection<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown;    // velocity
        m_trialUnkID = m_unknown;   // velocity
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setRANSsolution(gsField<T>& solution) { m_solution = solution;}

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMVisitors.hpp)
#endif