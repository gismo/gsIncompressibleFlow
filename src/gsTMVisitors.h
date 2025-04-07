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
    using Base::m_dofMappers;
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
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
    index_t m_unknown;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_dofMappers;
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

    gsTMVisitorTimeIterationSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    Base(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    { }


protected: // *** Member functions ***

    void evaluate(index_t testFunID);
    
    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void defineTerms()
    {
        m_terms.push_back( new gsTMTerm_VecCoeffGradVal<T>() );

        // production term going to rhs
        m_terms.push_back( new gsTMTerm_ProductionRhs<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
    
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown;
        m_trialUnkID = m_unknown;
    }

public: // *** Member functions *** 

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
    
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;
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
    using Base::m_dofMappers;
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

    gsTMVisitorNonlinearSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr, index_t unk) :
    Base(paramsPtr), m_TMModelPtr(TMModelPtr), m_unknown(unk)
    { 
                
    }

protected: // *** Member functions ***

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void defineTerms()
    {
        m_numLhsTerms = 3;
        m_numRhsTerms = 1;
        
        // further, define all lhs terms first and all rhs terms then
        // diffusion term with coefficient (m_konst1 * F1 + m_konst2 * (1 - F1)) * turbulent viscosity + m_konst3
        m_terms.push_back( new gsTMTerm_CoeffGradGrad<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        // nonlinear reaction term
        m_terms.push_back( new gsTMTerm_CoeffValVal<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        // blended term 2 * (1 - F1) * sigma0mega2 / omega * grad(k) * grad(omega) for omega equation
        m_terms.push_back( new gsTMTerm_BlendCoeff<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        
        m_terms.push_back( new gsTMTerm_BlendCoeffRhs<T>(m_paramsPtr, m_TMModelPtr, m_unknown) );
        
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = m_unknown; 
        m_trialUnkID = m_unknown;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void assemble();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setCurrentSolution(const gsMatrix<T>& solution) { m_solution = solution; }

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMVisitors.hpp)
#endif