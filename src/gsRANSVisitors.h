/** @file gsRANSVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsRANSTerms.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

namespace gismo
{

template <class T, int MatOrder>
class gsRANSVisitorUU : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsRANSVisitorUU> Ptr; 
    typedef memory::unique_ptr<gsRANSVisitorUU> uPtr;

protected:  // *** Class members ***

    gsMatrix<T> m_solution;
    real_t m_viscosity;
    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    gsVector<T> m_TurbulentViscosityVals;

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
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsRANSVisitorUU() {}

    gsRANSVisitorUU(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { 
        initMembers();
    }


public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    virtual void defineTerms()
    {
        m_terms.push_back( new gsRANSTerm_SymmetricGradient<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 0;   // velocity
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void resizeMatVec()
    {
        short_t dim = m_paramsPtr->getPde().dim();
        index_t nblocksE = dim * (dim + 1) / 2; // #diag + #upper-off-diag blocks
        m_locMatVec.resize(nblocksE + 1); // +1 for matrix A which is equal for all components

        // m_locMatVec = [A, E11, E22, (E33), E12, (E13), (E23)]
    }


public: // *** Getters/setters ***

    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    void setRANSsolution(gsMatrix<T> sol) { m_solution = sol;}

}; // gsRANSVisitorUU


template <class T, int MatOrder>
class gsRANSVisitorUU_full : public gsRANSVisitorUU<T, MatOrder>
{

public:
    typedef gsRANSVisitorUU<T, MatOrder> Base;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsRANSVisitorUU_full> Ptr; 
    typedef memory::unique_ptr<gsRANSVisitorUU_full> uPtr;

protected: // *** Base class members ***

    using Base::m_solution;
    using Base::m_viscosity;
    using Base::m_TMsolverPtr;
    using Base::m_TurbulentViscosityVals;
    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsRANSVisitorUU_full() {}

    gsRANSVisitorUU_full(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsRANSTerm_SymmetricGradient_full<T>() );
    }

    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    virtual void resizeMatVec()
    {
        short_t dim = m_paramsPtr->getPde().dim();
        m_locMatVec.resize(dim * dim + 1); // +1 for matrix A which is equal for all components

        // m_locMatVec = [A, E11, E12, (E13), E21. E22, (E23), (E31), (E32), (E33)]
    }

};


// ================================================================================================================
// For T-CSD stabilization
//

template <class T, int MatOrder>
class gsRANSVisitorTCSDStabilization_time : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;

public:
    gsField<T> m_solution;
    real_t m_viscosity;
    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    gsVector<T> m_TurbulentViscosityVals;
    gsMatrix<T> m_tauS;

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
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsRANSVisitorTCSDStabilization_time() {}

    gsRANSVisitorTCSDStabilization_time(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { 
        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();


    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TCSDStabilization_time<T>() );
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

public: // *** Getters/setters ***

    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    void setRANSsolution(gsField<T>& solution) { m_solution = solution;}

};

template <class T, int MatOrder>
class gsRANSVisitorTCSDStabilization_advection : public gsRANSVisitorTCSDStabilization_time<T, MatOrder>
{

public:
    typedef gsRANSVisitorTCSDStabilization_time<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_tauS;
    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_hasPeriodicBC;
    using Base::m_periodicTransformMat;
    

public: // *** Constructor/destructor ***

    gsRANSVisitorTCSDStabilization_advection() {}

    gsRANSVisitorTCSDStabilization_advection(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    {   }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TCSDStabilization_advection<T>() );
    }

public: // *** Member functions *** 

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSVisitors.hpp)
#endif