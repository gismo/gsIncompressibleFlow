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
class gsRANSVisitorUUSymmetricGradient : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
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
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradient() {}

    gsRANSVisitorUUSymmetricGradient(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { 
        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    /// @brief Evaluates turbulent viscosity.
    void evaluate(index_t testFunID);

    /// @brief Evaluates turbulent viscosity.
    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void defineTerms()
    {
        m_terms.push_back( new gsRANSTerm_SymmetricGradient<T>() );
    }

    virtual void defineTestTrialUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_trialUnkID = 0;   // velocity
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    void setRANSsolution(gsMatrix<T> sol) { m_solution = sol;}

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSVisitors.hpp)
#endif