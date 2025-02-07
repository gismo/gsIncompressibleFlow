/** @file gsINSVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsRANSTerms.h>

namespace gismo
{

template <class T, int MatOrder>
class gsRANSVisitorUUSymmetricGradientDiag : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradientDiag() {}

    gsRANSVisitorUUSymmetricGradientDiag(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsRANSTerm_SymmetricGradientDiag<T>(m_paramsPtr->getPde().viscosity()) );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

public: // *** Member functions ***    

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ==========================================================================================================

template <class T, int MatOrder>
class gsRANSVisitorUUSymmetricGradientOffdiag : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradientOffdiag() {}

    gsRANSVisitorUUSymmetricGradientOffdiag(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsRANSTerm_SymmetricGradientOffdiag<T>(m_paramsPtr->getPde().viscosity()) );

        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

public: // *** Member functions ***    

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSVisitors.hpp)
#endif