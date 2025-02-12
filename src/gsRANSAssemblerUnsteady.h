/** @file gsRANSAssemblerUnsteady.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsRANSVisitors.h>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>

#include <gsCore/gsField.h>

namespace gismo
{

template<class T, int MatOrder>
class gsRANSAssemblerUnsteady: public gsINSAssemblerUnsteady<T, MatOrder>
{

public:
    typedef gsINSAssemblerUnsteady<T, MatOrder> Base;

protected: // *** Class members ***

    gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder> m_visitorRANSsymgraddiag;
    gsRANSVisitorUUSymmetricGradientOffdiag<T, MatOrder> m_visitorRANSsymgradoffdiag;
    gsSparseMatrix<T, MatOrder> m_matRANSsymgraddiag;
    gsSparseMatrix<T, MatOrder> m_matRANSsymgradoffdiag;
    gsMatrix<T> m_rhsRANS;
    gsMatrix<T> m_rhsRANSsymgradoffdiag;
    gsTMSolverSST<T, MatOrder>* m_TMsolver;
    //gsTMSolverSST<T, MatOrder>* m_TMsolver;
    bool m_bComputeTMfirst;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_pshift;
    using Base::m_udofs;
    using Base::m_tarDim;
    using Base::m_nnzPerRowU;
    using Base::m_dofMappers;
    using Base::m_solution;
    using Base::m_baseMatrix;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_currentVelField;

public: // *** Constructor/destructor ***

    gsRANSAssemblerUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr, gsTMSolverSST<T, MatOrder>* TMsolver): 
    Base(paramsPtr)
    { 
        initMembers(TMsolver);
    }

    virtual ~gsRANSAssemblerUnsteady()
    { }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers(gsTMSolverSST<T, MatOrder>* TMsolver);

    /// @brief Assemble the linear part of the matrix.
    virtual void assembleLinearPart();

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();


//public: // *** Member functions ***

}; //gsRANSAssemblerUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSAssemblerUnsteady.hpp)
#endif