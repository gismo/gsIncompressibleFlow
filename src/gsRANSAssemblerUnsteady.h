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
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsRANSVisitors.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>

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

    gsRANSVisitorUU<T, MatOrder> m_visitorRANSsymgrad;
    //gsRANSVisitorUU_full<T, MatOrder> m_visitorRANSsymgrad; // TODO: choose automatically for RbR loop
    gsSparseMatrix<T, MatOrder> m_matRANSsymgrad;
    gsMatrix<T> m_rhsRANS;
    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    bool m_bComputeTMfirst;
    gsField<T> m_oldTimeFieldU, m_currentFieldU;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_pshift;
    using Base::m_udofs;
    using Base::m_tarDim;
    using Base::m_nnzPerOuterU;
    using Base::m_solution;
    using Base::m_baseMatrix;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_currentVelField;
    using Base::m_blockTimeDiscr;
    using Base::m_rhsTimeDiscr;

public: // *** Constructor/destructor ***

    gsRANSAssemblerUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr): 
    Base(paramsPtr)
    { 
        initMembers();
    }

    virtual ~gsRANSAssemblerUnsteady()
    { }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

    /// @brief Update sizes of members (when DOF numbers change after constructing the assembler).
    virtual void updateSizes();

    /// @brief Assemble the linear part of the matrix.
    virtual void assembleLinearPart();

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();

    /// @brief Update the assembler in new nonlinear iteration.
    /// @param solVector 
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);


//public: // *** Member functions ***

public: // Getter/setters

    /// @brief Set turbulence solver to the member function.
    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

}; //gsRANSAssemblerUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSAssemblerUnsteady.hpp)
#endif