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

    typename gsRANSVisitorUU<T, MatOrder>::uPtr m_visitorRANSsymgradPtr;
    gsRANSVisitorTCSDStabilization_time<T, MatOrder> m_visitorRANS_TCSD_time;
    gsRANSVisitorTCSDStabilization_advection<T, MatOrder> m_visitorRANS_TCSD_advection;
    gsRANSVisitorSUPGStabilization_diffusion<T, MatOrder> m_visitorRANS_SUPG_diffusion;
    gsRANSVisitorSUPGStabilization_presssure<T, MatOrder> m_visitorRANS_SUPG_pressure;
    
    gsSparseMatrix<T, MatOrder> m_matRANSsymgrad, m_matRANS_TCSD_time, m_matRANS_TCSD_advection, m_matRANS_SUPG_diffusion, m_matRANS_SUPG_pressure, m_matRANS_TCSD_velocity_full;
    gsMatrix<T> m_rhsRANS, m_rhsRANS_TCSD_time, m_rhsRANS_TCSD_advection, m_rhsRANS_SUPG_diffusion, m_rhsRANS_SUPG_pressure, m_rhsRANS_TCSD_velocity_full;
    
    typename gsTMSolverBase<T, MatOrder>::Ptr m_TMsolverPtr = NULL;
    typename gsTMModelData<T>::Ptr m_TMModelPtr = NULL;
    bool m_bComputeTMfirst;
    gsField<T> m_oldTimeFieldU, m_currentFieldU;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_pshift;
    using Base::m_udofs;
    using Base::m_pdofs;
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

    gsRANSAssemblerUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMSolverBase<T, MatOrder>::Ptr TMsolverPtr): 
    Base(paramsPtr), m_TMsolverPtr(TMsolverPtr)
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

    virtual void makeBlockUU(gsSparseMatrix<T, MatOrder>& result, bool linPartOnly = false);

    virtual void makeRhsU(gsMatrix<T>& result, bool linPartOnly = false);

    /// @brief Add the nonlinear part to the given matrix and right-hand side.
    virtual void fillSystem();

    /// @brief Update the assembler in new nonlinear iteration.
    /// @param solVector 
    /// @param[in] updateSol    true - save solVector into m_solution (false is used in the inner Picard iteration for unsteady problem)
    virtual void update(const gsMatrix<T> & solVector, bool updateSol = true);

}; //gsRANSAssemblerUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSAssemblerUnsteady.hpp)
#endif