/** @file gsRANSAssemblerUnsteady.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova), B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>

namespace gismo
{

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();
    
    m_TMModelPtr = m_TMsolverPtr->getTMModel();

    if (m_paramsPtr->options().getString("assemb.loop") == "EbE")
        m_visitorRANSsymgradPtr = std::make_unique< gsRANSVisitorUU<T, MatOrder> >(m_paramsPtr);
    else
        m_visitorRANSsymgradPtr = std::make_unique< gsRANSVisitorUU_full<T, MatOrder> >(m_paramsPtr);

    m_visitorRANSsymgradPtr->initialize();
    m_visitorRANSsymgradPtr->setTurbulenceModel(m_TMModelPtr);

    m_visitorUUnonlin.setTurbulenceModel(m_TMModelPtr);
    m_visitorUU_TCSD_time.setTurbulenceModel(m_TMModelPtr);
    //m_visitorUP_SUPG_pressure.setTurbulenceModel(m_TMModelPtr);
    //m_visitorPP_ResStab_continuity.setTurbulenceModel(m_TMModelPtr);

    updateSizes();

}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::updateSizes()
{
    Base::updateSizes();

    m_matRANSsymgrad.resize(m_pshift, m_pshift);
    m_rhsRANS.setZero(m_pshift, 1);

    m_currentFieldU = Base::constructSolution(m_solution, 0);
    m_paramsPtr->setVelocitySolution(m_currentFieldU);
    
    m_oldTimeFieldU = m_currentFieldU;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::makeBlockUU(gsSparseMatrix<T, MatOrder>& result, bool linPartOnly)
{
    Base::makeBlockUU(result, linPartOnly);
    result += m_matRANSsymgrad;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::makeRhsU(gsMatrix<T>& result, bool linPartOnly)
{
    Base::makeRhsU(result, linPartOnly);
    result += m_rhsRANS;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_matRANSsymgrad);
    m_rhs.topRows(m_pshift) += m_rhsRANS;

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
    {
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

        //m_visitorRANSsymgradPtr->setTurbulenceSolver(m_TMsolverPtr);
        m_visitorRANSsymgradPtr->setRANSsolution(solVector);
        m_matRANSsymgrad.resize(m_pshift, m_pshift);
        m_matRANSsymgrad.reserve(gsVector<index_t>::Constant(m_matRANSsymgrad.outerSize(), m_tarDim * m_nnzPerOuterU));
        m_rhsRANS.setZero(m_pshift, 1);
        this->assembleBlock(*m_visitorRANSsymgradPtr, 0, m_matRANSsymgrad, m_rhsRANS);
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
    {
        fillSystem();
    }
}

}