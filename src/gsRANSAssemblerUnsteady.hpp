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
    updateSizes();

    m_visitorRANSsymgrad = gsRANSVisitorUUSymmetricGradient<T, MatOrder>(m_paramsPtr);
    m_visitorRANSsymgrad.initialize();

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
void gsRANSAssemblerUnsteady<T, MatOrder>::updateDofMappers()
{
    Base::updateDofMappers();

    m_visitorRANSsymgrad.updateDofMappers(m_dofMappers);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_matRANSsymgrad);
    m_rhs.topRows(m_pshift) += m_rhsRANS;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
    {
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

        m_visitorRANSsymgrad.setTurbulenceSolver(m_TMsolverPtr);
        m_visitorRANSsymgrad.setRANSsolution(solVector);
        m_matRANSsymgrad.resize(m_pshift, m_pshift);
        m_matRANSsymgrad.reserve(gsVector<index_t>::Constant(m_matRANSsymgrad.outerSize(), m_tarDim * m_nnzPerOuterU));
        m_rhsRANS.setZero(m_pshift, 1);
        this->assembleBlock(m_visitorRANSsymgrad, 0, m_matRANSsymgrad, m_rhsRANS);
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
    {
        fillSystem();
    }
}

}