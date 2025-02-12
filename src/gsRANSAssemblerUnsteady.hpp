/** @file gsINSAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>

namespace gismo
{

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::initMembers(gsTMSolverSST<T, MatOrder>* TMsolver)
{
    Base::initMembers();
    Base::updateSizes();

    m_TMsolver = TMsolver;

    m_visitorRANSsymgraddiag = gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder>(m_paramsPtr, m_TMsolver);
    m_visitorRANSsymgraddiag.initialize();

    m_visitorRANSsymgradoffdiag = gsRANSVisitorUUSymmetricGradientOffdiag<T, MatOrder>(m_paramsPtr, m_TMsolver);
    m_visitorRANSsymgradoffdiag.initialize();
    
    // neco s turbulentni viskozitou nebo modelem?
    m_TMsolver = TMsolver;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaning
    m_matRANSsymgraddiag.resize(m_tarDim * m_udofs, m_tarDim * m_udofs);
    m_matRANSsymgraddiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgraddiag.outerSize(), m_tarDim * m_nnzPerRowU));

    this->assembleBlock(m_visitorRANSsymgraddiag, 0, m_matRANSsymgraddiag, m_rhsRANS);

    // matrix cleaning
    m_matRANSsymgradoffdiag.resize(m_pshift, m_pshift);
    m_matRANSsymgradoffdiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgradoffdiag.outerSize(), m_tarDim * m_nnzPerRowU));

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorRANSsymgradoffdiag, 0, m_matRANSsymgradoffdiag, dummyRhs);

    m_rhsRANSsymgradoffdiag = m_matRANSsymgradoffdiag * m_solution.topRows(m_pshift);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_matRANSsymgraddiag);
    m_rhs.topRows(m_pshift) += m_rhsRANS;
    m_rhs.topRows(m_pshift) += m_rhsRANSsymgradoffdiag;
}

}