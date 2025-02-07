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
void gsRANSAssemblerUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();
    updateSizes();

    m_visitorRANS = gsRANSVisitorUUlin<T, MatOrder>(m_paramsPtr);
    m_visitorRANS.initialize();
    
    // neco s turbulentni viskozitou nebo modelem?

}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaning
    m_blockRANS.resize(m_tarDim * m_udofs, m_tarDim * m_udofs);
    m_blockRANS.reserve(gsVector<index_t>::Constant(m_blockRANS.outerSize(), m_tarDim * m_nnzPerRowU));

    this->assembleBlock(m_visitorRANS, 0, m_blockRANS, m_rhsRANS);

}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_blockRANS);
    m_rhs.topRows(m_pshift) += m_rhsRANS;
}

}