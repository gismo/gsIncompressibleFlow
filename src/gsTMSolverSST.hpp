/** @file gsTMSolverSST.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

namespace gismo
{

// upravit
template<class T, int MatOrder>
void gsTMSolverSST<T, MatOrder>::initMembers()
{
    Base::initMembers();
}

// upravit
template<class T, int MatOrder>
void gsTMSolverSST<T, MatOrder>::evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
{
    m_TurbulentViscosityVals.setZero(quNodes.cols());
    if ((getAssembler()->isInitialized()) && (m_TMModelPtr->isInitialized()))
    {
        m_TMModelPtr->evalTurbulentViscosity(quNodes, numNodesPerElem, patchId);
        m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    }
    
}


}