/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

namespace gismo
{

template <class T, int MatOrder>
typename gsTMSolverBase<T, MatOrder>::tmPtr gsTMSolverBase<T, MatOrder>::make(std::string turbModel, typename gsFlowSolverParams<T>::Ptr paramsPtr)
{
    if (turbModel == "SST") 
    {
        return gsTMSolverSST<T, MatOrder>::make(paramsPtr);
    }
    //elseif (m_paramsPtr->options().getSwitch("TM.eval") == "SA") 
    //{ }
    else 
    {
        gsInfo << "Unknown identifier of a turbulent model entered! Using k-omega SST model." << std::endl;
        return gsTMSolverSST<T, MatOrder>::make(paramsPtr);
    }
}

} // namespace gismo