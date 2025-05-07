/** @file gsFlowBndEvaluators.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/


#pragma once
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>

namespace gismo
{

template<class T>
void gsFlowBndEvaluator<T>::evalOnPatchSide(index_t patchID, boxSide side, bool setZeroFirst)
{
    if (setZeroFirst)
        m_quantValue = 0.0;

    const gsBasis<T>* basis = &m_paramsPtr->getBasis(m_unkID).basis(patchID);
    short_t dim = m_paramsPtr->getPde().domain().targetDim();
    m_mapData.patchId = patchID;

    gsVector<int> numQuadNodes(dim);
    const int dir = side.direction();
    for (short_t i = 0; i < dim; ++i)
        numQuadNodes[i] = basis->degree(i) + 1; // TODO: check correctness
    numQuadNodes[dir] = 1;

    gsGaussRule<T> QuRule(numQuadNodes);

    gsMatrix<T> quNodes; 
    gsVector<T> quWeights; 

    typename gsBasis<T>::domainIter domIt = basis->domain()->beginBdr(side);
    typename gsBasis<T>::domainIter domItEnd = basis->domain()->endBdr(side);
    for (; domIt!=domItEnd; ++domIt)
    {
        QuRule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);
    
        m_mapData.points = quNodes;
        m_mapData.side = side;
        m_paramsPtr->getPde().patches().patch(patchID).computeMap(m_mapData);

        this->evalOnElement(patchID, side, quNodes, quWeights);
    }
}


// ===================================================================================================================


template<class T>
void gsFlowBndEvaluator_flowRate<T>::evalOnElement(index_t patchID, boxSide side, const gsMatrix<T>& quNodes, const gsVector<T>& quWeights)
{
    gsMatrix<T> solUVals = m_velocityField.value(quNodes, patchID);

    for (index_t k = 0; k < quWeights.rows(); ++k)
    {
        gsVector<T> normal = m_mapData.outNormal(k);
        
        // the normal norm is equal to integral measure
        m_quantValue += quWeights[k] * normal.dot(solUVals.col(k));
    }
}

} // namespace gismo
