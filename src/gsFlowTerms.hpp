/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowTerms.h>

namespace gismo
{

template<class T>
gsVector<T> gsFlowTerm<T>::getCoeffWeightsProduct(const gsVector<T>& quWeights)
{
    if (m_coeff.size() == 1)
        return m_coeff(0) * quWeights;
    else if (m_coeff.size() == quWeights.size())
        return m_coeff.array() * quWeights.array();
    else
        GISMO_ERROR("Wrong size of m_coeff in gsFlowTerm.");
}


template<class T>
void gsFlowTermValVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    this->evalCoeff(mapData);
    gsVector<T> wCoeff = this->getCoeffWeightsProduct(quWeights);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunVals = shapeFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = wCoeff(k) * mapData.measure(k);
        localMat += weight * (testFunVals.col(k) * shapeFunVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTermGradGrad<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    this->evalCoeff(mapData);
    gsVector<T> wCoeff = this->getCoeffWeightsProduct(quWeights);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = wCoeff(k) * mapData.measure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        localMat += weight * (testFunPhysGrad.transpose() * shapeFunPhysGrad);
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTermRhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    m_pRhsFun->eval_into(mapData.values[0], m_rhsVals);

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
    }
}


} // namespace gismo