/** @file gsINSTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{


template<class T>
void gsINSTermPvalUdiv<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    this->evalCoeff(mapData);
    gsVector<T> wCoeff = this->getCoeffWeightsProduct(quWeights);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunVals = shapeFunData[0];

    gsMatrix<T> testFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = wCoeff(k) * mapData.measure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);

        for (size_t i = 0; i != localMat.size(); ++i)
            localMat[i].noalias() += weight * (shapeFunVals.col(k) * testFunPhysGrad.row(i)).transpose();
    }
}

// ===================================================================================================================

template<class T>
void gsINSTermUdivPval<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    this->evalCoeff(mapData);
    gsVector<T> wCoeff = this->getCoeffWeightsProduct(quWeights);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    gsMatrix<T> shapeFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = wCoeff(k) * mapData.measure(k);

        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        for (size_t i = 0; i != localMat.size(); ++i)
            localMat[i].noalias() += weight * (shapeFunPhysGrad.row(i) * testFunVals(k) );
    }
}

// ===================================================================================================================

template<class T>
void gsINSTermUsolGradVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    this->computeCoeffSolU(mapData);
    this->evalCoeff(mapData);
    gsVector<T> wCoeff = this->getCoeffWeightsProduct(quWeights);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    gsMatrix<T> shapeFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = wCoeff(k) * mapData.measure(k);

        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        localMat.noalias() += weight * (testFunVals.col(k) * (m_solUVals.col(k).transpose() * shapeFunPhysGrad));
    }
}


} // namespace gismo