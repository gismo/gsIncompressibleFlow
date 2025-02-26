/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMTerms.h>

namespace gismo
{

template<class T>
void gsTMTerm_CoeffGradGrad<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        localMat.noalias() += weight * (m_k1 * m_turbViscosityVals(k) + m_k2) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);
    }
}

// ===================================================================================================================

template<class T>
void gsTMTerm_VecCoeffGradVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    this->computeCoeffSolU(mapData);
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    gsMatrix<T> shapeFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        localMat.noalias() += weight * (testFunVals.col(k) * (this->m_solUVals.col(k).transpose() * shapeFunPhysGrad));
    }
}

// ===================================================================================================================

template<class T>
void gsTMTerm_CoeffValVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunVals = shapeFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);
        localMat.noalias() += weight * m_konst * m_solVals(0,k) * (testFunVals.col(k) * shapeFunVals.col(k).transpose());
    }
}

/*
template<class T>
void gsFlowTerm_rhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    m_pRhsFun->eval_into(mapData.values[0], m_rhsVals);

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
    }
}

template<class T>
void gsFlowTerm_ValVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& shapeFunVals = shapeFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);
        localMat.noalias() += weight * (testFunVals.col(k) * shapeFunVals.col(k).transpose());
    }
}
*/

} // namespace gismo