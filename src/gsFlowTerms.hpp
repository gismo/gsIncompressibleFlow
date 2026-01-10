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
gsVector<T> gsFlowTerm<T>::getCoeffGeoMapProduct(const gsMapData<T>& mapData)
{
    evalCoeff(mapData);

    index_t npts = mapData.points.cols();
    gsVector<T> result(npts);

    for (index_t k = 0; k < npts; k++)
    {
        index_t coeffID = (m_coeff.size() == 1) ? 0 : k;
        result(k) = m_coeff(coeffID) * mapData.measure(k);
    }

    return result;
}


template<class T>
void gsFlowTerm_ValVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);
        localMat.noalias() += weight * (testFunVals.col(k) * trialFunVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTerm_GradGrad<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& trialFunGrads = trialFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, trialFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, trialFunGrads, trialFunPhysGrad);

        localMat.noalias() += weight * (testFunPhysGrad.transpose() * trialFunPhysGrad);
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTerm_TCSDStabilization_time<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    //this->computeCoeffSolU(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& trialFunVals = trialFunData[0];
    
    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, trialFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        
        localMat.noalias() += weight * m_tauS(0, k) * (trialFunVals.col(k) * (this->m_solUVals.col(k).transpose() * testFunPhysGrad));
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTerm_TCSDStabilization_advection<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    //this->computeCoeffSolU(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& trialFunGrads = trialFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, trialFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, trialFunGrads, trialFunPhysGrad);

        localMat.noalias() += weight * m_tauS(0, k) * ((this->m_solUVals.col(k).transpose() * trialFunPhysGrad).transpose() * (this->m_solUVals.col(k).transpose() * testFunPhysGrad));
    }
}

// ===================================================================================================================

template<class T>
void gsFlowTerm_rhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
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