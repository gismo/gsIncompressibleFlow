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
void gsINSTerm_PvalUdiv<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    gsMatrix<T> testFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);

        for (size_t i = 0; i != localMat.size(); ++i)
            localMat[i].noalias() += weight * (trialFunVals.col(k) * testFunPhysGrad.row(i)).transpose();
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_UdivPval<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunGrads = trialFunData[1];

    gsMatrix<T> trialFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, trialFunGrads, trialFunPhysGrad);

        for (size_t i = 0; i != localMat.size(); ++i)
            localMat[i].noalias() += weight * (trialFunPhysGrad.row(i) * testFunVals(k) );
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_UsolGradVal<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    this->computeCoeffSolU(mapData);
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunGrads = trialFunData[1];

    gsMatrix<T> trialFunPhysGrad;

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, trialFunGrads, trialFunPhysGrad);

        localMat.noalias() += weight * (testFunVals.col(k) * (this->m_solUVals.col(k).transpose() * trialFunPhysGrad));
    }
}


// ===================================================================================================================
// For weak imposition of Dirichlet boundary conditions
// ===================================================================================================================

template<class T>
void gsINSTerm_CoeffUvalUval_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_coeff.resize(nQuPoints);

    this->computeCoeffSolU(mapData);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        m_coeff(k) = normal.dot(m_solUVals.col(k));
    } 
}

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsUVal_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_rhsVals.setZero(mapData.dim.first, nQuPoints);

    this->computeCoeffSolU(mapData);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        m_rhsVals.col(k) = (normal.dot(m_solUVals.col(k))) * m_vals.col(k);
    } 
}

template<class T>
void gsINSTerm_RhsUVal_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);
    
    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);
        
        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_CoeffUvalUvalPenalty_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    short_t dim = mapData.points.rows();
    real_t gama1 = m_penalties[0];
    real_t gama2 = m_penalties[1];
    real_t h = (mapData.points.col(nQuPoints-1) - mapData.points.col(0)).norm();

    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        localMat[0].noalias() += weight * gama1/h * (testFunVals.col(k) * trialFunVals.col(k).transpose());
        for (index_t i = 0; i < dim; i++)
            for (index_t j = 0; j < dim; j++)
                localMat[j + i*dim].noalias() += weight * gama2/h * (normal(i) * normal(j) * testFunVals.col(k) * trialFunVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsUValPenalty_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    real_t gama1 = m_penalties[0];
    real_t gama2 = m_penalties[1];
    real_t h = (mapData.points.col(nQuPoints-1) - mapData.points.col(0)).norm();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();
        
        localMat.noalias() += weight * (gama1/h * (testFunData[0].col(k) *  m_vals.col(k).transpose()) + gama2/h * normal.dot(m_vals.col(k)) * (testFunData[0].col(k) * normal.asRowVector()));
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_CoeffUvalUdiv_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    gsMatrix<T> testFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();
        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);

        localMat.noalias() += weight * (trialFunVals.col(k) *  (normal.reshape(1, normal.rows()) * testFunPhysGrad));
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsUdiv_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_rhsVals.setZero(1, nQuPoints);

    for (index_t k = 0; k < nQuPoints; k++)
        m_rhsVals(0, k) = m_viscosity;
}

template<class T>
void gsINSTerm_RhsUdiv_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);

    const gsMatrix<T>& testFunGrads = testFunData[1];
    gsMatrix<T> testFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();
        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        
        gsMatrix<T> pom = normal.reshape(1, normal.rows()) * testFunPhysGrad;
        localMat.noalias() += weight * m_rhsVals(0, k) * (pom.transpose() * m_vals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_PvalUval_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);
        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        for (size_t i = 0; i < localMat.size()-1; ++i)
            localMat[i].noalias() += weight * (normal(i) * testFunVals.col(k) * trialFunVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsPvalU_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_rhsVals.setZero(1, nQuPoints);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        m_rhsVals(0, k) = normal.dot(m_vals.col(k));
    } 
}

template<class T>
void gsINSTerm_RhsPvalU_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals(0, k));
    }
}

// ===================================================================================================================

/*
template<class T>
void gsINSTerm_UvalPval_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);
        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        for (size_t i = 0; i < localMat.size()-1; ++i)
            localMat[i].noalias() += weight * (normal(i) * testFunVals.col(k) * trialFunVals.col(k).transpose());
    }
}
*/

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsUvalP_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_rhsVals.setZero(1, nQuPoints);

    for (index_t k = 0; k < nQuPoints; k++)
        m_rhsVals(0, k) = m_vals(k);
}

template<class T>
void gsINSTerm_RhsUvalP_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        gsVector<T> normal = mapData.outNormal(k);
        normal.normalize();

        localMat.noalias() += weight * (m_rhsVals(0, k) * testFunData[0].col(k) * normal.asRowVector());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_PvalPval_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const gsMatrix<T>& testFunVals = testFunData[0];
    const gsMatrix<T>& trialFunVals = trialFunData[0];

    const index_t nQuPoints = quWeights.rows();
    real_t gama3 = m_penalties[2];
    real_t h = (mapData.points.col(nQuPoints-1) - mapData.points.col(0)).norm();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k) * gama3 / h;

        localMat.noalias() += weight * (testFunVals.col(k) * trialFunVals.col(k).transpose());
    }
}

// ===================================================================================================================

template<class T>
void gsINSTerm_RhsPvalP_WeakDirichlet<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    m_rhsVals.setZero(1, nQuPoints);

    for (index_t k = 0; k < nQuPoints; k++)
        m_rhsVals(0, k) = m_vals(k);
}

template<class T>
void gsINSTerm_RhsPvalP_WeakDirichlet<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& trialFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);

    real_t gama3 = m_penalties[2];
    real_t h = (mapData.points.col(nQuPoints-1) - mapData.points.col(0)).norm();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k) * gama3 / h;

        localMat.noalias() += weight * (m_rhsVals(0, k) * testFunData[0].col(k));
    }
}

} // namespace gismo