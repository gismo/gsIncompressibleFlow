/** @file gsTMTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMTerms.h>

namespace gismo
{

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

        localMat.noalias() += weight * (testFunVals.col(k) * (m_solUVals.col(k).transpose() * shapeFunPhysGrad));
    }
}

// ===================================================================================================================

template<class T>
void gsTMTerm_CoeffGradGrad<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    gsVector<T> F1(nQuPoints);
    gsVector<T> turbViscosityVals(nQuPoints);
    m_coeff.resize(nQuPoints);
        
    real_t eps = m_TMModelPtr->get_eps();
    F1 = m_TMModelPtr->getF1Vals();
    turbViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();

    real_t k1, k2, k3;
    if (m_unknown == 2)
    {
        k1 = m_TMModelPtr->get_sigmaK1();
        k2 = m_TMModelPtr->get_sigmaK2();
        k3 = m_paramsPtr->getPde().viscosity();   
    }
    else
    {
        k1 = m_TMModelPtr->get_sigmaO1();
        k2 = m_TMModelPtr->get_sigmaO2();
        k3 = m_paramsPtr->getPde().viscosity();
    }
    
    // final coefficient
    m_coeff.resize(nQuPoints);    
    for (index_t k = 0; k < nQuPoints; k++)
    {
        m_coeff(k) = math::max((k1 * F1(k) + k2 * (1 - F1(k))) * turbViscosityVals(k) + k3, eps);
    }
}


// ===================================================================================================================

template<class T>
void gsTMTerm_CoeffValVal<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    real_t eps = m_TMModelPtr->get_eps();
    gsMatrix<T> OSolVals(1, nQuPoints);
    OSolVals = m_TMModelPtr->getOSolVals();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    
    m_coeff.resize(nQuPoints);
    if (m_unknown == 2)         // for k equation of SST model
    {
        for (index_t k = 0; k < nQuPoints; k++)
            m_coeff(k) = betaStar * math::max(OSolVals(0, k), eps);
    }
    else                        // for omega equation of SST model
    {
        gsVector<T> F1(nQuPoints);
        F1 = m_TMModelPtr->getF1Vals();
        for (index_t k = 0; k < nQuPoints; k++)
            m_coeff(k) = (beta1 * F1(k) + beta2 * (1 - F1(k))) * math::max(OSolVals(0, k), eps);
    }
}

// =====================================================================================================================

template<class T>
void gsTMTerm_BlendCoeff<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    gsVector<T> F1(nQuPoints);
    gsMatrix<T> KSolVals(1, nQuPoints);
    gsMatrix<T> OSolVals(1, nQuPoints);
    std::vector< gsMatrix<T> > KSolDers;
    std::vector< gsMatrix<T> > OSolDers;
        
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    F1 = m_TMModelPtr->getF1Vals();
    OSolVals = m_TMModelPtr->getOSolVals();
    KSolDers = m_TMModelPtr->getKSolDers();
    OSolDers = m_TMModelPtr->getOSolDers();
    
    m_coeff.setZero(nQuPoints);
    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsMatrix<T> u = KSolDers[k];
        gsMatrix<T> v = OSolDers[k];
        m_coeff(k) = math::max((-1) * math::max(2 * (1 - F1(k)) * sigmaO2, 0.) / math::max(math::pow(OSolVals(0, k), 2), math::pow(10, -15)) * (u.row(0).dot(v.row(0))), 0.);
    }
}


template<class T>
void gsTMTerm_BlendCoeff<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    if (m_unknown == 3)
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
    else
    {
        gsMatrix<T> dummyMat(testFunData[0].rows(), shapeFunData[0].rows());
        dummyMat.setZero();
        localMat.noalias() += dummyMat;
    }
}


// ============================================================================================================================

template<class T>
void gsTMTerm_BlendCoeffRhs<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    gsVector<T> F1(nQuPoints);
    gsMatrix<T> KSolVals(1, nQuPoints);
    gsMatrix<T> OSolVals(1, nQuPoints);
    std::vector< gsMatrix<T> > KSolDers;
    std::vector< gsMatrix<T> > OSolDers;
        
    real_t eps = m_TMModelPtr->get_eps();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    F1 = m_TMModelPtr->getF1Vals();
    OSolVals = m_TMModelPtr->getOSolVals();
    KSolDers = m_TMModelPtr->getKSolDers();
    OSolDers = m_TMModelPtr->getOSolDers();
    
    m_rhsVals.resize(1, nQuPoints);
    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsMatrix<T> u = KSolDers[k];
        gsMatrix<T> v = OSolDers[k];
        m_rhsVals(0, k) = math::max(2 * (1 - F1(k)) * sigmaO2 / math::max(OSolVals(0, k), eps) * (u.row(0).dot(v.row(0))), 0.);
    }
}

template<class T>
void gsTMTerm_BlendCoeffRhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    if (m_unknown == 3)
    {
        const index_t nQuPoints = quWeights.rows();
        evalCoeff(mapData);

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = quWeights(k) * mapData.measure(k);

            localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals(0, k));
        }
    }
    else
    {
        gsMatrix<T> dummyRhs(testFunData[0].rows(), 1);
        dummyRhs.setZero();
        localMat.noalias() += dummyRhs;
    }
}

// ============================================================================================================================

template<class T>
void gsTMTerm_ProductionRhs<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    index_t dim = mapData.dim.first;
    gsVector<T> turbViscosityVals;
    gsMatrix<T> KSolVals(1, nQuPoints);
    gsMatrix<T> OSolVals(1, nQuPoints);
    gsVector<T> F1(nQuPoints);
    gsVector<T> StrainRateMag;
    std::vector< gsMatrix<T> > StrainRateTensor;
    std::vector< gsMatrix<T> > USolDers;
        
    real_t eps = m_TMModelPtr->get_eps();
    real_t sigmaO1 = m_TMModelPtr->get_sigmaO1();
    real_t sigmaO2 = m_TMModelPtr->get_sigmaO2();
    real_t beta1 = m_TMModelPtr->get_beta1();
    real_t beta2 = m_TMModelPtr->get_beta2();
    real_t betaStar = m_TMModelPtr->get_betaStar();
    real_t kappa = m_TMModelPtr->get_kappa();
    KSolVals = m_TMModelPtr->getKSolVals();
    OSolVals = m_TMModelPtr->getOSolVals();
    USolDers = m_TMModelPtr->getUSolDers();
    F1 = m_TMModelPtr->getF1Vals();
    turbViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    StrainRateTensor = m_TMModelPtr->getStrainRateTensor();
    StrainRateMag = m_TMModelPtr->getStrainRateMagVals();
    
    m_rhsVals.resize(1, nQuPoints); 
    for (index_t k = 0; k < nQuPoints; k++)
    {
        //k-omega SST
        real_t StrainRateUDers = 0.0;
        for (index_t i = 0; i < dim; i++)
            for (index_t j = 0; j < dim; j++)
                StrainRateUDers += StrainRateTensor[k](i,j) * USolDers[k](i, j);
        m_rhsVals(0, k) = math::min(2 * turbViscosityVals(k) * StrainRateUDers, 20 * betaStar * math::max(KSolVals(0, k), eps) * math::max(OSolVals(0, k), eps));

        // k-omega SST menter 2009
        //m_rhsVals(0, k) = math::min(turbViscosityVals(k) * math::pow(StrainRateMag(k), 2), 10 * betaStar * KSolVals(0, k) * OSolVals(0, k));
    }

    if (m_unknown == 3)
    {
        real_t gamma1 = beta1/betaStar - (sigmaO1 * math::pow(kappa,2))/(math::sqrt(betaStar));
        real_t gamma2 = beta2/betaStar - (sigmaO2 * math::pow(kappa,2))/(math::sqrt(betaStar));

        for (index_t k = 0; k < nQuPoints; k++)
        {
            m_rhsVals(0, k) = ((gamma1 * F1(k) + gamma2 * (1 - F1(k))) / turbViscosityVals(k)) * m_rhsVals(0, k);
        }
    }
    
}

template<class T>
void gsTMTerm_ProductionRhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();
    evalCoeff(mapData);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals(0, k));
    }
}



} // namespace gismo