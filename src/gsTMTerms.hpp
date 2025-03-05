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
void gsTMTerm_CoeffGradGrad<T>::evalCoeff(const gsMapData<T>& mapData)
{
    real_t a1 = m_paramsPtr->getSSTModel().get_a1();
    real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
    real_t sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
    real_t visc = m_paramsPtr->getPde().viscosity();
    index_t nQuPoints = mapData.points.cols();
    index_t dim = mapData.dim.first;
    gsVector<T> F1(nQuPoints);
    std::vector< gsMatrix<T> > StrainRateTensor;

    if (!(m_paramsPtr->getSSTModel().isCurrent()))
    {
        gsField<T> USolField = m_paramsPtr->getVelocitySolution();
        gsField<T> KSolField = m_paramsPtr->getKSolution();
        gsField<T> OSolField = m_paramsPtr->getOmegaSolution();

        // evaluate k, omega
        gsVector<T> solKVals(dim, nQuPoints);
        solKVals = KSolField.value(mapData.points, mapData.patchId);
        gsVector<T> solOVals(dim, nQuPoints);
        solOVals = OSolField.value(mapData.points, mapData.patchId);

        //evaluate grad(k), grad(omega)
        //gsFunction<T> KSol;
        std::vector< gsMatrix<T> > KSolDers = KSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
        std::vector< gsMatrix<T> > OSolDers = OSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);

        // evaluate strainrate tensor S
        std::vector< gsMatrix<T> > USolDers = USolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
        gsVector<T> StrainRateMag(nQuPoints);
        StrainRateMag.setZero();
        real_t Sij;
        for (index_t k = 0; k < nQuPoints; k++)
        {
            gsMatrix<T> SS(dim, dim);
            SS.setZero();
            for (index_t i = 0; i < dim; i++)
                for (index_t j = 0; j < dim; j++)
                {
                    Sij = 0.5 * (USolDers[1](i * dim + j, k) + USolDers[1](j * dim + i, k));
                    SS(i, j) = Sij;
                    StrainRateMag(k) += 2 * Sij * Sij;
                }
            StrainRateTensor.push_back(SS);
        }

        // UPRAVIT !!! evaluate distance
        gsVector<T> Distance(nQuPoints);
        //Distance = DistanceField.value(mapData.points, mapData.patchId);
        for (index_t k = 0; k < nQuPoints; k++)
            Distance(k) = 1.0;

        // evaluate F2
        gsVector<T> F2(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            F2(k) = math::tanh(math::pow(math::max((2 * math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), 2));
            F2(k) = math::max(F2(k), 0.0);
            F2(k) = math::min(F2(k), 1.0);
        }

        // evaluate turbulent viscosity
        for (index_t k = 0; k < nQuPoints; k++)
            m_turbViscosityVals(k) = (a1 * solKVals(k)) / (math::max(a1 * solOVals(k), StrainRateMag(k) * F2(k)));

        //evaluate CDkomega
        real_t gradkdotgradomega;
        gsVector<T> CDkomega(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            gradkdotgradomega = 0.0;
            for (index_t i = 0; i < dim; i++)
                gradkdotgradomega += KSolDers[1](i, k) * OSolDers[1](i, k);
            CDkomega(k) = math::max(2 * sigmaO2 / solOVals(k) * gradkdotgradomega, math::pow(10, -10));
        }

        // evaluate F1
        for (index_t k = 0; k < nQuPoints; k++)
        {
            F1(k) = math::tanh(math::pow(math::min(math::max((math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), (4 * sigmaO2 * solKVals(k))/(CDkomega(k) * math::pow(Distance(k), 2))), 4));
            F1(k) = math::max(F1(k), 0.0);
            F1(k) = math::min(F1(k), 1.0);
        }

        m_paramsPtr->getSSTModel().setKSolVals(solKVals);
        m_paramsPtr->getSSTModel().setOSolVals(solOVals);
        m_paramsPtr->getSSTModel().setKSolDers(KSolDers);
        m_paramsPtr->getSSTModel().setKSolDers(OSolDers);
        m_paramsPtr->getSSTModel().setF1Vals(F1);
        m_paramsPtr->getSSTModel().setF2Vals(F2);
        m_paramsPtr->getSSTModel().setTurbulentViscosityVals(m_turbViscosityVals);
        m_paramsPtr->getSSTModel().setStrainRateMagVals(StrainRateMag);
        m_paramsPtr->getSSTModel().StrainRateTensor(StrainRateTensor);
        m_paramsPtr->getSSTModel().setCurrent();
    }
    else
    {
        F1 = m_paramsPtr->getSSTModel().getF1Vals();
        m_turbViscosityVals = m_paramsPtr->getSSTModel().getTurbulentViscosityVals();
    }
    // final coefficient
    m_coeff.resize(nQuPoints);    
    for (index_t k = 0; k < nQuPoints; k++)
    {
        m_coeff(k) = (m_k1 * F1(k) + m_k2 * (1 - F1(k))) * m_turbViscosityVals(k) + m_k3;
    }
}


// ===================================================================================================================

template<class T>
void gsTMTerm_CoeffValVal<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    index_t dim = mapData.dim.first;
    gsVector<T> m_solOVals;
    m_solOVals.resize(nQuPoints);
    if (!(m_paramsPtr->getSSTModel().isCurrent()))
    {
        gsField<T> OSolField = m_paramsPtr->getOmegaSolution();
        m_solOVals = OSolField.value(mapData.points, mapData.patchId);
    }
    else
    {
        m_solOVals = m_paramsPtr->getSSTModel().getOSolVals();
    }

    if (m_unknown == 0)
    {
        // constant coefficient betaStar
        real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
        for (index_t k = 0; k < nQuPoints; k++)
            m_coeff(k) = betaStar * m_solOVals(k);
    }
    else
    {
        // coefficient beta = beta1 * F1 + beta2 * (1 - F1)
        real_t sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
        real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
        real_t visc = m_paramsPtr->getPde().viscosity();
        gsVector<T> F1(nQuPoints);
        std::vector< gsMatrix<T> > StrainRateTensor;
        if (!(m_paramsPtr->getSSTModel().isCurrent()))
        {
            gsField<T> USolField = m_paramsPtr->getVelocitySolution();
            gsField<T> KSolField = m_paramsPtr->getKSolution();
            gsField<T> OSolField = m_paramsPtr->getOmegaSolution();

            // evaluate k, omega
            gsVector<T> solKVals(dim, nQuPoints);
            solKVals = KSolField.value(mapData.points, mapData.patchId);
            gsVector<T> solOVals(dim, nQuPoints);
            solOVals = OSolField.value(mapData.points, mapData.patchId);

            //evaluate grad(k), grad(omega)
            std::vector< gsMatrix<T> > KSolDers = KSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
            std::vector< gsMatrix<T> > OSolDers = OSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);

            // evaluate strainrate tensor S
            std::vector< gsMatrix<T> > USolDers = USolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
            gsVector<T> StrainRateMag(nQuPoints);
            StrainRateMag.setZero();
            real_t Sij;
            for (index_t k = 0; k < nQuPoints; k++)
            {
                gsMatrix<T> SS(dim, dim);
                SS.setZero();
                for (index_t i = 0; i < dim; i++)
                    for (index_t j = 0; j < dim; j++)
                    {
                        Sij = 0.5 * (USolDers[1](i * dim + j, k) + USolDers[1](j * dim + i, k));
                        SS(i, j) = Sij;
                        StrainRateMag(k) += 2 * Sij * Sij;
                    }
                StrainRateTensor.push_back(SS);
            }

            // UPRAVIT !!! evaluate distance
            gsVector<T> Distance(nQuPoints);
            //Distance = DistanceField.value(mapData.points, mapData.patchId);
            for (index_t k = 0; k < nQuPoints; k++)
                Distance(k) = 1.0;
            
            //evaluate CDkomega
            real_t gradkdotgradomega;
            gsVector<T> CDkomega(nQuPoints);
            for (index_t k = 0; k < nQuPoints; k++)
            {
                gradkdotgradomega = 0.0;
                for (index_t i = 0; i < dim; i++)
                    gradkdotgradomega += KSolDers[1](i, k) * OSolDers[1](i, k);
                CDkomega(k) = math::max(2 * sigmaO2 / solOVals(k) * gradkdotgradomega, math::pow(10, -10));
            }

            // evaluate F1
            for (index_t k = 0; k < nQuPoints; k++)
            {
                F1(k) = math::tanh(math::pow(math::min(math::max((math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), (4 * sigmaO2 * solKVals(k))/(CDkomega(k) * math::pow(Distance(k), 2))), 4));
                F1(k) = math::max(F1(k), 0.0);
                F1(k) = math::min(F1(k), 1.0);
            }
        }
        else
        {
            F1 = m_paramsPtr->getSSTModel().getF1Vals();
        }

        real_t beta1 = m_paramsPtr->getSSTModel().get_beta1();
        real_t beta2 = m_paramsPtr->getSSTModel().get_beta2();
        for (index_t k = 0; k < nQuPoints; k++)
            m_coeff(k) = (beta1 * F1(k) + beta2 * (1 - F1(k))) * m_solOVals(k);
    }
}

// =====================================================================================================================

template<class T>
void gsTMTerm_BlendCoeffRhs<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    index_t dim = mapData.dim.first;
    real_t sigmaO2;
    gsVector<T> F1;
    gsVector<T> OSolVals;
    std::vector< gsMatrix<T> > KSolDers;
    std::vector< gsMatrix<T> > OSolDers;
    
    if (!(m_paramsPtr->getSSTModel().isCurrent()))
    {
        real_t sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
        real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
        real_t visc = m_paramsPtr->getPde().viscosity();
        
        gsField<T> KSolField = m_paramsPtr->getKSolution();
        gsField<T> OSolField = m_paramsPtr->getOmegaSolution();

        // evaluate k, omega
        gsVector<T> solKVals(dim, nQuPoints);
        solKVals = KSolField.value(mapData.points, mapData.patchId);
        gsVector<T> solOVals(dim, nQuPoints);
        solOVals = OSolField.value(mapData.points, mapData.patchId);

        //evaluate grad(k), grad(omega)
        std::vector< gsMatrix<T> > KSolDers = KSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
        std::vector< gsMatrix<T> > OSolDers = OSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);

        // UPRAVIT !!! evaluate distance
        gsVector<T> Distance(nQuPoints);
        //Distance = DistanceField.value(mapData.points, mapData.patchId);
        for (index_t k = 0; k < nQuPoints; k++)
            Distance(k) = 1.0;

        //evaluate CDkomega
        real_t gradkdotgradomega;
        gsVector<T> CDkomega(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            gradkdotgradomega = 0.0;
            for (index_t i = 0; i < dim; i++)
                gradkdotgradomega += KSolDers[1](i, k) * OSolDers[1](i, k);
            CDkomega(k) = math::max(2 * sigmaO2 / solOVals(k) * gradkdotgradomega, math::pow(10, -10));
        }

        // evaluate F1
        for (index_t k = 0; k < nQuPoints; k++)
        {
            F1(k) = math::tanh(math::pow(math::min(math::max((math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), (4 * sigmaO2 * solKVals(k))/(CDkomega(k) * math::pow(Distance(k), 2))), 4));
            F1(k) = math::max(F1(k), 0.0);
            F1(k) = math::min(F1(k), 1.0);
        }
    }
    else
    {
        sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
        F1 = m_paramsPtr->getSSTModel().getF1Vals();
        OSolVals = m_paramsPtr->getSSTModel().getOSolVals();
        KSolDers = m_paramsPtr->getSSTModel().getKSolDers();
        OSolDers = m_paramsPtr->getSSTModel().getOSolDers();
    }

    m_rhsVals.resize(1, nQuPoints);    
    for (index_t k = 0; k < nQuPoints; k++)
    {
        m_rhsVals(0, k) = 2 * (1 - F1(k)) * sigmaO2 / OSolVals(k) * (KSolDers[1].col(k).dot(OSolDers[1].col(k)));
    }
}

template<class T>
void gsTMTerm_BlendCoeffRhs<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
{ 
    const index_t nQuPoints = quWeights.rows();

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
    }
}

// ============================================================================================================================

template<class T>
void gsTMTerm_ProductionRhs<T>::evalCoeff(const gsMapData<T>& mapData)
{
    index_t nQuPoints = mapData.points.cols();
    index_t dim = mapData.dim.first;
    real_t a1 = m_paramsPtr->getSSTModel().get_a1();
    real_t sigmaO1 = m_paramsPtr->getSSTModel().get_sigmaO1();
    real_t sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
    real_t beta1 = m_paramsPtr->getSSTModel().get_beta1();
    real_t beta2 = m_paramsPtr->getSSTModel().get_beta2();
    real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
    real_t kappa = m_paramsPtr->getSSTModel().get_kappa();
    real_t visc = m_paramsPtr->getPde().viscosity();
    gsVector<T> turbViscosityVals;
    gsVector<T> KSolVals;
    gsVector<T> OSolVals;
    gsVector<T> F1;
    std::vector< gsMatrix<T> > StrainRateTensor;
    std::vector< gsMatrix<T> > USolDers;
    
    if (!(m_paramsPtr->getSSTModel().isCurrent()))
    {
        gsField<T> USolField = m_paramsPtr->getVelocitySolution();
        gsField<T> KSolField = m_paramsPtr->getKSolution();
        gsField<T> OSolField = m_paramsPtr->getOmegaSolution();

        // evaluate k, omega
        gsVector<T> solKVals(dim, nQuPoints);
        solKVals = KSolField.value(mapData.points, mapData.patchId);
        gsVector<T> solOVals(dim, nQuPoints);
        solOVals = OSolField.value(mapData.points, mapData.patchId);

        //evaluate grad(k), grad(omega)
        std::vector< gsMatrix<T> > KSolDers = KSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
        std::vector< gsMatrix<T> > OSolDers = OSolField.function(mapData.patchId).evalAllDers(mapData.points, 1);

        // evaluate strainrate tensor S
        std::vector< gsMatrix<T> > USolDers = USolField.function(mapData.patchId).evalAllDers(mapData.points, 1);
        gsVector<T> StrainRateMag(nQuPoints);
        StrainRateMag.setZero();
        real_t Sij;
        for (index_t k = 0; k < nQuPoints; k++)
        {
            gsMatrix<T> SS(dim, dim);
            SS.setZero();
            for (index_t i = 0; i < dim; i++)
                for (index_t j = 0; j < dim; j++)
                {
                    Sij = 0.5 * (USolDers[1](i * dim + j, k) + USolDers[1](j * dim + i, k));
                    SS(i, j) = Sij;
                    StrainRateMag(k) += 2 * Sij * Sij;
                }
            StrainRateTensor.push_back(SS);
        }

        // UPRAVIT !!! evaluate distance
        gsVector<T> Distance(nQuPoints);
        //Distance = DistanceField.value(mapData.points, mapData.patchId);
        for (index_t k = 0; k < nQuPoints; k++)
            Distance(k) = 1.0;

        // evaluate F2
        gsVector<T> F2(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            F2(k) = math::tanh(math::pow(math::max((2 * math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), 2));
            F2(k) = math::max(F2(k), 0.0);
            F2(k) = math::min(F2(k), 1.0);
        }

        // evaluate turbulent viscosity
        for (index_t k = 0; k < nQuPoints; k++)
            turbViscosityVals(k) = (a1 * solKVals(k)) / (math::max(a1 * solOVals(k), StrainRateMag(k) * F2(k)));

        //evaluate CDkomega
        real_t gradkdotgradomega;
        gsVector<T> CDkomega(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            gradkdotgradomega = 0.0;
            for (index_t i = 0; i < dim; i++)
                gradkdotgradomega += KSolDers[1](i, k) * OSolDers[1](i, k);
            CDkomega(k) = math::max(2 * sigmaO2 / solOVals(k) * gradkdotgradomega, math::pow(10, -10));
        }

        // evaluate F1
        for (index_t k = 0; k < nQuPoints; k++)
        {
            F1(k) = math::tanh(math::pow(math::min(math::max((math::sqrt(solKVals(k)))/(betaStar * solOVals(k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(k))), (4 * sigmaO2 * solKVals(k))/(CDkomega(k) * math::pow(Distance(k), 2))), 4));
            F1(k) = math::max(F1(k), 0.0);
            F1(k) = math::min(F1(k), 1.0);
        }

        m_paramsPtr->getSSTModel().setKSolVals(solKVals);
        m_paramsPtr->getSSTModel().setOSolVals(solOVals);
        m_paramsPtr->getSSTModel().setKSolDers(KSolDers);
        m_paramsPtr->getSSTModel().setKSolDers(OSolDers);
        m_paramsPtr->getSSTModel().setF1Vals(F1);
        m_paramsPtr->getSSTModel().setF2Vals(F2);
        m_paramsPtr->getSSTModel().setTurbulentViscosityVals(turbViscosityVals);
        m_paramsPtr->getSSTModel().setStrainRateMagVals(StrainRateMag);
        m_paramsPtr->getSSTModel().StrainRateTensor(StrainRateTensor);
        m_paramsPtr->getSSTModel().setCurrent();
    }
    else
    {
        KSolVals = m_paramsPtr->getSSTModel().getKSolVals();
        OSolVals = m_paramsPtr->getSSTModel().getOSolVals();
        F1 = m_paramsPtr->getSSTModel().getF1Vals();
        turbViscosityVals = m_paramsPtr->getSSTModel().getTurbulentViscosityVals();
        USolDers = m_paramsPtr->getSSTModel().getUSolDers();
        StrainRateTensor = m_paramsPtr->getSSTModel().getStrainRateTensor();
    }

    m_rhsVals.resize(1, nQuPoints); 
    
    for (index_t k = 0; k < nQuPoints; k++)
    {
        real_t StrainRateUDers = 0.0;
        for (index_t i = 0; i < dim; i++)
            for (index_t j = 0; j < dim; j++)
                StrainRateUDers += StrainRateTensor[k](i,j) * USolDers[1](i*dim+j, k);
        m_rhsVals(0, k) = math::min(2 * turbViscosityVals(k) * StrainRateUDers, 10 * betaStar * KSolVals(k) * OSolVals(k));
    }

    if (m_unknown == 1)
    {
        real_t gamma1, gamma2;

        gamma1 = beta1/betaStar - (sigmaO1 * math::pow(kappa,2))/(math::sqrt(betaStar));
        gamma2 = beta2/betaStar - (sigmaO2 * math::pow(kappa,2))/(math::sqrt(betaStar));

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

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * mapData.measure(k);

        localMat.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
    }
}



} // namespace gismo