/** @file gsTMModels.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMModels.h>

namespace gismo
{

template <class T>
typename gsTMModelData<T>::Ptr gsTMModelData<T>::make(typename gsFlowSolverParams<T>::Ptr paramsPtr)
{
    std::string turbModel = paramsPtr->options().getString("TM");
    if (turbModel == "SST") 
    {
        return gsTMModelData_SST<T>::make(paramsPtr);
    }
    //elseif (m_paramsPtr->options().getSwitch("TM.eval") == "SA") 
    //{ }
    else 
    {
        paramsPtr->logger() << "Unknown identifier of a turbulent model entered! Using k-omega SST model.\n";
        return gsTMModelData_SST<T>::make(paramsPtr);
    }
}


template <class T>
void gsTMModelData<T>::plotTurbulentViscosity(typename gsFlowSolverParams<T>::Ptr paramsPtr, std::string str)
{
    gsMultiPatch<T> patches = paramsPtr->getPde().patches();    // multipatch representing the computational domain
    gsMultiBasis<T> basis = paramsPtr->getBasis(1);
    
    size_t np = patches.nPatches();
    gsMultiPatch<T>* turbViscMP = new gsMultiPatch<T>;
    for (size_t i = 0; i < np; i++)
    {
        index_t patchId = i;
        const gsBasis<T> & basisp = basis.piece(patchId);

        std::vector< gsVector<T> > rr;
        rr.reserve(patches.parDim());

        for (short_t j = 0; j < patches.parDim(); ++j)            // computing grid of point
        {
            rr.push_back(basisp.component(j).anchors().transpose());
        }
        gsMatrix<T> gridPts = gsPointGrid<T>(rr);

        evalTurbulentViscosity(gridPts, 1, patchId);
        gsVector<T> turbViscVals = getTurbulentViscosityVals();

        typename gsGeometry<T>::uPtr geo = basisp.interpolateAtAnchors(turbViscVals.transpose());    // interpolating distances at grid points 
        const gsMatrix<T> & turbViscCoeffs = geo->coefs();
        turbViscMP->addPatch(basisp.makeGeometry(turbViscCoeffs));
    }

    gsField<T> result = gsField<T>(paramsPtr->getPde().patches(), typename gsFunctionSet<T>::Ptr(turbViscMP), true);
    gsWriteParaview<T>(result, str, 10000);
}

// ============================================================================================================================

template <class T>
void gsTMModelData_SST<T>::evalDistance(gsMatrix<T>& quNodes, index_t patchId)
{
    index_t nQuPoints = quNodes.cols();
    gsField<T> distanceField; 
    distanceField = m_paramsPtr->getDistanceField();
    gsMatrix<T> Distance(1, nQuPoints);
    Distance = distanceField.value(quNodes, patchId);
    m_distance.setZero(nQuPoints);
    for (index_t k = 0; k < nQuPoints; k++)
        m_distance(k) = math::max(Distance(0, k), m_eps);
}

template <class T>
void gsTMModelData_SST<T>::evalVelocityQuantities(gsMatrix<T>& quNodes, index_t patchId)
{
    index_t nQuPoints = quNodes.cols();
    index_t dim = quNodes.rows();
    gsMultiBasis<T> basis = m_paramsPtr->getBasis(0);
    gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    
    gsMapData<T> mapData;
    unsigned geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    mapData.flags = geoFlags;
    mapData.patchId = patchId;
    mapData.points = quNodes;
    m_paramsPtr->getPde().patches().patch(patchId).computeMap(mapData);

    gsMatrix<index_t> actives;
    gsMatrix<T> parGrads, physGrad;
    basis.piece(patchId).deriv_into(quNodes, parGrads);
        
    gsMatrix<T> USolCoeffVec = USolField.coefficientVector(patchId);
    std::vector<gsMatrix<T> > USolDers(nQuPoints);
    for (index_t k = 0; k < nQuPoints; k++)
    {
        basis.piece(patchId).active_into(quNodes.col(k), actives);
        int numAct = actives.rows();
        gsMatrix<T> USolActCoeffs(dim, numAct);
        for (int j = 0; j < numAct; j++)
            USolActCoeffs.col(j) = USolCoeffVec.row(actives(j, 0)).transpose();

        transformGradients(mapData, k, parGrads, physGrad);
        USolDers[k].noalias() = USolActCoeffs * physGrad.transpose();
    }

    gsVector<T> StrainRateMag(nQuPoints);
    StrainRateMag.setZero();
    std::vector< gsMatrix<T> > StrainRateTensor;
    StrainRateTensor.resize(nQuPoints);
    for (index_t i = 0; i < nQuPoints; i++)
        StrainRateTensor[i].resize(dim, dim);
    real_t Sij;
    for (index_t k = 0; k < nQuPoints; k++)
    {
        gsMatrix<T> SS(dim, dim);
        SS.setZero();
        for (index_t i = 0; i < dim; i++)
            for (index_t j = 0; j < dim; j++)
            {
                Sij = 0.5 * (USolDers[k](i, j) + USolDers[k](j, i));
                SS(i, j) = Sij;
                StrainRateMag(k) += 2 * Sij * Sij;
            }
        StrainRateTensor[k] = SS;
        StrainRateMag(k) = math::sqrt(StrainRateMag(k));
    }
    m_USolDers = USolDers;
    m_StrainRateMag = StrainRateMag;
    m_StrainRateTensor = StrainRateTensor;
}

template <class T>
void gsTMModelData_SST<T>::evalKSol(gsMatrix<T>& quNodes, index_t patchId, index_t der)
{
    index_t nQuPoints = quNodes.cols();
    gsMultiBasis<T> basis = m_paramsPtr->getBasis(2);
    gsField<T> KSolField = m_paramsPtr->getKSolution();
    
    m_KSolVals.resize(1, nQuPoints);
    m_KSolVals = KSolField.function(patchId).eval(quNodes);
    
    if (der > 0)
    {
        gsMapData<T> mapData;
        unsigned geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        mapData.flags = geoFlags;
        mapData.patchId = patchId;
        mapData.points = quNodes;
        m_paramsPtr->getPde().patches().patch(patchId).computeMap(mapData);

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        basis.piece(patchId).deriv_into(quNodes, parGrads);
            
        gsMatrix<T> KSolCoeffVec = KSolField.coefficientVector(patchId);
        std::vector<gsMatrix<T> > KSolDers(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            basis.piece(patchId).active_into(quNodes.col(k), actives);
            int numAct = actives.rows();
            gsMatrix<T> KSolActCoeffs(1, numAct);
            for (int j = 0; j < numAct; j++)
                KSolActCoeffs(0, j) = KSolCoeffVec(actives(j, 0), 0);

            transformGradients(mapData, k, parGrads, physGrad);
            KSolDers[k].noalias() = KSolActCoeffs * physGrad.transpose();
        }

        m_KSolDers = KSolDers;
    }   
}

template <class T>
void gsTMModelData_SST<T>::evalOSol(gsMatrix<T>& quNodes, index_t patchId, index_t der)
{
    index_t nQuPoints = quNodes.cols();
    gsMultiBasis<T> basis = m_paramsPtr->getBasis(3);
    gsField<T> OSolField = m_paramsPtr->getOmegaSolution();
    
    m_OSolVals.resize(1, nQuPoints);
    m_OSolVals = OSolField.function(patchId).eval(quNodes);
    
    if (der > 0)
    {
        gsMapData<T> mapData;
        unsigned geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        mapData.flags = geoFlags;
        mapData.patchId = patchId;
        mapData.points = quNodes;
        m_paramsPtr->getPde().patches().patch(patchId).computeMap(mapData);

        gsMatrix<index_t> actives;
        gsMatrix<T> parGrads, physGrad;
        basis.piece(patchId).deriv_into(quNodes, parGrads);
            
        gsMatrix<T> OSolCoeffVec = OSolField.coefficientVector(patchId);
        std::vector<gsMatrix<T> > OSolDers(nQuPoints);
        for (index_t k = 0; k < nQuPoints; k++)
        {
            basis.piece(patchId).active_into(quNodes.col(k), actives);
            int numAct = actives.rows();
            gsMatrix<T> OSolActCoeffs(1, numAct);
            for (int j = 0; j < numAct; j++)
                OSolActCoeffs(0, j) = OSolCoeffVec(actives(j, 0), 0);

            transformGradients(mapData, k, parGrads, physGrad);
            OSolDers[k].noalias() = OSolActCoeffs * physGrad.transpose();
        }

        m_OSolDers = OSolDers;
    }   
}

template <class T>
void gsTMModelData_SST<T>::evalF1(gsMatrix<T>& quNodes, index_t patchId)
{
    index_t nQuPoints = quNodes.cols();
    index_t dim = quNodes.rows();
    real_t gradkdotgradomega;
    gsVector<T> CDkomega(nQuPoints);
    CDkomega.setZero();
    for (index_t k = 0; k < nQuPoints; k++)
    {
        gradkdotgradomega = 0.0;
        for (index_t i = 0; i < dim; i++)
            gradkdotgradomega += m_KSolDers[k](i) * m_OSolDers[k](i);
        CDkomega(k) = math::max(2 * m_sigmaO2 / math::max(m_OSolVals(0, k), m_eps) * gradkdotgradomega, math::pow(10, -10));
    }

    gsVector<T> F1(nQuPoints);
    F1.setZero();
    for (index_t k = 0; k < nQuPoints; k++)
    {
        F1(k) = math::tanh(math::pow(math::min(math::max((math::sqrt(math::max(m_KSolVals(0, k), m_eps)))/(m_betaStar * math::max(m_OSolVals(0, k), m_eps) * m_distance(k)), (500 * m_visc)/(math::pow(m_distance(k), 2) * math::max(m_OSolVals(0, k), m_eps))), (4 * m_sigmaO2 * m_KSolVals(0, k))/(CDkomega(k) * math::pow(m_distance(k), 2))), 4));
        F1(k) = math::max(F1(k), 0.0);
        F1(k) = math::min(F1(k), 1.0);
    }
    m_F1 = F1;
}

template <class T>
void gsTMModelData_SST<T>::evalF2(gsMatrix<T>& quNodes, index_t patchId)
{
    index_t nQuPoints = quNodes.cols();
    gsVector<T> F2(nQuPoints);
    F2.setZero();
    for (index_t k = 0; k < nQuPoints; k++)
    {
        F2(k) = math::tanh(math::pow(math::max((2 * math::sqrt(math::max(m_KSolVals(0, k), m_eps)))/(m_betaStar * math::max(m_OSolVals(0, k), m_eps) * m_distance(k)), (500 * m_visc)/(math::pow(m_distance(k), 2) * math::max(m_OSolVals(0, k), m_eps))), 2));
        F2(k) = math::max(F2(k), 0.0);
        F2(k) = math::min(F2(k), 1.0);
    }
    m_F2 = F2;
}

template <class T>
void gsTMModelData_SST<T>::evalTurbViscFromData(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
{
    index_t nQuPoints = quNodes.cols();
    gsVector<T> turbulentViscosityVals(nQuPoints);
    turbulentViscosityVals.setZero();
    for (index_t k = 0; k < nQuPoints; k++)
    {
        turbulentViscosityVals(k) = (m_a1 * m_KSolVals(0, k)) / (math::max(m_a1 * math::max(m_OSolVals(0, k), m_eps), m_StrainRateMag(k) * m_F2(k)));
        turbulentViscosityVals(k) = math::max(turbulentViscosityVals(k), m_eps);
    }

    if (m_average && numNodesPerElem != 1)
    {
        GISMO_ASSERT(nQuPoints % numNodesPerElem == 0, "Total number of quad nodes is not multiple of number of quad nodes per element!");

        index_t nElements = nQuPoints / numNodesPerElem;
        for (index_t e = 0; e < nElements; e++)
        {
            T avg = turbulentViscosityVals.middleRows(e * numNodesPerElem, numNodesPerElem).sum() / numNodesPerElem;
            turbulentViscosityVals.middleRows(e * numNodesPerElem, numNodesPerElem).setConstant(avg);
        }
    }

    m_turbulentViscosityVals = turbulentViscosityVals;
}

template <class T>
void gsTMModelData_SST<T>::updateModel(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
{
    // evaluate k, omega, grad(k), grad(omega)
    evalKSol(quNodes, patchId, 1);
    evalOSol(quNodes, patchId, 1);

    // evaluate strainrate tensor S
    evalVelocityQuantities(quNodes, patchId);

    // evaluate distance
    evalDistance(quNodes, patchId);

    // evaluate F2
    evalF2(quNodes, patchId);

    // evaluate turbulent viscosity
    evalTurbViscFromData(quNodes, numNodesPerElem, patchId);
    
    // evaluate F1
    evalF1(quNodes, patchId);
    
    m_isInitialized = true;
}

template <class T>
void gsTMModelData_SST<T>::evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
{
    // evaluate k, omega, grad(k), grad(omega)
    evalKSol(quNodes, patchId, 0);
    evalOSol(quNodes, patchId, 0);

    // evaluate strainrate tensor S
    evalVelocityQuantities(quNodes, patchId);

    // evaluate distance
    evalDistance(quNodes, patchId);

    // evaluate F2
    evalF2(quNodes, patchId);

    // evaluate turbulent viscosity
    evalTurbViscFromData(quNodes, numNodesPerElem, patchId);

}

} // namespace gismo