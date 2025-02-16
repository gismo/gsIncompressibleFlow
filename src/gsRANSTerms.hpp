/** @file gsFlowTerms.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsRANSTerms.h>

namespace gismo
{

template<class T>
void gsRANSTerm_SymmetricGradientDiag<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;
    
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    gsMatrix<T> lMat(nQuPoints, nQuPoints); // add one matrix to localMat for matrix A_nuT
    lMat.setZero();
    localMat.push_back(lMat);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        // Block A with turbulent viscosity, stored in localMat[0]
        localMat[0].noalias() += weight * (m_turbViscosityVals(k)) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);

        // NEMELA BY TADY A DALE U BLOKU EIJ BYT JEN TURBULENTNI VISKOZITA, BEZ KLASICKE VISKOZITY ???
        // Local blocks for symmetric gradient and corresponding diagonal blocks Eii, where Eii is stored in localMat[i], i=1,...,dim 
        for (index_t s = 0; s != dim; ++s)
            localMat[s+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(s).transpose() * shapeFunPhysGrad.row(s));

    }
}

// =====================================================================================================

template<class T>
void gsRANSTerm_SymmetricGradientOffdiag<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;
    
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    gsMatrix<T> lMat(nQuPoints, nQuPoints); // add one matrix to localMat for matrix A_nuT
    lMat.setZero();
    localMat.push_back(lMat);

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        // NEMELA BY TADY A DALE U BLOKU EIJ BYT JEN TURBULENTNI VISKOZITA, BEZ KLASICKE VISKOZITY ???
        // Local blocks for symmetric gradient and corresponding off-diagonal blocks Ekl, where E12 is stored in localMat[0], and, if dim=3, E13 in localMat[1] and E23 in localMat[2]
        localMat[0].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(1).transpose() * shapeFunPhysGrad.row(0)); 

        if (dim == 3)
        {
            localMat[1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(0));
            localMat[2].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(1));
        }
    }
}

} // namespace gismo