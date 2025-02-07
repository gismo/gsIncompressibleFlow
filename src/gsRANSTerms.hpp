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

// upravit
template<class T>
//void gsRANSTerm_Diffusion<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMatA, std::vector<gsMatrix<T> > localMatDiag)
void gsRANSTerm_Diffusion<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;
    
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        // localMat[i], where i=0 for matrix A, i=1,...,dim for matrices Eii, i=dim+1 for E12 (if dim=2) or i=dim+1,..,2*dim for E12,E13,E23 (if dim=3)
        localMat[0].template triangularView<gsEigen::Upper>() += weight * (m_viscosity + m_turbViscosityVals(k)) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);

        // NEMELA BY TADY A DALE U BLOKU EIJ BYT JEN TURBULENTNI VISKOZITA, BEZ KLASICKE VISKOZITY ???
        // Local block Ediag
        for (index_t s = 0; s != dim; ++s)
            localMat[s+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(s).transpose() * shapeFunPhysGrad.row(s));

        // Local blocks Enondiag
        localMat[dim+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(1).transpose() * shapeFunPhysGrad.row(0)); //dv1/dy * du2/dx

        if (dim == 3)
        {
            localMat[dim+2].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad[k].row(0)); //dv1/dz * du3/dx
            localMat[dim+3].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad[k].row(1)); //dv2/dz * du3/dy
        }
    }
}

} // namespace gismo