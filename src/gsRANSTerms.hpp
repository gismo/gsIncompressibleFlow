/** @file gsRANSTerms.hpp

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
void gsRANSTerm_SymmetricGradient<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;

    // localMat = [A, E11, E22, (E33), E12, (E13), (E23)]
    GISMO_ASSERT(localMat.size() == dim * (dim + 1) / 2 + 1, "gsRANSTerm_SymmetricGradient<T>::assemble(): Wrong size of localMat vector.");
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        // block A with turbulent viscosity, stored in localMat[0]
        localMat[0].noalias() += weight * (m_turbViscosityVals(k)) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);

        // local blocks for symmetric gradient and corresponding diagonal blocks Eii, where Eii is stored in localMat[i], i=1,...,dim 
        for (index_t s = 0; s != dim; ++s)
            localMat[s+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(s).transpose() * shapeFunPhysGrad.row(s));

        // local blocks for symmetric gradient and corresponding off-diagonal blocks Ekl, where E12 is stored in localMat[dim+1], and, 
        // if dim=3, E13 in localMat[dim+2] and E23 in localMat[dim+3]
        localMat[dim+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(1).transpose() * shapeFunPhysGrad.row(0)); 
        if (dim == 3)
        {
            localMat[dim+2].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(0));
            localMat[dim+3].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(1));
        }
    }
}    


template<class T>
void gsRANSTerm_SymmetricGradient_full<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;

    // localMat = [A, E11, E12, (E13), E21. E22, (E23), (E31), (E32), (E33)]
    GISMO_ASSERT(localMat.size() == dim * dim + 1, "gsRANSTerm_SymmetricGradient_full<T>::assemble(): Wrong size of localMat vector.");
    
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

        // block A with turbulent viscosity, stored in localMat[0]
        localMat[0].noalias() += weight * (m_turbViscosityVals(k)) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);

        index_t b = 1; // index of the current block in localMat

        // blocks Eij
        for (index_t i = 0; i != dim; ++i)
        {
            for (index_t j = 0; j != dim; ++j)
            {
                localMat[b].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(j).transpose() * shapeFunPhysGrad.row(i));
                ++b;
            }
        }
    }
}    

// ===================================================================================
// For SUPG stabilization

/*template<class T>
void gsRANSTerm_SG_SUPGstabilization_diffusion<T>::assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
{ 
    gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
    index_t dim = mapData.dim.second;

    // localMat = [A, E11, E22, (E33), E12, (E13), (E23)]
    GISMO_ASSERT(localMat.size() == dim * (dim + 1) / 2 + 1, "gsRANSTerm_SG_SUPGstabilization_diffusion<T>::assemble(): Wrong size of localMat vector.");
    const gsMatrix<T>& testFunGrads = testFunData[1];
    const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

    const index_t nQuPoints = quWeights.rows();
    gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad, turbViscPhysGrad;

    for (index_t k = 0; k < nQuPoints; k++)
    {
        const T weight = quWeights(k) * coeffMeasure(k);

        transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
        transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);
        //gsInfo << m_turbViscosityGrads[k] << std::endl;

        // (this->m_solUVals.col(k).transpose() * testFunPhysGrad) ... radek

        // block A with turbulent viscosity, stored in localMat[0]
        //localMat[0].noalias() += weight * (m_turbViscosityVals(k)) * (testFunPhysGrad.transpose() * shapeFunPhysGrad);
        localMat[0].noalias() -= weight * (((m_turbViscosityGrads[k] * shapeFunPhysGrad).transpose()) * (this->m_solUVals.col(k).transpose() * testFunPhysGrad));

        // local blocks for symmetric gradient and corresponding diagonal blocks Eii, where Eii is stored in localMat[i], i=1,...,dim 
        for (index_t s = 0; s != dim; ++s)
        //    localMat[s+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(s).transpose() * shapeFunPhysGrad.row(s));
            localMat[s+1].noalias() -= weight * (((m_turbViscosityGrads[k].col(s) * shapeFunPhysGrad.row(s)).transpose()) * (this->m_solUVals(s, k) * testFunPhysGrad.row(s)));

        // local blocks for symmetric gradient and corresponding off-diagonal blocks Ekl, where E12 is stored in localMat[dim+1], and, 
        // if dim=3, E13 in localMat[dim+2] and E23 in localMat[dim+3]
        //localMat[dim+1].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(1).transpose() * shapeFunPhysGrad.row(0));
        localMat[dim+1].noalias() -= weight * (((m_turbViscosityGrads[k].col(0) * shapeFunPhysGrad.row(0)).transpose()) * (this->m_solUVals(1, k) * testFunPhysGrad.row(1))); 
        if (dim == 3)
        {
        //    localMat[dim+2].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(0));
        //    localMat[dim+3].noalias() += weight * (m_turbViscosityVals(k) + m_viscosity) * (testFunPhysGrad.row(2).transpose() * shapeFunPhysGrad.row(1));
            localMat[dim+2].noalias() -= weight * (((m_turbViscosityGrads[k].col(0) * shapeFunPhysGrad.row(0)).transpose()) * (this->m_solUVals(2, k) * testFunPhysGrad.row(2)));
            localMat[dim+3].noalias() -= weight * (((m_turbViscosityGrads[k].col(1) * shapeFunPhysGrad.row(1)).transpose()) * (this->m_solUVals(2, k) * testFunPhysGrad.row(2)));
        }
    }
}*/

} // namespace gismo