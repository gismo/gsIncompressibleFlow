/** @file gsTMSolverSST.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

namespace gismo
{

// upravit
template<class T, int MatOrder>
void gsTMSolverSST<T, MatOrder>::initMembers()
{
    Base::initMembers();
}

// upravit
template<class T, int MatOrder>
void gsTMSolverSST<T, MatOrder>::evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t patchId)
{
    //for (index_t i = 0; i < quNodes.cols(); i++)
    //    m_TurbulentViscosityVals(i) = 0.01;

    m_TurbulentViscosityVals.setZero(quNodes.cols());
    if ((getAssembler()->isInitialized()) && (m_TMModelPtr->isInitialized()))
    {
        m_TMModelPtr->evalTurbulentViscosity(quNodes, patchId);
        m_TurbulentViscosityVals = m_TMModelPtr->getTurbulentViscosityVals();
    }
    
    
    //if ((getAssembler()->isInitialized()) && (getAssembler()->getSSTModelEvaluator()->isInitialized()))
    //    m_TurbulentViscosityVals = getAssembler()->getSSTModelEvaluator()->evalTurbulentViscosity(quNodes, patchId);
    
    // if ((getAssembler()->isInitialized()) && (m_solution.sum()))
    // {

    // real_t a1 = m_paramsPtr->getSSTModel().get_a1();
    // real_t betaStar = m_paramsPtr->getSSTModel().get_betaStar();
    // real_t visc = m_paramsPtr->getPde().viscosity();
    // index_t nQuPoints = quNodes.cols();
    // index_t dim = quNodes.rows();
    
    // m_TurbulentViscosityVals.setZero(quNodes.cols());
    
    // gsField<T> USolField = m_paramsPtr->getVelocitySolution();
    // gsField<T> KSolField = m_paramsPtr->getKSolution();
    // gsField<T> OSolField = m_paramsPtr->getOmegaSolution();
    
    // // evaluate k, omega
    // gsMatrix<T> solKVals(1, nQuPoints);
    // solKVals = KSolField.value(quNodes, patchId);
    // gsMatrix<T> solOVals(1, nQuPoints);
    // solOVals = OSolField.value(quNodes, patchId);
    
    //         //evaluate grad(k), grad(omega)
    //         //gsFunction<T> KSol = KSolField.function(mapData.patchId);
    //         //std::vector< gsMatrix<T> > KSolDers = KSol.evalAllDers(mapData.points, 1);
    //         //gsFunction<T> OSol = OSolField.function(mapData.patchId);
    //         //std::vector< gsMatrix<T> > OSolDers = OSol.evalAllDers(mapData.points, 1);
    
    // // evaluate strainrate tensor S
    // std::vector< gsMatrix<T> > USolDers = USolField.function(patchId).evalAllDers(quNodes, 1);
    // std::vector< gsMatrix<T> > StrainRateTensor;
    // gsVector<T> StrainRateMag(nQuPoints);
    // StrainRateMag.setZero();
    // real_t Sij;
    // for (index_t k = 0; k < nQuPoints; k++)
    // {
    //     gsMatrix<T> SS(dim, dim);
    //     SS.setZero();
    //     for (index_t i = 0; i < dim; i++)
    //         for (index_t j = 0; j < dim; j++)
    //         {
    //             Sij = 0.5 * (USolDers[1](i * dim + j, k) + USolDers[1](j * dim + i, k));
    //             SS(i, j) = Sij;
    //             StrainRateMag(k) += 2 * Sij * Sij;
    //         }
    //     StrainRateTensor.push_back(SS);
    // }
    
    // // UPRAVIT !!! evaluate distance
    // gsVector<T> Distance(nQuPoints);
    // //Distance = DistanceField.value(mapData.points, mapData.patchId);
    // for (index_t k = 0; k < nQuPoints; k++)
    //     Distance(k) = 1.0;
    
    // // evaluate F2
    // gsVector<T> F2(nQuPoints);
    // for (index_t k = 0; k < nQuPoints; k++)
    // {
    //     F2(k) = math::tanh(math::pow(math::max((2 * math::sqrt(solKVals(0, k)))/(betaStar * solOVals(0, k) * Distance(k)), (500 * visc)/(math::pow(Distance(k), 2) * solOVals(0, k))), 2));
    //     F2(k) = math::max(F2(k), 0.0);
    //     F2(k) = math::min(F2(k), 1.0);
    // }
    
    // // evaluate turbulent viscosity
    // for (index_t k = 0; k < nQuPoints; k++)
    //     m_TurbulentViscosityVals(k) = (a1 * solKVals(0, k)) / (math::max(a1 * solOVals(0, k), StrainRateMag(k) * F2(k)));

    // // TEST ONLY


    // }
    // else
    // {
    //     for (index_t k = 0; k < quNodes.cols(); k++)
    //         m_TurbulentViscosityVals(k) = 0.0;
    // }

    /*
    GISMO_ENSURE(m_pTMsolver != NULL, "uwbRANSBlockVisitor: No turbulent model set!");
        m_sTMEvaluatorType = m_pTMsolver->getTMEvaluator();
        GISMO_ASSERT(m_sTMEvaluatorType != "", "No evaluator type set.");

        int numTMvar = m_pTMsolver->getAssembler()->getNumVar();

        typename uwbTMEvaluator<T>::uPtr evaluatorTM = uwbTMEvaluator<T>::make(m_sTMEvaluatorType);

        evaluatorTM->initialize(m_viscosity, quNodes.cols());
        evaluatorTM->setKOmegaVariant(m_sTMEvaluatorType);

        gsField<T> solTMfield = m_pTMsolver->constructSolution();
        gsMatrix<T> solKOVals = solTMfield.value(quNodes, m_patchIndex);

        //--- KOGrads
        gsMatrix<T> solActKOCoeffs, physGradKO, bGradsKO;
        gsMatrix<index_t> activesKO;

        const gsMultiBasis<T> basisKO = m_pTMsolver->getAssembler()->getBlockAssembler().getSolBasis();
        basisKO.back().active_into(quNodes.col(0), activesKO);
        basisKO.back().deriv_into(quNodes, bGradsKO);

        const index_t numActKO = activesKO.rows();

        solActKOCoeffs.setZero(numTMvar, numActKO);
        for (int j = 0; j < numActKO; j++)
            solActKOCoeffs.col(j) = solTMfield.coefficientVector(m_patchIndex).row(activesKO(j)).transpose();

        index_t nQuPoints = quNodes.cols();
        std::vector<gsMatrix<T> > solKOGrads(nQuPoints);
        for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
        {
            geoEval.transformGradients(k, bGradsKO, physGradKO);
            solKOGrads[k].noalias() = solActKOCoeffs * physGradKO.transpose();
        }

        if (this->checkWallDistanceBasedTM())
        {
            gsField<T> solPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonSolution();
            gsMultiPatch<T> patchesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonPatches();
            gsMultiBasis<T> basesPoisson = m_pTMsolver->getAssembler()->getBlockAssembler().getPoissonBasis();

            unsigned evFlagsPoisson = NEED_MEASURE | NEED_GRAD_TRANSFORM;
            const gsGeometry<T>& geo = patchesPoisson.patch(m_patchIndex);
            typename gsGeometryEvaluator<T>::uPtr geoEvalPoisson(getEvaluator(evFlagsPoisson, geo));

            gsMatrix<index_t> activesPoisson;
            gsMatrix<T> basisGradsPoisson, physGradPoisson;
            basesPoisson.basis(m_patchIndex).active_into(quNodes.col(0), activesPoisson);
            const index_t numActPoisson = activesPoisson.rows();
            basesPoisson.basis(m_patchIndex).deriv_into(quNodes, basisGradsPoisson);

            geoEvalPoisson->evaluateAt(quNodes);

            gsMatrix<T> solPoissonVals = solPoisson.value(quNodes, m_patchIndex);
            gsMatrix<T> solActPoissonCoeffs;
            solActPoissonCoeffs.setZero(1, numActPoisson);
            for (int j = 0; j < numActPoisson; j++)
                solActPoissonCoeffs.col(j) = solPoisson.coefficientVector(m_patchIndex).row(activesPoisson(j)).transpose();

            std::vector<gsMatrix<T> > solPoissonGrads(nQuPoints);
            for (index_t k = 0; k < nQuPoints; ++k) // loop over quadrature nodes
            {
                geoEvalPoisson->transformGradients(k, basisGradsPoisson, physGradPoisson);
                solPoissonGrads[k].noalias() = solActPoissonCoeffs * physGradPoisson.transpose();
            }

            evaluatorTM->evalWallDistance(solPoissonVals, solPoissonGrads);
        }

        evaluatorTM->initAtElement(solUGrads, solKOVals, solKOGrads);
        evaluatorTM->evalQuantities_turbViscosity();
        m_turbViscosityVals = evaluatorTM->getTurbViscosity();

        for (int i = 0; i < quNodes.cols(); i++)
            m_diffusionCoeff(0, i) = m_viscosity + m_turbViscosityVals(i);

        m_bEffectiveViscSet = true;
    */
}


}