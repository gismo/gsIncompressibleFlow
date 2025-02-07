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

template<class T, int MatOrder>
void gsTMSolverSST<T, MatOrder>::evalTurbulentViscosity(/*std::vector<gsMatrix<T> >& solUGrads, */gsMatrix<T>& quNodes/*, gsGeometryEvaluator<T> & geoEval*/)
{
    m_TurbulentViscosityVals.setZero(m_quNodes.cols());
    
    for (index_t i = 0; i < m_quNodes.cols(); i++)
        m_TurbulentViscosityVals(i) = 0.01;   

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