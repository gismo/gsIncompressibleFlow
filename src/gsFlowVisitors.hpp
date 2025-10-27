/** @file gsFlowVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>

namespace gismo
{

template <class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::gatherEvalFlags()
{
    m_geoFlags = 0;
    m_testFunFlags = 0;
    m_trialFunFlags = 0;

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->updateEvalFlags(m_geoFlags, m_testFunFlags, m_trialFunFlags);
}    


template <class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::evalSingleFunData(const unsigned& basisFlags, const gsBasis<T>* basisPtr, const index_t funID, std::vector< gsMatrix<T> >& basisData)
{
    if(basisFlags & NEED_VALUE)
        basisPtr->evalSingle_into(funID, m_quNodes, basisData[0]);

    if(basisFlags & NEED_DERIV)
        basisPtr->derivSingle_into(funID, m_quNodes, basisData[1]);

    // currently not needed
    // if(basisFlags & NEED_DERIV2)
    //     basisPtr->deriv2Single_into(funID, m_quNodes, basisData[2]);
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::evalBasisData(const unsigned& basisFlags, const gsBasis<T>* basisPtr, gsMatrix<index_t>& activesUnique, std::vector< gsMatrix<T> >& basisData)
{
    gsMatrix<index_t> actives;
    basisPtr->active_into(m_quNodes, actives);
    activesUnique = createVectorOfUniqueIndices(actives);

    bool multipleElem = (activesUnique.rows() > actives.rows());
    std::unordered_map<int, int> activesUnique_val_to_ID;

    if (multipleElem)
    {
        for (int i = 0; i < activesUnique.size(); ++i)
            activesUnique_val_to_ID[activesUnique(i)] = i;
    }

    index_t dim = basisPtr->dim();
    index_t numAct = activesUnique.rows();

    if(basisFlags & NEED_VALUE)
    {
        gsMatrix<real_t> basisVals;
        basisPtr->eval_into(m_quNodes, basisVals);

        if (multipleElem)
        {
            basisData[0].setZero(numAct, m_quNodes.cols());

            for (index_t i = 0; i < actives.rows(); i++)
                for (index_t j = 0; j < actives.cols(); j++)
                    basisData[0](activesUnique_val_to_ID[actives(i, j)], j) = basisVals(i,j);
        }
        else
        {
            basisData[0] = basisVals;
        }
    }

    if(basisFlags & NEED_DERIV)
    {
        gsMatrix<real_t> basisDers;
        basisPtr->deriv_into(m_quNodes, basisDers);

        if (multipleElem)
        {
            basisData[1].setZero(dim*numAct, m_quNodes.cols());

            for (index_t i = 0; i < actives.rows(); i++)
            {
                for (index_t j = 0; j < actives.cols(); j++)
                {
                    index_t ii = activesUnique_val_to_ID[actives(i, j)];
                    basisData[1].block(dim*ii, j, dim, 1) = basisDers.block(dim*i, j, dim, 1);
                }
            }
        }
        else
        {
            basisData[1] = basisDers;
        }
    }

    // currently not needed
    // if(basisFlags & NEED_DERIV2)
    // {
    //     gsMatrix<real_t> basisDers2;
    //     basisPtr->deriv2_into(m_quNodes, basisDers2);

    //     if (multipleElem)
    //     {
    //         index_t dimSq = dim*dim;
    //         basisData[2].setZero(dimSq*numAct, m_quNodes.cols());

    //         for (index_t i = 0; i < actives.rows(); i++)
    //         {
    //             for (index_t j = 0; j < actives.cols(); j++)
    //             {
    //                 index_t ii = activesUnique_val_to_ID[actives(i, j)];
    //                 basisData[2].block(dimSq*ii, j, dimSq, 1) = basisDers2.block(dimSq*i, dimSq ,j, 1);
    //             }
    //         }
    //     }
    //     else
    //     {
    //         basisData[2] = basisDers2;
    //     }
    // }
    
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::initialize()
{
    defineTestTrialUnknowns();  
    deleteTerms();
    defineTerms();
    gatherEvalFlags();
    m_mapData.flags = m_geoFlags;
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::initOnPatch(index_t patchID)
{
    m_patchID = patchID;
    m_mapData.patchId = m_patchID;
    defineTestTrialBases();
    setupQuadrature();          

    m_testFunData.clear(); m_trialFunData.clear();
    m_testFunData.resize(2); // 0 - value, 1 - deriv (2nd derivative not needed at the moment)
    m_trialFunData.resize(2); // 0 - value, 1 - deriv
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::setCurrentSolution(std::vector<gsField<T> >& solutions)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsFlowTermNonlin<T>* termPtr = dynamic_cast< gsFlowTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solutions);
    }
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::setCurrentSolution(gsField<T>& solution)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsFlowTermNonlin<T>* termPtr = dynamic_cast< gsFlowTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solution);
    }
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::assemble()
{
    m_localMat.setZero(m_testFunActives.rows(), m_trialFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_localMat);
}

template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
{
    if(!m_hasPeriodicBC)
        localToGlobal_nonper(eliminatedDofs, globalMat, globalRhs);
    else
        localToGlobal_per(eliminatedDofs, globalMat, globalRhs);
}

template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::localToGlobal(gsMatrix<T>& globalRhs)
{
    if(!m_hasPeriodicBC)
        localToGlobal_nonper(globalRhs);
    else
        localToGlobal_per(globalRhs);
}

// ===================================================================================================================

template <class T, int MatOrder>
void gsFlowVisitorVectorValued<T, MatOrder>::assemble()
{
    resizeMatVec();

    for (size_t i = 0; i < m_locMatVec.size(); i++)
        m_locMatVec[i].setZero(m_testFunActives.rows(), m_trialFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_trialFunData, m_locMatVec);
}

// ===================================================================================================================
// ===================================================================================================================
// to be changed:

template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::setupQuadrature()
{
    gsVector<index_t> numQuadNodes(m_paramsPtr->getPde().dim()); 

    index_t maxDegTest = m_testBasisPtr->maxDegree();
    index_t maxDegTrial = m_trialBasisPtr->maxDegree();

    numQuadNodes.setConstant(math::max(maxDegTest, maxDegTrial)+1);

    m_quRule = gsGaussRule<T>(numQuadNodes);
}

template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::evaluate(index_t testFunID)
{
    // trial basis (on the whole support of testFunID)

    index_t dim = m_testBasisPtr->domainDim();
    gsMatrix<T> support = m_testBasisPtr->support(testFunID);

    gsMatrix<T> quNodes; // quad. nodes for the current element
    gsVector<T> quWeights; // weights for the current element
    std::vector< gsMatrix<T> > quNodesOnElem; // quad. nodes for all elements in support
    std::vector< gsVector<T> > quWeightsOnElem; // weights for all elements in support

    gsMatrix<T> lowerCorners, upperCorners; // bounding boxes of all elements in support
    elementsInSupport(*m_testBasisPtr, testFunID, lowerCorners, upperCorners);

    for (index_t i = 0; i < lowerCorners.cols(); i++)
    {
        m_quRule.mapTo((gsVector<T>)(lowerCorners.col(i)), (gsVector<T>)(upperCorners.col(i)), quNodes, quWeights);
        quNodesOnElem.push_back(quNodes);
        quWeightsOnElem.push_back(quWeights);
    }

    size_t numElemInSupport = quNodesOnElem.size();
    index_t numNodesInElem = quNodesOnElem[0].cols(); 
    index_t numNodesInSupport = numElemInSupport * numNodesInElem;
    m_quNodes.resize(dim, numNodesInSupport);
    m_quWeights.resize(numNodesInSupport);

    for (size_t e = 0; e < numElemInSupport; e++)
    {
        m_quNodes.middleCols(e*numNodesInElem, numNodesInElem) = quNodesOnElem[e];
        m_quWeights.middleRows(e*numNodesInElem, numNodesInElem) = quWeightsOnElem[e];
    }

    m_mapData.points = m_quNodes;
    m_paramsPtr->getPde().patches().patch(m_patchID).computeMap(m_mapData);

    evalBasisData(m_trialFunFlags, m_trialBasisPtr, m_trialFunActives, m_trialFunData);

    // test basis
    m_testFunActives.resize(1,1);
    m_testFunActives << testFunID;
    evalSingleFunData(m_testFunFlags, m_testBasisPtr, testFunID, m_testFunData);
}


template<class T, int MatOrder>
void gsFlowVisitor<T, MatOrder>::evaluate(const gsDomainIterator<T>* domIt)
{
    m_quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), m_quNodes, m_quWeights);
    m_mapData.points = m_quNodes;
    m_paramsPtr->getPde().patches().patch(m_patchID).computeMap(m_mapData);

    evalBasisData(m_testFunFlags, m_testBasisPtr, m_testFunActives, m_testFunData);

    if (m_trialUnkID == m_testUnkID)
        m_trialFunActives = m_testFunActives;
    else
        m_trialBasisPtr->active_into(m_quNodes.col(0), m_trialFunActives);

    if ( (m_trialUnkID == m_testUnkID) && (m_testFunFlags == m_trialFunFlags) )
        m_trialFunData = m_testFunData;
    else
        evalBasisData(m_trialFunFlags, m_trialBasisPtr, m_trialFunActives, m_trialFunData);
}

} // namespace gismo
