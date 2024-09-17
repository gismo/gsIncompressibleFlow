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

template<class T>
void gsFlowVisitor<T>::gatherEvalFlags()
{
    m_geoFlags = 0;
    m_testFunFlags = 0;
    m_shapeFunFlags = 0;

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->updateEvalFlags(m_geoFlags, m_testFunFlags, m_shapeFunFlags);
}    


template<class T>
void gsFlowVisitor<T>::evalBasisData(const unsigned& basisFlags, const gsBasis<T>* basisPtr, const gsMatrix<index_t>& basisActives, std::vector< gsMatrix<T> >& basisData)
{
    basisData.clear();
    basisData.resize(3); // 0 - value, 1 - deriv, 2 - deriv2

    index_t dim = m_params.getPde().dim();
    index_t numAct = basisActives.rows();

    if(basisFlags & NEED_VALUE)
    {
        basisData[0].setZero(numAct, m_quNodes.cols());

        gsMatrix<T> tmpData;
        
        for(index_t i = 0; i < numAct; i++)
        {
            basisPtr->evalSingle_into(basisActives(i), m_quNodes, tmpData);
            basisData[0].row(i) = tmpData;
        }
    }

    if(basisFlags & NEED_DERIV)
    {
        basisData[1].setZero(dim*numAct, m_quNodes.cols());

        gsMatrix<T> tmpData;

        for(index_t i = 0; i < numAct; i++)
        {
            basisPtr->derivSingle_into(basisActives(i), m_quNodes, tmpData);
            basisData[1].middleRows(dim*i, dim) = tmpData;
        }
    }

    if(basisFlags & NEED_DERIV2)
    {
        index_t dimSq = dim*dim;
        basisData[2].setZero(dimSq*numAct, m_quNodes.cols());

        gsMatrix<T> tmpData;
        
        for(index_t i = 0; i < numAct; i++)
        {
            basisPtr->deriv2Single_into(basisActives(i), m_quNodes, tmpData);
            basisData[2].middleRows(dimSq*i, dimSq) = tmpData;
        }
    }
}


template<class T>
void gsFlowVisitor<T>::initialize()
{
    defineTestShapeUnknowns();  
    m_params.createDofMappers(m_dofMappers);  

    deleteTerms();
    defineTerms();
    gatherEvalFlags();
    m_mapData.flags = m_geoFlags;
}


template<class T>
void gsFlowVisitor<T>::initOnPatch(index_t patchID)
{
    m_patchID = patchID;
    m_mapData.patchId = m_patchID;
    defineTestShapeBases();
    setupQuadrature();          
}


template<class T>
void gsFlowVisitor<T>::setCurrentSolution(std::vector<gsField<T> >& solutions)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsFlowTermNonlin<T>* termPtr = dynamic_cast< gsFlowTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solutions);
    }
}


template<class T>
void gsFlowVisitor<T>::setCurrentSolution(gsField<T>& solution)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsFlowTermNonlin<T>* termPtr = dynamic_cast< gsFlowTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solution);
    }
}


template<class T>
void gsFlowVisitor<T>::assemble()
{
    m_localMat.setZero(m_testFunActives.rows(), m_shapeFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}

// ===================================================================================================================

template <class T>
void gsFlowVisitorVectorValued<T>::assemble()
{
    m_locMatVec.resize(m_params.getPde().dim());

    for (size_t i = 0; i < m_locMatVec.size(); i++)
        m_locMatVec[i].setZero(m_testFunActives.rows(), m_shapeFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_locMatVec);
}

// ===================================================================================================================
// ===================================================================================================================
// to be changed:

template<class T>
void gsFlowVisitor<T>::setupQuadrature()
{
    gsVector<index_t> numQuadNodes(m_params.getPde().dim()); 

    index_t maxDegTest = m_testBasisPtr->maxDegree();
    index_t maxDegShape = m_shapeBasisPtr->maxDegree();

    numQuadNodes.setConstant(math::min(maxDegTest, maxDegShape)+1);

    m_quRule = gsGaussRule<T>(numQuadNodes);
}

template<class T>
void gsFlowVisitor<T>::evaluate(index_t testFunID)
{
    // shape basis (on the whole support of testFunID)

    index_t dim = m_params.getPde().dim();
    gsMatrix<T> support = m_testBasisPtr->support(testFunID);
    typename gsBasis<T>::domainIter domIt = m_params.getBases().front().piece(m_patchID).makeDomainIterator(boundary::none);

    gsMatrix<T> quNodes; // quad. nodes for the current element
    gsVector<T> quWeights; // weights for the current element
    std::vector< gsMatrix<T> > quNodesOnElem; // quad. nodes for all elements in support
    std::vector< gsVector<T> > quWeightsOnElem; // weights for all elements in support

    // loop over elements
    while(domIt->good())
    {
        bool inSupport = true; 

        // check if the current element lies in support of test function with testFunID
        for (index_t d = 0; d < dim; d++)
        {
            if ( (domIt->lowerCorner()[d] < support(d,0)) ||  (domIt->upperCorner()[d] > support(d,1)))
            {
                inSupport = false;
                break;
            }
        }

        // if so, compute and store the quadrature nodes and weights
        if (inSupport)
        {
            m_quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            quNodesOnElem.push_back(quNodes);
            quWeightsOnElem.push_back(quWeights);
        }

        domIt->next();
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
    m_params.getPde().patches().patch(m_patchID).computeMap(m_mapData);

    gsMatrix<index_t> allActives;
    m_shapeBasisPtr->active_into(m_quNodes, allActives);
    m_shapeFunActives = createVectorOfUniqueIndices(allActives);
    
    evalBasisData(m_shapeFunFlags, m_shapeBasisPtr, m_shapeFunActives, m_shapeFunData);

    // test basis
    m_testFunActives.resize(1,1);
    m_testFunActives << testFunID;
    evalBasisData(m_testFunFlags, m_testBasisPtr, m_testFunActives, m_testFunData);
}


template<class T>
void gsFlowVisitor<T>::evaluate(const gsDomainIterator<T>* domIt)
{
    m_quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), m_quNodes, m_quWeights);
    m_mapData.points = m_quNodes;
    m_params.getPde().patches().patch(m_patchID).computeMap(m_mapData);

    m_testBasisPtr->active_into(m_quNodes.col(0), m_testFunActives);
    m_shapeBasisPtr->active_into(m_quNodes.col(0), m_shapeFunActives);

    evalBasisData(m_testFunFlags, m_testBasisPtr, m_testFunActives, m_testFunData);
    evalBasisData(m_shapeFunFlags, m_shapeBasisPtr, m_shapeFunActives, m_shapeFunData);
}


} // namespace gismo