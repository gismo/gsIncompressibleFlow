/** @file gsINSVisitors.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSVisitors.h>

namespace gismo
{

template<class T>
void gsINSVisitor<T>::gatherEvalFlags()
{
    m_geoFlags = 0;
    m_testFunFlags = 0;
    m_shapeFunFlags = 0;

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->updateEvalFlags(m_geoFlags, m_testFunFlags, m_shapeFunFlags);
}    


template<class T>
void gsINSVisitor<T>::setupQuadrature()
{
    gsVector<index_t> numQuadNodes(m_params.getPde().dim()); 

    index_t maxDegTest = m_testBasisPtr->maxDegree();
    index_t maxDegShape = m_shapeBasisPtr->maxDegree();

    numQuadNodes.setConstant(math::min(maxDegTest, maxDegShape)+1);

    m_quRule = gsGaussRule<T>(numQuadNodes);
}


template<class T>
void gsINSVisitor<T>::initialize()
{
    defineTestShapeUnknowns();  
    m_params.createDofMappers(m_dofMappers);  

    deleteTerms();
    defineTerms();
    gatherEvalFlags();
    m_mapData.flags = m_geoFlags;
}


template<class T>
void gsINSVisitor<T>::initOnPatch(index_t patchID)
{
    m_patchID = patchID;
    m_mapData.patchId = m_patchID;
    defineTestShapeBases();
    setupQuadrature();          
}


template<class T>
void gsINSVisitor<T>::setCurrentSolution(std::vector<gsField<T> >& solutions)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsINSTermNonlin<T>* termPtr = dynamic_cast< gsINSTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solutions);
    }
}


template<class T>
void gsINSVisitor<T>::setCurrentSolution(gsField<T>& solution)
{ 
    for (size_t i = 0; i < m_terms.size(); i++)
    {
        gsINSTermNonlin<T>* termPtr = dynamic_cast< gsINSTermNonlin<T>* > (m_terms[i]);

        if (termPtr)
            termPtr->setCurrentSolution(solution);
    }
}


template<class T>
void gsINSVisitor<T>::evaluate(index_t testFunID)
{
    m_currentTestFunID = testFunID;

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
    index_t numAct = m_shapeFunActives.rows();

    // evaluate bases

    m_testFunData.clear();
    m_shapeFunData.clear();
    m_testFunData.resize(3); // 0 - value, 1 - deriv, 2 - deriv2
    m_shapeFunData.resize(3);

    if(m_testFunFlags & NEED_VALUE)
        m_testBasisPtr->evalSingle_into(testFunID, m_quNodes, m_testFunData[0]);

    if(m_testFunFlags & NEED_DERIV)
        m_testBasisPtr->derivSingle_into(testFunID, m_quNodes, m_testFunData[1]);

    if(m_testFunFlags & NEED_DERIV2)
        m_testBasisPtr->deriv2Single_into(testFunID, m_quNodes, m_testFunData[2]);


    if(m_shapeFunFlags & NEED_VALUE)
    {
        m_shapeFunData[0].setZero(numAct, numNodesInSupport);

        gsMatrix<T> tmpData;
        
        for(index_t i = 0; i < numAct; i++)
        {
            m_shapeBasisPtr->evalSingle_into(m_shapeFunActives(i), m_quNodes, tmpData);
            m_shapeFunData[0].row(i) = tmpData;
        }
    }

    if(m_shapeFunFlags & NEED_DERIV)
    {
        m_shapeFunData[1].setZero(dim*numAct, numNodesInSupport);

        gsMatrix<T> tmpData;

        for(index_t i = 0; i < numAct; i++)
        {
            m_shapeBasisPtr->derivSingle_into(m_shapeFunActives(i), m_quNodes, tmpData);
            m_shapeFunData[1].middleRows(dim*i, dim) = tmpData;
        }
    }
        

    if(m_shapeFunFlags & NEED_DERIV2)
    {
        index_t dimSq = dim*dim;
        m_shapeFunData[2].setZero(dimSq*numAct, numNodesInSupport);

        gsMatrix<T> tmpData;
        
        for(index_t i = 0; i < numAct; i++)
        {
            m_shapeBasisPtr->deriv2Single_into(m_shapeFunActives(i), m_quNodes, tmpData);
            m_shapeFunData[2].middleRows(dimSq*i, dimSq) = tmpData;
        }
    }

}


template<class T>
void gsINSVisitor<T>::assemble()
{
    m_localMat.setZero(1, m_shapeFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}

// ===================================================================================================================

template <class T>
void gsINSVisitorVectorValued<T>::assemble()
{
    m_locMatVec.resize(m_params.getPde().dim());

    for (size_t i = 0; i < m_locMatVec.size(); i++)
        m_locMatVec[i].setZero(1, m_shapeFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_locMatVec);
}

// ===================================================================================================================

template <class T>
void gsINSVisitorUU<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component
    index_t nComponents = globalMat.rows() / uCompSize;

    GISMO_ASSERT(nComponents == 1 || nComponents == dim, "Wrong matrix size in gsINSVisitorUU::localToGlobal.");

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t ii = testFunID(0);
    index_t numAct = m_shapeFunActives.rows();

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t j = 0; j < numAct; ++j)
        {
            const int jj = m_shapeFunActives(j);

            if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
            {
                for (index_t d = 0; d < nComponents; d++)
                    globalMat.coeffRef(ii + d*uCompSize, jj + d*uCompSize) += m_localMat(0, j);
            }
            else // is_boundary_index(jj)
            {
                const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                for (index_t d = 0; d < nComponents; d++)
                    globalRhs(ii + d*uCompSize, 0) -= m_localMat(0, j) * eliminatedDofs[m_shapeUnkID](bb, d);
            }
        }
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorPU<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal.");

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t ii = testFunID(0);
    index_t numAct = m_shapeFunActives.rows();

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t j = 0; j < numAct; ++j)
        {
            const int jj = m_shapeFunActives(j);

            if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
            {
                for (index_t d = 0; d < dim; d++)
                    globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](0, j);
            }
            else // is_boundary_index(jj)
            {
                const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                for (index_t d = 0; d < dim; d++)
                    globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](0, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
            }
        }
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorPU_withUPrhs<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[m_testUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.rows() == dim*uCompSize, "Wrong matrix size in gsINSVisitorPU::localToGlobal.");

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t ii = testFunID(0);
    index_t numAct = m_shapeFunActives.rows();

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t j = 0; j < numAct; ++j)
        {
            const int jj = m_shapeFunActives(j);

            if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
            {
                for (index_t d = 0; d < dim; d++)
                    globalMat.coeffRef(ii + d*uCompSize, jj) += m_locMatVec[d](0, j);
            }
            else // is_boundary_index(jj)
            {
                const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                for (index_t d = 0; d < dim; d++)
                    globalRhs(ii + d*uCompSize, 0) -= m_locMatVec[d](0, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
            }
        }
    }
    else // part arising from block B (assuming that the offdiag. blocks are symmetric)
    {
        const int bb = m_dofMappers[m_testUnkID].global_to_bindex(ii);
        for (index_t k = 0; k < numAct; k++)
        {
            const int kk = m_shapeFunActives(k);

            if (m_dofMappers[m_shapeUnkID].is_free_index(kk))
            {
                T tmp = 0;

                for (index_t d = 0; d < dim; d++)
                    tmp += m_locMatVec[d](0, k) * eliminatedDofs[m_testUnkID](bb, d);

                globalRhs(dim*uCompSize + kk, 0) += tmp;
            }
        }
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorUP<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[m_shapeUnkID].freeSize(); // number of dofs for one velocity component

    GISMO_ASSERT(globalMat.cols() == dim*uCompSize, "Wrong matrix size in gsINSVisitorUP::localToGlobal.");

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t ii = testFunID(0);
    index_t numAct = m_shapeFunActives.rows();

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t j = 0; j < numAct; ++j)
        {
            const int jj = m_shapeFunActives(j);

            if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
            {
                for (index_t d = 0; d < dim; d++)
                    globalMat.coeffRef(ii, jj + d*uCompSize) += m_locMatVec[d](0, j);
            }
            else // is_boundary_index(jj)
            {
                const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);
                
                T tmp = 0;

                for (index_t d = 0; d < dim; d++)
                    tmp -= m_locMatVec[d](0, j) * eliminatedDofs[m_shapeUnkID](bb, d);

                globalRhs(ii, 0) += tmp;
            }
        }
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorPP<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{
    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);
    m_dofMappers[m_shapeUnkID].localToGlobal(m_shapeFunActives, m_patchID, m_shapeFunActives);
    
    index_t ii = testFunID(0);
    index_t numAct = m_shapeFunActives.rows();

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t j = 0; j < numAct; ++j)
        {
            const int jj = m_shapeFunActives(j);

            if (m_dofMappers[m_shapeUnkID].is_free_index(jj))
            {
                globalMat.coeffRef(ii, jj) += m_localMat(0, j);
            }
            else // is_boundary_index(jj)
            {
                const int bb = m_dofMappers[m_shapeUnkID].global_to_bindex(jj);

                globalRhs(ii, 0) -= m_localMat(0, j) * eliminatedDofs[m_shapeUnkID](bb, 0);
            }
        }
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorRhsU<T>::assemble()
{
    m_localMat.setZero(1, m_params.getPde().dim());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}


template <class T>
void gsINSVisitorRhsU<T>::localToGlobal(gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[0].freeSize(); // number of dofs for one velocity component

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);

    index_t ii = testFunID(0);

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
    {
        for (index_t d = 0; d != dim; d++)
            globalRhs(ii + d*uCompSize, 0) += m_localMat(0, d);
    }
} 

// ===================================================================================================================

template <class T>
void gsINSVisitorRhsP<T>::assemble()
{
    m_localMat.setZero(1, 1);

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}


template <class T>
void gsINSVisitorRhsP<T>::localToGlobal(gsMatrix<T>& globalRhs)
{
    index_t dim = m_params.getPde().dim();
    const index_t uCompSize = m_dofMappers[0].freeSize(); // number of dofs for one velocity component

    gsMatrix<index_t> testFunID(1,1);
    testFunID << m_currentTestFunID;

    m_dofMappers[m_testUnkID].localToGlobal(testFunID, m_patchID, testFunID);

    index_t ii = testFunID(0);

    if (m_dofMappers[m_testUnkID].is_free_index(ii))
        globalRhs(dim*uCompSize + ii, 0) += m_localMat(0);
} 

} // namespace gismo