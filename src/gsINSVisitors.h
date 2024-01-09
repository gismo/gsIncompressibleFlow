/** @file gsINSVisitors.h
    
    @brief 
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova (Hornikova)
 */

#pragma once

#include <gsIncompressibleFlow/src/gsINSUtils.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSSolverParams.h>

namespace gismo
{

// ===================================================================================================================

// BASE CLASS
template <class T>
class gsINSVisitor
{

protected: // *** Class members ***

    // zvenku
    index_t m_patchID;
    gsINSSolverParams<T> m_params; // pde, bases, assemblerOptions, options, precOptions
    // pde members: dim, viscosity, rhs (f,g), gsBoundaryConditions, patches, unknownDim    // ulozit pointer/referenci?
 
    // definuje se tady, pak nemenne
    index_t m_testUnkID, m_shapeUnkID; // 0 - velocity, 1 - pressure, used, e.g., to reference the corresponding mapper
    gsQuadRule<T> m_quRule;
    std::vector< gsINSTerm<T>* > m_terms;
    unsigned m_geoFlags, m_testFunFlags, m_shapeFunFlags;
    const gsBasis<T>* m_testBasisPtr;
    const gsBasis<T>* m_shapeBasisPtr;
    std::vector< gsDofMapper > m_dofMappers;

    // aktualizuje se
    index_t m_currentTestFunID; // updated in evaluate()
    gsMatrix<T> m_localMat;
    gsMatrix<T> m_quNodes;
    gsVector<T> m_quWeights;
    gsMapData<T> m_mapData; // members: points, dim, patchID
    gsMatrix<index_t> m_shapeFunActives;
    std::vector< gsMatrix<T> > m_testFunData;
    std::vector< gsMatrix<T> > m_shapeFunData; 
    

public: // *** Constructor/destructor ***

    gsINSVisitor() {}

    gsINSVisitor(const gsINSSolverParams<T>& params) : m_params(params)
    { }

    ~gsINSVisitor()
    {
        deleteTerms();
    }

protected: // *** Member functions ***

    void deleteTerms()
    {
        for(size_t i = 0; i < m_terms.size(); i++)
            delete m_terms[i];

        m_terms.clear();
    }

    virtual void defineTerms()
    { GISMO_NO_IMPLEMENTATION }

    // is decided according to visitor type
    virtual void defineTestShapeUnknowns()
    { GISMO_NO_IMPLEMENTATION }

    void gatherEvalFlags()
    {
        m_geoFlags = 0;
        m_testFunFlags = 0;
        m_shapeFunFlags = 0;

        for (size_t i = 0; i < m_terms.size(); i++)
            m_terms[i]->updateEvalFlags(m_geoFlags, m_testFunFlags, m_shapeFunFlags);
    }    

    void defineTestShapeBases()
    { 
        m_testBasisPtr = &(m_params.getBases()[m_testUnkID].piece(m_patchID));
        m_shapeBasisPtr = &(m_params.getBases()[m_shapeUnkID].piece(m_patchID));
    }

    void setupQuadrature()
    {
        gsVector<index_t> numQuadNodes(m_params.getPde().dim()); 

        index_t maxDegTest = m_testBasisPtr->maxDegree();
        index_t maxDegShape = m_shapeBasisPtr->maxDegree();

        numQuadNodes.setConstant(math::min(maxDegTest, maxDegShape)+1);

        m_quRule = gsGaussRule<T>(numQuadNodes);
    }


public: // *** Member functions ***

    void initialize()
    {
        defineTestShapeUnknowns();  
        m_params.createDofMappers(m_dofMappers);  

        deleteTerms();
        defineTerms();
        gatherEvalFlags();
        m_mapData.flags = m_geoFlags;
    }

    void initOnPatch(index_t patchID)
    {
        m_patchID = patchID;
        m_mapData.patchId = m_patchID;
        defineTestShapeBases();
        setupQuadrature();          
    }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        for (size_t i = 0; i < m_terms.size(); i++)
        {
            gsINSTermNonlin<T>* termPtr = dynamic_cast< gsINSTermNonlin<T>* > (m_terms[i]);

            if (termPtr)
                termPtr->setCurrentSolution(solutions);
        }
    }

    void setCurrentSolution(gsField<T>& solution)
    { 
        for (size_t i = 0; i < m_terms.size(); i++)
        {
            gsINSTermNonlin<T>* termPtr = dynamic_cast< gsINSTermNonlin<T>* > (m_terms[i]);

            if (termPtr)
                termPtr->setCurrentSolution(solution);
        }
    }

    /// @brief 
    /// @param[in]  testFunID    the local test function index on the current patch
    void evaluate(index_t testFunID)
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
        m_shapeFunActives = getVectorOfUniqueIndices(allActives);
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


virtual void assemble()
{
    m_localMat.setZero(1, m_shapeFunActives.rows());

    for (size_t i = 0; i < m_terms.size(); i++)
        m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
}


virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
{ GISMO_NO_IMPLEMENTATION } 

virtual void localToGlobal(gsMatrix<T>& globalRhs)
{ GISMO_NO_IMPLEMENTATION } 

};

// ===================================================================================================================

template <class T>
class gsINSVisitorVectorValued : public gsINSVisitor<T>  // order: shape, test
{

public:
    typedef gsINSVisitor<T> Base;


protected: // *** Class members ***

    std::vector< gsMatrix<T> > m_locMatVec;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_mapData;
    using Base::m_shapeFunActives;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsINSVisitorVectorValued() {}

    gsINSVisitorVectorValued(const gsINSSolverParams<T>& params) : Base(params)
    { }


public: // *** Member functions ***

    virtual void assemble()
    {
        m_locMatVec.resize(m_params.getPde().dim());

        for (size_t i = 0; i < m_locMatVec.size(); i++)
            m_locMatVec[i].setZero(1, m_shapeFunActives.rows());

        for (size_t i = 0; i < m_terms.size(); i++)
            m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_locMatVec);
    }

};

// ===================================================================================================================
// ===================================================================================================================

// VELOCITY-VELOCITY VISITORS

template <class T>
class gsINSVisitorUU : public gsINSVisitor<T>
{

public:
    typedef gsINSVisitor<T> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_shapeFunActives;
    using Base::m_localMat;


public: // *** Constructor/destructor ***

    gsINSVisitorUU() {}

    gsINSVisitorUU(const gsINSSolverParams<T>& params) : Base(params)
    { }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
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

};

// ===================================================================================================================

template <class T>
class gsINSVisitorUUlin : public gsINSVisitorUU<T>
{

public:
    typedef gsINSVisitorUU<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUlin() {}

    gsINSVisitorUUlin(const gsINSSolverParams<T>& params) :
    Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermDiffusion<T>(m_params.getPde().viscosity()) );

        // if(m_params.options().getSwitch("unsteady"))
        //     m_terms.push_back( new gsINSTermTimeDiscr<T>(m_params.options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorUUnonlin : public gsINSVisitorUU<T>
{

public:
    typedef gsINSVisitorUU<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUnonlin() {}

    gsINSVisitorUUnonlin(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new typename gsINSTermNonlin<T>::ConvectionTerm() );

        // ... other terms, e.g. from stabilizations
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorUUtimeDiscr : public gsINSVisitorUUlin<T>
{

public:
    typedef gsINSVisitorUUlin<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUtimeDiscr() {}

    gsINSVisitorUUtimeDiscr(const gsINSSolverParams<T>& params) :
    Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermTimeDiscr<T>(m_params.options().getReal("timeStep")) );
    }

};


// ===================================================================================================================
// ===================================================================================================================

// VELOCITY-PRESSURE VISITORS
template <class T>
class gsINSVisitorPU : public gsINSVisitorVectorValued<T>  // order: shape, test
{

public:
    typedef gsINSVisitorVectorValued<T> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_patchID;
    using gsINSVisitor<T>::m_testUnkID;
    using gsINSVisitor<T>::m_shapeUnkID;
    using gsINSVisitor<T>::m_dofMappers;
    using gsINSVisitor<T>::m_currentTestFunID;
    using gsINSVisitor<T>::m_shapeFunActives;
    using gsINSVisitor<T>::m_terms;

public: // *** Constructor/destructor ***

    gsINSVisitorPU() {}

    gsINSVisitorPU(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermPvalUdiv<T>() );
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 1;   // pressure
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
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

                    globalRhs(dim*uCompSize + kk, 0) -= tmp;
                }
            }
        }
    } 

};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-PRESSURE VISITORS
template <class T>
class gsINSVisitorPP : public gsINSVisitor<T>
{

public:
    typedef gsINSVisitor<T> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_shapeFunActives;
    using Base::m_localMat;


public: // *** Constructor/destructor ***

    gsINSVisitorPP() {}

    gsINSVisitorPP(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 1;   // pressure
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
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

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPlin : public gsINSVisitorPP<T>
{

public:
    typedef gsINSVisitorPP<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlin() {}

    gsINSVisitorPPlin(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        // no default pressure-pressure terms in the Navier-Stokes eqns
        // optionally stabilization for inf-sup unstable discretizations
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPnonlin : public gsINSVisitorPP<T>
{

public:
    typedef gsINSVisitorPP<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPnonlin() {}

    gsINSVisitorPPnonlin(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        // no default pressure-pressure terms in the Navier-Stokes eqns
        // optionally stabilization for inf-sup unstable discretizations
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPmass : public gsINSVisitorPPlin<T>
{

public:
    typedef gsINSVisitorPPlin<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPmass() {}

    gsINSVisitorPPmass(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermValVal<T>() );
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPlaplace : public gsINSVisitorPPlin<T>
{

public:
    typedef gsINSVisitorPPlin<T> Base;


protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlaplace() {}

    gsINSVisitorPPlaplace(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermGradGrad<T>() );
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPconvection : public gsINSVisitorPPnonlin<T>
{

public:
    typedef gsINSVisitorPPnonlin<T> Base;

protected: // *** Base class members ***

    using gsINSVisitor<T>::m_params;
    using gsINSVisitor<T>::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPconvection() {}

    gsINSVisitorPPconvection(const gsINSSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermUsolGradVal<T>() );
    }

};

// ===================================================================================================================
// ===================================================================================================================

// RHS VISITORS

template <class T>
class gsINSVisitorRhsU : public gsINSVisitor<T>
{

public:
    typedef gsINSVisitor<T> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;


public: // *** Constructor/destructor ***

    gsINSVisitorRhsU() {}

    gsINSVisitorRhsU(const gsINSSolverParams<T>& params):
    Base(params), m_pRhsFun(params.getPde().force())
    {
        GISMO_ASSERT(m_pRhsFun->targetDim() == m_params.getPde().dim(), "Wrong RHS function passed into gsINSRhsU.");
    }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;    // velocity (not needed here)
    }

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermRhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble()
    {
        m_localMat.setZero(1, m_params.getPde().dim());

        for (size_t i = 0; i < m_terms.size(); i++)
            m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
    }


    virtual void localToGlobal(gsMatrix<T>& globalRhs)
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

};

// ===================================================================================================================

template <class T>
class gsINSVisitorRhsP : public gsINSVisitor<T>
{

public:
    typedef gsINSVisitor<T> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;


public: // *** Constructor/destructor ***

    gsINSVisitorRhsP() {}
    
    gsINSVisitorRhsP(const gsINSSolverParams<T>& params):
    Base(params), m_pRhsFun(params.getPde().source())
    {
        GISMO_ASSERT(m_pRhsFun == NULL || m_pRhsFun->targetDim() == 1, "Wrong RHS function passed into gsINSRhsP.");
    }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 1;    // pressure (not needed here)
    }

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermRhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble()
    {
        m_localMat.setZero(1, 1);

        for (size_t i = 0; i < m_terms.size(); i++)
            m_terms[i]->assemble(m_mapData, m_quWeights, m_testFunData, m_shapeFunData, m_localMat);
    }


    virtual void localToGlobal(gsMatrix<T>& globalRhs)
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

};



} // namespace gismo