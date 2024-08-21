/** @file gsFlowVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{

/// @brief      Base class for incompressible flow visitors.
/// @tparam T   real number type
template <class T>
class gsFlowVisitor
{

protected: // *** Type definitions ***

    typedef gsFlowTermValVal<T>     MassTerm;
    typedef gsINSTermPvalUdiv<T>    PressureGradTerm;
    typedef gsINSTermUsolGradVal<T> ConvectionTerm;


protected: // *** Class members ***

    // constant and same for all visitors
    index_t m_patchID;
    gsFlowSolverParams<T> m_params; // pde, bases, assemblerOptions, options, precOptions
    // pde members: dim, viscosity, rhs (f,g), gsBoundaryConditions, patches, unknownDim
 
    // constant for individual visitors
    index_t m_testUnkID, m_shapeUnkID; // used, e.g., to reference the corresponding mapper
    unsigned m_geoFlags, m_testFunFlags, m_shapeFunFlags;
    const gsBasis<T>* m_testBasisPtr;
    const gsBasis<T>* m_shapeBasisPtr;
    std::vector< gsDofMapper > m_dofMappers;
    std::vector< gsFlowTerm<T>* > m_terms;

    // updated repeatedly
    index_t m_currentTestFunID; // updated in evaluate()
    gsMatrix<T> m_localMat;
    gsMapData<T> m_mapData; // members: points, dim, patchID

    // will be changed:
    gsQuadRule<T> m_quRule;
    gsMatrix<T> m_quNodes;
    gsVector<T> m_quWeights;
    gsMatrix<index_t> m_shapeFunActives;
    std::vector< gsMatrix<T> > m_testFunData;
    std::vector< gsMatrix<T> > m_shapeFunData; 
    

public: // *** Constructor/destructor ***

    gsFlowVisitor() {}

    gsFlowVisitor(const gsFlowSolverParams<T>& params) : m_params(params)
    { }

    ~gsFlowVisitor()
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

    void gatherEvalFlags();

    void defineTestShapeBases()
    { 
        m_testBasisPtr = &(m_params.getBases()[m_testUnkID].piece(m_patchID));
        m_shapeBasisPtr = &(m_params.getBases()[m_shapeUnkID].piece(m_patchID));
    }

    void setupQuadrature();
    //{ GISMO_NO_IMPLEMENTATION }


public: // *** Member functions ***

    void initialize();

    void updateDofMappers(const std::vector<gsDofMapper>& mappers) { m_dofMappers = mappers; }

    void initOnPatch(index_t patchID);

    void setCurrentSolution(std::vector<gsField<T> >& solutions);

    void setCurrentSolution(gsField<T>& solution);

    /// @brief 
    /// @param[in]  testFunID    the local test function index on the current patch
    void evaluate(index_t testFunID);
    //{ GISMO_NO_IMPLEMENTATION }

    virtual void assemble();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

    virtual void localToGlobal(gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

};

// ===================================================================================================================

template <class T>
class gsFlowVisitorVectorValued : public gsFlowVisitor<T> 
{

public:
    typedef gsFlowVisitor<T> Base;


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

    gsFlowVisitorVectorValued() {}

    gsFlowVisitorVectorValued(const gsFlowSolverParams<T>& params) : Base(params)
    { }


public: // *** Member functions ***

    virtual void assemble();

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowVisitors.hpp)
#endif