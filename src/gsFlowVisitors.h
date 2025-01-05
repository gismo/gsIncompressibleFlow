/** @file gsFlowVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <unordered_map>

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsFlowPeriodicHelper.h>

namespace gismo
{

/// @brief              Base class for incompressible flow visitors.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template <class T, int MatOrder>
class gsFlowVisitor
{

protected: // *** Type definitions ***

    typedef gsFlowTerm_ValVal<T>     MassTerm;
    typedef gsINSTerm_PvalUdiv<T>    PressureGradTerm;
    typedef gsINSTerm_UsolGradVal<T> ConvectionTerm;


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
    gsMatrix<T> m_localMat;
    gsMatrix<index_t> m_testFunActives, m_shapeFunActives;
    gsMapData<T> m_mapData; // members: points, dim, patchID

    // will be changed:
    gsQuadRule<T> m_quRule;
    gsMatrix<T> m_quNodes;
    gsVector<T> m_quWeights;
    std::vector< gsMatrix<T> > m_testFunData;
    std::vector< gsMatrix<T> > m_shapeFunData; 

    // for periodicity in radially symmetric domains
    bool m_hasPeriodicBC;
    gsMatrix<T> m_periodicTransformMat;
    typename gsFlowPeriodicHelper<T>::Ptr m_testPeriodicHelperPtr, m_shapePeriodicHelperPtr;
    

public: // *** Constructor/destructor ***

    gsFlowVisitor() {}

    gsFlowVisitor(const gsFlowSolverParams<T>& params) : m_params(params)
    {
        m_hasPeriodicBC = false;
    }

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

    /// @brief Define terms to be evaluated according to visitor type and given parameters.
    virtual void defineTerms()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Set test and shape basis unknown IDs according to visitor type.
    virtual void defineTestShapeUnknowns()
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Gather evaluation flags from all terms.
    void gatherEvalFlags();

    /// @brief Set pointers to test and shape basis.
    virtual void defineTestShapeBases()
    { 
        m_testBasisPtr = &(m_params.getBases()[m_testUnkID].piece(m_patchID));
        m_shapeBasisPtr = &(m_params.getBases()[m_shapeUnkID].piece(m_patchID));
    }

    /// @brief Setup the quadrature rule.
    virtual void setupQuadrature();
    //{ GISMO_NO_IMPLEMENTATION }

    /// @brief Evaluate required data for the given basis function.
    /// @param[in]  basisFlags      evaluation flags
    /// @param[in]  basisPtr        the given basis
    /// @param[in ] funID           index of the function to be evaluated
    /// @param[out] basisData       resulting data
    void evalSingleFunData(const unsigned& basisFlags, const gsBasis<T>* basisPtr, const index_t funID, std::vector< gsMatrix<T> >& basisData);

    /// @brief Evaluate required data for the given basis.
    /// @param[in]  basisFlags      evaluation flags
    /// @param[in]  basisPtr        the given basis
    /// @param[out] activesUnique   vector of indices of basis functions that are nonzero in any of quNodes
    /// @param[out] basisData       resulting data
    void evalBasisData(const unsigned& basisFlags, const gsBasis<T>* basisPtr, gsMatrix<index_t>& activesUnique, std::vector< gsMatrix<T> >& basisData);
    
    /// @brief Map local matrix to the global matrix (with no radial periodic conditions).
    /// @param[in]  eliminatedDofs  coefficients of the eliminated Dirichlet DoFs
    /// @param[out] globalMat       resulting global matrix
    /// @param[out] globalRhs       resulting global rhs
    virtual void localToGlobal_nonper(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

    /// @brief Map local matrix to the global matrix (with radial periodic conditions).
    /// @param[in]  eliminatedDofs  coefficients of the eliminated Dirichlet DoFs
    /// @param[out] globalMat       resulting global matrix
    /// @param[out] globalRhs       resulting global rhs
    virtual void localToGlobal_per(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

    /// @brief Map local rhs vector to the global rhs vector (with no radial periodic conditions).
    /// @param[out] globalRhs resulting global rhs
    virtual void localToGlobal_nonper(gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

    /// @brief Map local rhs vector to the global rhs vector (with radial periodic conditions).
    /// @param[out] globalRhs resulting global rhs
    virtual void localToGlobal_per(gsMatrix<T>& globalRhs)
    { GISMO_NO_IMPLEMENTATION } 

public: // *** Member functions ***

    /// @brief Initialize the visitor.
    void initialize();

    /// @brief Update DoF mappers.
    /// @param[in] mappers new mappers
    void updateDofMappers(const std::vector<gsDofMapper>& mappers) { m_dofMappers = mappers; }

    /// @brief Initialize the visitor on the given patch.
    /// @param[in] patchID the patch number
    void initOnPatch(index_t patchID);

    /// @brief Set all needed current solution fields.
    /// @param[in] solutions vector of new solution fields
    void setCurrentSolution(std::vector<gsField<T> >& solutions);

    /// @brief Set the current solution field.
    /// @param[in] solution new solution field
    void setCurrentSolution(gsField<T>& solution);

    /// @brief Evaluate basis data on the support of a given test function (used for row-by-row assembly).
    /// @param[in] testFunID the local test function index on the current patch
    void evaluate(index_t testFunID);

    /// @brief Evaluate basis data on the current element (used for element-by-element assembly).
    /// @param[in] domIt domain iterator pointing to the current element
    void evaluate(const gsDomainIterator<T>* domIt);

    /// @brief Assemble the local matrix.
    virtual void assemble();

    /// @brief Map local matrix to the global matrix.
    /// @param[in]  eliminatedDofs  coefficients of the eliminated Dirichlet DoFs
    /// @param[out] globalMat       resulting global matrix
    /// @param[out] globalRhs       resulting global rhs
    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs); 

    /// @brief Map local rhs vector to the global rhs vector.
    /// @param[out] globalRhs resulting global rhs
    virtual void localToGlobal(gsMatrix<T>& globalRhs);

    void setPeriodicHelpers(typename gsFlowPeriodicHelper<T>::Ptr testPerHelperPtr, typename gsFlowPeriodicHelper<T>::Ptr shapePerHelperPtr)
    { 
        m_testPeriodicHelperPtr = testPerHelperPtr;
        m_shapePeriodicHelperPtr = shapePerHelperPtr;
        m_hasPeriodicBC = true;
        m_periodicTransformMat = m_params.getPde().bc().getTransformMatrix();
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsFlowVisitorVectorValued : public gsFlowVisitor<T, MatOrder> 
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    std::vector< gsMatrix<T> > m_locMatVec;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_mapData;
    using Base::m_testFunActives;
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

// ===================================================================================================================

// TODO

// /// @brief              Base class for incompressible flow boundary visitors.
// /// @tparam T           real number type
// /// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
// template <class T, int MatOrder>
// class gsFlowVisitorBnd : public gsFlowVisitor<T, MatOrder> 
// {

// public:
//     typedef gsFlowVisitor<T, MatOrder> Base;

// protected: // *** Class members ***

//     boxSide m_side;

// protected: // *** Base class members ***

//     //using Base::;

// protected: // *** Class members ***

    

// public: // *** Constructor/destructor ***

//     gsFlowVisitorBnd() {}

//     gsFlowVisitorBnd(const gsFlowSolverParams<T>& params) : Base(params)
//     { }
    

// protected: // *** Member functions ***

//      /// @brief Set pointers to test and shape basis.
//     virtual void defineTestShapeBases()
//     { 
//         m_testBasisPtr = &(m_params.getBases()[m_testUnkID].piece(m_patchID).boundaryBasis(m_side));
//         m_shapeBasisPtr = &(m_params.getBases()[m_shapeUnkID].piece(m_patchID).boundaryBasis(m_side));
//     }

//     /// @brief Setup the quadrature rule.
//     void setupQuadrature();


// public: // *** Member functions ***

//     /// @brief Initialize the visitor on the given patch.
//     /// @param[in] patchID the patch number
//     void initOnPatchSide(index_t patchID, boxSide side)
//     {
//         m_side = side;
//         Base::initOnPatch(patchID);
//     }

//     /// @brief Evaluate basis data on the support of a given test function (used for row-by-row assembly).
//     /// @param[in] testFunID the local test function index on the current patch
//     void evaluate(index_t testFunID);

//     /// @brief Evaluate basis data on the current element (used for element-by-element assembly).
//     /// @param[in] domIt domain iterator pointing to the current element
//     void evaluate(const gsDomainIterator<T>* domIt);
// };

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowVisitors.hpp)
#endif