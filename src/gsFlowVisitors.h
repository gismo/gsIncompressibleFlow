/** @file gsFlowVisitors.h
    
    @brief Base visitors for incompressible flow problems.

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

namespace gismo
{

/// @brief Base class for incompressible flow visitors.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
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
    typename gsFlowSolverParams<T>::Ptr m_paramsPtr; // pde, bases, assemblerOptions, options, precOptions
                                                     // pde members: dim, viscosity, rhs (f,g), gsBoundaryConditions, patches, unknownDim
 
    // constant for individual visitors
    index_t m_testUnkID, m_trialUnkID; // used, e.g., to reference the corresponding mapper
    unsigned m_geoFlags, m_testFunFlags, m_trialFunFlags;
    const gsBasis<T>* m_testBasisPtr;
    const gsBasis<T>* m_trialBasisPtr;
    std::vector< gsFlowTerm<T>* > m_terms;

    // updated repeatedly
    gsMatrix<T> m_localMat;
    gsMatrix<index_t> m_testFunActives, m_trialFunActives;
    gsMapData<T> m_mapData; // members: points, dim, patchID

    // will be changed:
    gsQuadRule<T> m_quRule;
    gsMatrix<T> m_quNodes;
    gsVector<T> m_quWeights;
    std::vector< gsMatrix<T> > m_testFunData;
    std::vector< gsMatrix<T> > m_trialFunData; 

    // for periodicity in radially symmetric domains
    bool m_hasPeriodicBC;
    gsMatrix<T> m_periodicTransformMat;
    
    gsMatrix<T> m_bcvals;

public: // *** Constructor/destructor ***

    gsFlowVisitor() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsFlowVisitor(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
    {
        if (m_paramsPtr->hasPeriodicBC())
        {
            m_hasPeriodicBC = true;
            m_periodicTransformMat = m_paramsPtr->getPde().bc().getTransformMatrix();
        }
        else 
            m_hasPeriodicBC = false;
    }


    /// @brief Copy constructor.
    gsFlowVisitor(const gsFlowVisitor<T, MatOrder>& other):
    gsFlowVisitor()
    {
        // shallow copy of all members
        *this = other;

        // deep copy of m_terms
        m_terms.clear();
        m_terms.reserve(other.m_terms.size());
        for (auto* t : other.m_terms)
            m_terms.push_back(t->clone().release());
    }


    ~gsFlowVisitor()
    {
        deleteTerms();
    }
    

protected: // *** Member functions ***

    /// Free the vector of terms created in this class.
    void deleteTerms()
    {
        for(size_t i = 0; i < m_terms.size(); i++)
            delete m_terms[i];

        m_terms.clear();
    }

    /// Define terms to be evaluated according to visitor type and given parameters.
    virtual void defineTerms()
    { GISMO_NO_IMPLEMENTATION }

    /// Set test and trial basis unknown IDs according to visitor type.
    virtual void defineTestTrialUnknowns()
    { GISMO_NO_IMPLEMENTATION }

    /// Gather evaluation flags from all terms.
    void gatherEvalFlags();

    /// Set pointers to test and trial basis.
    virtual void defineTestTrialBases()
    { 
        m_testBasisPtr = &(m_paramsPtr->getBasis(m_testUnkID).piece(m_patchID));
        m_trialBasisPtr = &(m_paramsPtr->getBasis(m_trialUnkID).piece(m_patchID));
    }

    /// Setup the quadrature rule.
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

    /// Initialize the visitor.
    void initialize();
    
    /// @brief Initialize the visitor on the given patch.
    /// @param[in] patchID the patch number
    void initOnPatch(index_t patchID);

    /// @brief Initialize the visitor on the given patch.
    /// @param[in] patchID  the patch number
    /// @param[in] side     the box side on the given patch
    virtual void initOnPatchSide(index_t patchID, boxSide side)
    { GISMO_NO_IMPLEMENTATION } 

    /// @brief Set all needed current solution fields.
    /// @param[in] solutions vector of new solution fields
    void setCurrentSolution(std::vector<gsField<T> >& solutions);

    /// @brief Set the current solution field.
    /// @param[in] solution new solution field
    void setCurrentSolution(gsField<T>& solution);

    /// @brief Evaluate basis data on the support of a given test function (used for row-by-row assembly).
    /// @param[in] testFunID the local test function index on the current patch
    virtual void evaluate(index_t testFunID);

    /// @brief Evaluate basis data on the current element (used for element-by-element assembly).
    /// @param[in] domIt domain iterator pointing to the current element
    virtual void evaluate(const gsDomainIterator<T>* domIt);

    /// Assemble the local matrix.
    virtual void assemble();

    /// @brief Add local contributions to the global sparse system.
    /// @param[in]  eliminatedDofs  coefficients of the eliminated Dirichlet DoFs
    /// @param[out] globalMat       reference to the global matrix block
    /// @param[out] globalRhs       reference to the global right-hand side
    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

    /// @brief Add local contributions to the global right-hand side vector.
    /// @param[out] globalRhs reference to the global right-hand side
    virtual void localToGlobal(gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

/**
 * @brief Base class for incompressible flow visitors with vector-valued functions.
 * 
 * Visitors giving several different matrix blocks depending on component of the vector-valued function.
 * 
 * @tparam T        real number type
 * @tparam MatOrder sparse matrix storage order (ColMajor/RowMajor)
 * @ingroup IncompressibleFlow
 */
template <class T, int MatOrder>
class gsFlowVisitorVectorValued : public gsFlowVisitor<T, MatOrder> 
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    std::vector< gsMatrix<T> > m_locMatVec;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_mapData;
    using Base::m_testFunActives;
    using Base::m_trialFunActives;
    using Base::m_quNodes;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_trialFunData;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsFlowVisitorVectorValued() {}

    /// @brief Constructor.
    /// @param[in] paramsPtr a shared pointer to the container of input parameters
    gsFlowVisitorVectorValued(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }

public: // *** Member functions ***

    virtual void assemble();

};

// ===================================================================================================================

/// @brief              Base class for incompressible flow boundary visitors.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsFlowVisitorBnd : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_mapData;
    using Base::m_testBasisPtr;
    using Base::m_trialBasisPtr;
    using Base::m_testUnkID;
    using Base::m_trialUnkID;
    using Base::m_patchID;
    using Base::m_quRule;
    using Base::m_quNodes;
    using Base::m_quWeights;
    using Base::m_bcvals;
    using Base::m_terms;

protected: // *** Class members ***

    boxSide m_side;
    //typename gsBasis<T>::uPtr testBasisPtr;
    //typename gsBasis<T>::uPtr trialBasisPtr;

public: // *** Constructor/destructor ***

    gsFlowVisitorBnd() {}

    gsFlowVisitorBnd(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }

protected: // *** Member functions ***

    /// @brief Setup the quadrature rule.
    void setupQuadrature();

public: // *** Member functions ***

    /// @brief Initialize the visitor on the given patch.
    /// @param[in] patchID  the patch number
    /// @param[in] side     the box side on the given patch
    void initOnPatchSide(index_t patchID, boxSide side)
    {
        m_side = side;
        m_mapData.side = side;
        Base::initOnPatch(patchID);
    }

    /// @brief Evaluate basis data on the current element (used for element-by-element assembly).
    /// @param[in] domIt domain iterator pointing to the current element
    virtual void evaluate(const gsDomainIterator<T>* domIt);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowVisitors.hpp)
#endif