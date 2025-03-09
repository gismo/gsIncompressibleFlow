/** @file gsTMVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsTMTerms.h>

namespace gismo
{

template <class T, int MatOrder>
class gsTMVisitorLinearSST : public gsFlowVisitor<T, MatOrder>
{
    
public:
    typedef gsFlowVisitor<T, MatOrder> Base;
    
protected: // *** Class members ***
    
    index_t m_unknown;
    
protected: // *** Base class members ***
    
    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;
    using Base::m_geoFlags;

    
public: // *** Constructor/destructor ***
    
    gsTMVisitorLinearSST() {}
    
    gsTMVisitorLinearSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    Base(paramsPtr), m_unknown(unk)
    { }
    
    
protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = m_unknown;
        m_shapeUnkID = m_unknown;
    }
    
    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );
    }

public: // *** Member functions *** 

    /// @brief Initialize the assembler.
    //virtual void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);
    
};

// ==============================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorTimeIterationSST : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;

public:
    gsField<T> m_USolField;
    index_t m_unknown;

protected: // *** Base class members ***

    using Base::m_localMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;
    using Base::m_geoFlags;

public: // *** Constructor/destructor ***

    gsTMVisitorTimeIterationSST() {}

    gsTMVisitorTimeIterationSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    Base(paramsPtr), m_unknown(unk)
    { }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    void evaluate(index_t testFunID);
    
    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        m_terms.push_back( new gsTMTerm_VecCoeffGradVal<T>() );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = m_unknown;
        m_shapeUnkID = m_unknown;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    //void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    //void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    //void setRANSsolution(gsField<T> sol) { m_velSolField = sol;}

};

// ====================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorNonlinearSST : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    //gsField<T> m_velSolField;
    index_t m_unknown;
    real_t m_konst1, m_konst2, m_konst3;
    real_t m_konst;
    gsField<T> m_currentSol;
    index_t m_numLhsTerms;
    index_t m_numRhsTerms;
    gsMatrix<T> m_solution;
    gsField<T> m_distanceField;
    //real_t m_viscosity;
    //typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    //gsVector<T> m_TurbulentViscosityVals;
    //gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;
    using Base::m_geoFlags;
    //using Base::m_TurbulentViscosityVals;
    //using Base::m_viscosity;
    //using Base::m_TMsolverPtr;

public: // *** Constructor/destructor ***

    gsTMVisitorNonlinearSST() {}

    gsTMVisitorNonlinearSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, real_t k1, real_t k2, real_t k3, index_t unk) :
    Base(paramsPtr), m_unknown(unk), m_konst1(k1), m_konst2(k2), m_konst3(k3)
    {  }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    //void evaluate(index_t testFunID);

    //void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        m_numLhsTerms = 2;
        m_numRhsTerms = 2;
        
        // diffusion term with coefficient (m_konst1 * F1 + m_konst2 * (1 - F1)) * turbulent viscosity + m_konst3
        m_terms.push_back( new gsTMTerm_CoeffGradGrad<T>(m_paramsPtr, m_konst1, m_konst2, m_konst3) );
        // nonlinear reaction term
        m_terms.push_back( new gsTMTerm_CoeffValVal<T>(m_paramsPtr, m_unknown) );

        // blended term 2 * (1 - F1) * sigma0mega2 / omega * grad(k) * grad(omega) going to rhs of omega equation
        //if (m_unknown == 3)
        m_terms.push_back( new gsTMTerm_BlendCoeffRhs<T>(m_paramsPtr, m_unknown) );
        // production term going to rhs
        m_terms.push_back( new gsTMTerm_ProductionRhs<T>(m_paramsPtr, m_unknown) );
        
        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = m_unknown; 
        m_shapeUnkID = m_unknown;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    //void initialize();

    virtual void assemble();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setCurrentSolution(const gsMatrix<T>& solution) { m_solution = solution; }

    //void setDistanceField(const gsField<T>& dfield)
    //{ 
    //    m_distanceField = dfield;
    //}

    //void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    //void setRANSsolution(gsField<T> sol) { m_velSolField = sol;}

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMVisitors.hpp)
#endif