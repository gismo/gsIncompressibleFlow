/** @file gsINSVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

//#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
//#include <gsIncompressibleFlow/src/gsINSAssembler.h>
//#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
//#include <gsIncompressibleFlow/src/gsINSVisitors.h>
//#include <gsIncompressibleFlow/src/gsINSTerms.h>
//#include <gsIncompressibleFlow/src/gsRANSTerms.h>
//#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
//#include <gsIncompressibleFlow/src/gsTMSolverSST.h>
#include <gsIncompressibleFlow/src/gsTMTerms.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>


namespace gismo
{

template <class T, int MatOrder>
class gsTMVisitorVelocityAdvection : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;

public:
    //gsField<T> m_velSolField;
    index_t m_unknown;
    //real_t m_viscosity;
    //typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    //gsVector<T> m_TurbulentViscosityVals;
    //gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

protected: // *** Base class members ***

    using Base::m_locMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    //using Base::m_TurbulentViscosityVals;
    //using Base::m_viscosity;
    //using Base::m_TMsolverPtr;

public: // *** Constructor/destructor ***

    gsTMVisitorVelocityAdvection() {}

    gsTMVisitorVelocityAdvection(typename gsFlowSolverParams<T>::Ptr paramsPtr, index_t unk) :
    Base(paramsPtr), m_unknown(unk)
    { 
//        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    //void evaluate(index_t testFunID);

    //void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsTMTerm_VecCoeffGradVal<T>() );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns(index_t unk)
    {
        m_testUnkID = unk;
        m_shapeUnkID = unk;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    //void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    //void setRANSsolution(gsField<T> sol) { m_velSolField = sol;}

};

// ====================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorDiffusion : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;

public:
    //gsField<T> m_velSolField;
    index_t m_unknown;
    real_t m_konst1, m_konst2;
    //real_t m_viscosity;
    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    gsVector<T> m_TurbulentViscosityVals;
    //gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

protected: // *** Base class members ***

    using Base::m_locMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    //using Base::m_TurbulentViscosityVals;
    //using Base::m_viscosity;
    //using Base::m_TMsolverPtr;

public: // *** Constructor/destructor ***

    gsTMVisitorDiffusion() {}

    gsTMVisitorDiffusion(typename gsFlowSolverParams<T>::Ptr paramsPtr, real_t k1, real_t k2, index_t unk) :
    Base(paramsPtr), m_unknown(unk), m_konst1(k1), m_konst2(k2)
    { 
//        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsTMTerm_CoeffGradGrad<T>(m_konst1, m_konst2) );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns(index_t unk)
    {
        m_testUnkID = unk; 
        m_shapeUnkID = unk;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    //void setRANSsolution(gsField<T> sol) { m_velSolField = sol;}

};

// ===========================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorReaction : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;

public:
    //gsField<T> m_velSolField;
    index_t m_unknown;
    real_t m_konst;
    gsField<T> m_currentSol;
    //real_t m_viscosity;
    //typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    //gsVector<T> m_TurbulentViscosityVals;
    //gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

protected: // *** Base class members ***

    using Base::m_locMat;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_terms;
    using Base::m_quNodes;
    //using Base::m_TurbulentViscosityVals;
    //using Base::m_viscosity;
    //using Base::m_TMsolverPtr;

public: // *** Constructor/destructor ***

    gsTMVisitorReaction() {}

    gsTMVisitorReaction(typename gsFlowSolverParams<T>::Ptr paramsPtr, real_t konst, index_t unk) :
    Base(paramsPtr), m_unknown(unk), m_konst(konst)
    { 
//        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsTMTerm_CoeffValVal<T>(m_konst) );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns(index_t unk)
    {
        m_testUnkID = unk; 
        m_shapeUnkID = unk;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setCurrentSolution(gsField<T>& solution)
    { 
        m_currentSol = solution;
    }  

    //void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    //void setRANSsolution(gsField<T> sol) { m_velSolField = sol;}

};

// ====================================================================================================================

template <class T, int MatOrder>
class gsTMVisitorProductuionRhsSST : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    //gsField<T> m_velSolField;
    //real_t m_viscosity;
    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    gsVector<T> m_TurbulentViscosityVals;
    //gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;
    gsMatrix<T> m_velSolCoeffs;
    gsMatrix<T> m_KSolCoeffs;
    gsMatrix<T> m_OmegaSolCoeffs;
    gsMatrix<T> m_RANSsolution;

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
    //using Base::m_TurbulentViscosityVals;
    //using Base::m_viscosity;
    //using Base::m_TMsolverPtr;

public: // *** Constructor/destructor ***

    gsTMVisitorProductuionRhsSST() {}

    gsTMVisitorProductuionRhsSST(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { 
//        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    //void initMembers();

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsTMTerm_CoeffGradGrad<T>(m_konst1, m_konst2) );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0; 
        m_shapeUnkID = 0;
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setTurbulenceSolver(typename gsTMSolverBase<T, MatOrder>::tmPtr TMsolver) { m_TMsolverPtr = TMsolver;}

    void setRANSsolution(gsMatrix<T> sol) { m_RANSsolution = sol;}

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMVisitors.hpp)
#endif