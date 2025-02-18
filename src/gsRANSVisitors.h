/** @file gsINSVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
//#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsRANSTerms.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>


namespace gismo
{

/*
template <class T, int MatOrder>
class gsRANSVisitorUUSymmetricGradient : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    real_t m_viscosity;
    gsTMSolverSST<T, MatOrder>* m_TMsolver = NULL;
    gsVector<T> m_TurbulentViscosityVals;

protected:   // *** Base class members ***

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

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradient() {}

    gsRANSVisitorUUSymmetricGradient(typename gsFlowSolverParams<T>::Ptr paramsPtr, gsTMSolverSST<T, MatOrder>* TMsolver) :
    Base(paramsPtr)
    { 
        initMembers(TMsolver);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers(gsTMSolverSST<T, MatOrder>* TMsolver);

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }
};
*/

// ========================================================================================================

template <class T, int MatOrder>
class gsRANSVisitorUUSymmetricGradientDiag : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    real_t m_viscosity;
    gsTMSolverBase<T, MatOrder>* m_TMsolver = NULL;
    gsVector<T> m_TurbulentViscosityVals;
    gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

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
    //using Base::m_TMsolver;

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradientDiag() {}

    gsRANSVisitorUUSymmetricGradientDiag(typename gsFlowSolverParams<T>::Ptr paramsPtr/*, gsTMSolverSST<T, MatOrder>* TMsolver*/) :
    Base(paramsPtr)
    { 
        initMembers(/*TMsolver*/);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers(/*gsTMSolverSST<T, MatOrder>* TMsolver*/);

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsRANSTerm_SymmetricGradientDiag<T>(/*m_viscosity, m_TurbulentViscosityVals*/) );
        
        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions *** 

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setTurbulenceSolver(gsTMSolverBase<T, MatOrder>* TMsolver) { m_TMsolver = TMsolver;}

    // @brief Returns a pointer to the assembler.
    //gsFlowAssemblerBase<T, MatOrder>* getAssembler() const
    //{
    //    //return dynamic_cast<gsRANSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    //    return m_assemblerPtr;
    //}

    //void setRANSassembler(gsINSAssemblerUnsteady<T, MatOrder>* assemblerPtr) { m_assemblerPtr = assemblerPtr; }

};

// ==========================================================================================================

template <class T, int MatOrder>
class gsRANSVisitorUUSymmetricGradientOffdiag : public gsFlowVisitorVectorValued<T, MatOrder>
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;

public:
    real_t m_viscosity;
    gsTMSolverBase<T, MatOrder>* m_TMsolver = NULL;
    gsVector<T> m_TurbulentViscosityVals;
    gsFlowAssemblerBase<T, MatOrder>* m_assemblerPtr = NULL;

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

public: // *** Constructor/destructor ***

    gsRANSVisitorUUSymmetricGradientOffdiag() {}

    gsRANSVisitorUUSymmetricGradientOffdiag(typename gsFlowSolverParams<T>::Ptr paramsPtr/*, gsTMSolverSST<T, MatOrder>* TMsolver*/) :
    Base(paramsPtr)
    { 
        initMembers(/*TMsolver*/);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers(/*gsTMSolverSST<T, MatOrder>* TMsolver*/);

    void evaluate(index_t testFunID);

    void evaluate(const gsDomainIterator<T>* domIt);

    // upravit pro RANS
    virtual void defineTerms()
    {
        // evaluate turbulent viscosity

        m_terms.push_back( new gsRANSTerm_SymmetricGradientOffdiag<T>(/*m_viscosity, m_TurbulentViscosityVals*/) );

        //if(m_paramsPtr->options().getSwitch("unsteady"))
        //    m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***    

    /// @brief Initialize the visitor.
    void initialize();

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

public: // Getter/setters

    void setTurbulenceSolver(gsTMSolverBase<T, MatOrder>* TMsolver) { m_TMsolver = TMsolver;}

    // @brief Returns a pointer to the assembler.
    //gsFlowAssemblerBase<T, MatOrder>* getAssembler() const
    //{
    //    //return dynamic_cast<gsRANSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    //    return m_assemblerPtr;
    //}

    //void setRANSassembler(gsINSAssemblerUnsteady<T, MatOrder>* assemblerPtr) { m_assemblerPtr = assemblerPtr; }
    
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSVisitors.hpp)
#endif