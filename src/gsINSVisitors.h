/** @file gsINSVisitors.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{

// ===================================================================================================================
// VELOCITY-VELOCITY VISITORS

template <class T, int MatOrder>
class gsINSVisitorUU : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_localMat;

public: // *** Constructor/destructor ***

    gsINSVisitorUU() {}

    gsINSVisitorUU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorUUlin : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUlin() {}

    gsINSVisitorUUlin(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_Diffusion<T>(m_paramsPtr->getPde().viscosity()) );

        // if(m_paramsPtr->options().getSwitch("unsteady"))
        //     m_terms.push_back( new gsFlowTermTimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );

        // ... other terms, e.g. from stabilizations
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorUUnonlin : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUnonlin() {}

    gsINSVisitorUUnonlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new typename gsFlowVisitor<T, MatOrder>::ConvectionTerm() );

        // ... other terms, e.g. from stabilizations
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorUUmass : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUmass() {}

    gsINSVisitorUUmass(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_ValVal<T>() );
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorUUtimeDiscr : public gsINSVisitorUU<T, MatOrder>
{

public:
    typedef gsINSVisitorUU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUtimeDiscr() {}

    gsINSVisitorUUtimeDiscr(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_TimeDiscr<T>(m_paramsPtr->options().getReal("timeStep")) );
    }

};


// ===================================================================================================================
// ===================================================================================================================

// VELOCITY-PRESSURE VISITORS
template <class T, int MatOrder>
class gsINSVisitorPU : public gsFlowVisitorVectorValued<T, MatOrder>  // order: shape, test
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;


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

public: // *** Constructor/destructor ***

    gsINSVisitorPU() {}

    gsINSVisitorPU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_PvalUdiv<T>() );
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 1;   // pressure
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorPU_withUPrhs : public gsINSVisitorPU<T, MatOrder>  // order: shape, test
{

public:
    typedef gsINSVisitorPU<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;

public: // *** Constructor/destructor ***

    gsINSVisitorPU_withUPrhs() {}

    gsINSVisitorPU_withUPrhs(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-VELOCITY VISITORS
template <class T, int MatOrder>
class gsINSVisitorUP : public gsFlowVisitorVectorValued<T, MatOrder>  // order: shape, test
{

public:
    typedef gsFlowVisitorVectorValued<T, MatOrder> Base;


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

public: // *** Constructor/destructor ***

    gsINSVisitorUP() {}

    gsINSVisitorUP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_UdivPval<T>() );
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-PRESSURE VISITORS
template <class T, int MatOrder>
class gsINSVisitorPP : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_shapeFunActives;
    using Base::m_localMat;


public: // *** Constructor/destructor ***

    gsINSVisitorPP() {}

    gsINSVisitorPP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 1;   // pressure
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, MatOrder>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPPlin : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;


// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPPlin() {}

//     gsINSVisitorPPlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // no default pressure-pressure terms in the Navier-Stokes eqns
//         // optionally stabilization for inf-sup unstable discretizations
//     }

// };

// // ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPPnonlin : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;


// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPPnonlin() {}

//     gsINSVisitorPPnonlin(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // no default pressure-pressure terms in the Navier-Stokes eqns
//         // optionally stabilization for inf-sup unstable discretizations
//     }

// };

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorPPmass : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPmass() {}

    gsINSVisitorPPmass(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_ValVal<T>() );
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorPPlaplace : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlaplace() {}

    gsINSVisitorPPlaplace(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_GradGrad<T>() );
    }

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorPPconvection : public gsINSVisitorPP<T, MatOrder>
{

public:
    typedef gsINSVisitorPP<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPconvection() {}

    gsINSVisitorPPconvection(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTerm_UsolGradVal<T>() );
    }

};

// ===================================================================================================================

// template <class T, int MatOrder>
// class gsINSVisitorPP_PCDrobinBC : public gsINSVisitorPP<T, MatOrder>
// {

// public:
//     typedef gsINSVisitorPP<T, MatOrder> Base;

// protected: // *** Base class members ***

//     using Base::m_paramsPtr;
//     using Base::m_terms;


// public: // *** Constructor/destructor ***

//     gsINSVisitorPP_PCDrobinBC() {}

//     gsINSVisitorPP_PCDrobinBC(typename gsFlowSolverParams<T>::Ptr paramsPtr):
//     Base(paramsPtr)
//     { }


// protected: // *** Member functions ***

//     virtual void defineTerms()
//     {
//         // TODO
//     }

// };

// ===================================================================================================================
// ===================================================================================================================

// RHS VISITORS

template <class T, int MatOrder>
class gsINSVisitorRhsU : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;


public: // *** Constructor/destructor ***

    gsINSVisitorRhsU() {}

    gsINSVisitorRhsU(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr), m_pRhsFun(paramsPtr->getPde().force())
    {
        GISMO_ASSERT(m_pRhsFun->targetDim() == m_paramsPtr->getPde().dim(), "Wrong RHS function passed into gsINSRhsU.");
    }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;    // velocity (not needed here)
    }

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTerm_rhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble();

    virtual void localToGlobal(gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T, int MatOrder>
class gsINSVisitorRhsP : public gsFlowVisitor<T, MatOrder>
{

public:
    typedef gsFlowVisitor<T, MatOrder> Base;


protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFun;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_terms;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_testFunActives;
    using Base::m_localMat;
    using Base::m_mapData;
    using Base::m_quWeights;
    using Base::m_testFunData;
    using Base::m_shapeFunData;


public: // *** Constructor/destructor ***

    gsINSVisitorRhsP() {}
    
    gsINSVisitorRhsP(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr), m_pRhsFun(paramsPtr->getPde().source())
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
        m_terms.push_back( new gsFlowTerm_rhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble();

    virtual void localToGlobal(gsMatrix<T>& globalRhs);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSVisitors.hpp)
#endif