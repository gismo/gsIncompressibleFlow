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

template <class T>
class gsINSVisitorUU : public gsFlowVisitor<T>
{

public:
    typedef gsFlowVisitor<T> Base;


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

    gsINSVisitorUU(const gsFlowSolverParams<T>& params) : Base(params)
    { }
        

protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 0;    // velocity
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T>
class gsINSVisitorUUlin : public gsINSVisitorUU<T>
{

public:
    typedef gsINSVisitorUU<T> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUlin() {}

    gsINSVisitorUUlin(const gsFlowSolverParams<T>& params) :
    Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTermDiffusion<T>(m_params.getPde().viscosity()) );

        // if(m_params.options().getSwitch("unsteady"))
        //     m_terms.push_back( new gsFlowTermTimeDiscr<T>(m_params.options().getReal("timeStep")) );

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

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUnonlin() {}

    gsINSVisitorUUnonlin(const gsFlowSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new typename gsFlowVisitor<T>::ConvectionTerm() );

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

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorUUtimeDiscr() {}

    gsINSVisitorUUtimeDiscr(const gsFlowSolverParams<T>& params) :
    Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTermTimeDiscr<T>(m_params.options().getReal("timeStep")) );
    }

};


// ===================================================================================================================
// ===================================================================================================================

// VELOCITY-PRESSURE VISITORS
template <class T>
class gsINSVisitorPU : public gsFlowVisitorVectorValued<T>  // order: shape, test
{

public:
    typedef gsFlowVisitorVectorValued<T> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_params;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_shapeFunActives;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsINSVisitorPU() {}

    gsINSVisitorPU(const gsFlowSolverParams<T>& params) : Base(params)
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

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPU_withUPrhs : public gsINSVisitorPU<T>  // order: shape, test
{

public:
    typedef gsINSVisitorPU<T> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_params;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_shapeFunActives;

public: // *** Constructor/destructor ***

    gsINSVisitorPU_withUPrhs() {}

    gsINSVisitorPU_withUPrhs(const gsFlowSolverParams<T>& params) : Base(params)
    { }


public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-VELOCITY VISITORS
template <class T>
class gsINSVisitorUP : public gsFlowVisitorVectorValued<T>  // order: shape, test
{

public:
    typedef gsFlowVisitorVectorValued<T> Base;


protected: // *** Base class members ***

    using Base::m_locMatVec;
    using Base::m_params;
    using Base::m_patchID;
    using Base::m_testUnkID;
    using Base::m_shapeUnkID;
    using Base::m_dofMappers;
    using Base::m_currentTestFunID;
    using Base::m_shapeFunActives;
    using Base::m_terms;

public: // *** Constructor/destructor ***

    gsINSVisitorUP() {}

    gsINSVisitorUP(const gsFlowSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsINSTermUdivPval<T>() );
    }

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 0;   // velocity
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================
// ===================================================================================================================

// PRESSURE-PRESSURE VISITORS
template <class T>
class gsINSVisitorPP : public gsFlowVisitor<T>
{

public:
    typedef gsFlowVisitor<T> Base;


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

    gsINSVisitorPP(const gsFlowSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTestShapeUnknowns()
    {
        m_testUnkID = 1;    // pressure
        m_shapeUnkID = 1;   // pressure
    }

public: // *** Member functions ***

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, gsSparseMatrix<T, RowMajor>& globalMat, gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPlin : public gsINSVisitorPP<T>
{

public:
    typedef gsINSVisitorPP<T> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlin() {}

    gsINSVisitorPPlin(const gsFlowSolverParams<T>& params) : Base(params)
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

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPnonlin() {}

    gsINSVisitorPPnonlin(const gsFlowSolverParams<T>& params) : Base(params)
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

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPmass() {}

    gsINSVisitorPPmass(const gsFlowSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTermValVal<T>() );
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPlaplace : public gsINSVisitorPPlin<T>
{

public:
    typedef gsINSVisitorPPlin<T> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPlaplace() {}

    gsINSVisitorPPlaplace(const gsFlowSolverParams<T>& params) : Base(params)
    { }


protected: // *** Member functions ***

    virtual void defineTerms()
    {
        m_terms.push_back( new gsFlowTermGradGrad<T>() );
    }

};

// ===================================================================================================================

template <class T>
class gsINSVisitorPPconvection : public gsINSVisitorPPnonlin<T>
{

public:
    typedef gsINSVisitorPPnonlin<T> Base;

protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_terms;


public: // *** Constructor/destructor ***

    gsINSVisitorPPconvection() {}

    gsINSVisitorPPconvection(const gsFlowSolverParams<T>& params) : Base(params)
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
class gsINSVisitorRhsU : public gsFlowVisitor<T>
{

public:
    typedef gsFlowVisitor<T> Base;


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

    gsINSVisitorRhsU(const gsFlowSolverParams<T>& params):
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
        m_terms.push_back( new gsFlowTermRhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble();

    virtual void localToGlobal(gsMatrix<T>& globalRhs);

};

// ===================================================================================================================

template <class T>
class gsINSVisitorRhsP : public gsFlowVisitor<T>
{

public:
    typedef gsFlowVisitor<T> Base;


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
    
    gsINSVisitorRhsP(const gsFlowSolverParams<T>& params):
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
        m_terms.push_back( new gsFlowTermRhs<T>(m_pRhsFun) );
    }

public: // *** Member functions ***

    virtual void assemble();

    virtual void localToGlobal(gsMatrix<T>& globalRhs);

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSVisitors.hpp)
#endif