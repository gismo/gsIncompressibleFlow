/** @file gsTMSolverSST.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for all flow solvers in gsIncompressibleFlow.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMSolverSST: public gsTMSolverBase<T, MatOrder>
{

public:
    typedef gsTMSolverBase<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_assemblerPtr;
    using Base::m_paramsPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;
    using Base::m_TurbulentViscosityVals;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMSolverSST(gsFlowSolverParams<T>& params):
    gsTMSolverSST(memory::make_shared_not_owned(&params))
    { }

    gsTMSolverSST(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }

    virtual ~gsTMSolverSST()
    { }


public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    //virtual void evalTurbulentViscosity(/*std::vector<gsMatrix<T> >& solUGrads, */gsMatrix<T>& quNodes/*, gsGeometryEvaluator<T> & geoEval*/);
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes);

    //virtual void plotTurbulentViscosity();

public: // *** Getters/setters ***

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsTMSolverSST"; }

    /// @brief Returns a pointer to the assembler.
    //virtual gsINSAssembler<T, MatOrder>* getAssembler() const
    //{ return dynamic_cast<gsTMAssembler<T, MatOrder>*>(m_assemblerPtr); }

};


} // namespace gismo