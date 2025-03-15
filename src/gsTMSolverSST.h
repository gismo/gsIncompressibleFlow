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
#include <gsIncompressibleFlow/src/gsTMAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsTMAssemblerSST.h>
//#include <gsIncompressibleFlow/src/gsTMModels.h>

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

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsTMSolverSST> tmPtr;

protected: // *** Class members ***

    //typename SSTModel<T>::Ptr m_SSTPtr;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr;     
    bool m_isSSTModelSet;

protected: // *** Base class members ***

    using Base::m_assemblerPtr;
    using Base::m_paramsPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;
    using Base::m_TurbulentViscosityVals;
    using Base::m_TMavgPicardIter;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMSolverSST(gsFlowSolverParams<T>& params, typename gsTMModelData<T>::tdPtr TMModelPtr):
    gsTMSolverSST(memory::make_shared_not_owned(&params), TMModelPtr)
    { }

    gsTMSolverSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr):
    Base(paramsPtr), m_TMModelPtr(TMModelPtr)
    { 
        m_assemblerPtr = new gsTMAssemblerSST<T, MatOrder>(paramsPtr, TMModelPtr);

        initMembers();
        m_paramsPtr->options().setSwitch("unsteady", true);

        m_isSSTModelSet = false;
    }

    virtual ~gsTMSolverSST()
    { }

protected: // *** Member functions ***

    /// @brief Initialize all members.
    void initMembers();

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static tmPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr)
    {
        return memory::make_shared_not_owned(new gsTMSolverSST<T, MatOrder>(paramsPtr, TMModelPtr));
    }

public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    //virtual void evalTurbulentViscosity(/*std::vector<gsMatrix<T> >& solUGrads, */gsMatrix<T>& quNodes/*, gsGeometryEvaluator<T> & geoEval*/);
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t patchId);

    //virtual void plotTurbulentViscosity();

public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    gsTMAssemblerSST<T, MatOrder>* getAssembler() const
    {
        return dynamic_cast<gsTMAssemblerSST<T, MatOrder>*>(m_assemblerPtr);
    }

    /// @brief Returns the average number of Picard iterations per time step.
    T getAvgPicardIterations() const { return m_TMavgPicardIter / m_iterationNumber; }

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsTMSolverSST"; }

    //SSTModel getModel() { return m_SSTModel; }

};

// ============================================================================================================================


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMSolverSST.hpp)
#endif