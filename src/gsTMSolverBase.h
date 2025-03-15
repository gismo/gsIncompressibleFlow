/** @file gsTMSolverBase.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <gsIncompressibleFlow/src/gsTMAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>

#include <gsIncompressibleFlow/src/gsTMModels.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A base class for all flow solvers in gsIncompressibleFlow.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMSolverBase: public gsFlowSolverBase<T, MatOrder>
{

public:
    typedef gsFlowSolverBase<T, MatOrder> Base;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsTMSolverBase> tmPtr;

public: // *** Class members ***

    typename gsTMModelData<T>::tdPtr m_TMModelPtr;    
    gsVector<T> m_TurbulentViscosityVals;
    T m_TMtime, m_TMtimeStepSize;
    T m_TMinnerIter, m_TMavgPicardIter;
    T m_TMinnerTol;   

protected: // *** Base class members ***

    using Base::m_assemblerPtr;
    using Base::m_paramsPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;
    using Base::m_outStream;


protected: // *** Base class function ***

    using Base::getAssembler;
    using Base::solutionChangeRelNorm;
    //using Base::initialize;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMSolverBase(gsFlowSolverParams<T>& params):
    Base(params)
    { }

    gsTMSolverBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }

    virtual ~gsTMSolverBase()
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance of the given preconditioner type.
    /// @param[in] precType the reqiured preconditioner type as a string
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner (assuming the following order: NS system matrix, mass matrix (velocity, pressure or both), other matrices)
    /// @param[in] opt a list of options for the preconditioner
    static tmPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr);


protected: // *** Member functions ***

    /// @brief Initialize members.
    void initMembers();


public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t patchId)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    //virtual void plotTurbulentViscosity();

public: // *** Getters/setters ***

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsTMSolverBase"; }

    gsVector<T> getTurbulentViscosity() 
    { 
        GISMO_ASSERT(m_TurbulentViscosityVals.rows() > 0, "Turbulent viscosity not evaluated yet.");    
        return m_TurbulentViscosityVals; 
    }

    bool isInitialized() { return getAssembler()->isInitialized(); }

    /// @brief Returns a pointer to the assembler.
    virtual gsTMAssemblerBase<T, MatOrder>* getAssembler() const
    { return dynamic_cast<gsTMAssemblerBase<T, MatOrder>*>(m_assemblerPtr); }

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMSolverBase.hpp)
#endif