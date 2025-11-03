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

/// @brief              A base class for all turbulent models solvers in gsIncompressibleFlow.
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


protected: // *** Base class function ***

    using Base::getAssembler;
    using Base::solutionChangeRelNorm;
    
public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMSolverBase(gsFlowSolverParams<T>& params):
    Base(params)
    { }

    /// @brief Constructor.
    gsTMSolverBase(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { }

    virtual ~gsTMSolverBase()
    { }

public: // *** Static functions ***

    /// @brief Returns a shared pointer to a newly created instance of the turbulence solver.
    /// @param[in] paramsPtr        a shared point to the instance of an object holding all parameters of the solver
    /// @param[in] TMModelPtr       a shared pointer to the chosen turbulence model
    static tmPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::tdPtr TMModelPtr);


protected: // *** Member functions ***

    /// @brief Initialize members.
    void initMembers();


public: // *** Member functions ***

    /// @brief Evaluates the turbulent viscosity.
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Perform next iteration step.
    virtual void nextIteration();


public: // *** Getters/setters ***

    /// @brief Returns the name of the class as a string.
    virtual std::string getName() { return "gsTMSolverBase"; }

    /// @brief Returns the turbulent viscosity values.
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