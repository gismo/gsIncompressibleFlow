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

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

namespace gismo
{

/// @brief              A class for a solver of the k-omega SST turbulence model.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsTMSolverSST: public gsTMSolverBase<T, MatOrder>
{

public:
    typedef gsTMSolverBase<T, MatOrder> Base;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsTMSolverSST> Ptr;

protected: // *** Class members ***
    
    bool m_isSSTModelSet;

protected: // *** Base class members ***

    using Base::m_TMModelPtr;
    using Base::m_assemblerPtr;
    using Base::m_paramsPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_TurbulentViscosityVals;
    using Base::m_TMavgPicardIter;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsTMSolverSST(gsFlowSolverParams<T>& params, typename gsTMModelData<T>::Ptr TMModelPtr):
    gsTMSolverSST(memory::make_shared_not_owned(&params), TMModelPtr)
    { }

    /// @brief Constructor.
    gsTMSolverSST(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr):
    Base(paramsPtr, TMModelPtr)
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

    /// @brief Returns a shared pointer to a newly created instance.
    /// @param[in] paramsPtr        a shared point to the instance of an object holding all parameters of the solver
    /// @param[in] TMModelPtr       a shared pointer to the chosen turbulence model
    static Ptr make(typename gsFlowSolverParams<T>::Ptr paramsPtr, typename gsTMModelData<T>::Ptr TMModelPtr)
    {
        return memory::make_shared_not_owned(new gsTMSolverSST<T, MatOrder>(paramsPtr, TMModelPtr));
    }

public: // *** Member functions ***

    /// @brief Evaluates the turbulent viscosity.
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId);

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

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTMSolverSST.hpp)
#endif