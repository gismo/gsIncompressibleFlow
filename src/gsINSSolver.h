/** @file gsINSSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

#include <gsCore/gsDebug.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief              A base class for incompressible Navier-Stokes solvers.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsINSSolver: public gsFlowSolverBase<T, MatOrder>
{

public:
    typedef gsFlowSolverBase<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolver(gsFlowSolverParams<T>& params):
    Base(params)
    { }

    virtual ~gsINSSolver()
    { }


public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    virtual void solveStokes();


public: // *** Getters/setters ***

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsINSSolver"; }

    /// @brief Returns a pointer to the assembler.
    virtual gsINSAssembler<T, MatOrder>* getAssembler() const
    { return dynamic_cast<gsINSAssembler<T, MatOrder>*>(m_assemblerPtr); }

}; // gsINSSolver

// ===================================================================================================================

/// @brief              The steady incompressible Navier-Stokes solver.
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T = real_t, int MatOrder = RowMajor>
class gsINSSolverSteady : public gsINSSolver<T, MatOrder>
{

public:
    typedef gsINSSolver<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_params;
    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;

    // m_initAssembT, m_assembT, m_solsetupT, m_solveT;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverSteady(gsFlowSolverParams<T>& params):
    Base(params)
    { 
        m_assemblerPtr = new gsINSAssemblerSteady<T, MatOrder>(m_params);

        Base::initMembers();
        m_params.options().setSwitch("unsteady", false);
    }


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration();


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    virtual gsINSAssemblerSteady<T, MatOrder>* getAssembler() const
    {
        return dynamic_cast<gsINSAssemblerSteady<T, MatOrder>*>(m_assemblerPtr);
    }

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsINSSolverSteady"; }


}; // gsINSSolverSteady

// ===================================================================================================================

/// @brief              The unsteady incompressible Navier-Stokes solver.
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T = real_t, int MatOrder = RowMajor>
class gsINSSolverUnsteady : public gsINSSolver<T, MatOrder>
{

public:
    typedef gsINSSolver<T, MatOrder> Base;


protected: // *** Class members ***

    T m_time, m_timeStepSize;
    T m_innerIter, m_avgPicardIter;
    T m_innerTol;

protected: // *** Base class members ***

    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_assemblerPtr;
    using Base::m_params;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverUnsteady(gsFlowSolverParams<T>& params): Base(params)
    { 
        m_assemblerPtr = new gsINSAssemblerUnsteady<T, MatOrder>(m_params);

        initMembers();
        m_params.options().setSwitch("unsteady", true);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    /// @brief Solve the generalized Stokes problem.
    virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    { GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    gsINSAssemblerUnsteady<T, MatOrder>* getAssembler() const
    {
        return dynamic_cast<gsINSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    }

    // @brief Returns the elapsed simulation time.
    T getSimulationTime() const { return m_time; }

    /// @brief Returns the average number of Picard iterations per time step.
    T getAvgPicardIterations() const { return m_avgPicardIter / m_iterationNumber; }

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsINSSolverUnsteady"; }


}; // gsINSSolverUnsteady

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSSolver.hpp)
#endif