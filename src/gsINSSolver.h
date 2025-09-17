/** @file gsINSSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
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
/// @ingroup IncompressibleFlow
template <class T, int MatOrder>
class gsINSSolver: public gsFlowSolverBase<T, MatOrder>
{

public:
    typedef gsFlowSolverBase<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolver(gsFlowSolverParams<T>& params):
    Base(params)
    {
        printOptions();
    }

    gsINSSolver(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    {
        printOptions();
    }

    virtual ~gsINSSolver()
    { }


public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    virtual void solveStokes();


protected: // *** Member functions ***

    /// @brief Print the options used by the solver to the logger.
    void printOptions()
    {
        std::stringstream sstr;
        sstr << "\n-----------------------------------\n";
        sstr << "Incompressible flow solver options:\n";
        sstr << m_paramsPtr->options() << "-----------------------------------\n";
        m_paramsPtr->logger().log(sstr.str(), true); // true = log to file only
    }


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
/// @ingroup IncompressibleFlow
template <class T = real_t, int MatOrder = RowMajor>
class gsINSSolverSteady : public gsINSSolver<T, MatOrder>
{

public:
    typedef gsINSSolver<T, MatOrder> Base;


protected: // *** Base class members ***

    using Base::m_paramsPtr;
    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;

    // m_initAssembT, m_assembT, m_solsetupT, m_solveT;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverSteady(const gsFlowSolverParams<T>& params):
    gsINSSolverSteady(memory::make_shared_not_owned(&params))
    { }

    /// @brief Constructor.
    gsINSSolverSteady(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { 
        m_assemblerPtr = new gsINSAssemblerSteady<T, MatOrder>(m_paramsPtr);

        Base::initMembers();
        m_paramsPtr->options().setSwitch("unsteady", false);
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
/// @ingroup IncompressibleFlow
template <class T = real_t, int MatOrder = RowMajor>
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
    using Base::m_paramsPtr;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSSolverUnsteady(gsFlowSolverParams<T>& params, bool createAssembler = true):
    gsINSSolverUnsteady(memory::make_shared_not_owned(&params), createAssembler)
    { }

    /// @brief Constructor.
    gsINSSolverUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr, bool createAssembler = true):
    Base(paramsPtr)
    { 
        if (createAssembler)
        {
            m_assemblerPtr = new gsINSAssemblerUnsteady<T, MatOrder>(m_paramsPtr);
            initMembers();
        }

        m_paramsPtr->options().setSwitch("unsteady", true);
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::string fileNamePrefix, unsigned plotPts);


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    void solveWithAnimation(const int totalIter, const int iterStep, std::string fileNamePrefix = "", const T epsilon = 1e-3, unsigned plotPts = 10000, const int minIterations = 1);

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