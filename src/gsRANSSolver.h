/** @file gsINSSolver.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
//#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

#include <gsCore/gsDebug.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief              A base class for incompressible Reynolds-Averaged Navier-Stokes solvers.
/// @tparam T           real number type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T, int MatOrder>
class gsRansSolver: public gsINSSolver<T, MatOrder>
{

public:
    typedef gsINSSolver<T, MatOrder> Base;


protected: // *** Base class members ***

    // ???

    using Base::m_assemblerPtr;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsRANSSolver(gsFlowSolverParams<T>& params):
    Base(params)
    { }

    virtual ~gsRANSSolver()
    { }


public: // *** Member functions ***

    /// @brief Compute the Stokes problem and save the solution into m_solution.
    //virtual void solveStokes();

    // funkce pro inicializaci RANS vypoctu a vypoctu turbulentniho modelu (linearizovane rovnic)???
    //
    //


public: // *** Getters/setters ***

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsRANSSolver"; }

    /// @brief Returns a pointer to the assembler.
    virtual gsRANSAssembler<T, MatOrder>* getAssembler() const
    { return dynamic_cast<gsRANSAssembler<T, MatOrder>*>(m_assemblerPtr); }

}; // gsINSSolver

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

    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::string fileNameSuffix, unsigned plotPts);


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    void solveWithAnimation(const int totalIter, const int iterStep, std::string fileNameSuffix = "", const T epsilon = 1e-3, unsigned plotPts = 10000, const int minIterations = 1);

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

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSSolver.hpp)
#endif