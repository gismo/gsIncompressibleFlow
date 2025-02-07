/** @file gsRANSSolverUnsteady.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova, B. Bastl
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.h>

#include <gsCore/gsDebug.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

/// @brief              The unsteady incompressible Reynolds-Averaged Navier-Stokes solver.
/// @tparam T           coefficient type
/// @tparam MatOrder    sparse matrix storage order (ColMajor/RowMajor)
template<class T = real_t, int MatOrder = RowMajor>
class gsRANSSolverUnsteady : public gsINSSolverUnsteady<T, MatOrder>
{

public:
    typedef gsINSSolverUnsteady<T, MatOrder> Base;


public: // *** Class members ***

    //T m_time, m_timeStepSize;
    //T m_innerIter, m_avgPicardIter;
    //T m_innerTol;

    // nove definovane zde
    gsTMSolverSST<T>* m_TMsolver;

protected: // *** Base class members ***

    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_assemblerPtr;
    using Base::m_params;
    using Base::m_outFile;
    using Base::m_fileOutput;
    using Base::m_dispOutput;
    using Base::m_time;
    using Base::m_timeStepSize;
    using Base::m_innerIter;
    using Base::m_avgPicardIter;
    using Base::m_innerTol;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsRANSSolverUnsteady(gsFlowSolverParams<T>& params):
    gsRANSSolverUnsteady(memory::make_shared_not_owned(&params))
    { }

    /// @brief Constructor.
    gsRANSSolverUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    { 
        // create turbulence solver
        //gsTMSolverSST<T> m_TMsolver;
        
        // create assembler
        m_assemblerPtr = new gsRANSAssemblerUnsteady<T, MatOrder>(m_params);
        
        initMembers();

        m_params.options().setSwitch("unsteady", true);
    }

    virtual ~gsRANSSolverUnsteady()
    {  }

protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::ofstream& fileN, std::ofstream& fileTM, std::string fileNameSuffix, unsigned plotPts);


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    void solveWithAnimation(const int totalIter, const int iterStep, std::string fileNameSuffix = "", const T epsilon = 1e-3, unsigned plotPts = 10000, bool plotTurb = false, const int minIterations = 1);

    /// @brief Solve the generalized Stokes problem.
    //virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1)
    //{ GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    gsRANSAssemblerUnsteady<T, MatOrder>* getAssembler() const
    {
        return dynamic_cast<gsRANSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    }

    // @brief Returns the elapsed simulation time.
    // T getSimulationTime() const { return m_time; } ... lze pouzit definici z Base

    /// @brief Returns the average number of Picard iterations per time step.
    //T getAvgPicardIterations() const { return m_avgPicardIter / m_iterationNumber; } .. lze pouzit definici z Base

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsRANSSolverUnsteady"; }


}; // gsINSSolverUnsteady

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSSolverUnsteady.hpp)
#endif