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
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>
#include <gsIncompressibleFlow/src/gsTMSolverBase.h>

#include <gsIncompressibleFlow/src/gsTMModels.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

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

    typename gsTMSolverBase<T, MatOrder>::tmPtr m_TMsolverPtr = NULL;
    typename gsTMModelData<T>::tdPtr m_TMModelPtr = NULL;
    bool m_bComputeTMfirst;
    real_t m_turbT;

protected: // *** Base class members ***

    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_assemblerPtr;
    using Base::m_paramsPtr;
    using Base::m_time;
    using Base::m_timeStepSize;
    using Base::m_innerIter;
    using Base::m_avgPicardIter;
    using Base::m_innerTol;
    using Base::m_clock;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsRANSSolverUnsteady(gsFlowSolverParams<T>& params):
    gsRANSSolverUnsteady(memory::make_shared_not_owned(&params))
    { }

    /// @brief Constructor.
    gsRANSSolverUnsteady(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr, false)
    { 
        // create turbulence model
        m_TMModelPtr = gsTMModelData<T>::make(paramsPtr);

        // create turbulence solver
        m_TMsolverPtr = gsTMSolverBase<T, MatOrder>::make(paramsPtr, m_TMModelPtr);
        
        // create assembler
        m_assemblerPtr = new gsRANSAssemblerUnsteady<T, MatOrder>(paramsPtr);
                        
        initMembers();

        m_paramsPtr->options().setSwitch("unsteady", true);
    }

    virtual ~gsRANSSolverUnsteady()
    {  }

protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers();

    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, std::ofstream& fileK, std::ofstream& fileO, std::ofstream& fileTV, std::string fileNamePrefix, unsigned plotPts);


public: // *** Member functions ***

    /// @brief Initialize the solver.
    virtual void initialize();

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    void solveWithAnimation(const int totalIter, const int iterStep, std::string fileNamePrefix = "", const T epsilon = 1e-3, unsigned plotPts = 10000, bool plotTurb = false, const int minIterations = 1);

    /// @brief Construct solution field for the unknown \a unk for the current solution vector.
    /// @param unk the considered unknown (0 - velocity, 1 - pressure, >=2 - turbulent variables)
    /// @return 
    gsField<T> constructSolution(int unk, bool customSwitch = false) const
    { 
        if (unk < 2)
            return Base::constructSolution(unk, customSwitch);
        else
            return m_TMsolverPtr->getAssembler()->constructSolution(unk, customSwitch); 
    }


public: // *** Getters/setters ***

    /// @brief Returns a pointer to the assembler.
    gsRANSAssemblerUnsteady<T, MatOrder>* getAssembler() const
    {
        return dynamic_cast<gsRANSAssemblerUnsteady<T, MatOrder>*>(m_assemblerPtr);
    }

    /// @brief Returns the total number of DOFs for turbulent model (the matrix size).
    int numDofsTM() { return (m_TMsolverPtr->getAssembler())->numDofs(); }

    /// @brief Retrurns the name of the class as a string.
    virtual std::string getName() { return "gsRANSSolverUnsteady"; }


}; // gsRANSSolverUnsteady

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRANSSolverUnsteady.hpp)
#endif