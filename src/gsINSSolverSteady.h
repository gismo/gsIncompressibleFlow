/** @file gsINSSolverSteady.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolverBase.h>
#include <gsIncompressibleFlow/src/gsINSAssemblerSteady.h>

namespace gismo
{

/// @brief The steady incompressible Navier-Stokes solver.
/// @tparam T coefficient type
template<class T>
class gsINSSolverSteady : public gsINSSolverBase<T>
{

public:
    typedef gsINSSolverBase<T> Base;

protected: // *** Base class members ***

    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;


public: // *** Constructor/destructor ***

    gsINSSolverSteady(gsINSSolverParams<T>& params)
    {
        m_pAssembler = new gsINSAssemblerSteady<T>(params);

        Base::initMembers();
    }

    virtual ~gsINSSolverSteady()
    {
    }


public: // *** Member functions ***

    /// @brief Perform next iteration step.
    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        this->updateAssembler();

        // in the first iteration the pattern will get analyzed
        if (!m_iterationNumber)
            this->initIteration();

        this->applySolver(m_solution);

        m_iterationNumber++;
    }

    /// @brief Returns a pointer to the assembler.
    virtual gsINSAssemblerSteady<T>* getAssembler() const
    {
        return dynamic_cast<gsINSAssemblerSteady<T>*>(m_pAssembler);
    }

}; //gsINSSolverSteady


// ===========================================================================


/// @brief The steady incompressible Navier-Stokes solver with iterative solvers for the linear system.
/// @tparam T           coefficient type
/// @tparam SolverType  type of the linear solver
template<class T, class SolverType>
class gsINSSolverSteadyIter : public gsINSSolverSteady<T>, public gsINSSolverBaseIter<T, SolverType>
{

public:

    typedef gsINSSolverSteady<T> Base;
    typedef gsINSSolverBaseIter<T, SolverType> BaseIter;

protected: // *** Base class members ***

    using gsINSSolverBase<T>::m_solution;
    using gsINSSolverBase<T>::m_iterationNumber;

public: // *** Constructor/destructor ***

    /// @brief Constructor of the object.
    /// @param params[in] list of optional parameters 
    gsINSSolverSteadyIter(gsINSSolverParams<T>& params) :
        Base(params), BaseIter(params, *this->getAssembler())
    { }

    virtual ~gsINSSolverSteadyIter()
    {
    }

public: // *** Member functions ***

    /// @brief Compute the Stokes problem.
    virtual void solveStokes()
    {
        m_solution = BaseIter::getStokesSolution();
    }

    /// @brief Prepare for the solution process.
    virtual void initIteration()
    {
        BaseIter::initIteration();
    }

    /// @brief Solve the linear system.
    /// @param[out] solution a reference to the vector, where the computed solution will be stored
    virtual void applySolver(gsMatrix<T>& solution)
    {
        BaseIter::applySolver(solution);
    }

    /// @brief Initialize the solver.
    virtual void initialize()
    {
        Base::initialize();
        BaseIter::initialize(m_solution);
    }

    /// @brief Eliminate given DOFs as homogeneous Dirichlet boundary.
    /// @param[in] boundaryDofs     indices of the given boundary DOFs
    /// @param[in] unk              the considered unknown (0 - velocity, 1 - pressure)
    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        Base::markDofsAsEliminatedZeros(boundaryDofs, unk);
        BaseIter::reinitMembers();
    }

    /// @brief Returns the total time spent on matrix assembly.
    const T getAssemblyTime() const 
    { 
        return (Base::getAssemblyTime() + BaseIter::getAssemblyTime());
    }

    /// @brief Returns the total time spent on linear solver setup.
    const T getSolverSetupTime() const 
    { 
        return (Base::getSolverSetupTime() + BaseIter::getSolverSetupTime());
    }

    /// @brief Returns the total time spent on solving of the linear systems.
    const T getSolveTime() const 
    { 
        return (Base::getSolveTime() + BaseIter::getSolveTime());
    }

}; //gsINSSolverSteadyIter

} //namespace gismo
