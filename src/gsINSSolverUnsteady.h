/** @file gsINSSolverUnsteady.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek, E. Turnerova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSSolverBase.h>
#include <gsIncompressibleFlow/src/gsINSAssemblerUnsteady.h>

#include <fstream>

namespace gismo
{

/// @brief The unsteady incompressible Navier-Stokes solver.
/// @tparam T coefficient type
template<class T>
class gsINSSolverUnsteady : public gsINSSolverBase<T>
{

public:
    typedef gsINSSolverBase<T> Base;

protected: // *** Class members ***

    T m_time, m_timeStepSize;
    T m_innerIter, m_avgPicardIter;
    T m_innerTol;

 protected: // *** Base class members ***

    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;
    using Base::m_clock;
    using Base::m_assembT;
    using Base::m_solsetupT;
    using Base::m_solveT;
    using Base::m_solver;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param params[in]       list of optional parameters
    /// @param createAssembler  create new assembler for this solver (true/false)
    gsINSSolverUnsteady(gsINSSolverParams<T>& params, bool createAssembler = true)
    {
        //create assembler
        m_pAssembler = new gsINSAssemblerUnsteady<T>(params);

        initMembers();

        m_timeStepSize = params.options().getReal("timeStep");
        m_innerIter = params.options().getInt("maxIt_picard");
        m_innerTol = params.options().getReal("tol_picard");
        m_avgPicardIter = 0;
    }

    virtual ~gsINSSolverUnsteady()
    {
    }


protected: // *** Member functions ***

    /// @brief Initialize all members.
    virtual void initMembers()
    {
        Base::initMembers();
        m_time = 0;
    }

    /// @brief Auxiliary function for plotting in Paraview.
    /// @param[out] file a reference to the output file stream
    void startAnimationFile(std::ofstream& file)
    {
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"Collection\" version=\"0.1\">";
        file << "<Collection>\n";
    }

    /// @brief Auxiliary function for plotting in Paraview.
    /// @param[out] file a reference to the output file stream
    void endAnimationFile(std::ofstream& file)
    {
        file << "</Collection>\n";
        file << "</VTKFile>\n";
        file.close();
    }

    /// @brief Plots the current velocity and pressure in Paraview.
    /// @param[out] fileU   a reference to the output file stream for velocity
    /// @param[out] fileP   a reference to the output file stream for pressure
    /// @param[in]  plotPts number of sample points for plotting
    void plotCurrentTimeStep(std::ofstream& fileU, std::ofstream& fileP, unsigned plotPts);


public: // *** Member functions ***

    /// @brief Update the assembler with current solution during Picard iteration (does not save the current solution in the assembler).
    /// @param[in] sol a const reference to the current solution vector
    virtual void updateAssemblerPicard(const gsMatrix<T>& sol)
    {
        m_clock.restart();
        getAssembler()->updatePicard(sol);
        m_assembT += m_clock.stop();
    }

    /// @brief Perform next iteration step.
    virtual void nextIteration();

    /// @brief Perform next time step (calls nextIteration()).
    void nextTimeStep() { nextIteration(); }

    /// @brief Perform several time steps.
    /// @param[in] numberOfIterations the number of time steps to be performed
    void nextTimeStep(const unsigned numberOfSteps)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        gsInfo << "Simulation ... \n";
        for (unsigned step = 0; step < numberOfSteps; step++)
        {
            nextIteration();
            gsInfo << "Step number " << m_iterationNumber << ", time " << m_time << " s, done. \n";
        }
    }

    /// @brief Solve the unsteady problem and create an animation in Paraview.
    /// @param[in] totalIter        number of time steps to be performed
    /// @param[in] iterStep         time steps interval for exporting solution snapshots in Paraview
    /// @param[in] epsilon          stopping tolerance
    /// @param[in] plotPts          number of sample points for plotting
    /// @param[in] minIterations    the minimum number of time steps
    virtual void solveWithAnimation(const int totalIter, const int iterStep, const T epsilon = 1e-3, unsigned plotPts = 10000, const int minIterations = 1);

protected: // *** Member functions ***

    /// @brief Initialize the generalized Stokes problem solution.
    /// @param[out] stokesMatrix    the generalized Stokes matrix
    void initGeneralizedStokesSolution(gsSparseMatrix<T>& stokesMatrix, gsMatrix<T>& stokesRhs);

public: // *** Getters/setters ***

    /// @brief Set a given solution vector as the initial solution and set the simulation time to zero.
    /// @param[in] solVector the given solution
    void setInitialCondition(const gsMatrix<T> & solVector)
    {
        this->setSolution(solVector);

        m_iterationNumber = 0;
        m_time = 0;
    }

    /// @brief Compute solution of the Stokes problem and set it as the initial solution.
    void setStokesInitialCondition()
    {
        gsInfo << "Setting Stokes solution as initial condition...\n";

        Base::solveStokes();

        m_iterationNumber = 0;
        m_time = 0;
    }

    /// @brief Solve the generalized Stokes problem.
    virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1);

    /// @brief Change the time step size.
    /// @param[in] timeStepSize the new time step size
    void changeTimeStepSize(const T timeStepSize)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->changeTimeStep(timeStepSize);
    }


    /// @brief Compute and display the relative norm of the solution change given the two successive solutions.
    /// @param[in] solOld the old solution
    /// @param[in] solNew the new solution
    void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     [u, p] Picard's solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << this->solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

        gsInfo << "\n";
    }

    /// @brief Returns the current iteration number.
    unsigned getTimeStepNumber() const { return m_iterationNumber; }

    /// @brief Returns the elapsed simulation time.
    T getSimulationTime() const { return m_time; }

    /// @brief Returns the average number of Picard iterations per time step.
    T getAvgPicardIterations() const { return m_avgPicardIter / m_iterationNumber; }

    /// @brief Returns a pointer to the assembler.
    virtual gsINSAssemblerUnsteady<T>* getAssembler() const
    {
        return dynamic_cast<gsINSAssemblerUnsteady<T>*>(m_pAssembler);
    }


}; // class gsINSSolverUnsteady


// ===========================================================================


/// @brief The unsteady incompressible Navier-Stokes solver with iterative solvers for the linear system.
/// @tparam T           coefficient type
/// @tparam SolverType  type of the linear solver
template<class T, class SolverType>
class gsINSSolverUnsteadyIter : public gsINSSolverUnsteady<T>, public gsINSSolverBaseIter<T, SolverType>
{

public:

    typedef gsINSSolverUnsteady<T> Base;
    typedef gsINSSolverBaseIter<T, SolverType> BaseIter;

protected: // *** Base class members ***

    using Base::m_avgPicardIter;
    using Base::m_time;
    using Base::m_timeStepSize;
    using BaseIter::m_matrices;
    using BaseIter::m_precType;
    using BaseIter::m_matGammaPart;
    using BaseIter::m_rhsGammaPart;
    using BaseIter::m_pPrec;
    using BaseIter::m_maxLinIter;
    using BaseIter::m_linTol;
    using gsINSSolverBase<T>::m_solution;
    using gsINSSolverBase<T>::m_iterationNumber;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param params[in]       list of optional parameters
    /// @param createAssembler  create new assembler for this solver (true/false)
    gsINSSolverUnsteadyIter(gsINSSolverParams<T>& params, bool createAssembler = true) :
        Base(params, createAssembler), BaseIter(params, *this->getAssembler())
    { }

    virtual ~gsINSSolverUnsteadyIter()
    {
    }


public: // *** Member functions ***

    /// @brief Compute the Stokes problem.
    virtual void solveStokes()
    {
        m_solution = BaseIter::getStokesSolution();
    }

    /// @brief Compute solution of the Stokes problem and set it as the initial solution.
    void setStokesInitialCondition()
    {
        gsInfo << "Setting Stokes solution as initial condition...\n";

        m_solution = BaseIter::getStokesSolution();

        m_iterationNumber = 0;
        m_time = 0;
    }

    /// @brief Solve the generalized Stokes problem.
    virtual void solveGeneralizedStokes(const int maxIterations, const T epsilon, const int minIterations = 1);

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

}; // class gsINSSolverUnsteadyIter

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSSolverUnsteady.hpp)
#endif
