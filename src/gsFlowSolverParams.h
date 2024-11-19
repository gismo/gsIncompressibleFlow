/** @file gsFlowSolverParams.h
 
    @brief A class that holds all parameters needed by the incompressible flow solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsNavStokesPde.h>
#include <gsIncompressibleFlow/src/gsINSPreconditioners.h>

namespace gismo
{

/** @brief
    A class that holds all parameters needed by the incompressible flow solver.

    - the INS PDE representation
    - discretization bases
    - list of parameters/options for the solver
    - list of assembler options
    - list of preconditioner options
 */

template<class T>
class gsFlowSolverParams
{

protected: // *** Class members ***

    gsNavStokesPde<T>               m_pde;
    std::vector<gsMultiBasis<T> >   m_bases;
    gsAssemblerOptions              m_assembOpt;
    gsOptionList                    m_opt;
    gsOptionList                    m_precOpt;

public: // *** Constructor/destructor ***

    gsFlowSolverParams() {}
    

    /// @brief Constructor of the object.
    /// @param pde an incompressible Navier-Stokes problem
    /// @param bases vector of discretization bases (velocity, pressure)
    gsFlowSolverParams(const gsNavStokesPde<T>& pde, const std::vector<gsMultiBasis<T> >& bases)
        : m_pde(pde), m_bases(bases)
    {
        m_assembOpt.dirStrategy = dirichlet::elimination;
        m_assembOpt.dirValues = dirichlet::interpolation;
        m_assembOpt.intStrategy = iFace::glue;

        m_opt = gsFlowSolverParams<T>::defaultOptions();
        m_precOpt = gsINSPreconditioner<T, RowMajor>::defaultOptions();
    }

    ~gsFlowSolverParams()
    {
    }

public: // *** Static functions ***

    /// @brief Returns a list of default options for the incompressible flow solver.
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;

        // nonlinear iteration
        opt.addInt("nonlin.maxIt", "Maximum number of Picard iterations in one time step", 10);
        opt.addReal("nonlin.tol", "Stopping tolerance for Picard iteration", 1e-5);
        
        // solving linear systems
        opt.addString("lin.solver", "The type of linear system solver (direct / iter / petsc)", "direct");
        opt.addString("lin.precType", "Preconditioner to be used with iterative linear solver", "PCDmod_FdiagEqual");
        opt.addInt("lin.maxIt", "Maximum number of iterations for linear solver (if iterative)", 200);
        opt.addReal("lin.tol", "Stopping tolerance for linear solver (if iterative)", 1e-6);

        // asssembly 
        //opt.addString("assemb.quad", "The numerical quadrature (Gauss/WQ)", "Gauss");
        opt.addString("assemb.loop", "EbE = element by element, RbR = row by row", "EbE");
        //opt.addSwitch("assemb.sumFact", "Use sum factorization for integration", false);
        opt.addSwitch("fillGlobalSyst", "Fill the global linear systems from blocks", true);
        
        // time-dependent problem
        opt.addSwitch("unsteady", "Assemble the velocity mass matrix", false);
        opt.addReal("timeStep", "Time step size", 0.1);

        // output
        opt.addSwitch("fileOutput", "Create an output file", false);
        opt.addSwitch("quiet", "Do not display output in terminal", false);
        opt.addString("outFile", "Name of the output file (or the full path to it)", "");

        // parallel 
        opt.addSwitch("parallel", "Currently running in parallel", false);

        // geometry jacobian evaluation
        opt.addInt("jac.npts", "Number of points along a patch side (in each direction) for geometry jacobian check", 100);
        opt.addReal("jac.dist", "Distance from boundary (in the parametric space) for geometry jacobian check", 1e-2);
        opt.addReal("jac.tol", "Critical value of geometry jacobian to throw warning", 1e-4);

        return opt;
    }


public: // *** Member functions ***

    /// @brief Creates DOF mappers for velocity and pressure.
    void createDofMappers(std::vector<gsDofMapper>& mappers)
    {
        mappers.resize(2);
    
        m_bases.front().getMapper(m_assembOpt.dirStrategy, m_assembOpt.intStrategy, m_pde.bc(), mappers.front(), 0);
        m_bases.back().getMapper(m_assembOpt.dirStrategy, m_assembOpt.intStrategy,  m_pde.bc(), mappers.back(), 1);
    }


public: // *** Getters/setters ***

    /// @brief Returns a const reference to the PDE.
    const gsNavStokesPde<T>& getPde() const { return m_pde; }

    /// @brief Returns a const reference to the boundary conditions.
    const gsBoundaryConditions<T>& getBCs() const { return m_pde.bc(); }

    /**
     * @brief Returns a reference to the discretization bases.
     *
     * There is also a const version returning a const reference.
     */
    std::vector<gsMultiBasis<T> >&       getBases() { return m_bases; }
    const std::vector<gsMultiBasis<T> >& getBases() const { return m_bases; }

    /**
     * @brief Returns a reference to the assembler option list.
     *
     * There is also a const version returning a const reference.
     */
    gsAssemblerOptions& assemblerOptions() { return m_assembOpt; }
    const gsAssemblerOptions& assemblerOptions() const { return m_assembOpt; }

    /// Set assembler options given in \a opt.
    void setAssemblerOptions(const gsAssemblerOptions& opt) { m_assembOpt = opt; }

    /**
     * @brief Returns a reference to the preconditioner option list.
     *
     * There is also a const version returning a const reference.
     */
    gsOptionList& precOptions() { return m_precOpt; }
    const gsOptionList& precOptions() const { return m_precOpt; }

    /// @brief Set preconditioner options given in \a opt.
    void setPrecOptions(const gsOptionList& opt) { m_precOpt = opt; }

    /**
     * @brief Returns a reference to the INS solver option list.
     *
     * There is also a const version returning a const reference.
     */
    gsOptionList& options() { return m_opt; }
    const gsOptionList& options() const { return m_opt; }

    /// @brief Set INS solver options given in \a opt.
    void setOptions(const gsOptionList& opt) { m_opt = opt; }


}; // class gsFlowSolverParams

} // namespace gismo