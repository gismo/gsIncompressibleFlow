/** @file gsINSSolverParams.h
 
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
class gsINSSolverParams
{

protected: // *** Class members ***

    gsNavStokesPde<T>               m_pde;
    std::vector<gsMultiBasis<T> >   m_bases;
    gsAssemblerOptions              m_assembOpt;
    gsOptionList                    m_opt;
    gsOptionList                    m_precOpt;

public: // *** Constructor/destructor ***

    /// @brief Constructor of the object.
    /// @param pde an incompressible Navier-Stokes problem
    /// @param bases vector of discretization bases (velocity, pressure)
    gsINSSolverParams(const gsNavStokesPde<T>& pde, const std::vector<gsMultiBasis<T> >& bases)
        : m_pde(pde), m_bases(bases)
    {
        m_assembOpt.dirStrategy = dirichlet::elimination;
        m_assembOpt.dirValues = dirichlet::interpolation;
        m_assembOpt.intStrategy = iFace::glue;

        m_opt = gsINSSolverParams<T>::defaultOptions();
        m_precOpt = gsINSPreconditioner<T>::defaultOptions();
    }

    /// @brief Constructor of the object with given options for the solver.
    /// @param pde an incompressible Navier-Stokes problem
    /// @param bases vector of discretization bases (velocity, pressure)
    /// @param opt list of options for the solver
    gsINSSolverParams(const gsNavStokesPde<T>& pde, const std::vector<gsMultiBasis<T> >& bases, gsOptionList opt)
        : m_opt(opt)
    { }

    ~gsINSSolverParams()
    {
    }

public: // *** Static functions ***

    /// @brief Returns a list of default options for the incompressible flow solver.
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;

        int nthreads =
        #ifdef _OPENMP
            std::round(omp_get_max_threads() / 2);
        #else
            1;
        #endif

        opt.addInt("numThreads", "Number of threads to be used if OPENMP enabled", nthreads);
        opt.addInt("maxIt_picard", "Maximum number of Picard iterations in one time step", 10);
        opt.addInt("maxIt_lin", "Maximum number of iterations for linear solver (if iterative)", 200);

        opt.addReal("timeStep", "Time step size", 0.1);
        opt.addReal("tol_picard", "Stopping tolerance for Picard iteration", 1-5);
        opt.addReal("tol_lin", "Stopping tolerance for linear solver (if iterative)", 1-6);

        opt.addString("precType", "Preconditioner to be used with iterative linear solver", "PCDmod_FdiagEqual");

        //opt.addSwitch("", "", true);

        return opt;
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
    std::vector<gsMultiBasis<T> > &       getBases() { return m_bases; }
    const std::vector<gsMultiBasis<T> > & getBases() const { return m_bases; }

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

}; // class gsINSSolverParams

} // namespace gismo