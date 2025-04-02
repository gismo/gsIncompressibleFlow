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

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowSolverParams> Ptr;
    typedef memory::unique_ptr<gsFlowSolverParams> uPtr;

protected: // *** Class members ***

    typename gsNavStokesPde<T>::Ptr m_pdePtr;
    std::vector<gsMultiBasis<T> >   m_bases;
    std::vector<gsMultiBasis<T> >   m_basesTM;
    gsBoundaryConditions<T>         m_BC;
    gsAssemblerOptions              m_assembOpt;
    gsOptionList                    m_opt;
    gsOptionList                    m_precOpt;

    bool m_isBndSet;
    std::vector<std::pair<int, boxSide> > m_bndIn, m_bndOut, m_bndWall;

    gsField<T> m_USolField;
    gsField<T> m_KSolField;
    gsField<T> m_OSolField;
    gsField<T> m_distanceField;
    
public: // *** Constructor/destructor ***

    gsFlowSolverParams() {}
    
    /// @brief Constructor of the object.
    /// @param pde an incompressible Navier-Stokes problem
    /// @param bases vector of discretization bases (velocity, pressure)
    gsFlowSolverParams(const gsNavStokesPde<T>& pde, const std::vector<gsMultiBasis<T> >& bases, const std::vector<gsMultiBasis<T> >& basesTM = {})
        : m_pdePtr(memory::make_shared_not_owned(&pde)), m_bases(bases), m_basesTM(basesTM)
    {
        m_assembOpt.dirStrategy = dirichlet::elimination;
        m_assembOpt.dirValues = dirichlet::interpolation;
        m_assembOpt.intStrategy = iFace::glue;

        m_opt = gsFlowSolverParams<T>::defaultOptions();
        m_precOpt = gsINSPreconditioner<T, RowMajor>::defaultOptions();

        m_BC = m_pdePtr->bc();
        m_isBndSet = false;

        if (!basesTM.empty())
        {
            for (size_t i = 0; i < basesTM.size(); i++)
                m_bases.push_back(basesTM[i]);
        }
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
        opt.addString("lin.krylov", "The Krylov subspace method from G+Smo/Eigen (for lin.solver = iter)", "gmres");
        opt.addString("lin.precType", "Preconditioner to be used with iterative linear solver", "MSIMPLER_FdiagEqual");
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

        // Turbulent model
        opt.addString("TM", "Chosen tubulence model identifier. Current choices: SST", "SST");
        opt.addInt("TM.maxIt", "Maximum number of Picard iterations in turbulent model in one time step", 10);
        opt.addInt("TM.maxItFirst", "Maximum number of Picard iterations in turbulent model in the first time step", 10);
        opt.addReal("TM.tol", "Stopping tolerance for Picard iteration in turbulent model", 1e-5);
        opt.addSwitch("TM.limitproduction","Using limiter for production term in turbulence model", false);
        opt.addReal("TM.uFreeStream", "Magnitude of a free-stream velocity", 1.0);
        opt.addReal("TM.turbIntensity", "Turbulent intensity", 0.05);
        opt.addReal("TM.viscosityRatio", "Specifies approximate ratio of turbulent viscosity to kinematic viscosity", 50.0);
        opt.addInt("TM.addRefsDF", "Number of additional uniform refinements of pressure basis for distance field computation", 2);

        return opt;
    }


public: // *** Member functions ***

    /// @brief Creates DOF mappers.
    void createDofMappers(std::vector<gsDofMapper>& mappers)
    {
        mappers.resize(m_bases.size());
   
        for (size_t i = 0; i < m_bases.size(); i++)
        {
            m_bases[i].getMapper(m_assembOpt.dirStrategy, m_assembOpt.intStrategy, m_BC, mappers[i], i);
        }
    }

    /// @brief Set boundary parts (vectors of pairs [patch, side]).
    /// @param[in] bndIn    inflow boundary part 
    /// @param[in] bndOut   outlfow boundary part
    /// @param[in] bndWall  solid wall boundary part
    void setBndParts(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall)
    {
        m_bndIn = bndIn;
        m_bndOut = bndOut;
        m_bndWall = bndWall;

        m_isBndSet = true;
    }


public: // *** Getters/setters ***

    /// @brief Returns a const reference to the PDE.
    const gsNavStokesPde<T>& getPde() const { return *m_pdePtr; }

    /// @brief Returns a const reference to the boundary conditions.
    const gsBoundaryConditions<T>& getBCs() const { return m_BC; }

    /// @brief Returns a const reference to the boundary conditions.
    void setBCs(gsBoundaryConditions<T>& bcs) { m_BC = bcs; }

    /**
     * @brief Returns a reference to the discretization bases for velocity and pressure
     *
     * There is also a const version returning a const reference.
     */
    std::vector<gsMultiBasis<T> >&       getBases() { return m_bases; }
    const std::vector<gsMultiBasis<T> >& getBases() const { return m_bases; }

    /**
     * @brief Returns a reference to the discretization bases for turbulence model
     *
     * There is also a const version returning a const reference.
     */
    std::vector<gsMultiBasis<T> >&       getBasesTM() { return m_basesTM; }
    const std::vector<gsMultiBasis<T> >& getBasesTM() const { return m_basesTM; }

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

    /// @brief Get vector of [patch, side] corresponding to the inflow boundary.
    std::vector<std::pair<int, boxSide> > getBndIn()
    {
        GISMO_ASSERT(m_isBndSet, "Boundary parts are not set in gsFlowSolverParams, call setBndParts(...).");
        return m_bndIn;
    }

    /// @brief Get vector of [patch, side] corresponding to the outflow boundary.
    std::vector<std::pair<int, boxSide> > getBndOut()
    {
        GISMO_ASSERT(m_isBndSet, "Boundary parts are not set in gsFlowSolverParams, call setBndParts(...).");
        return m_bndOut;
    }

    /// @brief Get vector of [patch, side] corresponding to the solid wall boundary.
    std::vector<std::pair<int, boxSide> > getBndWall()
    {
        GISMO_ASSERT(m_isBndSet, "Boundary parts are not set in gsFlowSolverParams, call setBndParts(...).");
        return m_bndWall;
    }

    gsField<T> getVelocitySolution() { return m_USolField; }
    gsField<T> getKSolution() { return m_KSolField; }
    gsField<T> getOmegaSolution() { return m_OSolField; }

    void setVelocitySolution(gsField<T> sol) { m_USolField = sol; }
    void setKSolution(gsField<T> sol) { m_KSolField = sol; }
    void setOmegaSolution(gsField<T> sol) { m_OSolField = sol; }

    void setDistanceField(gsField<T>& dfield) { m_distanceField = dfield; }

    gsField<T>& getDistanceField() { return m_distanceField; }

}; // class gsFlowSolverParams

// ======================================================================================================================================



} // namespace gismo