/** @file uwbINSSolverIterativeBase.h

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSAssemblerBase.h"
#include "uwbLinSolvers.h"
#include "uwbPreconditioners.h"

namespace gismo
{

template<class T, class SolverType>
class uwbINSSolverIterativeBase
{

public:
    uwbINSSolverIterativeBase()
    {}

    uwbINSSolverIterativeBase(uwbINSSolverParams<T>& params, uwbINSAssemblerBase<T>& assembler) :
        m_assemblerRef(assembler)
    {
        m_maxLinIter = params.settings().get(constantsINS::iter_maxIt);
        m_linTol = params.settings().get(constantsINS::iter_tol);
        m_precType = params.settings().getPrecondType();
        m_precOpt = params.getPrecOptions();
        
        initMembers();
    }

    virtual ~uwbINSSolverIterativeBase()
    {
    }

protected:
    virtual void initMembers()
    {
        m_assembT = 0;
        m_solsetupT = 0;
        m_solveT = 0;
        m_linIterVector.clear();

        m_precOpt.addInt("dim", "Problem dimension", m_assemblerRef.getTarDim());
        m_precOpt.addReal("visc", "Viscosity", m_assemblerRef.getViscosity());
        m_precOpt.addInt("udofs", "Number of velocity dofs", m_assemblerRef.getUdofs());
        m_precOpt.addInt("pdofs", "Number of pressure dofs", m_assemblerRef.getPdofs());
    }

    virtual void reinitMembers() { initMembers(); }

public:
    virtual void initialize(const gsMatrix<T>& solution)
    {
        GISMO_ASSERT(m_assemblerRef.isInitialized(), "Assembler must be initialized first, call initialize()");

        m_clock.restart();
        m_assemblerRef.getBlockAssembler().assembleMassMatrix(); // this is assembled 2 times in unsteady case (TODO?)
        m_assemblerRef.getBlockAssembler().assemblePressureMassMatrix();
        m_assemblerRef.update(solution);
        m_assembT += m_clock.stop();

        m_matrices.clear();
        m_matrices.insert(std::make_pair("matNS", m_assemblerRef.matrix()));
        m_matrices.insert(std::make_pair("matMu", m_assemblerRef.getVelocityMassMatrix()));
        m_matrices.insert(std::make_pair("matMp", m_assemblerRef.getPressureMassMatrix()));

        initPCDmat();

        m_clock.restart();
        m_pPrec = uwbINSPreconditioner<real_t>::make(m_precType, m_matrices, m_precOpt);
        m_solsetupT += m_clock.stop();
    }

    virtual void initIteration()
    {}

    virtual void applySolver(gsMatrix<T>& solution)
    {
        m_matrices.at("matNS") = m_assemblerRef.matrix();

        updatePCDmat();

        m_clock.restart();
        m_pPrec->update(m_matrices);
        
        SolverType solver(m_assemblerRef.matrix(), m_pPrec);
        solver.setMaxIterations(m_maxLinIter);
        solver.setTolerance(m_linTol);
        m_solsetupT += m_clock.stop();

        //solution.setZero(); // zero initial solution for iterative solver in each iteration

        m_clock.restart();
        solver.solve(m_assemblerRef.rhs(), solution);
        m_solveT += m_clock.stop();

        m_linIterVector.push_back(solver.iterations());
    }

    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
    {
        GISMO_ASSERT(solution.rows() > 0, "The solution in applySolver() is empty!");

        int usize = m_assemblerRef.getPshift();
        int pdofs = m_assemblerRef.getPdofs();

        m_matrices.at("matNS") = m_assemblerRef.matrix();

        updatePCDmat();

        m_clock.restart();
        m_pPrec->update(m_matrices);

        SolverType solver(m_assemblerRef.matrix(), m_pPrec);
        solver.setMaxIterations(m_maxLinIter);
        solver.setTolerance(m_linTol);
        m_solsetupT += m_clock.stop();

        gsMatrix<T> newsol(solution.rows(), solution.cols());
        //newsol.setZero();
        newsol = solution;

        m_clock.restart();
        solver.solve(m_assemblerRef.rhs(), newsol);
        solution.topRows(usize) = alpha_u * newsol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
        solution.bottomRows(pdofs) = alpha_p * newsol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
        m_solveT += m_clock.stop();

        m_linIterVector.push_back(solver.iterations());
    }

    virtual T residualRelNorm(const gsMatrix<T>& solution) const
    {
        gsMatrix<T> residual = m_assemblerRef.rhs() - m_assemblerRef.matrix() * solution;
        T resNorm = residual.norm() / m_assemblerRef.rhs().norm();

        return resNorm;
    }

    gsMatrix<T> getStokesSolution()
    {
        GISMO_ASSERT(m_assemblerRef.isInitialized(), "Assembler must be initialized first, call initialize()");

        gsInfo << "Setting Stokes initial condition...\n";

        gsSparseMatrix<T> stokesMatrix;
        gsMatrix<T> stokesRhs;

        m_assemblerRef.fillStokesSystem_into(stokesMatrix, stokesRhs);
        
        m_matrices.clear();
        m_matrices.insert(std::make_pair("matNS", m_assemblerRef.matrix()));
        m_matrices.insert(std::make_pair("matMu", m_assemblerRef.getVelocityMassMatrix()));
        m_matrices.insert(std::make_pair("matMp", m_assemblerRef.getPressureMassMatrix()));

        initPCDmat();

        typename uwbINSPreconditioner<T>::Ptr prec = uwbINSPreconditioner<real_t>::make(m_precType, m_matrices, m_precOpt);

        SolverType solver(stokesMatrix, prec);
        solver.setMaxIterations(m_maxLinIter);
        solver.setTolerance(m_linTol);

        gsMatrix<T> solution;
        solver.solve(stokesRhs, solution);

        return solution;
    }

    std::vector<index_t> getLinIterVector() const { return m_linIterVector; }

    T getAvgLinIterations() const
    {
        index_t linIterSum = 0;

        for (size_t i = 0; i < m_linIterVector.size(); i++)
            linIterSum += m_linIterVector[i];

        return (T)linIterSum / m_linIterVector.size();
    }

    const T getAssemblyTime() const { return m_assembT; }
    const T getSolverSetupTime() const { return m_solsetupT; }
    const T getSolveTime() const { return m_solveT; }

protected:
    void initPCDmat()
    {
        if (m_precType.substr(0, 3) == "PCD")
        {
            gsSparseMatrix<T> Ap, Fp;

            m_clock.restart();
            m_assemblerRef.fillPCDblocks(Ap, Fp, m_precOpt.getInt("pcd_bcType"), m_precOpt.getSwitch("pcd_assembAp"), m_precOpt.getSwitch("pcd_assembFp"), m_precOpt.getSwitch("lumpingM"));
            m_assembT += m_clock.stop();

            m_matrices.insert(std::make_pair("matFp", Fp));
            m_matrices.insert(std::make_pair("matAp", Ap));
        }
    }

    void updatePCDmat()
    {
        if (m_precType.substr(0, 3) == "PCD")
        {
            gsSparseMatrix<T> Ap, Fp;

            // TODO: rudundant assembly of Ap, modify!
            m_clock.restart();
            m_assemblerRef.fillPCDblocks(Ap, Fp, m_precOpt.getInt("pcd_bcType"), m_precOpt.getSwitch("pcd_assembAp"), m_precOpt.getSwitch("pcd_assembFp"), m_precOpt.getSwitch("lumpingM"));
            m_assembT += m_clock.stop();

            m_matrices.at("matFp") = Fp;
        }
    }

protected:

    std::string m_precType;
    gsOptionList m_precOpt;
    std::map<std::string, gsSparseMatrix<> > m_matrices;
    int m_maxLinIter;
    T m_linTol;
    typename uwbINSPreconditioner<T>::Ptr m_pPrec;
    uwbINSAssemblerBase<T>& m_assemblerRef;
    std::vector<index_t> m_linIterVector;

    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT;

}; //uwbINSSolverIterativeBase

} //namespace gismo
