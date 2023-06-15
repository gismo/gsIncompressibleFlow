/** @file gsINSSolverBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, E. Turnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSSolverBase.h>

namespace gismo
{

template<class T>
void gsINSSolverBase<T>::initMembers()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
    m_iterationNumber = 0;
    m_relNorm = std::numeric_limits<T>::infinity();
    m_assembT = 0;
    m_solsetupT = 0;
    m_solveT = 0;

    #ifdef GISMO_WITH_PARDISO
        pardisoSetup(m_solver);
    #endif
}


template<class T>
void gsINSSolverBase<T>::applySolver(gsMatrix<T>& solution)
{
    m_clock.restart();
    m_solver.factorize(getAssembler()->matrix());
    m_solsetupT += m_clock.stop();

    m_clock.restart();
    solution = m_solver.solve(getAssembler()->rhs());
    m_solveT += m_clock.stop();
}


template<class T>
void gsINSSolverBase<T>::applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
{
    GISMO_ASSERT(solution.rows() > 0,"The solution in applySolver() is empty!");

    int usize = getAssembler()->getPshift();
    int pdofs = getAssembler()->getPdofs();

    m_clock.restart();
    m_solver.factorize(getAssembler()->matrix());
    m_solsetupT += m_clock.stop();

    m_clock.restart();
    gsMatrix<T> newsol = m_solver.solve(getAssembler()->rhs());
    solution.topRows(usize) = alpha_u * newsol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
    solution.bottomRows(pdofs) = alpha_p * newsol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
    m_solveT += m_clock.stop();
}


template<class T>
void gsINSSolverBase<T>::solve(const int maxIterations, const T epsilon, const int minIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
    int iter = 0;
    m_relNorm = solutionChangeRelNorm();

    while ((iter < minIterations) || ((m_relNorm > epsilon) && (iter < maxIterations)))
    {
        gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

        nextIteration();
        m_relNorm = solutionChangeRelNorm();
        gsInfo << " Solution change relative norm: " << m_relNorm << "\n";

        iter++;
    }
}


template<class T>
void gsINSSolverBase<T>::solveStokes()
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    gsInfo << "Computing the steady Stokes problem...\n";

    gsSparseMatrix<T> stokesMatrix;
    gsMatrix<T> stokesRhs;

    getAssembler()->fillStokesSystem_into(stokesMatrix, stokesRhs);
    
    m_solver.analyzePattern(stokesMatrix);
    m_solver.factorize(stokesMatrix);
    m_solution = m_solver.solve(stokesRhs);
}


template<class T>
T gsINSSolverBase<T>::solutionChangeRelNorm() const
{
    T relNorm;

    if (m_iterationNumber)
    {
        gsMatrix<T> solChangeVector = getAssembler()->getSolution() - m_solution;
        relNorm = solChangeVector.norm() / m_solution.norm();

    }
    else
    {
        relNorm = std::numeric_limits<T>::infinity();
    }

    return relNorm;
}


template<class T>
T gsINSSolverBase<T>::solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsMatrix<T> solChangeVector = solOld - solNew;
    T relNorm = solChangeVector.norm() / solNew.norm();

    return relNorm;
}


template<class T>
void gsINSSolverBase<T>::dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsInfo << "     [u, p] solution change relative norm: ";

    for (int i = 0; i < solOld.cols(); i++)
        gsInfo << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

    gsInfo << "\n";
}


// ===========================================================================


template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::initMembers()
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


template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::initPrecMat()
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
    else if (m_precType.substr(0, 2) == "AL")
    {
        gsINSBlockPrecondAL<T>::fillALgammaPart_into(m_matGammaPart, m_rhsGammaPart, m_matrices, m_assemblerRef.rhs(), m_precOpt);
    }
}

template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::updatePrecMat()
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


template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::initialize(const gsMatrix<T>& solution)
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

    initPrecMat();

    if (m_precType.substr(0, 2) == "AL")
        m_matrices.at("matNS") += m_matGammaPart;

    m_clock.restart();
    m_pPrec = gsINSPreconditioner<real_t>::make(m_precType, m_matrices, m_precOpt);
    m_solsetupT += m_clock.stop();
}


template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::applySolver(gsMatrix<T>& solution)
{
    m_matrices.at("matNS") = m_assemblerRef.matrix();
    gsMatrix<T> rhs = m_assemblerRef.rhs();

    updatePrecMat();

    if (m_precType.substr(0, 2) == "AL")
    {
        m_matrices.at("matNS") += m_matGammaPart;
        rhs += m_rhsGammaPart;
    }

    m_clock.restart();
    m_pPrec->update(m_matrices);
    
    SolverType solver(m_matrices.at("matNS"), m_pPrec);
    solver.setMaxIterations(m_maxLinIter);
    solver.setTolerance(m_linTol);
    m_solsetupT += m_clock.stop();

    //solution.setZero(); // zero initial solution for iterative solver in each iteration

    m_clock.restart();
    solver.solve(rhs, solution);
    m_solveT += m_clock.stop();

    m_linIterVector.push_back(solver.iterations());
}


template<class T, class SolverType>
void gsINSSolverBaseIter<T,SolverType>::applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
{
    GISMO_ASSERT(solution.rows() > 0, "The solution in applySolver() is empty!");

    int usize = m_assemblerRef.getPshift();
    int pdofs = m_assemblerRef.getPdofs();

    m_matrices.at("matNS") = m_assemblerRef.matrix();
    gsMatrix<T> rhs = m_assemblerRef.rhs();

    updatePrecMat();

    if (m_precType.substr(0, 2) == "AL")
    {
        m_matrices.at("matNS") += m_matGammaPart;
        rhs += m_rhsGammaPart;
    }

    m_clock.restart();
    m_pPrec->update(m_matrices);

    SolverType solver(m_matrices.at("matNS"), m_pPrec);
    solver.setMaxIterations(m_maxLinIter);
    solver.setTolerance(m_linTol);
    m_solsetupT += m_clock.stop();

    gsMatrix<T> newsol = solution;
    //newsol.setZero(); // zero initial solution for iterative solver in each iteration

    m_clock.restart();
    solver.solve(rhs, newsol);
    solution.topRows(usize) = alpha_u * newsol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
    solution.bottomRows(pdofs) = alpha_p * newsol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
    m_solveT += m_clock.stop();

    m_linIterVector.push_back(solver.iterations());
}


template<class T, class SolverType>
gsMatrix<T> gsINSSolverBaseIter<T,SolverType>::getStokesSolution()
{
    GISMO_ASSERT(m_assemblerRef.isInitialized(), "Assembler must be initialized first, call initialize()");

    gsInfo << "Computing the steady Stokes problem...\n";

    gsSparseMatrix<T> stokesMatrix;
    gsMatrix<T> stokesRhs, stokesSol;

    m_assemblerRef.fillStokesSystem_into(stokesMatrix, stokesRhs);
    
    m_matrices.clear();
    m_matrices.insert(std::make_pair("matNS", stokesMatrix));
    m_matrices.insert(std::make_pair("matMu", m_assemblerRef.getVelocityMassMatrix()));
    m_matrices.insert(std::make_pair("matMp", m_assemblerRef.getPressureMassMatrix()));

    initPrecMat();

    if (m_precType.substr(0, 2) == "AL")
    {
        m_matrices.at("matNS") += m_matGammaPart;
        stokesRhs += m_rhsGammaPart;
    }

    typename gsINSPreconditioner<T>::Ptr prec = gsINSPreconditioner<real_t>::make(m_precType, m_matrices, m_precOpt);

    SolverType solver(m_matrices.at("matNS"), prec);
    solver.setMaxIterations(m_maxLinIter);
    solver.setTolerance(m_linTol);

    solver.solve(stokesRhs, stokesSol);

    return stokesSol;
}


template<class T, class SolverType>
T gsINSSolverBaseIter<T,SolverType>::getAvgLinIterations() const
{
    index_t linIterSum = 0;

    for (size_t i = 0; i < m_linIterVector.size(); i++)
        linIterSum += m_linIterVector[i];

    return (T)linIterSum / m_linIterVector.size();
}

} //namespace gismo
