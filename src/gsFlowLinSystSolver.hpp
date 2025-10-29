/** @file gsFlowLinSystSolver.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>

namespace gismo
{

template<class T, int MatOrder>
T gsFlowLinSystSolver<T, MatOrder>::stopwatchStart()
{

#ifdef GISMO_WITH_MPI

    MPI_Barrier(m_paramsPtr->getMpiComm());
    return MPI_Wtime();

#else

    m_clock.restart();
    return 0.0;

#endif

}


template<class T, int MatOrder>
T gsFlowLinSystSolver<T, MatOrder>::stopwatchStop()
{

#ifdef GISMO_WITH_MPI

    MPI_Barrier(m_paramsPtr->getMpiComm());
    return MPI_Wtime();

#else

    return m_clock.stop();

#endif

}


template<class T, int MatOrder>
void gsFlowLinSystSolver<T, MatOrder>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p, index_t usize, index_t pdofs)
{
    GISMO_ASSERT(solution.rows() > 0, "The solution in applySolver() is empty!");

    gsMatrix<T> newSol;
    this->applySolver(mat, rhs, newSol);

    real_t time0 = stopwatchStart();
    solution.topRows(usize) = alpha_u * newSol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
    solution.bottomRows(pdofs) = alpha_p * newSol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
    real_t time1 = stopwatchStop();

    m_solveT += time1 - time0;
}


// ===================================================================================================================
// ===================================================================================================================


template<class T, int MatOrder>
void gsFlowLinSystSolver_direct<T, MatOrder>::setupSolver(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    m_solver.analyzePattern(mat);
    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}


template<class T, int MatOrder>
void gsFlowLinSystSolver_direct<T, MatOrder>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    real_t time0 = stopwatchStart();
    m_solver.factorize(mat);

    real_t time1 = stopwatchStop();
    
    solution = m_solver.solve(rhs);
    real_t time2 = stopwatchStop();

    m_setupT += time1 - time0;
    m_solveT += time2 - time1;
}

// ===================================================================================================================

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iter<T, MatOrder, SolverType>::setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    m_precPtr = gsIncompleteLUOp< gsSparseMatrix<T, MatOrder> >::make(mat);
    //m_precPtr = gsGaussSeidelOp< gsSparseMatrix<T, MatOrder>, gsGaussSeidel::symmetric >::make(mat);
    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iter<T, MatOrder, SolverType>::applySolver(const gsSparseMatrix<T, MatOrder>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    real_t time0 = stopwatchStart();

    this->setupPreconditioner(mat);
    
    SolverType solver(mat, m_precPtr);
    solver.setMaxIterations(m_paramsPtr->options().getInt("lin.maxIt"));
    solver.setTolerance(m_paramsPtr->options().getReal("lin.tol"));

    real_t time1 = stopwatchStop();
    
    solver.solve(rhs, solution);

    real_t time2 = stopwatchStop();

    m_setupT += time1 - time0;
    m_solveT += time2 - time1;

    m_linIterVector.push_back(solver.iterations());
}

// ===================================================================================================================

template<class T, int MatOrder, class SolverType>
void gsFlowLinSystSolver_iterSP<T, MatOrder, SolverType>::setupPreconditioner(const gsSparseMatrix<T, MatOrder>& mat)
{
    real_t time0 = stopwatchStart();
    
    m_matrices.clear();
    m_matrices.insert(std::make_pair("matNS", mat));
    m_matrices.insert(std::make_pair("matMu", m_assemblerPtr->getMassMatrix(0)));
    m_matrices.insert(std::make_pair("matMp", m_assemblerPtr->getMassMatrix(1)));

    if (m_precType.substr(0, 3) == "PCD")
    {
        gsWarn << "PCD preconditioner not fully implemented yet, using LSC instead.";
        m_precType.replace(0, 3, "LSC");

        // gsSparseMatrix<T> Ap, Fp;
        // m_assemblerPtr->fillPCDblocks(Ap, Fp, m_precOpt.getInt("pcd_bcType"), m_precOpt.getSwitch("pcd_assembAp"), m_precOpt.getSwitch("pcd_assembFp"), m_precOpt.getSwitch("lumpingM"));
        // m_matrices.insert(std::make_pair("matFp", Fp));
        // m_matrices.insert(std::make_pair("matAp", Ap));
    }

    m_precPtr = gsINSPreconditioner<T, MatOrder>::make(m_precType, m_matrices, m_precOpt);

    real_t time1 = stopwatchStop();

    m_setupT += time1 - time0;
}

// ===================================================================================================================

#ifdef gsPetsc_ENABLED

template<class T>
gsOptionList gsFlowLinSystSolver_PETSc<T>::getDefaultOptions()
{
    index_t maxIt = m_paramsPtr->options().getInt("lin.maxIt");
    real_t linTol = m_paramsPtr->options().getReal("lin.tol");

    gsOptionList opt;

    opt.addString("-ksp_type", "", "fgmres");
    opt.addString("-ksp_initial_guess_nonzero", "", "true");
    opt.addString("-ksp_max_it", "", util::to_string(maxIt));
    opt.addString("-ksp_rtol", "", util::to_string(linTol));
    
    opt.addString("-pc_type", "", "asm");
    opt.addString("-pc_asm_overlap", "", "1");
    opt.addString("-sub_pc_type", "", "ilu");
    opt.addString("-sub_pc_factor_levels", "", "0");
    opt.addString("-sub_pc_factor_mat_ordering_type", "", "rcm");
    opt.addString("-sub_pc_factor_reuse_ordering", "", "");

    return opt;
}

template<class T>
void gsFlowLinSystSolver_PETSc<T>::setupSolver(const gsSparseMatrix<T, RowMajor>& mat)
{
    real_t time0 = stopwatchStart();

    index_t M = mat.rows();
    index_t N = mat.cols();
    MPI_Comm comm = m_paramsPtr->getMpiComm();

    petsc_computeMatLayout(M, m_rowLocInfo, comm);
    petsc_computeMatLayout(N, m_colLocInfo, comm);
    petsc_setupMatrix(m_petscMat, M, N, comm);
    PetscCallVoid( MatCreateVecs(m_petscMat, &m_petscSol, &m_petscRhs) );

    PetscCallVoid( KSPCreate(comm, &m_ksp) );
    PetscCallVoid( KSPGetPC(m_ksp, &m_pc) );

    this->applyOptions();

    real_t time1 = stopwatchStop();
    m_setupT += time1 - time0;
}

template<class T>
void gsFlowLinSystSolver_PETSc<T>::applySolver(const gsSparseMatrix<T, RowMajor>& mat, const gsMatrix<T>& rhs, gsMatrix<T>& solution)
{
    MPI_Comm comm = m_paramsPtr->getMpiComm();

    real_t time0 = stopwatchStart();
    petsc_copySparseMat(mat, m_petscMat, m_rowLocInfo, m_colLocInfo, comm);
    petsc_copyVec(rhs, m_petscRhs, comm);
    real_t time1 = stopwatchStop();
    PetscCallVoid( KSPSetOperators(m_ksp, m_petscMat, m_petscMat) );
    PetscCallVoid( KSPSolve(m_ksp, m_petscRhs, m_petscSol) );
    real_t time2 = stopwatchStop();

    int nIter = 0;
    PetscCallVoid( KSPGetIterationNumber(m_ksp, &nIter) );
    m_linIterVector.push_back(nIter);

    real_t time3 = stopwatchStop();
    petsc_copyVecToGismo(m_petscSol, solution, comm);
    real_t time4 = stopwatchStop();

    // Clear petsc vector
    PetscCallVoid( VecZeroEntries(m_petscRhs) );

    m_solveT += time2 - time1;
    m_dataCopyT += (time1 - time0) + (time4 - time3);
}

template<class T>
void gsFlowLinSystSolver_PETSc<T>::applyOptions()
{
    if (m_paramsPtr->options().getSwitch("petsc.optFromFile"))
    {
        gsOptionList petscOpt;
        gsFileData<> fd(m_paramsPtr->options().getString("petsc.optPath"));
        fd.getId(0, petscOpt);
        applyOptions(petscOpt);
    }
    else
        applyOptions(getDefaultOptions());
}

template<class T>
void gsFlowLinSystSolver_PETSc<T>::applyOptions(gsOptionList petscOpt)
{
    PetscCallVoid( PetscOptionsClear(NULL) );

    for ( auto & opt : petscOpt.getAllEntries() )
        PetscCallVoid( PetscOptionsSetValue(NULL, opt.label.c_str(), opt.val.c_str()) );

    PetscCallVoid( KSPSetFromOptions(m_ksp) );
    PetscCallVoid( PCSetFromOptions(m_pc) );
}

// ===========================================================

template<class T>
gsOptionList gsFlowLinSystSolver_PETSc_SP<T>::getDefaultOptions()
{
    short_t dim = m_paramsPtr->getPde().dim();
    index_t maxIt = m_paramsPtr->options().getInt("lin.maxIt");
    real_t linTol = m_paramsPtr->options().getReal("lin.tol");

    gsOptionList opt;

    opt.addString("-ksp_type", "", "fgmres");
    opt.addString("-ksp_initial_guess_nonzero", "", "true");
    opt.addString("-ksp_max_it", "", util::to_string(maxIt));
    opt.addString("-ksp_rtol", "", util::to_string(linTol));

    opt.addString("-pc_type", "", "fieldsplit"); // PCFIELDSPLIT preconditioner
    opt.addString("-pc_fieldsplit_detect_saddle_point", "", ""); // automatically detect saddle-point structure
    opt.addString("-pc_fieldsplit_type", "", "schur"); // Schur complement preconditioner
    opt.addString("-pc_fieldsplit_schur_fact_type", "", "upper"); // upper triangular factorization
    opt.addString("-pc_fieldsplit_schur_precondition", "", "self"); // matrix to generate Schur complement preconditioner
    
    // subsystem solvers:
    // fielssplit_0 - velocity block
    // fieldsplit_1 - pressure block (Schur complement)

    // universal settings:
    opt.addString("-fieldsplit_0_mat_block_size", "", util::to_string(dim)); // block size = dimension of velocity (not sure, if this actually does anything)
    opt.addString("-fieldsplit_1_pc_type", "", "lsc"); // LSC approximation of Schur complement
    opt.addString("-fieldsplit_1_pc_lsc_scale_diag", "", ""); // scale with diagonal of velocity mass matrix

    // direct solver (MUMPS):
    // opt.addString("-fieldsplit_0_ksp_type", "", "preonly"); // use only preconditioner
    // opt.addString("-fieldsplit_0_pc_type", "", "lu"); // direct solver for velocity block
    // opt.addString("-fieldsplit_0_pc_factor_mat_solver_type", "", "mumps"); // use MUMPS as direct solver
    // opt.addString("-fieldsplit_1_lsc_ksp_type", "", "preonly"); // use only preconditioner
    // opt.addString("-fieldsplit_1_lsc_pc_type", "", "lu"); // direct solver for Schur complement approximation
    // opt.addString("-fieldsplit_1_lsc_pc_factor_mat_solver_type", "", "mumps"); // use MUMPS as direct solver

    // Schwarz method for subsystems:
    opt.addString("-fieldsplit_0_ksp_type", "", "gmres"); 
    opt.addString("-fieldsplit_0_ksp_rtol", "", "1e-2"); 
    opt.addString("-fieldsplit_0_pc_type", "", "asm"); 
    opt.addString("-fieldsplit_0_pc_asm_overlap", "", "1"); 
    opt.addString("-fieldsplit_0_sub_pc_type", "", "ilu"); 
    opt.addString("-fieldsplit_0_sub_pc_factor_levels", "", "0"); 
    opt.addString("-fieldsplit_0_sub_pc_factor_mat_ordering_type", "", "rcm"); 
    opt.addString("-fieldsplit_0_sub_pc_factor_reuse_ordering", "", ""); 
    // opt.addString("-fieldsplit_1_lsc_ksp_type", "", "gmres");
    // opt.addString("-fieldsplit_1_lsc_ksp_rtol", "", "1e-2");
    // opt.addString("-fieldsplit_1_lsc_pc_type", "", "asm");
    // opt.addString("-fieldsplit_1_lsc_pc_asm_overlap", "", "1");
    // opt.addString("-fieldsplit_1_lsc_sub_pc_type", "", "ilu");
    // opt.addString("-fieldsplit_1_lsc_sub_pc_factor_levels", "", "0");
    // opt.addString("-fieldsplit_1_lsc_sub_pc_factor_mat_ordering_type", "", "rcm");
    // opt.addString("-fieldsplit_1_lsc_sub_pc_factor_reuse_ordering", "", "");

    // Jacobi for Schur complement:
    opt.addString("-fieldsplit_1_lsc_ksp_type", "", "cg");
    opt.addString("-fieldsplit_1_lsc_ksp_rtol", "", "1e-2");
    opt.addString("-fieldsplit_1_lsc_pc_type", "", "jacobi");

    return opt;
}

template<class T>
void gsFlowLinSystSolver_PETSc_SP<T>::applyOptions()
{
    if (m_paramsPtr->options().getSwitch("petsc.optFromFileSP"))
    {
        gsOptionList petscOpt;
        gsFileData<> fd(m_paramsPtr->options().getString("petsc.optPathSP"));
        fd.getId(0, petscOpt);
        Base::applyOptions(petscOpt);
    }
    else
        Base::applyOptions(getDefaultOptions());
}

template<class T>
void gsFlowLinSystSolver_PETSc_SP<T>::setupSolver(const std::vector< gsSparseMatrix<T, RowMajor> >& matBlocks)
{
    real_t time0 = stopwatchStart();

    MPI_Comm comm = m_paramsPtr->getMpiComm();

    short_t dim = m_paramsPtr->getPde().patches().dim();
    index_t usize = matBlocks[1].rows();
    index_t nu = usize / dim;
    index_t np = matBlocks[1].cols();

    petsc_computeMatLayout(nu, m_uLocInfo, comm);
    petsc_computeMatLayout(np, m_pLocInfo, comm);

    PetscCallVoid( KSPCreate(comm, &m_ksp) );
    PetscCallVoid( KSPGetPC(m_ksp, &m_pc) );

    this->applyOptions();

    real_t time1 = stopwatchStop();
    m_setupT += time1 - time0;
}

template<class T>
void gsFlowLinSystSolver_PETSc_SP<T>::applySolver(const std::vector< gsSparseMatrix<T, RowMajor> >& matBlocks, const std::vector< gsMatrix<T> >& rhsBlocks, gsMatrix<T>& solution)
{
    MPI_Comm comm = m_paramsPtr->getMpiComm();

    short_t dim = m_paramsPtr->getPde().patches().dim();
    index_t usize = matBlocks[1].rows();
    index_t nu = usize / dim;
    index_t np = matBlocks[1].cols();
    
    Mat pMatBlocks[4];
    Vec pRhsBlocks[2];

    real_t time0 = stopwatchStop();

    // velocity-velocity block
    petsc_setupMatrix(pMatBlocks[0], usize, usize, comm, dim, dim);
    petsc_copySparseMat(matBlocks[0], pMatBlocks[0], m_uLocInfo, m_uLocInfo, comm, dim, dim);

    // velocity-pressure block
    petsc_setupMatrix(pMatBlocks[1], usize, np, comm, dim, 1);
    petsc_copySparseMat(matBlocks[1], pMatBlocks[1], m_uLocInfo, m_pLocInfo, comm, dim, 1);

    // pressure-velocity block
    petsc_setupMatrix(pMatBlocks[2], np, usize, comm, 1, dim);
    petsc_copySparseMat(matBlocks[2], pMatBlocks[2], m_pLocInfo, m_uLocInfo, comm, 1, dim);

    // pressure-pressure block
    if (matBlocks.size() == 4)
    {
        petsc_setupMatrix(pMatBlocks[3], np, np, comm);
        petsc_copySparseMat(matBlocks[3], pMatBlocks[3], m_pLocInfo, m_pLocInfo, comm);
    }
    else
        pMatBlocks[3] = NULL;

    // rhs
    PetscCallVoid( MatCreateVecs(pMatBlocks[0], NULL, &pRhsBlocks[0]) ); // velocity rhs
    PetscCallVoid( MatCreateVecs(pMatBlocks[2], NULL, &pRhsBlocks[1]) ); // pressure rhs

    gsMatrix<T> locRhsU, locRhsP;
    locRhsU.resize(m_uLocInfo.first, dim);
    locRhsP = rhsBlocks[1].middleRows(m_pLocInfo.second, m_pLocInfo.first);

    for (index_t d = 0; d < dim; d++)
        locRhsU.col(d) = rhsBlocks[0].middleRows(m_uLocInfo.second + d*nu, m_uLocInfo.first);

    petsc_copyVec(locRhsU, pRhsBlocks[0], comm);
    petsc_copyVec(locRhsP, pRhsBlocks[1], comm);

    real_t time1 = stopwatchStop();

    PetscCallVoid( MatCreateNest(comm, 2, NULL, 2, NULL, pMatBlocks, &m_petscMat) );
    PetscCallVoid( MatNestSetVecType(m_petscMat, VECNEST) );
    PetscCallVoid( VecCreateNest(comm, 2, NULL, pRhsBlocks, &m_petscRhs) );
    PetscCallVoid( MatCreateVecs(m_petscMat, &m_petscSol, NULL) );  
    PetscCallVoid( KSPSetOperators(m_ksp, m_petscMat, m_petscMat) );

    real_t time2 = stopwatchStop();
    PetscCallVoid( KSPSolve(m_ksp, m_petscRhs, m_petscSol) );
    real_t time3 = stopwatchStop();

    int nIter = 0;
    PetscCallVoid( KSPGetIterationNumber(m_ksp, &nIter) );
    m_linIterVector.push_back(nIter);

    Vec pSolU, pSolP;
    PetscCallVoid( VecNestGetSubVec(m_petscSol, 0, &pSolU) );
    PetscCallVoid( VecNestGetSubVec(m_petscSol, 1, &pSolP) );

    real_t time4 = stopwatchStop();

    gsMatrix<T> solU, solP;
    petsc_copyVecToGismo(pSolU, solU, comm, dim);
    petsc_copyVecToGismo(pSolP, solP, comm, 1);

    solution.topRows(usize) = solU;
    solution.bottomRows(np) = solP;
    
    real_t time5 = stopwatchStop();

    // Clear petsc vector
    PetscCallVoid( VecZeroEntries(m_petscRhs) );

    m_solveT += time3 - time2;
    m_dataCopyT += (time1 - time0) + (time5 - time4);
}

#endif

}