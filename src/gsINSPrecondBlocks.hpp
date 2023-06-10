/** @file gsINSPrecondBlocks.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSPrecondBlocks.h>

namespace gismo
{

template <class T>
gsSparseSolver<T>* gsINSPrecondBlock<T>::createLinSolver()
{
    if (m_opt.getSwitch("iter"))
    {
        return new typename gsSparseSolver<T>::BiCGSTABILUT;
    }
    else
    {
        #ifdef GISMO_WITH_PARDISO
        return new typename gsSparseSolver<T>::PardisoLU;
        #else
        return new typename gsSparseSolver<T>::LU;
        #endif
    }
}


template <class T>
void gsINSPrecondBlock<T>::setupLinSolver(gsSparseSolver<T>& solver)
{
    if (m_opt.getSwitch("iter"))
    {
        typename gsSparseSolver<T>::BiCGSTABILUT* pSolver = dynamic_cast<typename gsSparseSolver<T>::BiCGSTABILUT*>(&solver);

        pSolver->preconditioner().setDroptol(m_opt.getReal("dropTol"));
        pSolver->preconditioner().setFillfactor(m_opt.getInt("fill"));
        pSolver->setMaxIterations(m_opt.getInt("maxIt"));
        pSolver->setTolerance(m_opt.getReal("tol"));
    }
    else
    {
        #ifdef GISMO_WITH_PARDISO
            typename gsSparseSolver<T>::PardisoLU* pSolver = dynamic_cast<typename gsSparseSolver<T>::PardisoLU*>(&solver);
            gsINSSolverBase<T>::pardisoSetup(*pSolver);
            #endif
    }
}


template <class T>
void gsINSPrecondBlockF<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
    x.resize(m_size, 1);

    int udofs = m_opt.getInt("udofs");
    for (int i = 0; i < m_opt.getInt("dim"); i++)
        x.middleRows(i*udofs, udofs) = m_pSolver->solve(input.middleRows(i*udofs, udofs));
}


template <class T>
void gsINSPrecondBlockFwhole<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
    x.resize(m_size, 1);

    x = m_pSolver->solve(input);
}


template <class T>
void gsINSPrecondBlockFdiag<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
    x.resize(m_size, 1);

    int udofs = m_opt.getInt("udofs");

    // solving system with lower block-triangular part of Agamma
    for (int i = 0; i < m_opt.getInt("dim"); i++)
        x.middleRows(i*udofs, udofs) = m_solvers[i]->solve(input.middleRows(i*udofs, udofs));

}


template <class T>
void gsINSPrecondBlockFmod<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
    x.resize(m_size, 1);

    int udofs = m_opt.getInt("udofs");

    // solving system with lower block-triangular part of Agamma
    for (int i = 0; i < m_opt.getInt("dim"); i++)
    {
        gsMatrix<T> rhsi = input.middleRows(i*udofs, udofs);

        for (int j = 0; j < i; j++)
            rhsi -= m_matRef.block(i*udofs, j*udofs, udofs, udofs) * x.middleRows(j*udofs, udofs);

        x.middleRows(i*udofs, udofs) = m_solvers[i]->solve(rhsi);
    }
}


template <class T>
void gsINSPrecondBlockBt<T>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    GISMO_ASSERT(input.rows() == this->cols(), "Wrong input size.");
    x.resize(this->rows(), 1);

    x.noalias() = m_mat * input;
}


template <class T>
void gsINSPrecondSchurLSC<T>::apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
{
    GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
    x.resize(this->rows(), 1);

    // solve -(mat1 * mat2^-1 * mat1) * x = input
    gsMatrix<T> tmp;
    tmp = m_pSolver->solve(-1 * input); // mat1 * tmp = -input, where tmp = mat2^-1 * mat1 * x
    x = m_pSolver->solve(m_mat2 * tmp); // mat1 * x = mat2 * tmp
}


template <class T>
void gsINSPrecondSchurPCD<T>::apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
{
    GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
    x.resize(this->rows(), 1);

    // solve - Ap * Fp^-1 * Mp * x = input
    gsMatrix<T> tmp1, tmp2;
    tmp1 = m_pSolver->solve(-1 * input); // Ap * y = - input
    tmp2 = m_matFp * tmp1; // z = Fp * y
    x = m_pMassSolver->solve(tmp2); // Mp * x = z
}


template <class T>
void gsINSPrecondSchurPCDmod<T>::apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
{
    GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
    x.resize(this->rows(), 1);

    // solve - Mp * Fp^-1 * Ap * x = input
    gsMatrix<T> tmp1, tmp2;
    tmp1 = m_pMassSolver->solve(-1 * input); // Mp * y = - input
    tmp2 = m_matFp * tmp1; // z = Fp * y
    x = m_pSolver->solve(tmp2);  // Ap * x = z
}

} // namespace gismo
