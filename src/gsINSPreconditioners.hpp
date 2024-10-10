/** @file gsINSPreconditioners.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSPreconditioners.h>

namespace gismo
{

template <class T, int MatOrder>
typename gsINSPreconditioner<T, MatOrder>::uPtr gsINSPreconditioner<T, MatOrder>::make(std::string precType, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
{
    if (precType == "LSC_FdiagEqual")
        return gsINSBlockPrecondLSC<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "LSC_Fdiag")
        return gsINSBlockPrecondLSC<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "LSC_Fwhole")
        return gsINSBlockPrecondLSC<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "LSC_Fmod")
        return gsINSBlockPrecondLSC<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCD_FdiagEqual")
        return gsINSBlockPrecondPCD<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCD_Fdiag")
        return gsINSBlockPrecondPCD<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCD_Fwhole")
        return gsINSBlockPrecondPCD<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCD_Fmod")
        return gsINSBlockPrecondPCD<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCDmod_FdiagEqual")
        return gsINSBlockPrecondPCDmod<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCDmod_Fdiag")
        return gsINSBlockPrecondPCDmod<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCDmod_Fwhole")
        return gsINSBlockPrecondPCDmod<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "PCDmod_Fmod")
        return gsINSBlockPrecondPCDmod<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "AL_Fwhole")
        return gsINSBlockPrecondAL<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "AL_Fmod")
        return gsINSBlockPrecondAL<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLE_FdiagEqual")
        return gsINSBlockPrecondSIMPLE<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLE_Fdiag")
        return gsINSBlockPrecondSIMPLE<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLE_Fwhole")
        return gsINSBlockPrecondSIMPLE<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLE_Fmod")
        return gsINSBlockPrecondSIMPLE<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLER_FdiagEqual")
        return gsINSBlockPrecondSIMPLER<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLER_Fdiag")
        return gsINSBlockPrecondSIMPLER<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLER_Fwhole")
        return gsINSBlockPrecondSIMPLER<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "SIMPLER_Fmod")
        return gsINSBlockPrecondSIMPLER<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "MSIMPLER_FdiagEqual")
        return gsINSBlockPrecondMSIMPLER<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "MSIMPLER_Fdiag")
        return gsINSBlockPrecondMSIMPLER<T, MatOrder, gsINSPrecondBlockFdiag<T, MatOrder> >::make(mat, opt);
    else if (precType == "MSIMPLER_Fwhole")
        return gsINSBlockPrecondMSIMPLER<T, MatOrder, gsINSPrecondBlockFwhole<T, MatOrder> >::make(mat, opt);
    else if (precType == "MSIMPLER_Fmod")
        return gsINSBlockPrecondMSIMPLER<T, MatOrder, gsINSPrecondBlockFmod<T, MatOrder> >::make(mat, opt);
    else if (precType == "StokesDiag_FdiagEqual")
        return gsBlockPrecondStokes<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    else if (precType == "StokesTriang_FdiagEqual")
        return gsBlockPrecondStokesTriang<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    /*else if (precType == "SchurEx_LSC_FdiagEqual")
        return uwbINSBlockPrecondSchurEx<T, gsINSPrecondBlockF<T>, gsINSPrecondSchurLSC<T> >::make(mat, opt);
    else if (precType == "SchurEx_SIMPLE_FdiagEqual")
        return uwbINSBlockPrecondSchurEx<T, gsINSPrecondBlockF<T>, uwbINSPrecondSchurSIMPLE<T> >::make(mat, opt);
    else if (precType == "SchurEx_MSIMPLER_FdiagEqual")
        return uwbINSBlockPrecondSchurEx<T, gsINSPrecondBlockF<T>, gsINSPrecondSchurMSIMPLER<T> >::make(mat, opt);*/
    else
    {
        gsInfo << "Invalid preconditioner type, using LSC_FdiagEqual.\n";
        return gsINSBlockPrecondLSC<T, MatOrder, gsINSPrecondBlockF<T, MatOrder> >::make(mat, opt);
    }
}


template <class T, int MatOrder>
gsOptionList gsINSPreconditioner<T, MatOrder>::defaultOptions()
{
    gsOptionList opt;

    opt.addInt("dim", "Problem dimension", 0);
    opt.addReal("visc", "Viscosity", 0.0);
    opt.addInt("udofs", "Number of velocity dofs", 0);
    opt.addInt("pdofs", "Number of pressure dofs", 0);

    opt.addInt("pcd_bcType", "Type of BCs for Fp, Ap in PCD preconditioner", 3);
    opt.addReal("gamma", "Parameter gamma for AL", 1.0);
    opt.addReal("alphaP", "Pressure relaxation parameter for SIMPLE-type", 1.0);
    opt.addSwitch("lumpingM", "Use lumped diagonal mass matrices in preconditioners", true);
    opt.addSwitch("lumpingA", "Use lumped diagonal blockA in SIMPLE-type preconditioners", false);
    opt.addSwitch("pcd_assembAp", "Assemble pressure Poisson matrix for definition of Ap in PCD", true);
    opt.addSwitch("pcd_assembFp", "Assemble pressure Poisson matrix for definition of Fp in PCD", true);

    opt.addSwitch("iter", "Subsystems iteratively", false);
    opt.addInt("maxIt", "Max number of iterations for inner solves", 10);
    opt.addReal("tol", "Stopping tolerance for inner solves", 1e-2);
    opt.addInt("fill", "Fill factor for ILUT precond (inner solves)", 2);
    opt.addReal("dropTol", "Drop tolerance for ILUT precond  (inner solves)", 1e-8);

    return opt;
}


template <class T, int MatOrder>
void gsINSBlockPrecondBase<T, MatOrder>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    int uSize = m_Finv->rows();
    int pSize = m_Sinv->rows();

    x.resize(uSize + pSize, 1);

    gsMatrix<T> tmpX;
    m_Sinv->apply(input.bottomRows(pSize), tmpX);
    x.bottomRows(pSize) = tmpX;

    gsMatrix<T> tmp;
    m_Bt->apply(x.bottomRows(pSize), tmp);
    tmp = input.topRows(uSize) - tmp;

    m_Finv->apply(tmp, tmpX);
    x.topRows(uSize) = tmpX;
}


template <class T, int MatOrder, class BlockFType>
void gsINSBlockPrecondAL<T, MatOrder, BlockFType>::fillALgammaPart_into(gsSparseMatrix<T, MatOrder>& matGammaPart, gsMatrix<T>& rhsGammaPart, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsMatrix<T>& rhs, const gsOptionList& opt)
{
    gsInfo << "Filling the extra part of Agamma matrix... ";

    int dim = opt.getInt("dim");
    int udofs = opt.getInt("udofs");
    int pdofs = opt.getInt("pdofs");
    int uSize = dim * udofs;
    int numDofs = uSize + pdofs;
    real_t gamma = opt.getReal("gamma");

    const gsSparseMatrix<T, MatOrder>& matNS = mat.at("matNS");
    const gsSparseMatrix<T, MatOrder>& presM = mat.at("matMp");

    gsSparseMatrix<T, MatOrder> presMinv(pdofs, pdofs); // approximation of pressure mass matrix inverse
    presMinv.setIdentity();
    for (int i = 0; i < pdofs; i++)
        presMinv.coeffRef(i, i) = 1 / presM.coeff(i, i);

    gsSparseMatrix<T, MatOrder> Mgamma = gamma * matNS.block(0, uSize, uSize, pdofs) * presMinv * matNS.block(uSize, 0, pdofs, uSize);

    gsVector<index_t> nonZerosVector;
    nonZerosVector.setZero(numDofs);
    nonZerosVector.topRows(uSize) = getNnzVectorPerOuter(Mgamma);    

    matGammaPart.resize(numDofs, numDofs);
    matGammaPart.reserve(nonZerosVector);

    for (int outer = 0; outer < Mgamma.outerSize(); ++outer)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(Mgamma, outer); it; ++it)
            matGammaPart.insert(it.row(), it.col()) = it.value();

    matGammaPart.makeCompressed();

    rhsGammaPart.setZero(numDofs, 1);
    rhsGammaPart.topRows(uSize) = gamma * matNS.block(0, uSize, uSize, pdofs) * presMinv * rhs.bottomRows(pdofs);

    gsInfo << "Done.\n";
}

template <class T, int MatOrder, class BlockFType>
void gsINSBlockPrecondAL<T, MatOrder, BlockFType>::fillALmodifSystem_into(gsSparseMatrix<T, MatOrder>& matGamma, gsMatrix<T>& rhsGamma, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsMatrix<T>& rhs, const gsOptionList& opt)
{
    gsInfo << "Filling the Agamma matrix... ";

    gsSparseMatrix<T, MatOrder> matGammaPart;
    gsMatrix<T> rhsGammaPart;

    fillALgammaPart_into(matGammaPart, rhsGammaPart, mat, rhs, opt);

    matGamma = mat.at("matNS") + matGammaPart;
    rhsGamma = rhs + rhsGammaPart;

    gsInfo << "Done.\n";
}


template <class T, int MatOrder, class BlockFType>
void gsINSBlockPrecondSIMPLE<T, MatOrder, BlockFType>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    int uSize = m_Finv->rows();
    int pSize = m_Sinv->rows();

    x.resize(uSize + pSize, 1);

    gsMatrix<T> uStar;
    m_Finv->apply(input.topRows(uSize), uStar);

    gsMatrix<T> tmp1;
    tmp1.noalias() = m_B * uStar;

    gsMatrix<T> dp;
    m_Sinv->apply(input.bottomRows(pSize) - tmp1, dp);

    gsMatrix<T> tmp2;
    m_Bt->apply(dp, tmp2);

    x.topRows(uSize).noalias() = uStar - m_Dinv * tmp2;
    x.bottomRows(pSize) = m_alphaP * dp;
}


template <class T, int MatOrder, class BlockFType>
void gsINSBlockPrecondSIMPLER<T, MatOrder, BlockFType>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    int uSize = m_Finv->rows();
    int pSize = m_Sinv->rows();

    x.resize(uSize + pSize, 1);

    // pressure estimate
    gsMatrix<T> tmp1;
    tmp1.noalias() = input.bottomRows(pSize) - m_BDinv * input.topRows(uSize);
    gsMatrix<T> pStar;
    m_Sinv->apply(tmp1, pStar);

    // velocity solve
    gsMatrix<T> BtpStar, uStar;
    m_Bt->apply(pStar, BtpStar);
    m_Finv->apply(input.topRows(uSize) - BtpStar, uStar);

    // pressure correction
    gsMatrix<T> dp;
    m_Sinv->apply(input.bottomRows(pSize) - m_B * uStar, dp);

    gsMatrix<T> tmp2;
    m_Bt->apply(dp, tmp2);

    x.topRows(uSize).noalias() = uStar - m_Dinv * tmp2;
    x.bottomRows(pSize) = pStar + m_alphaP * dp;
}


template <class T, int MatOrder, class BlockFType>
void gsINSBlockPrecondMSIMPLER<T, MatOrder, BlockFType>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    int uSize = m_Finv->rows();
    int pSize = m_Sinv->rows();

    x.resize(uSize + pSize, 1);

    // pressure estimate
    gsMatrix<T> tmp1;
    tmp1.noalias() = input.bottomRows(pSize) - m_BMinv * input.topRows(uSize);
    gsMatrix<T> pStar;
    m_Sinv->apply(tmp1, pStar);

    // velocity solve
    gsMatrix<T> BtpStar, uStar;
    m_Bt->apply(pStar, BtpStar);
    m_Finv->apply(input.topRows(uSize) - BtpStar, uStar);

    // pressure correction
    gsMatrix<T> dp;
    m_Sinv->apply(input.bottomRows(pSize) - m_B * uStar, dp);

    gsMatrix<T> tmp2;
    m_Bt->apply(dp, tmp2);

    x.topRows(uSize).noalias() = uStar - m_velMinv * tmp2;
    x.bottomRows(pSize) = pStar + m_alphaP * dp;
}


template <class T, int MatOrder, class BlockFType>
void gsBlockPrecondStokes<T, MatOrder, BlockFType>::apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
{
    int uSize = m_Finv->rows();
    int pSize = m_Sinv->rows();

    x.resize(uSize + pSize, 1);

    gsMatrix<T> tmpX;
    m_Sinv->apply(input.bottomRows(pSize), tmpX);
    x.bottomRows(pSize) = tmpX;

    m_Finv->apply(input.topRows(uSize), tmpX);
    x.topRows(uSize) = tmpX;
}

// ===============================================================================



} // namespace gismo