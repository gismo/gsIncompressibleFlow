/** @file uwbPreconditioners.h

Author(s): H. Hornikova
*/
#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

#include "uwbPreconditionerBlocks.h"

namespace gismo
{

// forward declaration
template <class T>
class uwbINSPreconditioner;

// ======================================== Block preconditioner ========================================

/** @brief
Base class for block preconditioners for linear systems arising from linearized Navier-Stokes of form
\f[
\left[ \begin{array}{cc}
A & B^T \\
B & 0
\end{array}
\right]
\left[ \begin{array}{c}
u \\
p
\end{array}
\right] =
\left[ \begin{array}{c}
f \\
g
\end{array}
\right],
\f].
*/

template <class T>
class uwbINSBlockPrecondBase : public uwbINSPreconditioner<T>
{
public:
    /// Shared pointer for uwbINSBlockPrecondBase
    typedef memory::shared_ptr<uwbINSBlockPrecondBase> Ptr;

    /// Unique pointer for uwbINSBlockPrecondBase
    typedef memory::unique_ptr<uwbINSBlockPrecondBase> uPtr;

    uwbINSBlockPrecondBase(uwbINSPrecondBlock<T>* Ainv, gsLinearOperator<T>* B, uwbINSPrecondBlock<T>* Sinv, const gsOptionList& opt) :
        m_Ainv(Ainv), m_Sinv(Sinv), m_Bt(B), m_opt(opt)
    {}

    /// Make function returning a smart pointer
    static uPtr make()
    {
        return memory::make_unique(new uwbINSBlockPrecondBase<T>());
    }

    virtual int rows() const
    {
        return (m_Ainv->rows() + m_Sinv->rows());
    }

    virtual int cols() const
    {
        return (m_Ainv->cols() + m_Bt->cols());
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
        int pSize = m_Sinv->rows();

        x.resize(uSize + pSize, 1);

        gsMatrix<T> tmpX;
        m_Sinv->apply(input.bottomRows(pSize), tmpX);
        x.bottomRows(pSize) = tmpX;

        gsMatrix<T> tmp;
        m_Bt->apply(x.bottomRows(pSize), tmp);
        tmp = input.topRows(uSize) - tmp;

        m_Ainv->apply(tmp, tmpX);
        x.topRows(uSize) = tmpX;
    }

protected:
    void update(uwbINSPrecondBlock<T>* Ainv)
    {
        m_Ainv.reset(Ainv);
    }

    void update(uwbINSPrecondBlock<T>* Ainv, uwbINSPrecondBlock<T>* Sinv)
    {
        m_Ainv.reset(Ainv);
        m_Sinv.reset(Sinv);
    }

protected:
    typename uwbINSPrecondBlock<T>::Ptr m_Ainv, m_Sinv;
    typename gsLinearOperator<T>::Ptr m_Bt;
    const gsOptionList m_opt;

}; // class uwbINSBlockPrecondBase

// ======================================== LSC precond ========================================

/** @brief
    Least-squares commutator preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondLSC : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondLSC
    typedef memory::shared_ptr<uwbINSBlockPrecondLSC> Ptr;

    /// Unique pointer for uwbINSBlockPrecondLSC
    typedef memory::unique_ptr<uwbINSBlockPrecondLSC> uPtr;

    uwbINSBlockPrecondLSC(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
             new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
             new uwbINSPrecondSchurLSC<T>(mat, opt), opt)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondLSC<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt), new uwbINSPrecondSchurLSC<T>(mat, m_opt));
    }

    virtual std::string getName()
    {
        std::string name = "LSC_" + m_Ainv->getName();
        return name;
    }

protected:

    using Base::m_Ainv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondLSC

// ======================================== PCD precond ========================================

/*
Pressure convection-diffusion preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondPCD : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondPCD
    typedef memory::shared_ptr<uwbINSBlockPrecondPCD> Ptr;

    /// Unique pointer for uwbINSBlockPrecondPCD
    typedef memory::unique_ptr<uwbINSBlockPrecondPCD> uPtr;

    uwbINSBlockPrecondPCD(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurPCD<T>(mat, opt), opt)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondPCD<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt), new uwbINSPrecondSchurPCD<T>(mat, m_opt));
    }

    virtual std::string getName()
    {
        std::string name = "PCD_" + m_Ainv->getName();
        return name;
    }

protected:

    using Base::m_Ainv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondPCD

// ======================================== PCDmod precond ========================================

/*
Modified pressure convection-diffusion preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondPCDmod : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondPCDmod
    typedef memory::shared_ptr<uwbINSBlockPrecondPCDmod> Ptr;

    /// Unique pointer for uwbINSBlockPrecondPCDmod
    typedef memory::unique_ptr<uwbINSBlockPrecondPCDmod> uPtr;

    uwbINSBlockPrecondPCDmod(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurPCDmod<T>(mat, opt), opt)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondPCDmod<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt), new uwbINSPrecondSchurPCDmod<T>(mat, m_opt));
    }

    virtual std::string getName()
    {
        std::string name = "PCDmod_" + m_Ainv->getName();
        return name;
    }

protected:

    using Base::m_Ainv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondPCDmod

// ======================================== AL precond ========================================

/** @brief
    Augmented Lagrangian preconditioner.
*/

template <class T, class BlockAType = uwbINSPrecondBlockAwhole<T> >
class uwbINSBlockPrecondAL : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondAL
    typedef memory::shared_ptr<uwbINSBlockPrecondAL> Ptr;

    /// Unique pointer for uwbINSBlockPrecondAL
    typedef memory::unique_ptr<uwbINSBlockPrecondAL> uPtr;

    uwbINSBlockPrecondAL(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurAL<T>(mat, opt), opt)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondAL<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt));
    }

    virtual std::string getName()
    {
        std::string name = "AL_" + m_Ainv->getName();
        return name;
    }

    static void fillALmodifSystem_into(gsSparseMatrix<T>& matGamma, gsMatrix<T>& rhsGamma, const std::map<std::string, gsSparseMatrix<T> >& mat, const gsMatrix<T>& rhs, const gsOptionList& opt)
    {
        gsInfo << "Filling the Agamma matrix... ";

        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int pdofs = opt.getInt("pdofs");
        int uSize = dim * udofs;
        int numDofs = uSize + pdofs;
        real_t gamma = opt.getReal("gamma");

        const gsSparseMatrix<T>& matNS = mat.at("matNS");
        const gsSparseMatrix<T>& presM = mat.at("matMp");

        gsSparseMatrix<> presMinv(pdofs, pdofs); // approximation of pressure mass matrix inverse
        presMinv.setIdentity();
        for (int i = 0; i < pdofs; i++)
            presMinv.coeffRef(i, i) = 1 / presM.coeff(i, i);

        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(numDofs);

        gsSparseMatrix<> Mgamma = gamma * matNS.block(0, uSize, uSize, pdofs) * presMinv * matNS.block(uSize, 0, pdofs, uSize);

        for (int i = 0; i < uSize; i++)
            nonZerosPerColumnVector(i) += Mgamma.col(i).nonZeros();

        matGamma = matNS;
        matGamma.reserve(nonZerosPerColumnVector);

        for (int col = 0; col < uSize; ++col)
            for (typename gsSparseMatrix<>::InnerIterator it(Mgamma, col); it; ++it)
                matGamma.coeffRef(it.row(), col) += Mgamma(it.row(), col);

        matGamma.makeCompressed();

        rhsGamma = rhs;
        rhsGamma.topRows(uSize) += gamma * matNS.block(0, uSize, uSize, pdofs) * presMinv * rhs.bottomRows(pdofs);

        gsInfo << "Done.\n";
    }

protected:
    using Base::m_Ainv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondAL

// ======================================== SIMPLE precond ========================================

/** @brief
    SIMPLE preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondSIMPLE : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondSIMPLE
    typedef memory::shared_ptr<uwbINSBlockPrecondSIMPLE> Ptr;

    /// Unique pointer for uwbINSBlockPrecondSIMPLE
    typedef memory::unique_ptr<uwbINSBlockPrecondSIMPLE> uPtr;

    uwbINSBlockPrecondSIMPLE(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurSIMPLE<T>(mat, opt), opt)
    {
        int uSize = opt.getInt("dim") * opt.getInt("udofs");
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T>& matNS = mat.at("matNS");

        gsSparseMatrix<T> blockA = matNS.block(0, 0, uSize, uSize);
        uwbINSPreconditioner<T>::diagInvMatrix_into(blockA, m_Dinv, 1, opt.getSwitch("lumpingA"));

        m_B = matNS.block(uSize, 0, pdofs, uSize);

        m_alphaP = opt.askReal("alphaP", 1.0);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondSIMPLE<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt), new uwbINSPrecondSchurSIMPLE<T>(mat, m_opt));
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
        int pSize = m_Sinv->rows();

        x.resize(uSize + pSize, 1);

        gsMatrix<T> uStar;
        m_Ainv->apply(input.topRows(uSize), uStar);

        gsMatrix<T> tmp1;
        tmp1.noalias() = m_B * uStar;

        gsMatrix<T> dp;
        m_Sinv->apply(input.bottomRows(pSize) - tmp1, dp);

        gsMatrix<T> tmp2;
        m_Bt->apply(dp, tmp2);

        x.topRows(uSize).noalias() = uStar - m_Dinv * tmp2;
        x.bottomRows(pSize) = m_alphaP * dp;
    }

    virtual std::string getName()
    {
        std::string name = "SIMPLE_" + m_Ainv->getName();
        return name;
    }


protected:
    gsSparseMatrix<T> m_Dinv, m_B;
    real_t m_alphaP;

    using Base::m_Ainv;
    using Base::m_Bt;
    using Base::m_Sinv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondSIMPLE

// ======================================== SIMPLER precond ========================================

/** @brief
    SIMPLER preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondSIMPLER : public uwbINSBlockPrecondSIMPLE<T, BlockAType>
{
public:
    typedef uwbINSBlockPrecondSIMPLE<T, BlockAType> Base;

    /// Shared pointer for uwbINSBlockPrecondSIMPLER
    typedef memory::shared_ptr<uwbINSBlockPrecondSIMPLER> Ptr;

    /// Unique pointer for uwbINSBlockPrecondSIMPLER
    typedef memory::unique_ptr<uwbINSBlockPrecondSIMPLER> uPtr;

    uwbINSBlockPrecondSIMPLER(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(mat, opt)
    {
        m_BDinv = m_B * m_Dinv;
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondSIMPLER<T, BlockAType>(mat, opt));
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
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
        m_Ainv->apply(input.topRows(uSize) - BtpStar, uStar);

        // pressure correction
        gsMatrix<T> dp;
        m_Sinv->apply(input.bottomRows(pSize) - m_B * uStar, dp);

        gsMatrix<T> tmp2;
        m_Bt->apply(dp, tmp2);

        x.topRows(uSize).noalias() = uStar - m_Dinv * tmp2;
        x.bottomRows(pSize) = pStar + m_alphaP * dp;
    }

    virtual std::string getName()
    {
        std::string name = "SIMPLER_" + m_Ainv->getName();
        return name;
    }


protected:
    gsSparseMatrix<T> m_BDinv;

    using Base::m_Dinv;
    using Base::m_B;
    using Base::m_alphaP;
    using uwbINSBlockPrecondBase<T>::m_Ainv;
    using uwbINSBlockPrecondBase<T>::m_Bt;
    using uwbINSBlockPrecondBase<T>::m_Sinv;

}; // class uwbINSBlockPrecondSIMPLER

// ======================================== MSIMPLER precond ========================================

/** @brief
    MSIMPLER preconditioner.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSBlockPrecondMSIMPLER : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondMSIMPLER
    typedef memory::shared_ptr<uwbINSBlockPrecondMSIMPLER> Ptr;

    /// Unique pointer for uwbINSBlockPrecondMSIMPLER
    typedef memory::unique_ptr<uwbINSBlockPrecondMSIMPLER> uPtr;

    uwbINSBlockPrecondMSIMPLER(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurMSIMPLER<T>(mat, opt), opt)
    {
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int uSize = dim * udofs;
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T>& matNS = mat.at("matNS");
        const gsSparseMatrix<T>& velM = mat.at("matMu");

       // approximation of velocity mass matrix inverse
        uwbINSPreconditioner<T>::diagInvMatrix_into(velM, m_velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));

        m_B = matNS.block(uSize, 0, pdofs, uSize);
        m_BMinv = m_B * m_velMinv;

        m_alphaP = opt.askReal("alphaP", 1.0);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondMSIMPLER<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt));
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
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
        m_Ainv->apply(input.topRows(uSize) - BtpStar, uStar);

        // pressure correction
        gsMatrix<T> dp;
        m_Sinv->apply(input.bottomRows(pSize) - m_B * uStar, dp);

        gsMatrix<T> tmp2;
        m_Bt->apply(dp, tmp2);

        x.topRows(uSize).noalias() = uStar - m_velMinv * tmp2;
        x.bottomRows(pSize) = pStar + m_alphaP * dp;
    }

    virtual std::string getName()
    {
        std::string name = "MSIMPLER_" + m_Ainv->getName();
        return name;
    }

protected:
    gsSparseMatrix<T> m_velMinv, m_B, m_BMinv;
    real_t m_alphaP;

    using Base::m_Ainv;
    using Base::m_Bt;
    using Base::m_Sinv;
    using Base::m_opt;

}; // class uwbINSBlockPrecondMSIMPLER

// ======================================== Stokes diagonal precond ========================================

/*
Block diagonal preconditioner for Stokes.
*/
template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbBlockPrecondStokes : public uwbINSPreconditioner<T>
{
public:
    /// Shared pointer for uwbBlockPrecondStokes
    typedef memory::shared_ptr<uwbBlockPrecondStokes> Ptr;

    /// Unique pointer for uwbBlockPrecondStokes
    typedef memory::unique_ptr<uwbBlockPrecondStokes> uPtr;

    uwbBlockPrecondStokes(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        m_Ainv(new BlockAType(mat.at("matNS"), opt)), m_Sinv(new uwbINSPrecondSchurStokes<T>(mat, opt)), m_opt(opt)
    {}

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbBlockPrecondStokes<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        m_Ainv.reset(new BlockAType(mat.at("matNS"), m_opt));
    }

    virtual int rows() const
    {
        return (m_Ainv->rows() + m_Sinv->rows());
    }

    virtual int cols() const
    {
        return (m_Ainv->cols() + m_Sinv->cols());
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
        int pSize = m_Sinv->rows();

        x.resize(uSize + pSize, 1);

        gsMatrix<T> tmpX;
        m_Sinv->apply(input.bottomRows(pSize), tmpX);
        x.bottomRows(pSize) = tmpX;

        m_Ainv->apply(input.topRows(uSize), tmpX);
        x.topRows(uSize) = tmpX;
    }

    virtual std::string getName()
    {
        std::string name = "StokesDiag_" + m_Ainv->getName();
        return name;
    }

protected:
    typename uwbINSPrecondBlock<T>::Ptr m_Ainv, m_Sinv;
    const gsOptionList m_opt;

}; // class uwbBlockPrecondStokes

// ======================================== Stokes triangular precond ========================================

template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbBlockPrecondStokesTriang : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbBlockPrecondStokesTriang
    typedef memory::shared_ptr<uwbBlockPrecondStokesTriang> Ptr;

    /// Unique pointer for uwbBlockPrecondStokesTriang
    typedef memory::unique_ptr<uwbBlockPrecondStokesTriang> uPtr;

    uwbBlockPrecondStokesTriang(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurStokes<T>(mat, opt), opt)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbBlockPrecondStokesTriang<T, BlockAType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt));
    }

    virtual std::string getName()
    {
        std::string name = "StokesTriang_" + m_Ainv->getName();
        return name;
    }

protected:

    using Base::m_Ainv;
    using Base::m_opt;

}; // class uwbBlockPrecondStokesTriang


// ======================================== SchurEx precond ========================================

/*
A block triangular preconditioner using the exact Schur complement, solving the systems with S using GMRES preconditioned by S approximation.
*/
/*template <class T, class BlockAType = uwbINSPrecondBlockA<T>, class PrecSType = uwbINSPrecondSchurLSC<T> >
class uwbINSBlockPrecondSchurEx : public uwbINSBlockPrecondBase<T>
{
public:
    typedef uwbINSBlockPrecondBase<T> Base;

    /// Shared pointer for uwbINSBlockPrecondSchurEx
    typedef memory::shared_ptr<uwbINSBlockPrecondSchurEx> Ptr;

    /// Unique pointer for uwbINSBlockPrecondSchurEx
    typedef memory::unique_ptr<uwbINSBlockPrecondSchurEx> uPtr;

    uwbINSBlockPrecondSchurEx(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(new BlockAType(mat.at("matNS"), opt),
            new uwbINSPrecondBlockBt<T>(mat.at("matNS"), opt),
            new uwbINSPrecondSchurExact<T>(mat, opt), opt),
        m_precS(new PrecSType(mat, opt))
    {
        m_maxItS = opt.askInt("maxItS", 20);
        m_tolS = opt.askReal("tolS", 1e-2);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSBlockPrecondSchurEx<T, BlockAType, PrecSType>(mat, opt));
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        Base::update(new BlockAType(mat.at("matNS"), m_opt), new uwbINSPrecondSchurExact<T>(mat, m_opt));
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        int uSize = m_Ainv->rows();
        int pSize = m_Sinv->rows();

        x.resize(uSize + pSize, 1);

        gsMatrix<T> tmpX;
        gsGMRes<T> solverS(m_Sinv, m_precS);
        solverS.setMaxIterations(m_maxItS);
        solverS.setTolerance(m_tolS);
        solverS.solve(input.bottomRows(pSize), tmpX);
        x.bottomRows(pSize) = tmpX;

        gsInfo << "PrecS iterations: " << solverS.iterations() << "\n";

        gsMatrix<T> tmp;
        m_Bt->apply(x.bottomRows(pSize), tmp);
        tmp = input.topRows(uSize) - tmp;

        m_Ainv->apply(tmp, tmpX);
        x.topRows(uSize) = tmpX;
    }

    virtual std::string getName()
    {
        std::string name = "SchurEx_" + m_Ainv->getName();
        return name;
    }

protected:
    typename uwbINSPrecondBlock<T>::Ptr m_precS;
    int m_maxItS;
    real_t m_tolS;

    using Base::m_Ainv;
    using Base::m_Bt;
    using Base::m_Sinv;
    using Base::m_opt;


};*/ // class uwbINSBlockPrecondSchurEx

// ======================================== Base class ========================================

template <class T>
class uwbINSPreconditioner : public gsLinearOperator<T>
{
public:
    /// Shared pointer for uwbINSPreconditioner
    typedef memory::shared_ptr<uwbINSPreconditioner> Ptr;

    /// Unique pointer for uwbINSPreconditioner
    typedef memory::unique_ptr<uwbINSPreconditioner> uPtr;

    uwbINSPreconditioner() {}

public:
    virtual int rows() const
    { 
        GISMO_NO_IMPLEMENTATION
    }

    virtual int cols() const
    { 
        GISMO_NO_IMPLEMENTATION 
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual void update(const std::map<std::string, gsSparseMatrix<T> >& mat)
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual std::string getName()
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Make function returning a smart pointer
    static uPtr make()
    {
        return memory::make_unique(new uwbINSPreconditioner<T>());
    }

    // mat is a vector of required matrices
    //     assuming the following order: NS system matrix, mass matrix (velocity, pressure or both), other matrices

    static uPtr make(std::string precType, const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        if (precType == "LSC_AdiagEqual")
            return uwbINSBlockPrecondLSC<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "LSC_Adiag")
            return uwbINSBlockPrecondLSC<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "LSC_Awhole")
            return uwbINSBlockPrecondLSC<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "LSC_Amod")
            return uwbINSBlockPrecondLSC<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "PCD_AdiagEqual")
            return uwbINSBlockPrecondPCD<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "PCD_Adiag")
            return uwbINSBlockPrecondPCD<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "PCD_Awhole")
            return uwbINSBlockPrecondPCD<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "PCD_Amod")
            return uwbINSBlockPrecondPCD<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "PCDmod_AdiagEqual")
            return uwbINSBlockPrecondPCDmod<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "PCDmod_Adiag")
            return uwbINSBlockPrecondPCDmod<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "PCDmod_Awhole")
            return uwbINSBlockPrecondPCDmod<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "PCDmod_Amod")
            return uwbINSBlockPrecondPCDmod<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "AL_Awhole")
            return uwbINSBlockPrecondAL<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "AL_Amod")
            return uwbINSBlockPrecondAL<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "SIMPLE_AdiagEqual")
            return uwbINSBlockPrecondSIMPLE<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "SIMPLE_Adiag")
            return uwbINSBlockPrecondSIMPLE<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "SIMPLE_Awhole")
            return uwbINSBlockPrecondSIMPLE<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "SIMPLE_Amod")
            return uwbINSBlockPrecondSIMPLE<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "SIMPLER_AdiagEqual")
            return uwbINSBlockPrecondSIMPLER<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "SIMPLER_Adiag")
            return uwbINSBlockPrecondSIMPLER<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "SIMPLER_Awhole")
            return uwbINSBlockPrecondSIMPLER<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "SIMPLER_Amod")
            return uwbINSBlockPrecondSIMPLER<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "MSIMPLER_AdiagEqual")
            return uwbINSBlockPrecondMSIMPLER<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "MSIMPLER_Adiag")
            return uwbINSBlockPrecondMSIMPLER<T, uwbINSPrecondBlockAdiag<T> >::make(mat, opt);
        else if (precType == "MSIMPLER_Awhole")
            return uwbINSBlockPrecondMSIMPLER<T, uwbINSPrecondBlockAwhole<T> >::make(mat, opt);
        else if (precType == "MSIMPLER_Amod")
            return uwbINSBlockPrecondMSIMPLER<T, uwbINSPrecondBlockAmod<T> >::make(mat, opt);
        else if (precType == "StokesDiag_AdiagEqual")
            return uwbBlockPrecondStokes<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        else if (precType == "StokesTriang_AdiagEqual")
            return uwbBlockPrecondStokesTriang<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        /*else if (precType == "SchurEx_LSC_AdiagEqual")
            return uwbINSBlockPrecondSchurEx<T, uwbINSPrecondBlockA<T>, uwbINSPrecondSchurLSC<T> >::make(mat, opt);
        else if (precType == "SchurEx_SIMPLE_AdiagEqual")
            return uwbINSBlockPrecondSchurEx<T, uwbINSPrecondBlockA<T>, uwbINSPrecondSchurSIMPLE<T> >::make(mat, opt);
        else if (precType == "SchurEx_MSIMPLER_AdiagEqual")
            return uwbINSBlockPrecondSchurEx<T, uwbINSPrecondBlockA<T>, uwbINSPrecondSchurMSIMPLER<T> >::make(mat, opt);*/
        else
        {
            gsInfo << "Invalid preconditioner type, using LSC_AdiagEqual.\n";
            return uwbINSBlockPrecondLSC<T, uwbINSPrecondBlockA<T> >::make(mat, opt);
        }
    }

    // default preconditioner options in a gsOptionList object
    static gsOptionList defaultOptions()
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

    // repeat - number of the diagonal block repetition (e.g. for velocity components)
    static void diagInvMatrix_into(const gsSparseMatrix<T>& mat, gsSparseMatrix<T>& diagInv, int repeat, bool lumping = false)
    {
        GISMO_ENSURE(mat.nonZeros() != 0, "diagInvMatrix_into(): The matrix is empty!");

        int varDofs = mat.rows();

        diagInv.resize(repeat*varDofs, repeat*varDofs);
        diagInv.reserve(gsVector<int>::Constant(diagInv.cols(), 1));

        const gsSparseMatrix<T>* matPtr = &mat;
        gsSparseMatrix<T> lumped(varDofs, varDofs);

        if (lumping)
        {
            lumped.reserve(gsVector<int>::Constant(varDofs, 1));

            for (int j = 0; j < varDofs; j++)
                lumped.insert(j, j) = math::abs(mat.at(j, j)); // abs value because of "lumping" diag(A) in SIMPLE-type prec., does not change lumped mass matrix in IgA

            for (int j = 0; j < varDofs; j++)
            {
                for (typename gsSparseMatrix<T>::InnerIterator it(mat, j); it; ++it)
                {
                    int i = it.row();

                    if (i != j)
                        lumped.coeffRef(i, i) += math::abs(it.value());
                }
            }

            matPtr = &lumped;
        }

        for (int i = 0; i < varDofs; i++)
        {
            T tmp = 1 / matPtr->coeff(i, i);

            for (int s = 0; s < repeat; s++)
                diagInv.coeffRef(i + s * varDofs, i + s * varDofs) = tmp;
        }
    }

}; // class uwbINSPreconditioner

} // namespace gismo

//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(uwbPreconditioners.hpp)
//#endif
