/** @file gsINSPreconditioners.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/
#pragma once

#include <gsIncompressibleFlow/src/gsINSPrecondBlocks.h>
#include <gsIncompressibleFlow/src/gsFlowUtils.h>

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

// === Base class for INS preconditioners === //

template <class T, int MatOrder>
class gsINSPreconditioner : public gsLinearOperator<T>
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPreconditioner> Ptr;
    typedef memory::unique_ptr<gsINSPreconditioner> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSPreconditioner() {}

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    static uPtr make()
    { return memory::make_unique(new gsINSPreconditioner<T, MatOrder>()); }

    /// @brief Returns a unique pointer to a newly created instance of the given preconditioner type.
    /// @param[in] precType the reqiured preconditioner type as a string
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner (assuming the following order: NS system matrix, mass matrix (velocity, pressure or both), other matrices)
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(std::string precType, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt);

    /// @brief Returns default preconditioner options as a gsOptionList object.
    static gsOptionList defaultOptions();


public: // *** Member functions ***

    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    { GISMO_NO_IMPLEMENTATION }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the preconditioner.
    virtual int rows() const
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Returns the number of columns of the preconditioner.
    virtual int cols() const
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    { GISMO_NO_IMPLEMENTATION }

}; // class gsINSPreconditioner


// === Block preconditioner === //

/** @brief
Base class for block preconditioners for linear systems arising from linearized Navier-Stokes of form
\f[
\left[ \begin{array}{cc}
F & B^T \\
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
template <class T, int MatOrder>
class gsINSBlockPrecondBase : public gsINSPreconditioner<T, MatOrder>
{

protected: // *** Class members ***

    typename gsINSPrecondBlock<T>::Ptr m_Finv, m_Sinv;
    typename gsLinearOperator<T>::Ptr m_Bt;
    const gsOptionList m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondBase> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondBase> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] Finv a pointer to the block implementing multiplication by the matrix \fF^{-1}\f
    /// @param[in] B    a pointer to the block implementing multiplication by the matrix \fB^T\f
    /// @param[in] Sinv a pointer to the block implementing multiplication by the matrix \fS^{-1}\f
    /// @param[in] opt  a list of options for the preconditioner
    gsINSBlockPrecondBase(gsINSPrecondBlock<T>* Finv, gsLinearOperator<T>* B, gsINSPrecondBlock<T>* Sinv, const gsOptionList& opt) :
        m_Finv(Finv), m_Sinv(Sinv), m_Bt(B), m_opt(opt)
    {}

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    static uPtr make(gsINSPrecondBlock<T>* Finv, gsLinearOperator<T>* B, gsINSPrecondBlock<T>* Sinv, const gsOptionList& opt)
    { return memory::make_unique(new gsINSBlockPrecondBase<T, MatOrder>(Finv, B, Sinv, opt)); }

public: // *** Member functions ***

    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

protected: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// To be used when only the \fF\f block changes.
    /// @param[in] Finv a pointer to the new \fF\f block
    void update(gsINSPrecondBlock<T>* Finv)
    {
        m_Finv.reset(Finv);
    }

    /// @brief Update the preconditioner (new linearization or time step).
    /// To be used when both \fF\f and \fS\f block change.
    /// @param[in] Finv a pointer to the new \fF\f block
    /// @param[in] Sinv a pointer to the new \fS\f block
    void update(gsINSPrecondBlock<T>* Finv, gsINSPrecondBlock<T>* Sinv)
    {
        m_Finv.reset(Finv);
        m_Sinv.reset(Sinv);
    }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the preconditioner.
    virtual int rows() const
    {
        return (m_Finv->rows() + m_Sinv->rows());
    }

    /// @brief Returns the number of columns of the preconditioner.
    virtual int cols() const
    {
        return (m_Finv->cols() + m_Bt->cols());
    }

}; // class gsINSBlockPrecondBase


// === LSC preconditioner === //

/// @brief Least-squares commutator preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondLSC : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Class members ***

    using Base::m_Finv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondLSC> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondLSC> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondLSC(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
             new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
             new gsINSPrecondSchurLSC<T, MatOrder>(mat, opt), opt)
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondLSC<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt), new gsINSPrecondSchurLSC<T, MatOrder>(mat, m_opt));
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "LSC_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondLSC


// === PCD preconditioner === //

/// @brief Pressure convection-diffusion preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondPCD : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Class members ***

    using Base::m_Finv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondPCD> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondPCD> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondPCD(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurPCD<T, MatOrder>(mat, opt), opt)
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondPCD<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt), new gsINSPrecondSchurPCD<T, MatOrder>(mat, m_opt));
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "PCD_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondPCD


// === PCDmod preconditioner === //

/// @brief Modified pressure convection-diffusion preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondPCDmod : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Class members ***

    using Base::m_Finv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondPCDmod> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondPCDmod> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondPCDmod(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurPCDmod<T, MatOrder>(mat, opt), opt)
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondPCDmod<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt), new gsINSPrecondSchurPCDmod<T, MatOrder>(mat, m_opt));
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "PCDmod_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondPCDmod


// === AL preconditioner === //

/// @brief Augmented Lagrangian preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockFwhole<T, MatOrder> >
class gsINSBlockPrecondAL : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

 protected: // *** Class members ***

    using Base::m_Finv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondAL> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondAL> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondAL(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurAL<T, MatOrder>(mat, opt), opt)
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondAL<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Static functions ***

    /// @brief Fill the extra part of the Augmented Lagrangian linear system.
    /// @param[out] matGammaPart    a reference to the resulting matrix
    /// @param[out] rhsGammaPart    a reference to the resulting right-hand side
    /// @param[in]  mat         a const reference to std::map of labeled blocks of the original system
    /// @param[in]  rhs         the original right-hand side
    /// @param[in]  opt         a list of options for the preconditioner
    static void fillALgammaPart_into(gsSparseMatrix<T, MatOrder>& matGammaPart, gsMatrix<T>& rhsGammaPart, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsMatrix<T>& rhs, const gsOptionList& opt);

    /// @brief Fill the Augmented Lagrangian linear system.
    /// @param[out] matGamma    a reference to the resulting matrix
    /// @param[out] rhsGamma    a reference to the resulting right-hand side
    /// @param[in]  mat         a const reference to std::map of labeled blocks of the original system
    /// @param[in]  rhs         the original right-hand side
    /// @param[in]  opt         a list of options for the preconditioner
    static void fillALmodifSystem_into(gsSparseMatrix<T, MatOrder>& matGamma, gsMatrix<T>& rhsGamma, const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsMatrix<T>& rhs, const gsOptionList& opt);

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt));
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "AL_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondAL


// === SIMPLE preconditioner === //

/// @brief SIMPLE preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondSIMPLE : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Class members ***

    gsSparseMatrix<T, MatOrder> m_Dinv, m_B;
    real_t m_alphaP;

protected: // *** Base class members ***

    using Base::m_Finv;
    using Base::m_Bt;
    using Base::m_Sinv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondSIMPLE> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondSIMPLE> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondSIMPLE(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurSIMPLE<T, MatOrder>(mat, opt), opt)
    {
        int uSize = opt.getInt("dim") * opt.getInt("udofs");
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T, MatOrder>& matNS = mat.at("matNS");

        gsSparseMatrix<T, MatOrder> blockA = matNS.block(0, 0, uSize, uSize);
        diagInvMatrix_into(blockA, m_Dinv, 1, opt.getSwitch("lumpingA"));

        m_B = matNS.block(uSize, 0, pdofs, uSize);

        m_alphaP = opt.askReal("alphaP", 1.0);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondSIMPLE<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt), new gsINSPrecondSchurSIMPLE<T, MatOrder>(mat, m_opt));
    }

    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "SIMPLE_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondSIMPLE


// === SIMPLER preconditioner === //

/// @brief SIMPLER preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondSIMPLER : public gsINSBlockPrecondSIMPLE<T, MatOrder, BlockFType>
{

public:
    typedef gsINSBlockPrecondSIMPLE<T, MatOrder, BlockFType> Base;

protected: // *** Class members ***

    gsSparseMatrix<T, MatOrder> m_BDinv;

protected: // *** Base class members ***

    using Base::m_Dinv;
    using Base::m_B;
    using Base::m_alphaP;
    using gsINSBlockPrecondBase<T, MatOrder>::m_Finv;
    using gsINSBlockPrecondBase<T, MatOrder>::m_Bt;
    using gsINSBlockPrecondBase<T, MatOrder>::m_Sinv;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondSIMPLER> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondSIMPLER> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondSIMPLER(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(mat, opt)
    {
        m_BDinv = m_B * m_Dinv;
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondSIMPLER<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "SIMPLER_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondSIMPLER


// === MSIMPLER preconditioner === //

/// @brief MSIMPLER preconditioner.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsINSBlockPrecondMSIMPLER : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Class members ***

    gsSparseMatrix<T, MatOrder> m_velMinv, m_B, m_BMinv;
    real_t m_alphaP;

protected: // *** Base class members ***

    using Base::m_Finv;
    using Base::m_Bt;
    using Base::m_Sinv;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSBlockPrecondMSIMPLER> Ptr;
    typedef memory::unique_ptr<gsINSBlockPrecondMSIMPLER> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsINSBlockPrecondMSIMPLER(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurMSIMPLER<T, MatOrder>(mat, opt), opt)
    {
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int uSize = dim * udofs;
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T, MatOrder>& matNS = mat.at("matNS");
        const gsSparseMatrix<T, MatOrder>& velM = mat.at("matMu");

       // approximation of velocity mass matrix inverse
        diagInvMatrix_into(velM, m_velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));

        m_B = matNS.block(uSize, 0, pdofs, uSize);
        m_BMinv = m_B * m_velMinv;

        m_alphaP = opt.askReal("alphaP", 1.0);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSBlockPrecondMSIMPLER<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt));
    }

    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "MSIMPLER_" + m_Finv->getName();
        return name;
    }

}; // class gsINSBlockPrecondMSIMPLER


// === Stokes diagonal preconditioner === //

/// @brief Block diagonal preconditioner for the Stokes problem.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsBlockPrecondStokes : public gsINSPreconditioner<T, MatOrder>
{

protected: // *** Class members ***

    typename gsINSPrecondBlock<T>::Ptr m_Finv, m_Sinv;
    const gsOptionList m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsBlockPrecondStokes> Ptr;
    typedef memory::unique_ptr<gsBlockPrecondStokes> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsBlockPrecondStokes(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        m_Finv(new BlockFType(mat.at("matNS"), opt)), m_Sinv(new gsINSPrecondSchurStokes<T, MatOrder>(mat, opt)), m_opt(opt)
    {}

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsBlockPrecondStokes<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        m_Finv.reset(new BlockFType(mat.at("matNS"), m_opt));
    }


    /// @brief Apply the preconditioner.
    /// Computes the vector \fx = P^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the preconditioner.
    virtual int rows() const
    {
        return (m_Finv->rows() + m_Sinv->rows());
    }

    /// @brief Returns the number of columns of the preconditioner.
    virtual int cols() const
    {
        return (m_Finv->cols() + m_Sinv->cols());
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "StokesDiag_" + m_Finv->getName();
        return name;
    }

}; // class gsBlockPrecondStokes


// === Stokes triangular preconditioner === //

/// @brief Block triangular preconditioner for the Stokes problem.
/// @tparam T           coefficient type
/// @tparam BlockFType  type of block \fF\f
template <class T, int MatOrder, class BlockFType = gsINSPrecondBlockF<T, MatOrder> >
class gsBlockPrecondStokesTriang : public gsINSBlockPrecondBase<T, MatOrder>
{

public:
    typedef gsINSBlockPrecondBase<T, MatOrder> Base;

protected: // *** Base class members ***

    using Base::m_Finv;
    using Base::m_opt;

public: // *** Smart pointers ***
    
    typedef memory::shared_ptr<gsBlockPrecondStokesTriang> Ptr;
    typedef memory::unique_ptr<gsBlockPrecondStokesTriang> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a reference to the std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    gsBlockPrecondStokesTriang(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt) :
        Base(new BlockFType(mat.at("matNS"), opt),
            new gsINSPrecondBlockBt<T, MatOrder>(mat.at("matNS"), opt),
            new gsINSPrecondSchurStokes<T, MatOrder>(mat, opt), opt)
    { }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsBlockPrecondStokesTriang<T, MatOrder, BlockFType>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Update the preconditioner (new linearization or time step).
    /// @param[in] mat a const reference to std::map of updated matrices
    virtual void update(const std::map<std::string, gsSparseMatrix<T, MatOrder> >& mat)
    {
        Base::update(new BlockFType(mat.at("matNS"), m_opt));
    }

    /// @brief Returns the preconditioner name as a string.
    virtual std::string getName()
    {
        std::string name = "StokesTriang_" + m_Finv->getName();
        return name;
    }

}; // class gsBlockPrecondStokesTriang

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSPreconditioners.hpp)
#endif
