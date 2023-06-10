/** @file gsINSPrecondBlocks.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova
*/

#pragma once

#include <gsIncompressibleFlow/src/gsINSUtils.h>

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

// === BASE BLOCKS ==== //

/// @brief A base class for individual blocks of block preconditioners.
/// @tparam T coefficient type
template <class T>
class gsINSPrecondBlock : public gsLinearOperator<T>
{

public:
    typedef gsLinearOperator<T> Base;

protected: // *** Class members ***

    gsSparseSolver<T>* m_pSolver;
    gsOptionList m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlock> Ptr; 
    typedef memory::unique_ptr<gsINSPrecondBlock> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    gsINSPrecondBlock()
    {
        m_pSolver = NULL;
    }

    /// @brief Constructor.
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondBlock(const gsOptionList& opt) : gsINSPrecondBlock()
    {
        m_opt = opt;

        m_pSolver = createLinSolver();
        setupLinSolver(*m_pSolver);
    }

    virtual ~gsINSPrecondBlock()
    {
        if (m_pSolver)
        {
            delete m_pSolver;
            m_pSolver = NULL;
        }
    }

public: // *** Member functions ***

    /// @brief Returns a pointer to new linear solver (direct or iterative).
    gsSparseSolver<T>* createLinSolver();

    /// @brief Set up the linear solver for the block.
    /// @param[out] solver a reference to the solver object
    void setupLinSolver(gsSparseSolver<T>& solver);

    virtual std::string getName() const
    { GISMO_NO_IMPLEMENTATION }

}; // class gsINSPrecondBlock


/** @brief
    Base class for block F of (modified) block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    F & B^T \\
    0 & S
    \end{array}
    \right],
    \f]

    where several systems have to be solved.
*/
template <class T>
class gsINSPrecondBlockMod : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    std::vector<gsSparseSolver<T>* > m_solvers;

protected: // *** Base class members ***

    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockMod> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockMod> uPtr;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondBlockMod(const gsOptionList& opt) : Base()
    {
        m_opt = opt;
        int dim = m_opt.getInt("dim");

        m_solvers.resize(dim);
        
        for (int i = 0; i < dim; i++)
        {
            m_solvers[i] = this->createLinSolver();
            this->setupLinSolver(*m_solvers[i]);
        }
    }

    virtual ~gsINSPrecondBlockMod()
    {
        for (int i = 0; i < (int)m_solvers.size(); i++)
        {
            if (m_solvers[i])
            {
                delete m_solvers[i];
                m_solvers[i] = NULL;
            }
        }
    }

}; // class gsINSPrecondBlockMod


// === BLOCK F ==== //

/** @brief
    Class for block F of block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    F & B^T \\
    0 & S
    \end{array}
    \right],
    \f].

    Implements action of \fF^{-1}\f. Asssumes that F is block-diagonal with identical diagonal subblocks, solves several linear systems with the diagonal subblock.
*/
template <class T>
class gsINSPrecondBlockF : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    const gsSparseMatrix<T> m_blockF1;
    int m_size;

protected: // *** Base class members ***

    using Base::m_pSolver;
    using Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockF> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockF> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    gsINSPrecondBlockF(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt), m_blockF1(matNS.block(0, 0, opt.getInt("udofs"), opt.getInt("udofs")))
    {
        m_size = opt.getInt("dim") * opt.getInt("udofs");

        m_pSolver->compute(m_blockF1);
    }
    
public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondBlockF<T>(matNS, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = F^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;


public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_size; }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_size; }

    /// @brief Returns the block name as a string.
    virtual std::string getName() const
    {
        return "FdiagEqual";
    }

}; // class gsINSPrecondBlockF


// === BLOCK Fwhole ==== //

/** @brief
Class for block F of block preconditioner of the form
\f[
\left[ \begin{array}{cc}
F & B^T \\
0 & S
\end{array}
\right],
\f].

Implements action of \fF^{-1}\f, solves linear system with the whole block F.
*/
template <class T>
class gsINSPrecondBlockFwhole : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    const gsSparseMatrix<T> m_blockA;
    int m_size;

protected: // *** Base class members ***

    using Base::m_pSolver;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockFwhole> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockFwhole> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    gsINSPrecondBlockFwhole(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt), m_blockA(matNS.block(0, 0, opt.getInt("dim")*opt.getInt("udofs"), opt.getInt("dim")*opt.getInt("udofs")))
    {
        m_size = opt.getInt("dim") * opt.getInt("udofs");

        m_pSolver->compute(m_blockA);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondBlockFwhole<T>(matNS, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = F^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_size; }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_size; }

    /// @brief Returns the block name as a string.
    virtual std::string getName() const
    {
        return "Fwhole";
    }

}; // class gsINSPrecondBlockFwhole


// === BLOCK Fdiag ==== //

/** @brief
Base class for block F of block preconditioner of the form
\f[
\left[ \begin{array}{cc}
F & B^T \\
0 & S
\end{array}
\right],
\f].

Implements action of \fF^{-1}\f. Neglects the off-diagonal blocks, solves subsystems with the diagonal blocks and assumes that they are not equal.
*/
template <class T>
class gsINSPrecondBlockFdiag : public gsINSPrecondBlockMod<T>
{

public:
    typedef gsINSPrecondBlockMod<T> Base;

protected: // *** Class members ***
    std::vector<gsSparseMatrix<T> > m_diagBlocks;
    int m_size;

protected: // *** Base class members ***

    using Base::m_solvers;
    using Base::Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockFdiag> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockFdiag> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    gsINSPrecondBlockFdiag(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt)
    {
        int dim = opt.getInt("dim");
        m_size = dim * opt.getInt("udofs");

        m_diagBlocks.resize(dim);

        for (int i = 0; i < dim; i++)
        {
            int udofs = opt.getInt("udofs");
            m_diagBlocks[i] = matNS.block(i*udofs, i*udofs, udofs, udofs);
            m_solvers[i]->compute(m_diagBlocks[i]);
        }
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondBlockFdiag<T>(matNS, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = F^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_size; }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_size; }

    /// @brief Returns the block name as a string.
    virtual std::string getName() const
    {
        return "Fdiag";
    }

}; // class gsINSPrecondBlockFdiag


// === BLOCK Fmod ==== //

/** @brief
    Base class for block F of (modified) block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    F & B^T \\
    0 & S
    \end{array}
    \right],
    \f].

    Implements action of \fF^{-1}\f. Assumes that the block F is not block-diagonal, solves a system with lower block-triangular submatrix of F (neglects the upper triangle).
*/
template <class T>
class gsINSPrecondBlockFmod : public gsINSPrecondBlockMod<T>
{

public:
    typedef gsINSPrecondBlockMod<T> Base;

protected: // *** Class members ***

    const gsSparseMatrix<T>& m_matRef;
    std::vector<gsSparseMatrix<T> > m_diagBlocks;
    int m_size;

protected: // *** Base class members ***

    using Base::m_solvers;
    using Base::Base::m_opt;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockFmod> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockFmod> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    gsINSPrecondBlockFmod(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt), m_matRef(matNS)
    {
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        m_size = dim * udofs;

       m_diagBlocks.resize(dim);

        for (int i = 0; i < dim; i++)
        {
            m_diagBlocks[i] = m_matRef.block(i*udofs, i*udofs, udofs, udofs);
            m_solvers[i]->compute(m_diagBlocks[i]);
        }
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondBlockFmod<T>(matNS, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = F^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_size; }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_size; }

    /// @brief Returns the block name as a string.
    virtual std::string getName() const
    {
        return "Fmod";
    }

}; // class gsINSPrecondBlockFmod


// === BLOCK Bt ==== //

template <class T>
class gsINSPrecondBlockBt : public gsLinearOperator<T>
{

public:
    typedef gsLinearOperator<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_mat;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondBlockBt> Ptr;
    typedef memory::unique_ptr<gsINSPrecondBlockBt> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    gsINSPrecondBlockBt(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        m_mat = matNS.block(0, dim * udofs, dim * udofs, matNS.cols() - dim * udofs);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] matNS a const reference to the saddle-point system matrix
    /// @param[in] opt   a list of options for the preconditioner
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondBlockBt<T>(matNS, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = B^T y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat.cols(); }

}; // class gsINSPrecondBlockBt


// === Schur LSC ==== //

template <class T>
class gsINSPrecondSchurLSC : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_mat1, m_mat2;

protected: // *** Base class members ***

    using Base::m_pSolver;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurLSC> Ptr;
    typedef memory::unique_ptr<gsINSPrecondSchurLSC> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurLSC(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(opt)
    {
        // S = - (B * velM^-1 * BT) * (B * velM^-1 * F * velM^-1 * BT)^-1 * (B * velM^-1 * BT)
        // denote S = - mat1 * mat2^-1 * mat1

        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int uSize = dim * udofs; // number of all velocity dofs (all components)
        int pdofs = opt.getInt("pdofs"); // number of pressure dofs

        const gsSparseMatrix<T>& matNS = mat.at("matNS");
        const gsSparseMatrix<T>& velM = mat.at("matMu");

        gsSparseMatrix<T> velMinv(uSize, uSize); // approximation of velocity mass matrix inverse

        diagInvMatrix_into(velM, velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));
        
        // TODO: could be more efficient if the matrices not involving block F were not updated (they do not change)
        gsSparseMatrix<T> BMinv, MinvBt;
        BMinv = matNS.block(uSize, 0, pdofs, uSize) * velMinv; // B * velM^-1
        MinvBt = velMinv * matNS.block(0, uSize, uSize, pdofs); // velM^-1 * BT

        m_mat1 = BMinv * matNS.block(0, uSize, uSize, pdofs); // B * velM^-1 * BT
        m_mat2 = BMinv * matNS.block(0, 0, uSize, uSize) * MinvBt; // B * velM^-1 * F * velM^-1 * BT

        m_pSolver->compute(m_mat1);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurLSC<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat2.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat2.cols(); }

}; // class gsINSPrecondSchurLSC


// === Schur PCD ==== //

template <class T>
class gsINSPrecondSchurPCD : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_matAp, m_matFp, m_matMp;
    gsSparseSolver<T>* m_pMassSolver;

protected: // *** Base class members ***

    using Base::m_pSolver;
    
public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurPCD> Ptr; 
    typedef memory::unique_ptr<gsINSPrecondSchurPCD> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurPCD(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(opt)
    {
        // S = - Ap * Fp^-1 * Mp

        m_matAp = mat.at("matAp");
        m_matFp = mat.at("matFp");
        m_matMp = mat.at("matMp");
        
        m_pMassSolver = this->createLinSolver();
        this->setupLinSolver(*m_pMassSolver);

        m_pSolver->compute(m_matAp);
        m_pMassSolver->compute(m_matMp);
    }

    ~gsINSPrecondSchurPCD()
    {
        if (m_pMassSolver)
        {
            delete m_pMassSolver;
            m_pMassSolver = NULL;
        }
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurPCD<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const;

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_matAp.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_matAp.cols(); }

}; // class gsINSPrecondSchurPCD


// === Schur PCDmod ==== // 

template <class T>
class gsINSPrecondSchurPCDmod : public gsINSPrecondSchurPCD<T>
{

public:
    typedef gsINSPrecondSchurPCD<T> Base;

protected: // *** Base class members ***

    using Base::m_matAp;
    using Base::m_matFp;
    using Base::m_matMp;
    using Base::m_pMassSolver;
    using Base::Base::m_pSolver;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurPCDmod> Ptr;
    typedef memory::unique_ptr<gsINSPrecondSchurPCDmod> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurPCDmod(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(mat, opt)
    {
        // S = - Mp * Fp^-1 * Ap
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurPCDmod<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const;

}; // class gsINSPrecondSchurPCDmod


// === Schur AL ==== // 

template <class T>
class gsINSPrecondSchurAL : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    real_t m_const;
    gsSparseMatrix<T> m_mat;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurAL> Ptr;  
    typedef memory::unique_ptr<gsINSPrecondSchurAL> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurAL(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        // S = - (gamma + nu) Mp

        const gsSparseMatrix<T>& presM = mat.at("matMp");

        m_const = opt.getReal("visc") + opt.getReal("gamma");

        diagInvMatrix_into(presM, m_mat, 1, opt.getSwitch("lumpingM"));
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurAL<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = - m_const * m_mat * input;
    }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat.cols(); }

}; // class gsINSPrecondSchurAL


// === Schur SIMPLE ==== // 

template <class T>
class gsINSPrecondSchurSIMPLE : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_mat;

protected: // *** Base class members ***

    using Base::m_pSolver;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurSIMPLE> Ptr;  
    typedef memory::unique_ptr<gsINSPrecondSchurSIMPLE> uPtr;


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurSIMPLE(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(opt)
    {
        // S = - B Dinv Bt
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int uSize = dim * udofs;
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T>& matNS = mat.at("matNS");

        gsSparseMatrix<T> blockA = matNS.block(0, 0, uSize, uSize);
        gsSparseMatrix<T> Dinv(uSize, uSize); // approximation of A inverse
        diagInvMatrix_into(blockA, Dinv, 1, opt.getSwitch("lumpingA"));

        m_mat = - matNS.block(uSize, 0, pdofs, uSize) * Dinv * matNS.block(0, uSize, uSize, pdofs);

        m_pSolver->compute(m_mat);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurSIMPLE<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = m_pSolver->solve(input);
    }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat.cols(); }

}; // class gsINSPrecondSchurSIMPLE


// === Schur MSIMPLER ==== // 

template <class T>
class gsINSPrecondSchurMSIMPLER : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_mat;

protected: // *** Base class members ***

    using Base::m_pSolver;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurMSIMPLER> Ptr;
    typedef memory::unique_ptr<gsINSPrecondSchurMSIMPLER> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurMSIMPLER(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(opt)
    {
        // S = - B velMinv Bt
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        int uSize = dim * udofs;
        int pdofs = opt.getInt("pdofs");

        const gsSparseMatrix<T>& matNS = mat.at("matNS");
        const gsSparseMatrix<T>& velM = mat.at("matMu");

        gsSparseMatrix<T> velMinv(uSize, uSize); // approximation of velocity mass matrix inverse
        diagInvMatrix_into(velM, velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));

        m_mat = - matNS.block(uSize, 0, pdofs, uSize) * velMinv * matNS.block(0, uSize, uSize, pdofs);

        m_pSolver->compute(m_mat);
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurMSIMPLER<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = m_pSolver->solve(input);
    }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat.cols(); }

}; // class gsINSPrecondSchurMSIMPLER


// === Schur Stokes ==== // 

template <class T>
class gsINSPrecondSchurStokes : public gsINSPrecondBlock<T>
{

public:
    typedef gsINSPrecondBlock<T> Base;

protected: // *** Class members ***

    gsSparseMatrix<T> m_mat;
    T m_visc;

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsINSPrecondSchurStokes> Ptr;
    typedef memory::unique_ptr<gsINSPrecondSchurStokes> uPtr;

public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    gsINSPrecondSchurStokes(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        // S = diag(Mp)

        m_visc = opt.getReal("visc");

        const gsSparseMatrix<T>& presM = mat.at("matMp");
        diagInvMatrix_into(presM, m_mat, 1, opt.getSwitch("lumpingM"));
    }

public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the block
    /// @param[in] opt a list of options for the preconditioner
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new gsINSPrecondSchurStokes<T>(mat, opt));
    }

public: // *** Member functions ***

    /// @brief Apply the block.
    /// Computes the vector \fx = S^{-1} y\f.
    /// @param[in]  input   a const reference the vector \f y \f
    /// @param[out] x       a reference to the resulting vector \f x \f
    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = - m_visc * m_mat * input;
    }

public: // *** Getters/setters ***

    /// @brief Returns the number of rows of the block.
    int rows() const { return m_mat.rows(); }

    /// @brief Returns the number of columns of the block.
    int cols() const { return m_mat.cols(); }

}; // class gsINSPrecondSchurStokes


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSPrecondBlocks.hpp)
#endif