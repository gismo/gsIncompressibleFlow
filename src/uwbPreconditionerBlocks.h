/** @file uwbPreconditionerBlocks.h

Author(s): H. Hornikova
*/
#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

// forward declaration
template <class T>
class uwbINSPreconditioner;

// ======================================== base blocks ========================================

/** @brief
    Base class for a block of a block preconditioner, where a single linear system needs to be solved.
*/
template <class T>
class uwbINSPrecondBlock : public gsLinearOperator<T>
{
public:
    typedef gsLinearOperator<T> Base;

    /// Shared pointer for uwbINSPrecondBlock
    typedef memory::shared_ptr<uwbINSPrecondBlock> Ptr;

    /// Unique pointer for uwbINSPrecondBlock   
    typedef memory::unique_ptr<uwbINSPrecondBlock> uPtr;

    uwbINSPrecondBlock()
    {
        m_pSolver = NULL;
    }

    uwbINSPrecondBlock(const gsOptionList& opt) : uwbINSPrecondBlock()
    {
        m_opt = opt;

        m_pSolver = createLinSolver();
        setupLinSolver(*m_pSolver);
    }

    virtual ~uwbINSPrecondBlock()
    {
        if (m_pSolver)
        {
            delete m_pSolver;
            m_pSolver = NULL;
        }
    }

    gsSparseSolver<T>* createLinSolver()
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

    void setupLinSolver(gsSparseSolver<T>& solver)
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
                uwbINSSolverBase<T>::pardisoSetup(*pSolver);
                #endif
	    }
    }

    virtual std::string getName() const
    { GISMO_NO_IMPLEMENTATION }

protected:
    gsSparseSolver<T>* m_pSolver;
    gsOptionList m_opt;

}; // class uwbINSPrecondBlock

/** @brief
    Base class for block A of (modified) block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    A & B^T \\
    0 & S
    \end{array}
    \right],
    \f]

    where several systems have to be solved.
*/
template <class T>
class uwbINSPrecondBlockMod : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondBlockMod
    typedef memory::shared_ptr<uwbINSPrecondBlockMod> Ptr;

    /// Unique pointer for uwbINSPrecondBlockMod
    typedef memory::unique_ptr<uwbINSPrecondBlockMod> uPtr;

    uwbINSPrecondBlockMod(const gsOptionList& opt) : Base()
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

    virtual ~uwbINSPrecondBlockMod()
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

protected:
    std::vector<gsSparseSolver<T>* > m_solvers;

    using Base::m_opt;

}; // class uwbINSPrecondBlockMod

// ======================================== blockA ========================================

/** @brief
    Class for block A of block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    A & B^T \\
    0 & S
    \end{array}
    \right],
    \f].

    Implements action of \fA^{-1}\f. Asssumes that A is block-diagonal with identical diagonal subblocks, solves several linear systems with the diagonal subblock.
*/
template <class T>
class uwbINSPrecondBlockA : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondBlockA
    typedef memory::shared_ptr<uwbINSPrecondBlockA> Ptr;

    /// Unique pointer for uwbINSPrecondBlockA
    typedef memory::unique_ptr<uwbINSPrecondBlockA> uPtr;

    uwbINSPrecondBlockA(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt), m_blockA1(matNS.block(0, 0, opt.getInt("udofs"), opt.getInt("udofs")))
    {
        m_size = opt.getInt("dim") * opt.getInt("udofs");

        m_pSolver->compute(m_blockA1);
    }

    /// Make function returning a smart pointer
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondBlockA<T>(matNS, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
        x.resize(m_size, 1);

        int udofs = m_opt.getInt("udofs");
        for (int i = 0; i < m_opt.getInt("dim"); i++)
            x.middleRows(i*udofs, udofs) = m_pSolver->solve(input.middleRows(i*udofs, udofs));
    }

    int rows() const { return m_size; }

    int cols() const { return m_size; }

    virtual std::string getName() const
    {
        return "AdiagEqual";
    }

protected:
    const gsSparseMatrix<T> m_blockA1;
    int m_size;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;
    using Base::m_opt;

}; // class uwbINSPrecondBlockA

// ======================================== blockAwhole ========================================

/** @brief
Class for block A of block preconditioner of the form
\f[
\left[ \begin{array}{cc}
A & B^T \\
0 & S
\end{array}
\right],
\f].

Implements action of \fA^{-1}\f, solves linear system with the whole block A.
*/
template <class T>
class uwbINSPrecondBlockAwhole : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondBlockAwhole
    typedef memory::shared_ptr<uwbINSPrecondBlockAwhole> Ptr;

    /// Unique pointer for uwbINSPrecondBlockAwhole   
    typedef memory::unique_ptr<uwbINSPrecondBlockAwhole> uPtr;

    uwbINSPrecondBlockAwhole(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
        Base(opt), m_blockA(matNS.block(0, 0, opt.getInt("dim")*opt.getInt("udofs"), opt.getInt("dim")*opt.getInt("udofs")))
    {
        m_size = opt.getInt("dim") * opt.getInt("udofs");

        m_pSolver->compute(m_blockA);
    }

    /// Make function returning a smart pointer
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondBlockAwhole<T>(matNS, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
        x.resize(m_size, 1);

        x = m_pSolver->solve(input);
    }

    int rows() const { return m_size; }

    int cols() const { return m_size; }

    virtual std::string getName() const
    {
        return "Awhole";
    }

protected:
    const gsSparseMatrix<T> m_blockA;
    int m_size;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;

}; // class uwbINSPrecondBlockAwhole

// ======================================== blockAdiag ========================================

/** @brief
Base class for block A of block preconditioner of the form
\f[
\left[ \begin{array}{cc}
A & B^T \\
0 & S
\end{array}
\right],
\f].

Implements action of \fA^{-1}\f. Neglects the off-diagonal blocks, solves subsystems with the diagonal blocks and assumes that they are not equal.
*/
template <class T>
class uwbINSPrecondBlockAdiag : public uwbINSPrecondBlockMod<T>
{
public:
    typedef uwbINSPrecondBlockMod<T> Base;

    /// Shared pointer for uwbINSPrecondBlockAdiag
    typedef memory::shared_ptr<uwbINSPrecondBlockAdiag> Ptr;

    /// Unique pointer for uwbINSPrecondBlockAdiag
    typedef memory::unique_ptr<uwbINSPrecondBlockAdiag> uPtr;

    uwbINSPrecondBlockAdiag(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
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

    /// Make function returning a smart pointer
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondBlockAdiag<T>(matNS, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == m_size, "Wrong input size.");
        x.resize(m_size, 1);

        int udofs = m_opt.getInt("udofs");

        // solving system with lower block-triangular part of Agamma
        for (int i = 0; i < m_opt.getInt("dim"); i++)
            x.middleRows(i*udofs, udofs) = m_solvers[i]->solve(input.middleRows(i*udofs, udofs));

    }

    int rows() const { return m_size; }

    int cols() const { return m_size; }

    virtual std::string getName() const
    {
        return "Adiag";
    }

protected:
    std::vector<gsSparseMatrix<T> > m_diagBlocks;
    int m_size;

    // members from uwbINSPrecondBlockMod
    using Base::m_solvers;

    // members from uwbINSPrecondBlock
    using Base::Base::m_opt;

}; // class uwbINSPrecondBlockAdiag


// ======================================== blockAmod ========================================

/** @brief
    Base class for block A of (modified) block preconditioner of the form
    \f[
    \left[ \begin{array}{cc}
    A & B^T \\
    0 & S
    \end{array}
    \right],
    \f].

    Implements action of \fA^{-1}\f. Assumes that the block A is not block-diagonal, solves a system with lower block-triangular submatrix of A (neglects the upper triangle).
*/
template <class T>
class uwbINSPrecondBlockAmod : public uwbINSPrecondBlockMod<T>
{
public:
    typedef uwbINSPrecondBlockMod<T> Base;

    /// Shared pointer for uwbINSPrecondBlockAmod
    typedef memory::shared_ptr<uwbINSPrecondBlockAmod> Ptr;

    /// Unique pointer for uwbINSPrecondBlockAmod
    typedef memory::unique_ptr<uwbINSPrecondBlockAmod> uPtr;

    uwbINSPrecondBlockAmod(const gsSparseMatrix<T>& matNS, const gsOptionList& opt) :
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

    /// Make function returning a smart pointer
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondBlockAmod<T>(matNS, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
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

    int rows() const { return m_size; }

    int cols() const { return m_size; }

    virtual std::string getName() const
    {
        return "Amod";
    }

protected:
    const gsSparseMatrix<T>& m_matRef;
    std::vector<gsSparseMatrix<T> > m_diagBlocks;
    int m_size;

    // members from uwbINSPrecondBlockMod
    using Base::m_solvers;

    // members from uwbINSPrecondBlock
    using Base::Base::m_opt;

}; // class uwbINSPrecondBlockAmod

// ======================================== block BT ========================================

template <class T>
class uwbINSPrecondBlockBt : public gsLinearOperator<T>
{
public:
    typedef gsLinearOperator<T> Base;

    /// Shared pointer for uwbINSPrecondBlockA
    typedef memory::shared_ptr<uwbINSPrecondBlockBt> Ptr;

    /// Unique pointer for uwbINSPrecondBlockA   
    typedef memory::unique_ptr<uwbINSPrecondBlockBt> uPtr;

    uwbINSPrecondBlockBt(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        int dim = opt.getInt("dim");
        int udofs = opt.getInt("udofs");
        m_mat = matNS.block(0, dim * udofs, dim * udofs, matNS.cols() - dim * udofs);
    }

    /// Make function returning a smart pointer
    static uPtr make(const gsSparseMatrix<T>& matNS, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondBlockBt<T>(matNS, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->cols(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x.noalias() = m_mat * input;
    }

    int rows() const { return m_mat.rows(); }

    int cols() const { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;

}; // class uwbINSPrecondBlockBt

// ======================================== Schur LSC ========================================

template <class T>
class uwbINSPrecondSchurLSC : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurLSC
    typedef memory::shared_ptr<uwbINSPrecondSchurLSC> Ptr;

    /// Unique pointer for uwbINSPrecondSchurLSC   
    typedef memory::unique_ptr<uwbINSPrecondSchurLSC> uPtr;

    uwbINSPrecondSchurLSC(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
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

        uwbINSPreconditioner<T>::diagInvMatrix_into(velM, velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));
        
        // TODO: could be more efficient if the matrices not involving block F were not updated (they do not change)
        gsSparseMatrix<T> BMinv, MinvBt;
        BMinv = matNS.block(uSize, 0, pdofs, uSize) * velMinv; // B * velM^-1
        MinvBt = velMinv * matNS.block(0, uSize, uSize, pdofs); // velM^-1 * BT

        m_mat1 = BMinv * matNS.block(0, uSize, uSize, pdofs); // B * velM^-1 * BT
        m_mat2 = BMinv * matNS.block(0, 0, uSize, uSize) * MinvBt; // B * velM^-1 * F * velM^-1 * BT

        m_pSolver->compute(m_mat1);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurLSC<T>(mat, opt));
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        // solve -(mat1 * mat2^-1 * mat1) * x = input
        gsMatrix<T> tmp;
        tmp = m_pSolver->solve(-1 * input); // mat1 * tmp = -input, where tmp = mat2^-1 * mat1 * x
        x = m_pSolver->solve(m_mat2 * tmp); // mat1 * x = mat2 * tmp
    }

    int rows() const { return m_mat2.rows(); }

    int cols() const { return m_mat2.cols(); }

private:
    gsSparseMatrix<T> m_mat1, m_mat2;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;

}; // class uwbINSPrecondSchurLSC

// ======================================== Schur PCD ========================================

template <class T>
class uwbINSPrecondSchurPCD : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurPCD
    typedef memory::shared_ptr<uwbINSPrecondSchurPCD> Ptr;

    /// Unique pointer for uwbINSPrecondSchurPCD   
    typedef memory::unique_ptr<uwbINSPrecondSchurPCD> uPtr;

    uwbINSPrecondSchurPCD(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
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

    ~uwbINSPrecondSchurPCD()
    {
        if (m_pMassSolver)
        {
            delete m_pMassSolver;
            m_pMassSolver = NULL;
        }
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurPCD<T>(mat, opt));
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        // solve - Ap * Fp^-1 * Mp * x = input
        gsMatrix<T> tmp1, tmp2;
        tmp1 = m_pSolver->solve(-1 * input); // Ap * y = - input
        tmp2 = m_matFp * tmp1; // z = Fp * y
        x = m_pMassSolver->solve(tmp2); // Mp * x = z
    }

    int rows() const { return m_matAp.rows(); }

    int cols() const { return m_matAp.cols(); }

protected:
    gsSparseMatrix<T> m_matAp, m_matFp, m_matMp;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;
    gsSparseSolver<T>* m_pMassSolver;

}; // class uwbINSPrecondSchurPCD

   // ======================================== Schur PCDmod ========================================

template <class T>
class uwbINSPrecondSchurPCDmod : public uwbINSPrecondSchurPCD<T>
{
public:
    typedef uwbINSPrecondSchurPCD<T> Base;

    /// Shared pointer for uwbINSPrecondSchurPCDmod
    typedef memory::shared_ptr<uwbINSPrecondSchurPCDmod> Ptr;

    /// Unique pointer for uwbINSPrecondSchurPCDmod   
    typedef memory::unique_ptr<uwbINSPrecondSchurPCDmod> uPtr;

    uwbINSPrecondSchurPCDmod(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(mat, opt)
    {
        // S = - Mp * Fp^-1 * Ap
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurPCDmod<T>(mat, opt));
    }

    void apply(const gsMatrix<T>& input, gsMatrix<T>& x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        // solve - Mp * Fp^-1 * Ap * x = input
        gsMatrix<T> tmp1, tmp2;
        tmp1 = m_pMassSolver->solve(-1 * input); // Mp * y = - input
        tmp2 = m_matFp * tmp1; // z = Fp * y
        x = m_pSolver->solve(tmp2);  // Ap * x = z
    }

private:

    using Base::m_matAp;
    using Base::m_matFp;
    using Base::m_matMp;
    using Base::m_pMassSolver;

    // members from uwbINSPrecondBlock
    using Base::Base::m_pSolver;

}; // class uwbINSPrecondSchurPCDmod

// ======================================== Schur AL ========================================

template <class T>
class uwbINSPrecondSchurAL : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurAL
    typedef memory::shared_ptr<uwbINSPrecondSchurAL> Ptr;

    /// Unique pointer for uwbINSPrecondSchurAL   
    typedef memory::unique_ptr<uwbINSPrecondSchurAL> uPtr;

    uwbINSPrecondSchurAL(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        // S = - (gamma + nu) Mp

        const gsSparseMatrix<T>& presM = mat.at("matMp");

        m_const = opt.getReal("visc") + opt.getReal("gamma");

        uwbINSPreconditioner<T>::diagInvMatrix_into(presM, m_mat, 1, opt.getSwitch("lumpingM"));
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurAL<T>(mat, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = - m_const * m_mat * input;
    }

    int rows() const { return m_mat.rows(); }

    int cols() const { return m_mat.cols(); }

private:
    real_t m_const;
    gsSparseMatrix<T> m_mat;

}; // class uwbINSPrecondSchurAL


// ======================================== Schur SIMPLE ========================================

template <class T>
class uwbINSPrecondSchurSIMPLE : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurSIMPLE
    typedef memory::shared_ptr<uwbINSPrecondSchurSIMPLE> Ptr;

    /// Unique pointer for uwbINSPrecondSchurSIMPLE   
    typedef memory::unique_ptr<uwbINSPrecondSchurSIMPLE> uPtr;

    uwbINSPrecondSchurSIMPLE(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
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
        uwbINSPreconditioner<T>::diagInvMatrix_into(blockA, Dinv, 1, opt.getSwitch("lumpingA"));

        m_mat = - matNS.block(uSize, 0, pdofs, uSize) * Dinv * matNS.block(0, uSize, uSize, pdofs);

        m_pSolver->compute(m_mat);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurSIMPLE<T>(mat, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = m_pSolver->solve(input);
    }

    int rows() const { return m_mat.rows(); }

    int cols() const { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;

}; // class uwbINSPrecondSchurSIMPLE

// ======================================== Schur MSIMPLER ========================================

template <class T>
class uwbINSPrecondSchurMSIMPLER : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurMSIMPLER
    typedef memory::shared_ptr<uwbINSPrecondSchurMSIMPLER> Ptr;

    /// Unique pointer for uwbINSPrecondSchurMSIMPLER   
    typedef memory::unique_ptr<uwbINSPrecondSchurMSIMPLER> uPtr;

    uwbINSPrecondSchurMSIMPLER(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
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
        uwbINSPreconditioner<T>::diagInvMatrix_into(velM, velMinv, uSize / velM.rows(), opt.getSwitch("lumpingM"));

        m_mat = - matNS.block(uSize, 0, pdofs, uSize) * velMinv * matNS.block(0, uSize, uSize, pdofs);

        m_pSolver->compute(m_mat);
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurMSIMPLER<T>(mat, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = m_pSolver->solve(input);
    }

    int rows() const { return m_mat.rows(); }

    int cols() const { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;

    // members from uwbINSPrecondBlock
    using Base::m_pSolver;

}; // class uwbINSPrecondSchurMSIMPLER

   // ======================================== Schur Stokes ========================================

template <class T>
class uwbINSPrecondSchurStokes : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurStokes
    typedef memory::shared_ptr<uwbINSPrecondSchurStokes> Ptr;

    /// Unique pointer for uwbINSPrecondSchurStokes   
    typedef memory::unique_ptr<uwbINSPrecondSchurStokes> uPtr;

    uwbINSPrecondSchurStokes(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        // S = diag(Mp)

        m_visc = opt.getReal("visc");

        const gsSparseMatrix<T>& presM = mat.at("matMp");
        uwbINSPreconditioner<T>::diagInvMatrix_into(presM, m_mat, 1, opt.getSwitch("lumpingM"));
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurStokes<T>(mat, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        x = - m_visc * m_mat * input;
    }

    int rows() const { return m_mat.rows(); }

    int cols() const { return m_mat.cols(); }

private:
    gsSparseMatrix<T> m_mat;
    T m_visc;

}; // class uwbINSPrecondSchurStokes

// ======================================== Schur exact ========================================

/*
Class performing multiplication by the exact Schur complement matrix.
*/
/*template <class T, class BlockAType = uwbINSPrecondBlockA<T> >
class uwbINSPrecondSchurExact : public uwbINSPrecondBlock<T>
{
public:
    typedef uwbINSPrecondBlock<T> Base;

    /// Shared pointer for uwbINSPrecondSchurExact
    typedef memory::shared_ptr<uwbINSPrecondSchurExact> Ptr;

    /// Unique pointer for uwbINSPrecondSchurExact   
    typedef memory::unique_ptr<uwbINSPrecondSchurExact> uPtr;

    uwbINSPrecondSchurExact(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt) :
        Base(opt), m_mat(mat.at("matNS")), m_Ainv(new BlockAType(mat.at("matNS"), opt))
    {
        // S = - B Ainv Bt

        m_size = opt.getInt("pdofs");
        m_uSize = opt.getInt("dim") * opt.getInt("udofs");
    }

    /// Make function returning a smart pointer
    static uPtr make(const std::map<std::string, gsSparseMatrix<T> >& mat, const gsOptionList& opt)
    {
        return memory::make_unique(new uwbINSPrecondSchurExact<T, BlockAType>(mat, opt));
    }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        GISMO_ASSERT(input.rows() == this->rows(), "Wrong input size.");
        x.resize(this->rows(), 1);

        // A w = - Bt * input
        // B * w = S * input
        gsMatrix<T> tmp = - m_mat.block(0, m_uSize, m_uSize, m_size) * input;
        gsMatrix<T> w;
        m_Ainv->apply(tmp, w);
        x = m_mat.block(m_uSize, 0, m_size, m_uSize) * w;
    }

    int rows() const { return m_size; }

    int cols() const { return m_size; }

private:
    typename uwbINSPrecondBlock<T>::Ptr m_Ainv;
    gsSparseMatrix<T> m_mat;
    int m_size, m_uSize;

};*/ // class uwbINSPrecondSchurExact

} // namespace gismo
