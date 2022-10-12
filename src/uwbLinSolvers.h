/** @file uwbLinSolvers.h

    Author(s): H. Hornikova
*/
#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

// Preconditioned iterative solver using the biconjugate gradient stabilized method.
template <class T>
class uwbBiCGStab : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;
    
    typedef gsMatrix<T>  VectorType;
    
    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<uwbBiCGStab> Ptr;
    typedef memory::unique_ptr<uwbBiCGStab> uPtr;

    template< typename OperatorType >
    explicit uwbBiCGStab( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr() )
    : Base(mat, precond) 
    {
        eps = std::numeric_limits<T>::epsilon() *std::numeric_limits<T>::epsilon();
    }

    bool initIteration(const VectorType& rhs, VectorType& x)
    {
        if (Base::initIteration(rhs, x))
            return true;

        m_rhs = rhs;
        index_t n = m_rhs.rows();

        m_mat->apply(x, tmp);
        res = m_rhs - tmp;
        r0 = res;
        r0_sqNorm = r0.squaredNorm();

        m_error = res.norm() / m_rhs_norm;
        if (m_error < m_tol)
            return true;
        
        m_restarts = 0;
        rho = 1;
        alpha = 1;
        omega = 1;

        v = VectorType::Zero(n,1);
        p = VectorType::Zero(n,1);

        return false;
    }

    bool step(VectorType& x)
    {
        rho_prev = rho;
        rho = res.asVector().dot(r0.asVector());

        if (math::abs(rho) < eps * r0_sqNorm)
        {
            m_mat->apply(x, tmp);
            res = m_rhs - tmp;
            r0 = res;
            r0_sqNorm = r0.squaredNorm();
            rho = r0_sqNorm;
            
            if (m_restarts++ == 0)
                m_num_iter = -1;
        }

        beta = (alpha / omega) * (rho / rho_prev);
        p = res + beta * (p - omega * v);

        m_precond->apply(p, q);
        m_mat->apply(q, v);

        alpha = rho / r0.asVector().dot(v.asVector());
        s = res - alpha * v;

        m_precond->apply(s, z);
        m_mat->apply(z, t);

        t_norm = t.squaredNorm();
        if (t_norm > T(0))
            omega = s.asVector().dot(t.asVector()) / t_norm;
        else
            omega = T(0);

        x += alpha * q + omega * z;
        res = s - omega * t;

        m_error = res.norm() / m_rhs_norm;
        if (m_error < m_tol)
            return true;

        return false;
    }

    std::string detail() const
    {
        std::ostringstream os;
        print(os);
        os << " Tolerance            : " << this->tolerance() << "\n";
        os << " Solver error         : " << this->error() << "\n";
        os << " Number of restarts   : " << m_restarts << "\n";
        os << " Number of iterations : " << this->iterations() << " (max=" << m_max_iters << ")\n";
        return os.str();
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "uwbBiCGStab\n";
        return os;
    }

    int restarts() const { return m_restarts; }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    int m_restarts;
    T rho_prev, rho, alpha, beta, omega, t_norm, eps, r0_sqNorm;
    VectorType m_rhs;
    VectorType res, r0;
    VectorType p, v, q, s, t, z, tmp;
};

template <class T>
class uwbGCR : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T> VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<uwbGCR> Ptr;
    typedef memory::unique_ptr<uwbGCR> uPtr;

    template< typename OperatorType >
    explicit uwbGCR(const OperatorType& mat, const LinOpPtr & precond = LinOpPtr())
        : Base(mat, precond)
    { }

    bool initIteration(const VectorType& rhs, VectorType& x)
    {
        if (Base::initIteration(rhs, x))
            return true;

        p_hist.clear();
        q_hist.clear();

        m_mat->apply(x, tmp);
        res = rhs - tmp;

        m_error = res.norm() / m_rhs_norm;
        if (m_error < m_tol)
            return true;

        m_precond->apply(res, z);
        p = z;
        p_hist.push_back(p);

        m_mat->apply(p, q);
        q_hist.push_back(q);

        return false;
    }

    bool step(VectorType& x)
    {
        alpha = (res.asVector().dot(q.asVector())) / q.squaredNorm();

        x = x + alpha * p;
        res = res - alpha * q;

        m_error = res.norm() / m_rhs_norm;
        if (m_error < m_tol)
            return true;

        m_precond->apply(res, z);
        m_mat->apply(z, Az);

        beta.clear();
        for (size_t i = 0; i < q_hist.size(); i++)
        {
            T tmpBeta = -(Az.asVector().dot(q_hist[i].asVector())) / q_hist[i].squaredNorm();
            beta.push_back(tmpBeta);
        }
            
        p = z;
        q = Az;
        for (size_t i = 0; i < q_hist.size(); i++)
        {
            p += beta[i] * p_hist[i];
            q += beta[i] * q_hist[i];
        }
        p_hist.push_back(p);
        q_hist.push_back(q);

        return false;
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "uwbGCR\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    T alpha;
    VectorType res, z, Az, tmp, p, q;
    std::vector<T> beta;
    std::vector<VectorType> p_hist, q_hist;
};

template<class T = real_t>
class uwbGMResRight : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;
    typedef gsMatrix<T> VectorType;
    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<uwbGMResRight> Ptr;
    typedef memory::unique_ptr<uwbGMResRight> uPtr;

    template< typename OperatorType >
    explicit uwbGMResRight(const OperatorType& mat, const LinOpPtr& precond = LinOpPtr())
        : Base(mat, precond) {}

    template< typename OperatorType >
    static uPtr make(const OperatorType& mat, const LinOpPtr& precond = LinOpPtr())
    { return uPtr(new uwbGMResRight(mat, precond)); }

    bool initIteration(const typename Base::VectorType& rhs,
        typename Base::VectorType& x)
    {
        if (Base::initIteration(rhs, x))
            return true;

        m_mat->apply(x, tmp);
        residual = rhs - tmp;
        beta = residual.norm(); // This is  ||r||

        m_error = beta / m_rhs_norm;
        if (m_error < m_tol)
            return true;

        v.push_back(residual / beta);
        g.setZero(2, 1);
        g(0, 0) = beta;
        Omega = gsMatrix<T>::Identity(2, 2);
        Omega_prev = gsMatrix<T>::Identity(2, 2);

        return false;
    }

    void finalizeIteration(typename uwbGMResRight<T>::VectorType& x)
    {
        //Remove last row of H and g
        H.resize(m_num_iter, m_num_iter);
        H = H_prev.block(0, 0, m_num_iter, m_num_iter);
        g_tmp.resize(m_num_iter, 1);
        g_tmp = g.block(0, 0, m_num_iter, 1);

        //Solve H*y = g;
        solveUpperTriangular(H, g_tmp);

        //Create the matrix from the column matrix in v.
        gsMatrix<T> V(m_mat->rows(), m_num_iter);
        for (index_t k = 0; k< m_num_iter; ++k)
        {
            V.col(k) = v[k];
        }
        //Update solution
        gsMatrix<T> pVy;
        m_precond->apply(V*y, pVy);
        x += pVy;

        // cleanup temporaries
        tmp.clear();
        g.clear();
        g_tmp.clear();
        h_tmp.clear();
        y.clear();
        w.clear();
        residual.clear();
        H_prev.clear();
        H.clear();
        Omega.clear();
        Omega_prev.clear();
        Omega_tmp.clear();
        Omega_prev_tmp.clear();
        v.clear();
    }

    bool step(typename uwbGMResRight<T>::VectorType& x)
    {
        GISMO_UNUSED(x); // The iterate x is never updated! Use finalizeIteration to obtain x.
        const index_t k = m_num_iter - 1;
        H.setZero(k + 2, k + 1);
        h_tmp.setZero(k + 2, 1);

        if (k != 0)
        {
            H.block(0, 0, k + 1, k) = H_prev;
        }

        Omega = gsMatrix<T>::Identity(k + 2, k + 2);
        m_precond->apply(v[k], tmp);
        m_mat->apply(tmp, w);

        for (index_t i = 0; i< k + 1; ++i)
        {
            h_tmp(i, 0) = (w.transpose()*v[i]).value(); //Typo h_l,k
            w = w - h_tmp(i, 0)*v[i];
        }
        h_tmp(k + 1, 0) = w.norm();

        if (math::abs(h_tmp(k + 1, 0)) < 1e-16) //If exact solution
            return true;

        v.push_back(w / h_tmp(k + 1, 0));

        h_tmp = Omega_prev*h_tmp;
        H.block(0, k, k + 2, 1) = h_tmp;

        //Find coef in rotation matrix
        T sk = H(k + 1, k) / (math::sqrt(H(k, k)*H(k, k) + H(k + 1, k)*H(k + 1, k)));
        T ck = H(k, k) / (math::sqrt(H(k, k)*H(k, k) + H(k + 1, k)*H(k + 1, k)));
        Omega(k, k) = ck; Omega(k, k + 1) = sk;
        Omega(k + 1, k) = -sk; Omega(k + 1, k + 1) = ck;

        //Rotate H and g
        H = Omega*H;
        H_prev = H;
        g_tmp.setZero(k + 2, 1);
        if (k != 0)
            g_tmp.block(0, 0, k + 1, 1) = g;
        else
            g_tmp = g;
        g.noalias() = Omega * g_tmp;

        T residualNorm2 = g(k + 1, 0)*g(k + 1, 0);
        m_error = math::sqrt(residualNorm2) / m_rhs_norm;
        if (m_error < m_tol)
            return true;

        //Resize rotation product
        Omega_prev_tmp.setZero(k + 3, k + 3);
        Omega_prev_tmp.block(0, 0, k + 2, k + 2) = Omega_prev;
        Omega_prev_tmp(k + 2, k + 2) = 1;

        Omega_tmp.setZero(k + 3, k + 3);
        Omega_tmp.block(0, 0, k + 2, k + 2) = Omega;
        Omega_tmp(k + 2, k + 2) = 1;

        Omega_prev = Omega_tmp*Omega_prev_tmp;
        return false;
    }

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const VectorType& R, const VectorType& gg)
    {
        y = R.template triangularView<Eigen::Upper>().solve(gg);
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "uwbGMResRight\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;


    gsMatrix<T> tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<T> residual;
    gsMatrix<T> H_prev, H, Omega, Omega_prev, Omega_tmp, Omega_prev_tmp;
    std::vector< gsMatrix<T> > v;
    T beta; 
};

} // namespace gismo
