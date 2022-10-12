/** @file uwbADREvaluators.h

    Author(s): E. Turnerova

*/

#pragma once

namespace gismo
{

template <class T>
class uwbADREvaluator
{
public:
    uwbADREvaluator(index_t nPoints, index_t dim) : m_nPoints(nPoints), m_dim(dim)
    { }

    virtual ~uwbADREvaluator()
    {  }

    virtual void initialize()
    {
        initMembers();
    }

    virtual void initMembers()
    {
        m_diffusionCoeff.setZero(m_nPoints);
        m_advectionCoeff.resize(m_nPoints);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s].setZero(m_dim);
        m_reactionCoeff.setZero(m_nPoints);
        m_rhs.setZero(m_nPoints);
        m_deg = 0;
        m_timeStep = 0.;
        m_elemDiam = 0.;
    }

    void initAtElement(T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff)
    {
        m_diffusionCoeff.setConstant(m_nPoints, diffusionCoeff);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionCoeff;
        m_reactionCoeff.setConstant(m_nPoints, reactionCoeff);
    }

    void initAtElement(gsMatrix<T>& advectionSolVals, gsMatrix<T>& diffusionSolVals)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff = diffusionSolVals.row(0);
        m_reactionCoeff.setZero(m_nPoints);
    }

    void initAtElement(gsMatrix<T>& advectionSolVals, gsMatrix<T>& diffusionSolVals, gsMatrix<T>& reactionSolVals)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff = diffusionSolVals.row(0);
        m_reactionCoeff = reactionSolVals.row(0);
    }

    //virtual void evalCoefficients() { GISMO_NO_IMPLEMENTATION }
    //virtual void evalRhs() { GISMO_NO_IMPLEMENTATION }

public:
    void setSUPGvars(index_t tauStabType, index_t deg, T timeStep, T elemDiam)
    {
        m_tauStabType = tauStabType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_elemDiam = elemDiam;
    }

    void setTauType(index_t tauType) { m_tauStabType = tauType; }

    void setCrosswindVars(index_t crosswindType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads, gsVector<T>& residual, T elemDiam)
    {
        m_crosswindType = crosswindType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_elemDiam = elemDiam;
    }

    void setIsoADVars(gsVector<T>& residual, T elemDiam)
    {
        m_residual = residual;
        m_elemDiam = elemDiam;
    }

    void setADVars(index_t tauType, index_t deg, T timeStep, T elemDiam)
    {
        m_tauStabType = tauType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_elemDiam = elemDiam;
    }

    //    void setADRtype(std::string ADRtype = "") { m_ADRtype = ADRtype; }

    const gsVector<T> & getRhs() const { return m_rhs; }
    const T getRhs(const int i) const { return m_rhs(i); }
    const gsVector<T> & getDiffusionCoefficient() const { return m_diffusionCoeff; }
    const T getDiffusionCoefficient(const int i) const { return m_diffusionCoeff(i); }
    const std::vector<gsVector<T>> & getAdvectionCoefficient() const { return m_advectionCoeff; }
    const gsVector<T> getAdvectionCoefficient(const int i) const { return m_advectionCoeff[i]; }
    const gsVector<T> & getReactionCoefficient() const { return m_reactionCoeff; }
    const T getReactionCoefficient(const int i) const { return m_reactionCoeff(i); }

    T getTauS(index_t k)
    {
        if (m_elemDiam <= 0.)
            GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the evaluator!");
        T tau_s = 0.;
        T normAdvectionCoeff = m_advectionCoeff[k].norm();
        T P = 0.;

        switch (m_tauStabType) {
          case 1:

            tau_s = m_elemDiam / (2 * normAdvectionCoeff);

            if (normAdvectionCoeff == 0)
                tau_s = 0.;

            break;

        case 2:
            tau_s = 1. / (math::sqrt(math::pow(2*normAdvectionCoeff / m_elemDiam, 2) + 9 * math::pow(4*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2)));
            break;

        case 3:
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2/m_timeStep, 2) + math::pow(2*normAdvectionCoeff / m_elemDiam, 2) + 9 * math::pow(4*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2)));
            break;

        case 4:

            P = m_elemDiam * normAdvectionCoeff / (2 * m_diffusionCoeff(k));
            tau_s = (m_elemDiam / (2 * normAdvectionCoeff)) * (1 / tanh(P) - 1 / P);

            if (P == 0 || normAdvectionCoeff == 0)
                tau_s = 0.;

          break;

        case 5:

            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            tau_s = m_elemDiam / (2 * m_deg * normAdvectionCoeff);

            if (normAdvectionCoeff == 0)
                tau_s = 0.;

          break;

        case 6:
            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            tau_s = 1. / (math::sqrt(math::pow(2*m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9 * math::pow(4*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2)));
            break;

        case 7:
            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2./m_timeStep, 2) +
                    math::pow(2.*m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4.*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2)));
            break;

        default: //default tau

            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            P = m_elemDiam * normAdvectionCoeff / (2 * m_diffusionCoeff(k));
            tau_s = (m_elemDiam / (2 * m_deg * normAdvectionCoeff)) * (1 / tanh(P) - 1 / P);

            if (P == 0 || normAdvectionCoeff == 0)
                tau_s = 0.;
        }
        return tau_s;
    }

    T getCrosswindStabParam(index_t k)
    {
        T cwStabParam = 0.;

        T normAdvection = m_advectionCoeff[k].norm();
        T denom;

        switch (m_crosswindType) {
          case 0: // John and Knobloch 2005
            denom = normAdvection * m_solGrads[k].norm() + math::abs(m_residual(k));
            cwStabParam = getTauS(k) * math::pow(normAdvection, 2) * math::abs(m_residual(k)) /
                          denom;
            if (denom == 0)
                cwStabParam = 0.;

            break;

        default:
            GISMO_ERROR("Wrong crosswind type set in ADR evaluator.");
        }

        return cwStabParam;
    }

    gsMatrix<T> getCrosswindProjection(index_t k)
    {
        gsMatrix<T> mIdentity;
        mIdentity.setIdentity(m_dim, m_dim);
        gsMatrix<T> projMatrix;
        if (getAdvectionCoefficient(k).norm() == 0)
            projMatrix.setZero(m_dim, m_dim);
        else
            projMatrix = mIdentity -
                         (getAdvectionCoefficient(k) * (getAdvectionCoefficient(k).transpose())) / getAdvectionCoefficient(k).norm();

        return projMatrix;
    }

    T getIsoADStabParam(index_t k)
    {
        T alpha = 3./2.;
        T nu = 1.9;

        // Johnson 1990: alpha and nu from (3/2, 2) such that Johnson suggested nu close to 2
        T isoStabParam = math::max(0., alpha * math::pow(m_elemDiam, nu) * math::abs(m_residual(k)) - getDiffusionCoefficient(k));
        return isoStabParam;
    }

protected:
    index_t     m_nPoints;
    index_t     m_dim;
    index_t     m_deg;
    index_t     m_tauStabType;
    index_t     m_crosswindType;
    T           m_timeStep;
    T           m_elemDiam;
    std::vector< gsMatrix<T> > m_solGrads;
    gsVector<T> m_residual;
    gsVector<T> m_rhs;
    gsVector<T> m_diffusionCoeff;
    std::vector<gsVector<T> > m_advectionCoeff;
    gsVector<T> m_reactionCoeff;
};

//=====================================================================================================================================

} // namespace gismo

