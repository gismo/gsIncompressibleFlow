/** @file uwbStabilizationEvaluator.h

    Author(s): E. Turnerova

*/

#pragma once

namespace gismo
{

template <class T>
class uwbStabilizationEvaluator
{
public:
    uwbStabilizationEvaluator()
    { }

    ~uwbStabilizationEvaluator()
    { }

    void initialize(index_t nPoints, index_t dim)
    {
        m_nPoints = nPoints;
        m_dim = dim;
        initMembers();
    }

    void initMembers()
    {
        m_diffusionCoeff.setZero(m_nPoints);
        m_advectionCoeff.resize(m_nPoints);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s].setZero(m_dim);
        m_reactionCoeff.setZero(m_nPoints);
        m_deg = 0;
        m_timeStep = 0.;
        m_elemDiam = 0.;
        m_r = -1;
        m_advElemDiam = 0.;
        m_diffElemDiam = 0.;
        m_alpha = 1.;
        m_bAdvDiffElemDiam = false;
    }

    void initAtElement(gsVector<T> advectionCoeff, T diffusionCoeff)
    {
        m_diffusionCoeff.setConstant(m_nPoints, diffusionCoeff);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionCoeff;
        m_reactionCoeff.setZero(m_nPoints);
    }

    void initAtElement(gsVector<T> advectionCoeff, T diffusionCoeff, T reactionCoeff)
    {
        m_diffusionCoeff.setConstant(m_nPoints, diffusionCoeff);
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionCoeff;
        m_reactionCoeff.setConstant(m_nPoints, reactionCoeff);
    }

    void initAtElement(gsMatrix<T>& advectionSolVals)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff.setZero(m_nPoints);
        m_reactionCoeff.setZero(m_nPoints);
    }

    void initAtElement(gsMatrix<T>& advectionSolVals, gsMatrix<T>& diffusionSolVals)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff = diffusionSolVals.row(0);
        m_reactionCoeff.setZero(m_nPoints);
    }

    /*void initAtElement(gsMatrix<T>& advectionSolVals, gsMatrix<T>& diffusionSolVals, gsMatrix<T>& reactionSolVals)
    {
        for (int s = 0; s < m_nPoints; ++s)
            m_advectionCoeff[s] = advectionSolVals.col(s);
        m_diffusionCoeff = diffusionSolVals.row(0);
        m_reactionCoeff = reactionSolVals.row(0);
    }*/

public:
    void setSUPGvars(index_t tauStabType, index_t deg, T timeStep, T elemDiam, index_t r = 2,
                     bool advDiffElemDiam = false, T advElemDiam = 0., T diffElemDiam = 0.)
    {
        m_tauStabType = tauStabType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_elemDiam = elemDiam;
        m_r = r;
        m_bAdvDiffElemDiam = advDiffElemDiam;
        m_advElemDiam = advElemDiam;
        m_diffElemDiam = diffElemDiam;
    }

    void setTauType(index_t tauType, index_t r = 2) { m_tauStabType = tauType; m_r = r; }

    /*void setCrosswindVars(index_t crosswindType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T elemDiam)
    {
        m_crosswindType = crosswindType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_elemDiam = elemDiam;
    }*/

    void setCrosswindVars(index_t crosswindType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T elemDiam)
    {
        m_crosswindType = crosswindType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_elemDiam = elemDiam;
    }

    void setTanhCSDvars(index_t csdType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T elemDiam, T alpha = 1.)
    {
        m_csdType = csdType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_elemDiam = elemDiam;
        m_alpha = alpha;
    }

    void setTanhCSDvars(index_t csdType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, std::vector< gsVector<T> >& residual_adv, T elemDiam, T alpha = 1.)
    {
        m_csdType = csdType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_residual_adv = residual_adv;
        m_elemDiam = elemDiam;
        m_alpha = alpha;
    }

    void setIsoADvars(index_t isoType, index_t deg, T timeStep, std::vector< gsMatrix<T> >& solGrads,
                          std::vector< gsVector<T> >& residual, T elemDiam, T alpha = 1.)
    {
        m_isoType = isoType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_solGrads = solGrads;
        m_residual = residual;
        m_elemDiam = elemDiam;
        m_alpha = alpha;
    }

    void setIsoADVars(std::vector< gsVector<T> >& residual, std::vector< gsMatrix<T> >& solGrads, T elemDiam)
    {
        m_residual = residual;
        m_solGrads = solGrads;
        m_elemDiam = elemDiam;
    }

    void setADVars(index_t tauType, index_t deg, T timeStep, T elemDiam, index_t r = 2)
    {
        m_tauStabType = tauType;
        m_deg = deg;
        m_timeStep = timeStep;
        m_elemDiam = elemDiam;
        m_r = r;
    }

    const gsVector<T> & getDiffusionCoefficient() const { return m_diffusionCoeff; }
    const T getDiffusionCoefficient(const int i) const { return m_diffusionCoeff(i); }
    const std::vector<gsVector<T> > & getAdvectionCoefficient() const { return m_advectionCoeff; }
    const gsVector<T> getAdvectionCoefficient(const int i) const { return m_advectionCoeff[i]; }
    const gsVector<T> & getReactionCoefficient() const { return m_reactionCoeff; }
    const T getReactionCoefficient(const int i) const { return m_reactionCoeff(i); }

    T getTauS(index_t k, T diffCoeff = -1., T reactCoeff = 0.)
    {
        if (m_elemDiam <= 0.)
            GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);

        //if (reactCoeff < 0.)
        //    m_reactionCoeff(k) = math::max(m_reactionCoeff(k), 0.);
        //else
            m_reactionCoeff(k) = math::max(reactCoeff, 0.);

        T tau_s = 0.;
        T normAdvectionCoeff = m_advectionCoeff[k].norm();
        T P = 0.;

        switch (m_tauStabType)
        {
        case 0:

            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            P = m_elemDiam * normAdvectionCoeff / (2. * m_diffusionCoeff(k));
            tau_s = (m_elemDiam / (2. * m_deg * normAdvectionCoeff)) * (1. / tanh(P) - 1. / P);

            if (P == 0. || normAdvectionCoeff == 0.)
                tau_s = 0.;

            break;

        case 1:

            tau_s = m_elemDiam / (2. * normAdvectionCoeff);

            if (normAdvectionCoeff < 1e-6)
                tau_s = 0.;

            break;

        case 2:
            tau_s = 1. / (math::sqrt(math::pow(2. * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, 2), 2)));
            break;

        case 3:
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, 2), 2)));
            break;

        case 4:

            P = m_elemDiam * normAdvectionCoeff / (2. * m_diffusionCoeff(k));
            tau_s = (m_elemDiam / (2. * normAdvectionCoeff)) * (1. / tanh(P) - 1. / P);

            if (P == 0. || normAdvectionCoeff == 0.)
                tau_s = 0.;

          break;

        case 5:

            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            tau_s = m_elemDiam / (2. * m_deg * normAdvectionCoeff);

            if (normAdvectionCoeff == 0.)
                tau_s = 0.;

          break;

        case 6:
            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            tau_s = 1. / (math::sqrt(math::pow(2. * m_deg * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, 2), 2)));
            break;

        case 7:
            GISMO_ASSERT(m_deg > 0, "'m_deg equals zero. It is not set or constant basis functions are used'.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * m_deg * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, 2), 2)));
            break;

        case 8:
            GISMO_ASSERT(m_r > 0, "diameter power is not positive. It is not set or constant basis functions are used'.");
            tau_s = 1. / (math::sqrt(math::pow(2. * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, m_r), 2)));
            break;

        case 9:
            GISMO_ASSERT(m_r > 0, "diameter power is not positive. It is not set or constant basis functions are used'.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, m_r), 2)));
            break;

        case 10:
            GISMO_ASSERT(m_deg > 0 && m_r > 0, "'m_deg equals zero or diameter power is not positive. It is not set or constant basis functions are used'.");
            tau_s = 1. / (math::sqrt(math::pow(2. * m_deg * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, m_r), 2)));
            break;

        case 11:
            GISMO_ASSERT(m_deg > 0 && m_r > 0, "'m_deg equals zero or diameter power is not positive. It is not set or constant basis functions are used'.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * m_deg * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, m_r), 2)));
            break;

        case 12:
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient not set in the stabilization evaluator.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep + m_reactionCoeff(k), 2) +
                    math::pow(2. * m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k)/math::pow(m_elemDiam, 2), 2)));
            break;

        case 13:
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient not set in the stabilization evaluator.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2./m_timeStep, 2) +
                    math::pow(2.*m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4.*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2) +
                    math::abs(m_reactionCoeff(k))));
            break;

        case 14:
            //GISMO_ASSERT(m_deg > 0 && reactCoeff >= 0., "'m_deg equals zero or reaction coefficient is not set.");
            tau_s = 1. / (math::sqrt(math::pow(2./m_timeStep, 2) +
                    math::pow(2.*m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4.*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2) +
                    math::abs(m_reactionCoeff(k))));
            break;

        case 15:
            //GISMO_ASSERT(m_deg > 0 && reactCoeff >= 0., "'m_deg equals zero or reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2./m_timeStep, 2) +
                    math::pow(2.*m_deg*normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4.*m_diffusionCoeff(k)/math::pow(m_elemDiam,2), 2) +
                    math::abs(m_reactionCoeff(k))));
            break;

        case 16: //Codina
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (2.*normAdvectionCoeff / m_elemDiam +
                    4.*math::abs(m_diffusionCoeff(k))/math::pow(m_elemDiam,2) +
                    math::abs(1./m_timeStep + m_reactionCoeff(k)));
            break;

        case 17: //Codina+deg
            //GISMO_ASSERT(m_deg > 0 && reactCoeff >= 0., "'m_deg equals zero or reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (2.*m_deg*normAdvectionCoeff / m_elemDiam +
                    4.*math::abs(m_diffusionCoeff(k))/math::pow(m_elemDiam,2) +
                    math::abs(1./m_timeStep + m_reactionCoeff(k)));
            break;

        case 18: //KLR
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = math::min(math::min(m_elemDiam/(2.*normAdvectionCoeff),
                                        1./math::abs(1./m_timeStep + m_reactionCoeff(k))),
                              math::pow(m_elemDiam,2)/math::abs(m_diffusionCoeff(k)));
            break;

        case 19: //KLR+deg
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient is not set.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = math::min(math::min(m_elemDiam/(2.*m_deg*normAdvectionCoeff),
                                        1./math::abs(1./m_timeStep + m_reactionCoeff(k))),
                              math::pow(m_elemDiam,2)/(math::pow(m_deg, 4)*math::abs(m_diffusionCoeff(k))));
            break;

        case 20: // = 3 with John+Schmeyer
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = (1. / m_timeStep) * (1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * normAdvectionCoeff / m_elemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_elemDiam, 2), 2))));
            break;

        case 21:
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient not set in the stabilization evaluator.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep + m_reactionCoeff(k), 2) +
                    math::pow(2. * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k)/math::pow(m_elemDiam, 2), 2)));
            break;

        case 22:
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient not set in the stabilization evaluator.");
            tau_s = 1. / (math::sqrt(math::pow(m_reactionCoeff(k), 2) +
                    math::pow(2. * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k)/math::pow(m_elemDiam, 2), 2)));
            break;

        case 23:
            //GISMO_ASSERT(reactCoeff >= 0., "reaction coefficient not set in the stabilization evaluator.");
            tau_s = 1. / (math::sqrt(math::pow(m_reactionCoeff(k), 2) +
                    math::pow(2. * m_deg * normAdvectionCoeff / m_elemDiam, 2) +
                    9. * math::pow(4. * m_diffusionCoeff(k)/math::pow(m_elemDiam, 2), 2)));
            break;

        case 222:
            GISMO_ASSERT(m_bAdvDiffElemDiam, "characteristic element diameter not set for advection and diffusion part.");
            tau_s = 1. / (math::sqrt(math::pow(2. * normAdvectionCoeff / m_advElemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemDiam, 2), 2)));
            break;

        case 333:
            GISMO_ASSERT(m_bAdvDiffElemDiam, "characteristic element diameter not set for advection and diffusion part.");
            if (m_timeStep == 0.)
                GISMO_ERROR("'m_timeStep' equals zero. It is not set or stationary problem is solved.");
            tau_s = 1. / (math::sqrt(math::pow(2. / m_timeStep, 2)
                  + math::pow(2. * normAdvectionCoeff / m_advElemDiam, 2)
                  + 9. * math::pow(4. * m_diffusionCoeff(k) / math::pow(m_diffElemDiam, 2), 2)));
            break;

        default: //default tau            
            gsInfo << "Wrong or no type of the stabilization parameter set in the TM visitor. Default tau is used!\n";

            tau_s = m_elemDiam / (2. * normAdvectionCoeff);

            if (normAdvectionCoeff == 0.)
                tau_s = 0.;
        }

        return tau_s;
    }

    T getCrosswindStabParam(index_t k, index_t variable = 0, T diffCoeff = -1., T reactCoeff = 0.)
    {
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);

        T cwStabParam = 0.;

        T normAdvection = m_advectionCoeff[k].norm();
        T denom;
        T C;

        switch (m_crosswindType) {
          case 0: // John and Knobloch 2005
            denom = normAdvection * (m_solGrads[k].row(variable)).norm() + math::abs(m_residual[variable](k));
            cwStabParam = getTauS(k, m_diffusionCoeff(k), reactCoeff) * math::pow(normAdvection, 2) * math::abs(m_residual[variable](k)) /
                          denom;
            if (denom == 0.)
                cwStabParam = 0.;

            break;

          case 1:
            if (m_elemDiam <= 0.)
                GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");
            C = 0.6;
            cwStabParam = 0.5 * math::max(0., C - (2 * m_diffusionCoeff(k) * (m_solGrads[k].row(variable)).norm())
                        / (m_elemDiam * math::abs(m_residual[variable](k))))
                        * m_elemDiam * math::abs(m_residual[variable](k)) / (m_solGrads[k].row(variable)).norm();

            break;

        case 2:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          //T elemdiam = m_elemDiam;
          //if (m_bAdvDiffElemDiam)
          //    elemdiam = m_diffElemDiam;

          cwStabParam = (m_elemDiam / 2.) * math::abs(m_residual[variable](k)) / m_solGrads[k].row(variable).norm();
          if (m_solGrads[k].row(variable).norm() == 0.)
              cwStabParam = 0.;

          break;

        case 3:
          if (m_elemDiam <= 0.)
              GISMO_ERROR("'m_elemDiam' equals zero. It is not set in the stabilization evaluator!");

          //T elemdiam = m_elemDiam;
          //if (m_bAdvDiffElemDiam)
          //    elemdiam = m_diffElemDiam;

          cwStabParam = math::pow(m_elemDiam / 2., 2) *
                  math::abs(m_residual[variable](k)) / math::abs(m_advectionCoeff[k](variable));
          if (m_advectionCoeff[k](variable) == 0.)
              cwStabParam = 0.;

          break;

        default:
            GISMO_ERROR("Wrong crosswind type set in ADR evaluator.");
            break;
        }

        return cwStabParam;
    }

    gsMatrix<T> getCrosswindProjection(index_t k)
    {
        gsMatrix<T> mIdentity;
        gsVector<T> advection = m_advectionCoeff[k];
        gsMatrix<T> tensorProduct(m_dim, m_dim);
        if (m_dim == 2)
            tensorProduct << advection(0)*advection(0), advection(0)*advection(1),
                             advection(1)*advection(0), advection(1)*advection(1);
        else if (m_dim == 3)
            tensorProduct << advection(0)*advection(0), advection(0)*advection(1), advection(0)*advection(2),
                             advection(1)*advection(0), advection(1)*advection(1), advection(1)*advection(2),
                             advection(2)*advection(0), advection(2)*advection(1), advection(2)*advection(2);
        else
            GISMO_ERROR("Wrong dimension in ADR evaluator. Computation possible only in 2D or 3D.");

        mIdentity.setIdentity(m_dim, m_dim);
        gsMatrix<T> projMatrix;
        if (advection.norm() == 0)
            projMatrix.setZero(m_dim, m_dim);
        else
            projMatrix = mIdentity - tensorProduct/math::pow(advection.norm(), 2);

        return projMatrix;
    }

    T getCSDstabParam(index_t k, index_t variable = 0, T diffCoeff = -1., T reactCoeff = 0., index_t r = 2)
    {
        m_r = r;
        /*if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);
        if (reactCoeff < 0.)
            m_reactionCoeff(k) = math::max(m_reactionCoeff(k), 0.);
        else
            m_reactionCoeff(k) = math::max(reactCoeff, 0.);*/

        //T normAdvectionCoeff = m_advectionCoeff[k].norm();

        T csdStabParam = 0.;
        switch (m_csdType)
        {
          case 0: //tanh-CSD
            csdStabParam = getTauS(k, diffCoeff, reactCoeff) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;

          case 1: //tanh-CSD h^alpha
            csdStabParam = math::pow(m_elemDiam, m_alpha) * getTauS(k, diffCoeff, reactCoeff) * math::pow(math::tanh(m_residual[variable](k)), 2);

          break;

          case 2: //tanh-tanh-CSD
            csdStabParam = getTauS(k, diffCoeff, reactCoeff) * math::pow(math::tanh(m_residual[variable](k)), 2)
                                      * math::pow(math::tanh(m_residual_adv[variable](k)), 2);

          break;

          default:
              GISMO_ERROR("Wrong CSDtype set in stabilization evaluator.");
              break;
          }

        return csdStabParam;
    }

    //------- original
    T getIsoADStabParamOriginal(index_t k, index_t variable = 0, T diffCoeff = -1.)
    {
        /*T alpha = 3./2.;
        T nu = 1.999;
        if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);*/

        // Johnson 1990: alpha and nu from (3/2, 2) such that Johnson suggested nu close to 2
        //T isoStabParam = math::max(0., alpha * math::pow(m_elemDiam, nu) * math::abs(m_residual[variable](k)) - m_diffusionCoeff(k));

        T denom = (m_solGrads[k].row(variable)).norm();
        //T isoStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(m_residual[variable](k)), 2) / denom;
        //T isoStabParam = m_elemDiam * math::abs(m_residual[variable](k)) / denom;
        T isoStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(math::abs(m_residual[variable](k)) / denom), 2);

        if (denom == 0.)
            isoStabParam = 0.;

        return isoStabParam;
    }

    T getRansIsoADStabParam(index_t k, index_t variable = 0)
    {
        T isoStabParam;

        switch (m_isoType)
        {
          case 3: //Nazarov h^alpha |Res|
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::abs(m_residual[variable](k));

            break;
          case 4: //Nazarov h^alpha tanh(Res)
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;
          default:
              GISMO_ERROR("Wrong isoADType set in stabilization evaluator.");
              break;
          }

        return isoStabParam;
    }

    T getIsoADStabParam(index_t k, index_t variable = 0, T diffCoeff = -1., index_t isoADtype = 6, index_t tauStabType = 2, T timeStep = -1.,
                        index_t deg = 0, T reactCoeff = 0., index_t r = 2)
    {
        m_r = r;
        m_tauStabType = tauStabType;
        m_timeStep = timeStep;
        m_deg = deg;
        /*if (diffCoeff < 0.)
            m_diffusionCoeff(k) = math::max(m_diffusionCoeff(k), 1e-15);
        else
            m_diffusionCoeff(k) = math::max(diffCoeff, 1e-15);
        if (reactCoeff < 0.)
            m_reactionCoeff(k) = math::max(m_reactionCoeff(k), 0.);
        else
            m_reactionCoeff(k) = math::max(reactCoeff, 0.);*/

        T normAdvectionCoeff = m_advectionCoeff[k].norm();

        T isoStabParam;

        switch (isoADtype)
        {
          case 1:
            isoStabParam = (1./getTauS(k, diffCoeff, reactCoeff)) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;

          case 2:
            isoStabParam = getTauS(k, diffCoeff, reactCoeff) * math::pow(math::tanh(m_residual[variable](k)), 2) * math::pow(normAdvectionCoeff, 2);

            break;

          case 3: //Nazarov h^alpha |Res|
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::abs(m_residual[variable](k));

            break;
          case 4: //Nazarov h^alpha tanh(Res)
            isoStabParam = math::pow(m_elemDiam, m_alpha) * math::pow(math::tanh(m_residual[variable](k)), 2);

            break;
          default:
              GISMO_ERROR("Wrong isoADType set in stabilization evaluator.");
              break;
          }

        return isoStabParam;
    }
    //----------

    /*T getIsoADStabParamTest(index_t k, index_t variable = 0)
    {
        T normAdvectionCoeff = m_advectionCoeff[k].norm();
        //T isoStabParam = m_elemDiam * math::pow(math::tanh(m_residual[variable](k)), 2) * normAdvectionCoeff / 2.;
        T isoStabParam = m_elemDiam * normAdvectionCoeff / 2.;
        //T isoStabParam = m_elemDiam * math::pow(math::tanh(m_residual[variable](k)), 2) / 2.;
        //T isoStabParam = (m_elemDiam / 2.) * math::pow(math::tanh(m_residual[variable](k)), 2) *
        //                 math::max(1., normAdvectionCoeff);

        return isoStabParam;
    }*/

protected:
    bool        m_bAdvDiffElemDiam;
    index_t     m_nPoints;
    index_t     m_dim;
    index_t     m_deg;
    index_t     m_r;
    index_t     m_tauStabType;
    index_t     m_crosswindType;
    index_t     m_csdType;
    index_t     m_isoType;
    T           m_timeStep;
    T           m_elemDiam;
    T           m_advElemDiam;
    T           m_diffElemDiam;
    T           m_alpha;
    std::vector< gsMatrix<T> > m_solGrads;
    std::vector< gsVector<T> > m_residual, m_residual_adv;
    gsVector<T> m_diffusionCoeff;
    std::vector<gsVector<T> > m_advectionCoeff;
    gsVector<T> m_reactionCoeff;
};

//=====================================================================================================================================

} // namespace gismo

