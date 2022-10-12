/** @file uwbTMEvaluators.h

    Author(s): E. Turnerova

	k-omega Wilcox 2006 - LRN
	https://turbmodels.larc.nasa.gov/wilcox.html
    k-omega SST Menter
    https://turbmodels.larc.nasa.gov/sst.html
    k-omega SST Menter 2009
    k-omega SST SAS
    k-omega SST SAS SS
    k-omega SST SAS SO
    k-omega SST SAS OO

*/

#pragma once

namespace gismo
{
template <class T>
class uwbTMEvaluator;
//-----------------------------------------------------------------------------------------------------------------
template <class T>
class uwbTMEvaluatorKOmega : public uwbTMEvaluator<T>
{

public:
    typedef uwbTMEvaluator<T> Base;

    uwbTMEvaluatorKOmega() 
    {
        m_betaStar0 = 0.09;
        m_KOtype = "";
    }

public:
    virtual void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals, const std::vector<gsMatrix<T> >& solKOmegaGrads)
    {
        initAtElement(solUGrads, solKOVals);
        m_solKOmegaGrads = solKOmegaGrads;
    }

    virtual void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals)
    {
        Base::initAtElement(solUGrads);
        m_solKOVals = solKOVals;
    }

protected:

    virtual void initMembers()
    {
        Base::initMembers();
        m_beta.setZero(m_nPoints);
        m_betaStar.setZero(m_nPoints);
        m_kDiffusionCoeff.setZero(m_nPoints);
        m_oDiffusionCoeff.setZero(m_nPoints);
        m_gamma.setZero(m_nPoints);
        m_blendCoeff.setZero(m_nPoints);
        m_rhsK.setZero(m_nPoints);
        m_rhsO.setZero(m_nPoints);
    }

public:
    const gsVector<T> & getBeta() const { return m_beta; }
    const T getBeta(const int i) const { return m_beta(i); }
    const gsVector<T> & getBetaStar() const { return m_betaStar; }
    const T getBetaStar(const int i) const { return m_betaStar(i); }
    const gsVector<T> & getKDiffusionCoefficient() const { return m_kDiffusionCoeff; }
    const T getKDiffusionCoefficient(const int i) const { return m_kDiffusionCoeff(i); }
    const gsVector<T> & getOmegaDiffusionCoefficient() const { return m_oDiffusionCoeff; }
    const T getOmegaDiffusionCoefficient(const int i) const { return m_oDiffusionCoeff(i); }
    const gsVector<T> & getRhsK() const { return m_rhsK; }
    const T getRhsK(const int i) const { return m_rhsK(i); }
    const gsVector<T> & getRhsOmega() const { return m_rhsO; }
    const T getRhsOmega(const int i) const { return m_rhsO(i); }
    const gsVector<T> & getBlendCoeff() const { return m_blendCoeff; }
    const T getBlendCoeff(const int i) const { return m_blendCoeff(i); }

    virtual const gsVector<T> & getF1() const { GISMO_NO_IMPLEMENTATION }
    virtual const T getF1(const int i) const { GISMO_NO_IMPLEMENTATION }
    virtual const gsVector<T> & getSourceQSAS() const { GISMO_NO_IMPLEMENTATION }
    virtual const T getSourceQSAS(const int i) const { GISMO_NO_IMPLEMENTATION }
    virtual const gsVector<T> & getReactionAtRhs() const { GISMO_NO_IMPLEMENTATION }
    virtual const T getReactionAtRhs(const int i) const { GISMO_NO_IMPLEMENTATION }

    virtual void evalQuantities_blendCoeff() { GISMO_NO_IMPLEMENTATION }
    virtual void evalTurbVariables() { GISMO_NO_IMPLEMENTATION }
    virtual void evalKOmegaDiffusionCoefficient() { GISMO_NO_IMPLEMENTATION }

    void setULaplacian(std::vector<gsMatrix<T> > solULaplaces){ m_solULaplaces = solULaplaces; }

protected:
    T m_betaStar0;

    gsMatrix<T> m_solKOVals;
    std::vector<gsMatrix<T> > m_solKOmegaGrads;

    std::vector<gsMatrix<T> > m_solULaplaces;

    gsVector<T> m_beta;
    gsVector<T> m_betaStar;

    gsVector<T> m_kDiffusionCoeff;
    gsVector<T> m_oDiffusionCoeff;
    gsVector<T> m_gamma;
    gsVector<T> m_blendCoeff;
    gsVector<T> m_rhsK;
    gsVector<T> m_rhsO;

    //std::string m_KOtype;

    using Base::m_nPoints;
    using Base::m_KOtype;
};

//==================================================================================================================
template <class T>
class uwbTMEvaluatorKOmegaWilcoxLRN : public uwbTMEvaluatorKOmega<T>
{

public:
    /// Unique pointer for uwbTMEvaluatorKOmegaWilcoxLRN
    typedef memory::unique_ptr<uwbTMEvaluatorKOmegaWilcoxLRN> uPtr;

    typedef uwbTMEvaluatorKOmega<T> Base;

    uwbTMEvaluatorKOmegaWilcoxLRN()
    {
        m_sigmaK = 0.6;
        m_sigmaO = 0.5;
        m_Clim = 7. / 8;
        m_Rbeta = 8;
        m_Rk = 6;
        m_Ro = 2.61;
        m_alpha0 = 1. / 9;
        m_beta0 = 0.0708;
        m_alphaStar0 = m_beta0 / 3;
        m_KOtype = "koWilcoxLRN";
    }

    /// Make function returning a smart pointer
    static uPtr make() { return memory::make_unique(new uwbTMEvaluatorKOmegaWilcoxLRN<T>()); }

protected:

    void initMembers()
    {
        Base::initMembers();

        m_alphaStar.setZero(m_nPoints);
        m_chiNumerator.setZero(m_nPoints);
        m_ReT.setZero(m_nPoints);
        m_kOmegaFraction.setZero(m_nPoints);
        m_omegaHat.setZero(m_nPoints);
    }

public:

    void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals, const std::vector<gsMatrix<T> >& solKOmegaGrads)
    {
        // add check of nPoints, if doesn't match, call initMembers with new value?

        Base::initAtElement(solUGrads, solKOVals, solKOmegaGrads);

        m_chiNumerator.setZero(m_nPoints);
    }

    void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals)
    {
        Base::initAtElement(solUGrads, solKOVals);

        m_chiNumerator.setZero(m_nPoints);
    }

    void evalQuantities_turbViscosity() //evaluate quantities needed to evaluate turbulent viscosity
    {
        evalCommonPart();
        evalOmegaHat();

        evalTurbulentViscosity_into(m_turbViscosity);
    }

    void evalQuantities_nonlinearBlocksPart() //evaluate quantities needed to evaluate nonlinearBlocks in TMvisitors
    {
        evalCommonPart();
        evalKOmegaDiffusionCoefficient();        
        evalBetaStar();
        evalBeta();
        evalBlendCoefficient();
    }

    void evalQuantities_rhsPart() //evaluate quantities needed to evaluate source terms of TM equations
    {
        evalCommonPart();
        evalOmegaHat();
        evalGamma();

        evalTurbulentViscosity_into(m_turbViscosity);

        evalKOmegaRhs();
    }

    void evalQuantities_diffusionCoeff() //evaluate quantities needed to evaluate nonlinearBlocks in TMvisitors
    {
        evalCommonPart();
        evalKOmegaDiffusionCoefficient();
    }

    void evalQuantities_reactionCoeff()
    {
        evalCommonPart();
        evalBetaStar();
        evalBeta();
    }

    void evalQuantities_blendCoeff()
    {
        evalBlendCoefficient();
    }
    
    void evalCommonPart()
    {
        uwbTMEvaluator<T>::evalStrainRateTimesUGrad();
        uwbTMEvaluator<T>::evalVorticityStrainRateMagnitude();

        for (int k = 0; k < m_nPoints; k++)
        {
            // eval k/omega 
            m_kOmegaFraction(k) = m_solKOVals.coeff(0, k) / math::max(m_solKOVals.coeff(1, k), 1e-15);
            
            // eval ReT 
            m_ReT(k) = math::max((1 / m_viscosity) * m_kOmegaFraction(k), 0.0); // set to zero if negative

            // eval alphaStar 
            m_alphaStar(k) = (m_alphaStar0 + m_ReT(k) / m_Rk) / (1 + m_ReT(k) / m_Rk);

            // eval strainRate quantities
            T vorticity_ij;
            T vorticity_jk;
            T strainRate_ki;
            index_t dim = m_solUGrads[0].cols(); // space dimension

            for (index_t i = 0; i != dim; ++i)
            {
                for (index_t j = 0; j != dim; ++j)
                {
                    vorticity_ij = 0.5 * (m_solUGrads[k].coeff(i, j) - m_solUGrads[k].coeff(j, i));
                    for (index_t kk = 0; kk != dim; ++kk)
                    {
                        vorticity_jk = 0.5 * (m_solUGrads[k].coeff(j, kk) - m_solUGrads[k].coeff(kk, j));
                        strainRate_ki = 0.5 * (m_solUGrads[k].coeff(kk, i) + m_solUGrads[k].coeff(i, kk));
                        m_chiNumerator(k) += vorticity_ij * vorticity_jk * strainRate_ki;
                    }
                }
            }
        }
    }

    void evalOmegaHat()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            T tmp = math::max(m_solKOVals.coeff(1, k), 1e-15);
            m_omegaHat(k) = math::max(tmp, m_Clim * math::sqrt(m_alphaStar(k) * math::pow(m_strainRateMagnitude(k), 2) / m_betaStar0));
        }
   }

    void evalTurbulentViscosity_into(gsVector<T> & turbViscosityVals)
    {
        turbViscosityVals.resize(m_nPoints);

        for (int k = 0; k < m_nPoints; k++)
            turbViscosityVals(k) = math::max(m_alphaStar(k) * (m_solKOVals.coeff(0, k) / m_omegaHat(k)), 1e-15);

        if (m_average)
            this->computeAverage(turbViscosityVals);
    }

    void evalBeta()
    {
        T f_beta;
        T chi_omega;
        for (int k = 0; k < m_nPoints; k++)
        {
            chi_omega = math::abs(m_chiNumerator(k) / math::pow(m_betaStar(k) * math::max(m_solKOVals.coeff(1, k), 1e-15), 3));
            f_beta = (1 + 85 * chi_omega) / (1 + 100 * chi_omega);
            m_beta(k) = m_beta0 * f_beta;
        }
    }

    void evalBetaStar()
    {
        for (int k = 0; k < m_nPoints; k++)
            m_betaStar(k) = m_betaStar0 * (((100. / 27) * m_beta0 + math::pow(m_ReT(k) / m_Rbeta, 4)) / (1 + math::pow(m_ReT(k) / m_Rbeta, 4)));

    }

    void evalKOmegaDiffusionCoefficient()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            m_kDiffusionCoeff(k) = math::max(m_viscosity + m_sigmaK * m_alphaStar(k) * m_kOmegaFraction(k), 1e-15);
            m_oDiffusionCoeff(k) = math::max(m_viscosity + m_sigmaO * m_alphaStar(k) * m_kOmegaFraction(k), 1e-15);
        }
    }

    void evalKOmegaRhs()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            m_rhsK(k) = 2 * m_turbViscosity(k) * m_strainRateTimesUGrad(k);
            m_rhsO(k) = 2 * m_gamma(k) * (m_solKOVals.coeff(1, k) / m_omegaHat(k)) * m_strainRateTimesUGrad(k); // m_alphaStar divided 
        }
        /*if (m_average)
        {
            this->computeAverage(m_rhsK);
            this->computeAverage(m_rhsO);
        }*/
    }

    void evalGamma()
    {
        for (int k = 0; k < m_nPoints; k++)
            m_gamma(k) = (13. / 25) * (m_alpha0 + (m_ReT(k) / m_Ro)) / (1 + (m_ReT(k) / m_Ro)); // without alphaStar -> divided by the turbulence viscosity quantity in the source term of the omega eqn.

    }

    void evalBlendCoefficient()
    {
        //sigma_d
        for (int k = 0; k < m_nPoints; k++)
        {
            m_blendCoeff(k) = 0;
            if (m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1)) > 0)
                m_blendCoeff(k) = 1 / 8;
        }
    }

protected:
    T m_sigmaK;
    T m_sigmaO;
    T m_Clim;
    T m_Rbeta;
    T m_Rk;
    T m_Ro;
    T m_alpha0;
    T m_beta0;
    T m_alphaStar0;

    gsVector<T> m_alphaStar;
    gsVector<T> m_chiNumerator;
    gsVector<T> m_ReT;
    gsVector<T> m_kOmegaFraction;
    gsVector<T> m_omegaHat;

    //memberes from uwbTMEvaluatorKOmega
    using Base::m_solKOVals;
    using Base::m_solKOmegaGrads;
    using Base::m_beta;
    using Base::m_betaStar;
    using Base::m_betaStar0;
    using Base::m_kDiffusionCoeff;
    using Base::m_oDiffusionCoeff;
    using Base::m_gamma;
    using Base::m_blendCoeff;
    using Base::m_rhsK;
    using Base::m_rhsO;
    using Base::m_KOtype;
    
    //memberes from uwbTMEvaluator
    using uwbTMEvaluator<T>::m_viscosity;
    using uwbTMEvaluator<T>::m_turbViscosity;
    using uwbTMEvaluator<T>::m_nPoints;
    using uwbTMEvaluator<T>::m_solUGrads;
    using uwbTMEvaluator<T>::m_strainRateTimesUGrad;
    using uwbTMEvaluator<T>::m_strainRateMagnitude;
    using uwbTMEvaluator<T>::m_average;
};

//=====================================================================================================================================
template <class T>
class uwbTMEvaluatorKOmegaSST : public uwbTMEvaluatorKOmega<T>
{

public:
    /// Unique pointer for uwbTMEvaluatorKOmegaSST
    typedef memory::unique_ptr<uwbTMEvaluatorKOmegaSST> uPtr;

    typedef uwbTMEvaluatorKOmega<T> Base;

    uwbTMEvaluatorKOmegaSST()
    {
        m_sigmaK1 = 0.85;
        m_sigmaK2 = 1.0;
        m_sigmaO1 = 0.5;
        m_sigmaO2 = 0.856;
        m_beta1 = 0.075;
        m_beta2 = 0.0828;
        m_kappa = 0.41;
        m_a1 = 0.31;
        m_gamma1 = m_beta1 / m_betaStar0 - m_sigmaO1 * math::pow(m_kappa, 2) / math::sqrt(m_betaStar0);
        m_gamma2 = m_beta2 / m_betaStar0 - m_sigmaO2 * math::pow(m_kappa, 2) / math::sqrt(m_betaStar0);
    }

    /// Make function returning a smart pointer
    static uPtr make() { return memory::make_unique(new uwbTMEvaluatorKOmegaSST<T>()); }

protected:

    void initMembers()
    {
        Base::initMembers();

        m_betaStar.setConstant(m_nPoints, m_betaStar0);
        m_F2.setZero(m_nPoints);
        m_argMax1.setZero(m_nPoints);
        m_argMax2.setZero(m_nPoints);
        m_sigmaK.setZero(m_nPoints);
        m_sigmaO.setZero(m_nPoints);
        m_F1.setZero(m_nPoints);
        m_F4.setZero(m_nPoints);
        m_production.setZero(m_nPoints);
        m_sourceQ_SAS.setZero(m_nPoints);
        m_reactionAtRhsO.setZero(m_nPoints);
    }

public:

    void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals, const std::vector<gsMatrix<T> >& solKOmegaGrads)
    {
        // add check of nPoints, if doesn't match, call initMembers with new value?

        Base::initAtElement(solUGrads, solKOVals, solKOmegaGrads);
    }

    void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals)
    {
        // add check of nPoints, if doesn't match, call initMembers with new value?

        Base::initAtElement(solUGrads, solKOVals);
    }

    void evalQuantities_turbViscosity() //evaluate quantities needed to evaluate turbulent visity
    {
        evalVorticityStrainRateMagnitude();
        evalF2();

        evalTurbulentViscosity_into(m_turbViscosity);
    }

    void evalQuantities_nonlinearBlocksPart() //evaluate quantities needed to evaluate nonlinearBlocks in TMvisitors
    {
        evalQuantities_turbViscosity();
        evalConstantsBlending();

        evalKOmegaDiffusionCoefficient();
        evalBlendCoefficient();
    }

    void evalQuantities_rhsPart() //evaluate quantities needed to evaluate source terms of TM equations
    {
        uwbTMEvaluator<T>::evalStrainRateTimesUGrad();
        evalVorticityStrainRateMagnitude();
        evalQuantities_turbViscosity();
        evalConstantsBlending();

        evalKOmegaRhs();
    }

    void evalTurbulentViscosity_into(gsVector<T> & turbViscosityVals)
    {
        turbViscosityVals.resize(m_nPoints);

        GISMO_ASSERT(checkTMEvaluator(), "TMEvaluator is not set.");
        if (m_KOtype == "koSSTMenter2009" || m_KOtype == "koSAS" || m_KOtype == "koSAS_SS" ||
            m_KOtype == "koSAS_SO" || m_KOtype == "koSAS_OO")
        {
            for (int k = 0; k < m_nPoints; k++)
            {
                turbViscosityVals(k) =  m_a1 * m_solKOVals.coeff(0, k) / math::max(m_a1 * math::max(m_solKOVals.coeff(1, k), 1e-15), m_strainRateMagnitude(k) * m_F2(k));
                turbViscosityVals(k) = math::max(turbViscosityVals(k), 1e-15);
            }
        }
        else
        {
            for (int k = 0; k < m_nPoints; k++)
            {
                turbViscosityVals(k) =  m_a1 * m_solKOVals.coeff(0, k) / math::max(m_a1 * math::max(m_solKOVals.coeff(1, k), 1e-15), m_vorticityMagnitude(k) * m_F2(k));
                turbViscosityVals(k) = math::max(turbViscosityVals(k), 1e-15);
            }
        }

        if (m_average)
            this->computeAverage(turbViscosityVals);
    }

    void evalQuantities_diffusionCoeff()
    {
        evalQuantities_turbViscosity();
        evalConstantsBlending();
        evalKOmegaDiffusionCoefficient();
    }

    void evalQuantities_reactionCoeff()
    {
        evalConstantsBlending();
    }

    void evalQuantities_blendCoeff()
    {
        evalF1();
        evalBlendCoefficient();
    }

    void evalTurbVariables()
    {
        evalVorticityStrainRateMagnitude();
        if (checkSASTypeEvaluator())
            evalSourceSASTerm();
        else
        {
            evalF1();
            uwbTMEvaluator<T>::evalStrainRateTimesUGrad();
        }
    }

    void evalF1()
    {
        evalMaxArguments();
        for (int k = 0; k < m_nPoints; k++)
        {
            T CD_ko;
            T tmp1 = math::max(m_solKOVals.coeff(1, k), 1e-15);
            T tmp2 = math::max(m_wallDistanceY(k), 1e-15);
            GISMO_ASSERT(checkTMEvaluator(), "TMEvaluator is not set.");
            // rho devided and 10^(-20)/rho \aprox 10^(-20)/rho as constant close to zero
            if (m_KOtype == "koSSTMenter2009" || m_KOtype == "koSAS" || m_KOtype == "koSAS_SS" ||
                m_KOtype == "koSAS_SO" || m_KOtype == "koSAS_OO")
                CD_ko = math::max(2 * m_sigmaO2 * 1 / tmp1 * m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1)), math::pow(10, -10));
            else
                CD_ko = math::max(2 * m_sigmaO2 * 1 / tmp1 * m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1)), math::pow(10, -20));
            T arg = math::min(math::max(m_argMax1(k), m_argMax2(k)), 4 * m_sigmaO2 * m_solKOVals.coeff(0, k) / (CD_ko * math::pow(tmp2, 2)));
            m_F1(k) = math::tanh(math::pow(arg, 4));
            m_F1(k) = math::max(m_F1(k), 0.);
            m_F1(k) = math::min(m_F1(k), 1.);
        }
    }

    void evalF2()
    {
        evalMaxArguments();
        for (int k = 0; k < m_nPoints; k++)
        {
            T arg = math::max(2 * m_argMax1(k), m_argMax2(k));
            m_F2(k) = math::tanh(math::pow(arg, 2));
            m_F2(k) = math::max(m_F2(k), 0.);
            m_F2(k) = math::min(m_F2(k), 1.);
        }
    }

    void evalMaxArguments()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            T tmp1 = math::max(m_solKOVals.coeff(1, k), 1e-15);
            T tmp2 = math::max(m_wallDistanceY(k), 1e-15);
            T tmp3 = math::max(m_solKOVals.coeff(0, k), 1e-15);
            m_argMax1(k) = math::sqrt(tmp3) / (m_betaStar0 * tmp1 * tmp2);
            m_argMax2(k) = 500 * m_viscosity / (math::pow(tmp2, 2) * tmp1);
        }
    }

    void evalF4()
    {
        evalVorticityStrainRateMagnitude();
        T C_RC = 1.4;
        T Ri;
        for (int k = 0; k < m_nPoints; k++)
        {
            Ri = m_vorticityMagnitude(k) / m_strainRateMagnitude(k) * (m_vorticityMagnitude(k) / m_strainRateMagnitude(k) - 1);
            m_F4(k) = 1 / (1 + C_RC * Ri);
        }
    }

    void evalConstantsBlending()
    {
        evalF1();
        //evalF4();
        for (int k = 0; k < m_nPoints; k++)
        {
            m_sigmaK(k) = m_F1(k) * m_sigmaK1 + (1 - m_F1(k)) * m_sigmaK2;
            m_sigmaO(k) = m_F1(k) * m_sigmaO1 + (1 - m_F1(k)) * m_sigmaO2;
            m_beta(k) = m_F1(k) * m_beta1 + (1 - m_F1(k)) * m_beta2;
            //m_beta(k) = m_F4(k) * (m_F1(k) * m_beta1 + (1 - m_F1(k)) * m_beta2);
            m_gamma(k) = m_F1(k) * m_gamma1 + (1 - m_F1(k)) * m_gamma2;
        }
    }

    void evalKOmegaDiffusionCoefficient()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            m_kDiffusionCoeff(k) = math::max(m_viscosity + m_sigmaK(k) * m_turbViscosity(k), 1e-15);
            m_oDiffusionCoeff(k) = math::max(m_viscosity + m_sigmaO(k) * m_turbViscosity(k), 1e-15);
        }
    }

    void evalKOmegaRhs()
    {
        evalProduction();
        for (int k = 0; k < m_nPoints; k++)
        {
            //========= strain rate approach ===============================
            m_rhsK(k) = m_production(k);
            m_rhsO(k) = m_gamma(k) / math::max(m_turbViscosity(k), 1e-15) * m_production(k);
            //========= vorticity approach ================================
            /*m_rhsK(k) = m_turbViscosity(k) * math::pow(m_vorticityMagnitude(k), 2);
            //production limiter
            m_rhsK(k) = math::min(m_rhsK(k), 20 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            m_rhsO(k) = m_gamma(k) * math::pow(m_vorticityMagnitude(k), 2);*/
        }

        /*if (m_average)
        {
            this->computeAverage(m_rhsK);
            this->computeAverage(m_rhsO);
        }*/

        evalReactionAtRhsO();
    }

    void evalReactionAtRhsO()
    {
        evalBlendCoefficient();
        for (int k = 0; k < m_nPoints; k++)
            m_reactionAtRhsO(k) = math::max(m_blendCoeff(k) * (m_solKOmegaGrads[k].row(0).dot(m_solKOmegaGrads[k].row(1))) / math::max(m_solKOVals.coeff(1, k), 1e-15), 0.);
    }

    void evalProduction()
    {
        GISMO_ASSERT(checkTMEvaluator(), "TMEvaluator is not set.");
        for (int k = 0; k < m_nPoints; k++)
        {
            if (m_KOtype == "koSSTMenter2009") // SST article Liu, Guan, Xu: A production limiter study of SST-SAS turbulence model for bluff body flows
            {
                m_production(k) = m_turbViscosity(k) * math::pow(m_strainRateMagnitude(k), 2);
                m_production(k) = math::min(m_production(k), 10 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            }
            else if (m_KOtype == "koSAS")
            {
                m_production(k) = m_turbViscosity(k) * math::pow(m_strainRateMagnitude(k), 2);
            }
            else if (m_KOtype == "koSAS_SS")
            {
                m_production(k) = m_turbViscosity(k) * math::pow(m_strainRateMagnitude(k), 2);
                m_production(k) = math::min(m_production(k), 10 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            }
            else if (m_KOtype == "koSAS_SO")
            {
                m_production(k) = m_turbViscosity(k) * m_strainRateMagnitude(k) * m_vorticityMagnitude(k);
                m_production(k) = math::min(m_production(k), 10 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            }
            else if (m_KOtype == "koSAS_OO")
            {
                m_production(k) = m_turbViscosity(k) * m_vorticityMagnitude(k) * m_vorticityMagnitude(k);
                m_production(k) = math::min(m_production(k), 10 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            }
            else //koSST
            {
                m_production(k) = 2 * m_turbViscosity(k) * m_strainRateTimesUGrad(k);
                m_production(k) = math::min(m_production(k), 20 * m_betaStar(k) * m_solKOVals.coeff(0, k) * m_solKOVals.coeff(1, k));
            }
        }
        if (checkSASTypeEvaluator())
            evalSourceSASTerm();
    }

    void evalSourceSASTerm()
    {
        T ksi = 3.51;
        T C = 2.;
        T sigma_phi = 0.67;
        T tmpTerm1, tmpTerm2, L_l, L_vk;
        for (int k = 0; k < m_nPoints; k++)
        {
            T tmp = math::max(m_solKOVals.coeff(1, k), 1e-15);
            T tmp2 = math::max(m_solKOVals.coeff(0, k), 1e-15);
            L_l = math::sqrt(tmp2) / (math::pow(m_betaStar0, 0.25) * tmp);
            L_vk = math::max(m_kappa * m_strainRateMagnitude(k) / math::max(m_solULaplaces[k].norm(), 1e-15), 1e-15);
            tmpTerm1 = ksi * m_kappa * math::pow(m_strainRateMagnitude(k), 2) * math::pow(L_l / L_vk, 2);
            tmpTerm2 = C * 2 * m_solKOVals.coeff(0, k) / sigma_phi *
                    math::max(math::pow((m_solKOmegaGrads[k].row(1)).norm(), 2) / math::pow(tmp, 2),
                              math::pow((m_solKOmegaGrads[k].row(0)).norm(), 2) / math::pow(tmp2, 2));
            m_sourceQ_SAS(k) = math::max(tmpTerm1 - tmpTerm2, 0.);
        }
        //if (m_average)
        //    this->computeAverage(m_sourceQ_SAS);
    }

    void evalBlendCoefficient()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            m_blendCoeff(k) = 2 * (1 - m_F1(k)) * m_sigmaO2;
            m_blendCoeff(k) = math::max(m_blendCoeff(k), 0.);
        }
    }

    bool checkTMEvaluator() const
    {
        return (m_KOtype == "koSST" || m_KOtype == "koSSTMenter2009" || m_KOtype == "koSAS"
             || m_KOtype == "koSAS_SS" || m_KOtype == "koSAS_SO" || m_KOtype == "koSAS_OO");
    }

    bool checkSASTypeEvaluator()
    {
        return (m_KOtype == "koSAS" || m_KOtype == "koSAS_SS" ||
                m_KOtype == "koSAS_SO" || m_KOtype == "koSAS_OO");
    }

    const gsVector<T> & getF1() const { return m_F1; }
    const T getF1(const int i) const { return m_F1(i); }
    const gsVector<T> & getSourceQSAS() const { return m_sourceQ_SAS; }
    const T getSourceQSAS(const int i) const { return m_sourceQ_SAS(i); }
    const gsVector<T> & getReactionAtRhs() const { return m_reactionAtRhsO; }
    const T getReactionAtRhs(const int i) const { return m_reactionAtRhsO(i); }

protected:
    T m_sigmaK1;
    T m_sigmaK2;
    T m_sigmaO1;
    T m_sigmaO2;
    T m_beta1;
    T m_beta2;
    T m_kappa;
    T m_a1;
    T m_gamma1;
    T m_gamma2;

    gsVector<T> m_sigmaK;
    gsVector<T> m_sigmaO;

    gsVector<T> m_F2;
    gsVector<T> m_argMax1;
    gsVector<T> m_argMax2;
    gsVector<T> m_F1;
    gsVector<T> m_F4;
    gsVector<T> m_production;
    gsVector<T> m_sourceQ_SAS;
    gsVector<T> m_reactionAtRhsO;

    //memberes from uwbTMEvaluatorKOmega
    using Base::m_solKOVals;
    using Base::m_solKOmegaGrads;
    using Base::m_beta;
    using Base::m_betaStar;
    using Base::m_betaStar0; //denoted as in the LRN k-omega model as the constant betaStarSST = betaStar0
    using Base::m_kDiffusionCoeff;
    using Base::m_oDiffusionCoeff;
    using Base::m_gamma;
    using Base::m_blendCoeff;
    using Base::m_rhsK;
    using Base::m_rhsO;

    using Base::m_KOtype;
    using Base::m_solULaplaces;

    //memberes from uwbTMEvaluator
    using uwbTMEvaluator<T>::m_viscosity;
    using uwbTMEvaluator<T>::m_turbViscosity;
    using uwbTMEvaluator<T>::m_nPoints;
    using uwbTMEvaluator<T>::m_solUGrads;
    using uwbTMEvaluator<T>::m_strainRateTimesUGrad;
    using uwbTMEvaluator<T>::m_wallDistanceY;
    using uwbTMEvaluator<T>::m_vorticityMagnitude;
    using uwbTMEvaluator<T>::m_strainRateMagnitude;
    using uwbTMEvaluator<T>::m_average;

    //functions from uwbTMEvaluator
    using uwbTMEvaluator<T>::evalVorticityStrainRateMagnitude;
};

//-----------------------------------------------------------------------------------------------------------------------------

template <class T>
class uwbTMEvaluator
{
public:
    /// Unique pointer for uwbTMEvaluator
    typedef memory::unique_ptr<uwbTMEvaluator> uPtr;

    uwbTMEvaluator()
    { }

    virtual ~uwbTMEvaluator()
    { }

    void initialize(T viscosity, index_t nPoints)
    {
        m_viscosity = viscosity;
        m_nPoints = nPoints;
        m_turbViscosity.setZero(m_nPoints);
        m_strainRateTimesUGrad.setZero(m_nPoints);
        m_wallDistanceY.setZero(nPoints);
        m_average = true;

        initMembers();
    }

    virtual void initMembers()
    {
        m_vorticityMagnitude.setZero(m_nPoints);
        m_strainRateMagnitude.setZero(m_nPoints);
    }

    virtual void initAtElement(const std::vector<gsMatrix<T> >& solUGrads)
    {
        m_solUGrads = solUGrads;
        m_strainRateTimesUGrad.setZero(m_nPoints);
        m_vorticityMagnitude.setZero(m_nPoints);
        m_strainRateMagnitude.setZero(m_nPoints);
    }

    virtual void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals, const std::vector<gsMatrix<T> >& solKOmegaGrads)
    { GISMO_NO_IMPLEMENTATION }

    virtual void initAtElement(const std::vector<gsMatrix<T> >& solUGrads, const gsMatrix<T>& solKOVals) { GISMO_NO_IMPLEMENTATION }

    virtual void evalTurbulentViscosity_into(gsVector<T> & turbViscosityVals) { GISMO_NO_IMPLEMENTATION }

    void evalStrainRateTimesUGrad()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            // eval strainRate
            T strainRate_ij;
            index_t dim = m_solUGrads[0].cols(); // space dimension

            for (index_t i = 0; i != dim; ++i)
            {
                for (index_t j = 0; j != dim; ++j)
                {
                    strainRate_ij = 0.5 * (m_solUGrads[k].coeff(i, j) + m_solUGrads[k].coeff(j, i));
                    m_strainRateTimesUGrad(k) += strainRate_ij * m_solUGrads[k].coeff(i, j);
                }
            }
        }
    }

    void evalVorticityStrainRateMagnitude()
    {
        for (int k = 0; k < m_nPoints; k++)
        {
            T vorticity_ij, strainRate_ij;
            index_t dim = m_solUGrads[0].cols(); // space dimension

            for (index_t i = 0; i != dim; ++i)
            {
                for (index_t j = 0; j != dim; ++j)
                {
                    vorticity_ij = 0.5 * (m_solUGrads[k].coeff(i, j) - m_solUGrads[k].coeff(j, i));
                    strainRate_ij = 0.5 * (m_solUGrads[k].coeff(i, j) + m_solUGrads[k].coeff(j, i));
                    m_vorticityMagnitude(k) += vorticity_ij * vorticity_ij;
                    m_strainRateMagnitude(k) += strainRate_ij * strainRate_ij;
                }
            }
            m_vorticityMagnitude(k) = math::sqrt(2 * m_vorticityMagnitude(k));
            m_strainRateMagnitude(k) = math::sqrt(2 * m_strainRateMagnitude(k));
        }
    }

    void computeAverage(gsVector<T>& turbFn)
    {
        T turbSumTmp = 0.;
        for (int k = 0; k < m_nPoints; k++)
            turbSumTmp += turbFn(k);
        turbSumTmp = turbSumTmp / m_nPoints;
        turbFn.setConstant(m_nPoints, turbSumTmp);
    }

    void setAveraging(bool average = true) { m_average = average; }
    void setKOmegaVariant(std::string KOtype = "koSST") { m_KOtype = KOtype; }

    const gsVector<T> & getTurbViscosity() const { return m_turbViscosity; }
    const T getTurbViscosity(const int i) const { return m_turbViscosity(i); }
    const gsVector<T> & getVorticityMagnitude() const { return m_vorticityMagnitude; }
    const T getVorticityMagnitude(const int i) const { return m_vorticityMagnitude(i); }
    const gsVector<T> & getStrainRateMagnitude() const { return m_strainRateMagnitude; }
    const T getStrainRateMagnitude(const int i) const { return m_strainRateMagnitude(i); }
    const gsVector<T> & getStrainRateUGrad() const { return m_strainRateTimesUGrad; }
    const T getStrainRateUGrad(const int i) const { return m_strainRateTimesUGrad(i); }

    virtual void evalQuantities_turbViscosity() { GISMO_NO_IMPLEMENTATION }

    virtual void evalQuantities_nonlinearBlocksPart() { GISMO_NO_IMPLEMENTATION }

    virtual void evalQuantities_rhsPart() { GISMO_NO_IMPLEMENTATION }

    virtual void evalQuantities_diffusionCoeff() { GISMO_NO_IMPLEMENTATION }

    virtual void evalQuantities_reactionCoeff() { GISMO_NO_IMPLEMENTATION }

    void evalWallDistance(const gsMatrix<T>& solPoissonVals, const std::vector<gsMatrix<T> >& solPoissonGrads)
    {
        for (int k = 0; k < m_nPoints; k++)
            m_wallDistanceY(k) = -solPoissonGrads[k].norm() + math::sqrt(math::pow(solPoissonGrads[k].norm(), 2) + 2 * solPoissonVals(k));
    }

    const gsVector<T> & getWallDistance() const { return m_wallDistanceY; }
    const T getWallDistance(const int i) const { return m_wallDistanceY(i); }
    virtual const gsVector<T> & getKDiffusionCoefficient() const { GISMO_NO_IMPLEMENTATION }
    virtual const T getKDiffusionCoefficient(const int i) const { GISMO_NO_IMPLEMENTATION }
    virtual const gsVector<T> & getOmegaDiffusionCoefficient() const { GISMO_NO_IMPLEMENTATION }
    virtual const T getOmegaDiffusionCoefficient(const int i) const { GISMO_NO_IMPLEMENTATION }

public:
    /// Make function returning a smart pointer
    static uPtr make()
    {
        return memory::make_unique(new uwbTMEvaluator<T>());
    }

    static uPtr make(std::string evaluatorType)
    {
        if (evaluatorType == "koWilcoxLRN")
            return uwbTMEvaluatorKOmegaWilcoxLRN<real_t>::make();
        else if (evaluatorType == "koSST" || evaluatorType == "koSSTMenter2009" || evaluatorType == "koSAS" ||
                 evaluatorType == "koSAS_SS" || evaluatorType == "koSAS_SO" || evaluatorType == "koSAS_OO")
            return uwbTMEvaluatorKOmegaSST<real_t>::make();
        else
        {
            gsInfo << "Invalid TM evaluator type, using koWilcoxLRN.\n";
            return uwbTMEvaluatorKOmegaWilcoxLRN<real_t>::make();
        }
    }

protected:
    T           m_viscosity;
    index_t     m_nPoints;
    gsVector<T> m_turbViscosity;
    std::vector<gsMatrix<T> > m_solUGrads;
    gsVector<T> m_strainRateTimesUGrad;
    gsVector<T> m_wallDistanceY;
    gsVector<T> m_vorticityMagnitude;
    gsVector<T> m_strainRateMagnitude;
    bool        m_average;
    std::string m_KOtype;
};

} // namespace gismo

