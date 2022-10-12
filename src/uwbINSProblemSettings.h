/** @file uwbINSProblemSettings.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once
#include <gismo.h>

namespace gismo
{

template<class T>
class uwbTMSolverBaseUnsteady;

struct constantsINS
{
    enum intConst { dec_innerIt, unst_innerIt, turb_innerIt, turb_innerFirstIt, iter_maxIt,
                    tauStabTypeSUPG, tauStabTypeCW, tauStabTypeAD, tauStabTypePSPG, tauStabTypeTanhCSD, crosswindType,
                    diamPowerSUPG, diamPowerCW, diamPowerAD, CWresidualType, isoADtype, tanhCSDtype, tanhCSDresidualType,
                    isoADresidualType,
                    int_LAST };
    enum realConst { timeStep, alpha_u, alpha_p, dec_innerTol, unst_innerTol, turb_innerTol, iter_tol, omega, theta,
                     productionXPoint, tanhCSDalpha, isoADalpha,
                     real_LAST };

    enum boolConst { SUPG, TCSD, CROSSWIND, PSPG, RANScrosswind, RANSad, RANSisoAD, tanhCSD,
                     TMsupg, TMafc, TMafcHO, TMcrosswind, TMad, TMisoAD, TMfct_lowOrder,
                     timeDerTerm, limitTMProduction,
                     bool_LAST };
};

struct decoupled
{
    enum method { iterative, coupled, none };
    enum projection { proj1, proj2 };
};

template<class T>
class uwbINSProblemSettings
{

public:
    uwbINSProblemSettings()
    {
        m_intConst.setZero(constantsINS::int_LAST);
        m_realConst.setZero(constantsINS::real_LAST);
        m_boolConst.setZero(constantsINS::bool_LAST);
        m_intConst[constantsINS::dec_innerIt] = 20;
        m_intConst[constantsINS::unst_innerIt] = 3;
        m_intConst[constantsINS::turb_innerIt] = 5;
        m_intConst[constantsINS::turb_innerFirstIt] = 5;
        m_intConst[constantsINS::iter_maxIt] = 500;
        m_intConst[constantsINS::tauStabTypeSUPG] = 0;
        m_intConst[constantsINS::tauStabTypeCW] = 0;
        m_intConst[constantsINS::tauStabTypeAD] = 0;
        m_intConst[constantsINS::tauStabTypePSPG] = 0;
        m_intConst[constantsINS::tauStabTypeTanhCSD] = 0;
        m_intConst[constantsINS::crosswindType] = 0;
        m_intConst[constantsINS::diamPowerSUPG] = 2;
        m_intConst[constantsINS::diamPowerCW] = 2;
        m_intConst[constantsINS::diamPowerAD] = 2;
        m_intConst[constantsINS::CWresidualType] = 1;
        m_intConst[constantsINS::isoADtype] = 6;
        m_intConst[constantsINS::tanhCSDtype] = 0;
        m_intConst[constantsINS::tanhCSDresidualType] = 1;
        m_intConst[constantsINS::isoADresidualType] = 1;
        m_realConst[constantsINS::timeStep] = 0.1;
        m_realConst[constantsINS::alpha_u] = 1;
        m_realConst[constantsINS::alpha_p] = 1;
        m_realConst[constantsINS::dec_innerTol] = 1e-3;
        m_realConst[constantsINS::unst_innerTol] = 1e-4;
        m_realConst[constantsINS::turb_innerTol] = 1e-4;
        m_realConst[constantsINS::iter_tol] = 1e-4;
        m_realConst[constantsINS::theta] = 1.0;
        m_realConst[constantsINS::productionXPoint] = 0.0;
        m_realConst[constantsINS::tanhCSDalpha] = 1.0;
        m_boolConst[constantsINS::SUPG] = false;
        m_boolConst[constantsINS::TCSD] = false;
        m_boolConst[constantsINS::CROSSWIND] = false;
        m_boolConst[constantsINS::PSPG] = false;
        m_boolConst[constantsINS::RANScrosswind] = false;
        m_boolConst[constantsINS::RANSad] = false;
        m_boolConst[constantsINS::RANSisoAD] = false;
        m_boolConst[constantsINS::tanhCSD] = false;
        m_boolConst[constantsINS::TMsupg] = false;
        m_boolConst[constantsINS::TMafc] = false;
        m_boolConst[constantsINS::TMafcHO] = false;
        m_boolConst[constantsINS::TMcrosswind] = false;
        m_boolConst[constantsINS::TMad] = false;
        m_boolConst[constantsINS::TMisoAD] = false;
        m_boolConst[constantsINS::TMfct_lowOrder] = false;
        m_boolConst[constantsINS::timeDerTerm] = true;
        m_boolConst[constantsINS::limitTMProduction] = false;
        m_decMethod = decoupled::none;
        m_projVersion = decoupled::proj1;
        m_pIgaDirGeometry = NULL;
        m_pTurbulenceSolver = NULL;
        m_precondType = "LSC_AdiagEqual";
        m_tmEvaluator = "koWilcoxLRN";
    }

    ~uwbINSProblemSettings()
    {
    }

public:
    void set(constantsINS::intConst name, int value)
    { m_intConst[name] = value; }

    void set(constantsINS::realConst name, T value)
    {
        if (name == constantsINS::theta)
            GISMO_ASSERT(value >= 0.0 && value <= 1.0, "Invalid value of parameter theta.");

        m_realConst[name] = value;
    }

    void set(constantsINS::boolConst name, bool value)
    {
        if ((name == constantsINS::CROSSWIND || name == constantsINS::RANScrosswind || name == constantsINS::PSPG) &&
            value == true)
        {
            GISMO_ASSERT(value == false, "Implementation of chosen stabilization is not finished. Select different type of stabilization.");
            gsWarn << "\nImplementation of chosen stabilization is not finished. Your choice is set to false. Select different type of stabilization.\n";
            m_boolConst[name] = false;
        }
        else
        m_boolConst[name] = value;
    }

    int get(constantsINS::intConst name) const
    { return m_intConst[name]; }

    T get(constantsINS::realConst name) const
    { return m_realConst[name]; }

    bool get(constantsINS::boolConst name) const
    { return m_boolConst[name]; }

    void setPrecondType(std::string name)
    { m_precondType = name; }


    std::string getPrecondType()
    { return m_precondType; }

    bool isRotation() const { return (m_realConst[constantsINS::omega] != 0); }

    void setDecoupledMethod(decoupled::method m) { m_decMethod = m; }
    decoupled::method getDecoupledMethod() { return m_decMethod; }

    void setProjVersion(decoupled::projection p) { m_projVersion = p; }
    decoupled::projection getProjVersion() { return m_projVersion; }

    bool isIgaDirichletGeometry() const { return (m_pIgaDirGeometry != NULL); }

    void setIgaDirichletGeometry(gsGeometry<T> * geom)
    {
        m_pIgaDirGeometry = geom;
    }

    gsGeometry<T>& getIgaDirichletGeometry() const
    {
        GISMO_ASSERT(isIgaDirichletGeometry(), "IgA Dirichlet conditions not set.");
        return *m_pIgaDirGeometry;
    }

    bool isTurbulenceSolver() const { return (m_pTurbulenceSolver != NULL); }

    void setTurbulenceSolver(uwbTMSolverBaseUnsteady<T>* Tsolver)
    {
        m_pTurbulenceSolver = Tsolver;
    }

    void setTurbulenceEvaluator(std::string tmEvaluator = "koWilcoxLRN")
    {
        m_tmEvaluator = tmEvaluator;
    }

    uwbTMSolverBaseUnsteady<T>* getTurbulenceSolver() const
    {
        GISMO_ASSERT(isTurbulenceSolver(), "Turbulence solver not set.");
        return m_pTurbulenceSolver;
    }

    std::string getTMEvaluator() const
    {
        GISMO_ASSERT(checkTMEvaluator(), "TMEvaluator not set. Set koWilcoxLRN or koSST/koSSTMenter2009/koSAS/koSAS_SS/koSAS_SO/koSAS_OO.");
        return m_tmEvaluator;
    }

    bool checkTMEvaluator() const
    {
        return (m_tmEvaluator == "koWilcoxLRN" || m_tmEvaluator == "koSST"
             || m_tmEvaluator == "koSSTMenter2009" || m_tmEvaluator == "koSAS"
             || m_tmEvaluator == "koSAS_SS" || m_tmEvaluator == "koSAS_SO"
             || m_tmEvaluator == "koSAS_OO");
    }

protected:

    // vector of constants
    gsVector<int> m_intConst;
    gsVector<T> m_realConst;
    gsVector<bool> m_boolConst;
    decoupled::method m_decMethod;
    decoupled::projection m_projVersion;
    gsGeometry<T>* m_pIgaDirGeometry;
    uwbTMSolverBaseUnsteady<T>* m_pTurbulenceSolver;
    std::string m_tmEvaluator;
    std::string m_precondType;
}; //class uwbINSProblemSettings

} //namespace gismo
