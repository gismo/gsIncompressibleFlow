/** @file uwbADRProblemSettings.h

Author(s): E. Turnerova
*/

#pragma once
#include <gismo.h>

namespace gismo
{

struct constantsADR
{
    enum intConst { dec_innerIt, unst_innerIt, iter_maxIt, prec_iter, tauStabType, crosswindType, int_LAST };
    enum realConst { timeStep, unst_innerTol, iter_tol, theta, precAL_gamma, real_LAST };
    enum boolConst { unsteady, SUPG, AFC, CROSSWIND, artificialDiffusion, isoArtificialDiffusion, FCT_lowOrder, bool_LAST };
};

struct decoupledMethod
{
    enum method { iterative, coupled, none, proj1, proj2 };
};

template<class T>
class uwbADRProblemSettings
{

public:
    uwbADRProblemSettings()
    {
        m_intConst.setZero(constantsADR::int_LAST);
        m_realConst.setZero(constantsADR::real_LAST);
        m_boolConst.setZero(constantsADR::bool_LAST);
        m_intConst[constantsADR::dec_innerIt] = 20;
        m_intConst[constantsADR::unst_innerIt] = 3;
        m_intConst[constantsADR::iter_maxIt] = 500;
        m_intConst[constantsADR::tauStabType] = 0;
        m_intConst[constantsADR::crosswindType] = 0;
        m_realConst[constantsADR::timeStep] = 0.1;
        m_realConst[constantsADR::unst_innerTol] = 1e-4;
        m_realConst[constantsADR::iter_tol] = 1e-4;
        m_realConst[constantsADR::theta] = 1.0;
        m_realConst[constantsADR::precAL_gamma] = 1.0;
        m_boolConst[constantsADR::unsteady] = false;
        m_boolConst[constantsADR::SUPG] = false;
        m_boolConst[constantsADR::AFC] = false;
        m_boolConst[constantsADR::artificialDiffusion] = false;
        m_boolConst[constantsADR::isoArtificialDiffusion] = false;
        m_boolConst[constantsADR::CROSSWIND] = false;
        m_boolConst[constantsADR::FCT_lowOrder] = false;
        m_pIgaDirGeometry = NULL;
        m_precondType = "MSIMPLER_AdiagEqual";
        m_adrEvaluator = "";
    }

    ~uwbADRProblemSettings() { }

public:
    void set(constantsADR::intConst name, int value)
    { m_intConst[name] = value; }

    void set(constantsADR::realConst name, T value)
    {
        if (name == constantsADR::theta)
            GISMO_ASSERT(value >= 0.0 && value <= 1.0, "Invalid value of parameter theta.");
        if (name == constantsADR::timeStep)
        {
            GISMO_ASSERT(value > 0.0, "Invalid value of parameter timeStep.");
            set(constantsADR::unsteady, true);
        }

        m_realConst[name] = value;
    }

    void set(constantsADR::boolConst name, bool value)
    { m_boolConst[name] = value; }

    int get(constantsADR::intConst name) const
    { return m_intConst[name]; }

    T get(constantsADR::realConst name) const
    { return m_realConst[name]; }

    bool get(constantsADR::boolConst name) const
    { return m_boolConst[name]; }

    void setPrecondType(std::string name)
    { m_precondType = name; }

    std::string getPrecondType()
    { return m_precondType; }

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

    void setADREvaluator(std::string adrEvaluator)
    {
        m_adrEvaluator = adrEvaluator;
    }

    std::string getADREvaluator() const
    {
        GISMO_ASSERT(checkADREvaluator(), "adrEvaluator not set. Set linConstCoeffs/linNonConstCoeffs/nonlinCoeffsField/nonlinCoeffsBurgers.");
        return m_adrEvaluator;
    }

    bool checkADREvaluator() const
    {
        GISMO_ASSERT(m_adrEvaluator != "linNonConstCoeffs", "ADR solver with non-constant coefficients not finished. Choose 'linConstCoeffs'/'nonlinCoeffsField'/'nonlinCoeffsBurgers'.");
        return (m_adrEvaluator == "linConstCoeffs" || m_adrEvaluator == "linNonConstCoeffs" || m_adrEvaluator == "nonlinCoeffsField" || m_adrEvaluator == "nonlinCoeffsBurgers");
    }

protected:

    // vector of constants
    gsVector<int> m_intConst;
    gsVector<T> m_realConst;
    gsVector<bool> m_boolConst;
    gsGeometry<T>* m_pIgaDirGeometry;
    std::string m_precondType;
    std::string m_adrEvaluator;
}; //class uwbADRProblemSettings

} //namespace gismo
