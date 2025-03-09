

#pragma once

namespace gismo
{

template <class T>
class SSTModel
{
    
protected: // *** Class members ***
    
    // constants
    real_t m_a1 = 0.31;
    real_t m_betaStar = 0.09;
    real_t m_sigmaK1 = 0.85;
    real_t m_sigmaK2 = 1.0;
    real_t m_sigmaO1 = 0.5;
    real_t m_sigmaO2 = 0.856;
    real_t m_beta1 = 3.0/40.0;
    real_t m_beta2 = 0.0828;
    real_t m_kappa = 0.41;
    
    gsMatrix<T> m_KSolVals;
    gsMatrix<T> m_OSolVals;
    std::vector< gsMatrix<T> > m_KSolDers;
    std::vector< gsMatrix<T> > m_OSolDers;
    std::vector< gsMatrix<T> > m_USolDers;
    gsVector<T> m_F1;
    gsVector<T> m_F2;
    gsVector<T> m_TurbulentViscosityVals;
    gsVector<T> m_StrainRateMag;
    std::vector< gsMatrix<T> > m_StrainRateTensor;
    
    bool m_isCurrent = false;
    
public: // *** Constructor/destructor ***
    
    SSTModel() {}
    
    ~SSTModel() {}
        
public: // *** Getters/setters ***
    
    real_t get_a1() { return m_a1; }
    real_t get_betaStar() { return m_betaStar; }
    real_t get_sigmaK1() { return m_sigmaK1; }
    real_t get_sigmaK2() { return m_sigmaK2; }
    real_t get_sigmaO1() { return m_sigmaO1; }
    real_t get_sigmaO2() { return m_sigmaO2; }
    real_t get_beta1() { return m_beta1; }
    real_t get_beta2() { return m_beta2; }
    real_t get_kappa() { return m_kappa; }
    
    gsMatrix<T> getKSolVals() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_KSolVals; 
    }
    gsMatrix<T> getOSolVals() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_OSolVals; 
    }
    std::vector< gsMatrix<T> > getKSolDers() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_KSolDers; 
    }
    std::vector< gsMatrix<T> > getOSolDers() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_OSolDers; 
    }
    std::vector< gsMatrix<T> > getUSolDers() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_USolDers; 
    }
    gsVector<T> getF1Vals() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_F1; 
    }
    gsVector<T> getF2Vals() 
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_F2;
    }
    gsVector<T> getTurbulentViscosityVals()
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_TurbulentViscosityVals;
    }
    gsVector<T> getStrainRateMagVals()
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_StrainRateMag;
    }
    std::vector< gsMatrix<T> > getStrainRateTensor()
    {
        GISMO_ASSERT(m_isCurrent, "Turbulent quiantities not evaluated yet.");
        return m_StrainRateTensor;
    }
    
    void setKSolVals(gsMatrix<T> vals) { m_KSolVals = vals; }
    void setOSolVals(gsMatrix<T> vals) { m_OSolVals = vals; }
    void setKSolDers(std::vector< gsMatrix<T> > vals) { m_KSolDers = vals; }
    void setOSolDers(std::vector< gsMatrix<T> > vals) { m_OSolDers = vals; }
    void setUSolDers(std::vector< gsMatrix<T> > vals) { m_USolDers = vals; }
    void setF1Vals(gsVector<T> vals) { m_F1 = vals; }
    void setF2Vals(gsVector<T> vals) { m_F2 = vals; }
    void setTurbulentViscosityVals(gsVector<T> vals) { m_TurbulentViscosityVals = vals; }
    void setStrainRateMagVals(gsVector<T> vals) { m_StrainRateMag = vals; }
    void StrainRateTensor(std::vector< gsMatrix<T> > vals) { m_StrainRateTensor = vals; }
    
    bool isCurrent() { return m_isCurrent; }
    void setCurrent() { m_isCurrent = true; }
    void setNotCurrent() { m_isCurrent = false; }
    
};

// ========================================================================================================



} // namespace gismo