/** @file gsFlowTerms.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

namespace gismo
{
template <class T>
class gsTMModelData
{
    
public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsTMModelData> tdPtr;

protected: // *** Class members ***

    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;

    real_t m_eps;
    real_t m_a1;
    real_t m_betaStar;
    real_t m_sigmaK1;
    real_t m_sigmaK2;
    real_t m_sigmaO1;
    real_t m_sigmaO2;
    real_t m_beta1;
    real_t m_beta2;
    real_t m_kappa;
    real_t m_visc;
    
    gsMatrix<T> m_KSolVals;
    gsMatrix<T> m_OSolVals;
    std::vector< gsMatrix<T> > m_KSolDers;
    std::vector< gsMatrix<T> > m_OSolDers;
    std::vector< gsMatrix<T> > m_USolDers;
    gsVector<T> m_F1;
    gsVector<T> m_F2;
    gsVector<T> m_StrainRateMag;
    std::vector< gsMatrix<T> > m_StrainRateTensor;
    gsVector<T> m_turbulentViscosityVals;

    bool m_isInitialized = false;


public: // *** Constructor/destructor ***
    
    gsTMModelData(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    m_paramsPtr(paramsPtr)
    { 
        m_visc = m_paramsPtr->getPde().viscosity();
    }
    
    ~gsTMModelData() {}


public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance of the given preconditioner type.
    /// @param[in] precType the reqiured preconditioner type as a string
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner (assuming the following order: NS system matrix, mass matrix (velocity, pressure or both), other matrices)
    /// @param[in] opt a list of options for the preconditioner
    static tdPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr);


public: // *** Class functions ***

    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t patchId)
    { GISMO_NO_IMPLEMENTATION }

    virtual void updateModel(gsMatrix<T>& quNodes, index_t patchId)
    { GISMO_NO_IMPLEMENTATION }


public: // *** Getters/setters ***

    real_t get_eps() { return m_eps; }
    real_t get_a1() { return m_a1; }
    real_t get_betaStar() { return m_betaStar; }
    real_t get_sigmaK1() { return m_sigmaK1; }
    real_t get_sigmaK2() { return m_sigmaK2; }
    real_t get_sigmaO1() { return m_sigmaO1; }
    real_t get_sigmaO2() { return m_sigmaO2; }
    real_t get_beta1() { return m_beta1; }
    real_t get_beta2() { return m_beta2; }
    real_t get_kappa() { return m_kappa; }

    gsMatrix<T> getKSolVals() { return m_KSolVals; }
    gsMatrix<T> getOSolVals() { return m_OSolVals; }
    std::vector< gsMatrix<T> > getKSolDers() { return m_KSolDers; }
    std::vector< gsMatrix<T> > getOSolDers() { return m_OSolDers; }
    std::vector< gsMatrix<T> > getUSolDers() { return m_USolDers; }
    gsVector<T> getF1Vals() { return m_F1; }
    gsVector<T> getF2Vals() { return m_F2; }
    gsVector<T> getStrainRateMagVals() { return m_StrainRateMag; }
    std::vector< gsMatrix<T> > getStrainRateTensor() { return m_StrainRateTensor; }    
    gsVector<T> getTurbulentViscosityVals() { return m_turbulentViscosityVals; }
    bool isInitialized() { return m_isInitialized; }

};


// ===============================================================================================================

template <class T>
class gsTMModelData_SST : public gsTMModelData<T>
{

public:

    typedef gsTMModelData<T> Base;


public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsTMModelData_SST> tdPtr;    

protected: // *** Class members ***

    gsVector<T> m_distance;

protected: // *** Base class members ***    

    using Base::m_paramsPtr;    

    using Base::m_eps;
    using Base::m_a1;
    using Base::m_betaStar;
    using Base::m_sigmaK1;
    using Base::m_sigmaK2;
    using Base::m_sigmaO1;
    using Base::m_sigmaO2;
    using Base::m_beta1;
    using Base::m_beta2;
    using Base::m_kappa;
    using Base::m_visc;

    using Base::m_KSolVals;
    using Base::m_OSolVals;
    using Base::m_KSolDers;
    using Base::m_OSolDers;
    using Base::m_USolDers;
    using Base::m_F1;
    using Base::m_F2;
    using Base::m_StrainRateMag;
    using Base::m_StrainRateTensor;    
    using Base::m_turbulentViscosityVals;
    using Base::m_isInitialized;


public: // *** Constructor/destructor ***
    
    gsTMModelData_SST(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    Base(paramsPtr)
    {
        m_eps = math::pow(10, -15);
        m_a1 = 0.31;
        m_betaStar = 0.09;
        m_sigmaK1 = 0.85;
        m_sigmaK2 = 1.0;
        m_sigmaO1 = 0.5;
        m_sigmaO2 = 0.856;
        m_beta1 = 3.0/40.0;
        m_beta2 = 0.0828;
        m_kappa = 0.41;

        m_isInitialized = false;
    }
    
    ~gsTMModelData_SST() {}


public: // *** Static functions ***

    /// @brief Returns a unique pointer to a newly created instance.
    /// @param[in] mat a const reference to std::map of labeled matrices needed for construction of the preconditioner
    /// @param[in] opt a list of options for the preconditioner
    static tdPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr)
    {
        return memory::make_shared_not_owned(new gsTMModelData_SST<T>(paramsPtr));
    }

protected: // *** Class functions ***    

    void evalDistance(gsMatrix<T>& quNodes, index_t patchId);
    
    void evalVelocityQuantities(gsMatrix<T>& quNodes, index_t patchId);
    
    void evalKSol(gsMatrix<T>& quNodes, index_t patchId, index_t der);
   
    void evalOSol(gsMatrix<T>& quNodes, index_t patchId, index_t der);
   
    void evalF1(gsMatrix<T>& quNodes, index_t patchId);
    
    void evalF2(gsMatrix<T>& quNodes, index_t patchId);
    
    void evalTurbViscFromData(gsMatrix<T>& quNodes, index_t patchId);
    

public: // *** Class functions ***

    void updateModel(gsMatrix<T>& quNodes, index_t patchId);
    
    void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t patchId);
    
        
public: // *** Getters/setters ***
    
    // void setKSolVals(gsMatrix<T> vals) { m_KSolVals = vals; }
    // void setOSolVals(gsMatrix<T> vals) { m_OSolVals = vals; }
    // void setKSolDers(std::vector< gsMatrix<T> > vals) { m_KSolDers = vals; }
    // void setOSolDers(std::vector< gsMatrix<T> > vals) { m_OSolDers = vals; }
    // void setUSolDers(std::vector< gsMatrix<T> > vals) { m_USolDers = vals; }
    // void setF1Vals(gsVector<T> vals) { m_F1 = vals; }
    // void setF2Vals(gsVector<T> vals) { m_F2 = vals; }
    // void setTurbulentViscosityVals(gsVector<T> vals) { m_turbulentViscosityVals = vals; }
    // void setStrainRateMagVals(gsVector<T> vals) { m_StrainRateMag = vals; }
    // void StrainRateTensor(std::vector< gsMatrix<T> > vals) { m_StrainRateTensor = vals; }
        
    //void setCurrent() { m_isCurrent = true; }
    //void setNotCurrent() { m_isCurrent = false; }
    
};

// ========================================================================================================



} // namespace gismo