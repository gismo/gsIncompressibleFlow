/** @file gsTMModels.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova, B. Bastl
 */

#pragma once

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

/// @brief              A base class for turbulence models.
/// @tparam T           real number type
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

    bool m_average = true;
    bool m_isInitialized = false;


public: // *** Constructor/destructor ***
    
    /// @brief Constructor.
    gsTMModelData(typename gsFlowSolverParams<T>::Ptr paramsPtr) :
    m_paramsPtr(paramsPtr)
    { 
        m_visc = m_paramsPtr->getPde().viscosity();
    }
    
    ~gsTMModelData() {}


public: // *** Static functions ***

    /// @brief Returns a sharedpointer to a newly created instance.
    /// @param[in] paramsPtr a shared point to the instance of an object holding all parameters of the solver
    static tdPtr make(typename gsFlowSolverParams<T>::Ptr paramsPtr);


public: // *** Class functions ***

    /// @brief Evaluates the turbulent viscosity.
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch   
    virtual void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Plots the turbulent viscosity.
    virtual void plotTurbulentViscosity(typename gsFlowSolverParams<T>::Ptr paramsPtr, std::string str = "turbVisc");

    /// @brief Update the current turbuelnce model quantities for the given quNodes
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch 
    virtual void updateModel(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId)
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

/// @brief              A class for the k-omega SST turbulnce model.
/// @tparam T           real number type
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
    using Base::m_average;


public: // *** Constructor/destructor ***
    
    /// @brief Constructor.
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
        m_beta1 = 0.075;
        m_beta2 = 0.0828;
        m_kappa = 0.41;

        m_isInitialized = false;
    }
    
    ~gsTMModelData_SST() {}


public: // *** Static functions ***

    /// @brief Returns a shared pointer to a newly created instance.
    /// @param[in] paramsPtr a shared point to the instance of an object holding all parameters of the solver
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
    
    void evalTurbViscFromData(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId);
    

public: // *** Class functions ***

    /// @brief Update the current turbuelnce model quantities for the given quNodes
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch 
    void updateModel(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId);
    
    /// @brief Evaluates the turbulent viscosity.
    /// @param[in] quNodes          a matrix holding evaluation points
    /// @param[in] numNodesPerElem  number of evaluation points per element
    /// @param[in] patchId          an index of the patch
    void evalTurbulentViscosity(gsMatrix<T>& quNodes, index_t numNodesPerElem, index_t patchId);
    
        
public: // *** Getters/setters ***
    
    
};

// ========================================================================================================



} // namespace gismo