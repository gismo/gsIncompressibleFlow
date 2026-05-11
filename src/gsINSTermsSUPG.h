/** @file gsINSTermsSUPG.h

    @brief SUPG stabilization terms for incompressible Navier-Stokes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{

/// @brief SUPG streamline diffusion term for incompressible Navier-Stokes:
///        tau_SUPG * (u . grad phi_i) * (u . grad phi_j)
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGConvection : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    /// @param viscosity kinematic viscosity
    /// @param timeStep time step size
    /// @param scaleFactor scaling factor for tau (default 1.0)
    gsINSTerm_SUPGConvection(T viscosity, T timeStep, T scaleFactor = T(1))
    : Base(), m_viscosity(viscosity), m_timeStep(timeStep), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    /// @brief Assemble the SUPG convection term
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (trialFunData.size() < 2 || trialFunData[1].size() == 0)
            return;

        this->computeCoeffSolU(mapData);
        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            // Advection velocity magnitude
            const T aNorm = m_solUVals.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            transformGradients(mapData, k, trialFunData[1], trialPhysGrad);

            // aGradTest  = u^T * grad(phi_i)  -> (1 x nTest)
            // aGradTrial = u^T * grad(phi_j)  -> (1 x nTrial)
            gsMatrix<T> aGradTest  = m_solUVals.col(k).transpose() * testPhysGrad;
            gsMatrix<T> aGradTrial = m_solUVals.col(k).transpose() * trialPhysGrad;

            localMat.noalias() += (weight * tau) * aGradTest.transpose() * aGradTrial;
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGConvection)
};

// ===================================================================================================================

/// @brief SUPG streamline diffusion term for ALE incompressible Navier-Stokes:
///        tau_SUPG * ((u-w) . grad phi_i) * ((u-w) . grad phi_j)
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGConvectionALE : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    index_t m_tarDim;
    gsMatrix<T> m_meshVelVals;
    const gsField<T>* m_meshVelField;
    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    /// @param meshVelField mesh velocity field
    /// @param viscosity kinematic viscosity
    /// @param timeStep time step size
    /// @param tarDim target dimension
    /// @param scaleFactor scaling factor for tau (default 1.0)
    gsINSTerm_SUPGConvectionALE(const gsField<T>* meshVelField,
                                T viscosity,
                                T timeStep,
                                index_t tarDim = 2,
                                T scaleFactor = T(1))
    : Base(), m_meshVelField(meshVelField), m_viscosity(viscosity),
      m_timeStep(timeStep), m_tarDim(tarDim), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    /// @brief Set mesh velocity field
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
    }

    /// @brief Evaluate mesh velocity at quadrature points
    void computeMeshVelocity(const gsMapData<T>& mapData)
    {
        if (m_meshVelField != nullptr)
            m_meshVelVals = m_meshVelField->value(mapData.points, mapData.patchId);
        else
            m_meshVelVals.setZero(m_tarDim, mapData.points.cols());
    }

    /// @brief Assemble the ALE SUPG convection term
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (trialFunData.size() < 2 || trialFunData[1].size() == 0)
            return;

        computeMeshVelocity(mapData);
        this->computeCoeffSolU(mapData);

        // Relative velocity (u - w)
        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;

        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            // Advection velocity magnitude (relative)
            const T aNorm = relativeVel.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            transformGradients(mapData, k, trialFunData[1], trialPhysGrad);

            // aGradTest  = (u-w)^T * grad(phi_i)  -> (1 x nTest)
            // aGradTrial = (u-w)^T * grad(phi_j)  -> (1 x nTrial)
            gsMatrix<T> aGradTest  = relativeVel.col(k).transpose() * testPhysGrad;
            gsMatrix<T> aGradTrial = relativeVel.col(k).transpose() * trialPhysGrad;

            localMat.noalias() += (weight * tau) * aGradTest.transpose() * aGradTrial;
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGConvectionALE)
};

// ===================================================================================================================

/// @brief SUPG pressure gradient term for ALE incompressible Navier-Stokes:
///        tau_SUPG * ((u-w) . grad phi_i) * (d psi_j / d x_d)
///        Uses vector-of-matrices assemble (one per velocity component d).
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGPressureGradALE : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    index_t m_tarDim;
    gsMatrix<T> m_meshVelVals;
    const gsField<T>* m_meshVelField;
    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    gsINSTerm_SUPGPressureGradALE(const gsField<T>* meshVelField,
                                   T viscosity,
                                   T timeStep,
                                   index_t tarDim = 2,
                                   T scaleFactor = T(1))
    : Base(), m_meshVelField(meshVelField), m_viscosity(viscosity),
      m_timeStep(timeStep), m_tarDim(tarDim), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    void setMeshVelocityField(const gsField<T>* meshVelField)
    { m_meshVelField = meshVelField; }

    void computeMeshVelocity(const gsMapData<T>& mapData)
    {
        if (m_meshVelField != nullptr)
            m_meshVelVals = m_meshVelField->value(mapData.points, mapData.patchId);
        else
            m_meshVelVals.setZero(m_tarDim, mapData.points.cols());
    }

    /// @brief Assemble SUPG pressure gradient (vector-of-matrices: one per velocity component)
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         std::vector< gsMatrix<T> >& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;
        if (trialFunData.size() < 2 || trialFunData[1].size() == 0)
            return;

        computeMeshVelocity(mapData);
        this->computeCoeffSolU(mapData);

        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;

        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = relativeVel.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            // testPhysGrad: (dim x nTest) — velocity test function gradients
            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            // trialPhysGrad: (dim x nTrial) — pressure trial function gradients
            transformGradients(mapData, k, trialFunData[1], trialPhysGrad);

            // aGradTest = (u-w)^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = relativeVel.col(k).transpose() * testPhysGrad;

            // localMat[d] += tau * aGradTest^T * (d psi_j / d x_d)
            for (size_t d = 0; d < localMat.size(); ++d)
                localMat[d].noalias() += (weight * tau) * aGradTest.transpose() * trialPhysGrad.row(d);
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGPressureGradALE)
};

// ===================================================================================================================

/// @brief SUPG time discretization term for ALE incompressible Navier-Stokes:
///        (1/dt) * tau_SUPG * ((u-w) . grad phi_i) * phi_j
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGTimeDiscrALE : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    index_t m_tarDim;
    gsMatrix<T> m_meshVelVals;
    const gsField<T>* m_meshVelField;
    T m_viscosity;
    T m_scaleFactor;
    T m_timeStep;

public:
    /// @brief Constructor
    gsINSTerm_SUPGTimeDiscrALE(const gsField<T>* meshVelField,
                                T viscosity,
                                T timeStep,
                                index_t tarDim = 2,
                                T scaleFactor = T(1))
    : Base(), m_meshVelField(meshVelField), m_viscosity(viscosity),
      m_timeStep(timeStep), m_tarDim(tarDim), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_VALUE;
    }

    void setMeshVelocityField(const gsField<T>* meshVelField)
    { m_meshVelField = meshVelField; }

    void computeMeshVelocity(const gsMapData<T>& mapData)
    {
        if (m_meshVelField != nullptr)
            m_meshVelVals = m_meshVelField->value(mapData.points, mapData.patchId);
        else
            m_meshVelVals.setZero(m_tarDim, mapData.points.cols());
    }

    /// @brief Assemble SUPG time discretization (single matrix, UU block)
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;

        computeMeshVelocity(mapData);
        this->computeCoeffSolU(mapData);

        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;

        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = relativeVel.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);

            // aGradTest = (u-w)^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = relativeVel.col(k).transpose() * testPhysGrad;

            // localMat += (1/dt) * tau * aGradTest^T * phi_j^T
            // trialFunData[0].col(k) = phi_j values -> (nTrial x 1)
            localMat.noalias() += (weight * tau / m_timeStep)
                * aGradTest.transpose() * trialFunData[0].col(k).transpose();
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGTimeDiscrALE)
};

// ===================================================================================================================

/// @brief SUPG pressure gradient term for (non-ALE) incompressible Navier-Stokes:
///        tau_SUPG * (u . grad phi_i) * (d psi_j / d x_d)
///        Uses vector-of-matrices assemble (one per velocity component d).
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGPressureGrad : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    gsINSTerm_SUPGPressureGrad(T viscosity, T timeStep, T scaleFactor = T(1))
    : Base(), m_viscosity(viscosity), m_timeStep(timeStep), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV;
    }

    /// @brief Assemble SUPG pressure gradient (vector-of-matrices: one per velocity component)
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         std::vector< gsMatrix<T> >& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;
        if (trialFunData.size() < 2 || trialFunData[1].size() == 0)
            return;

        this->computeCoeffSolU(mapData);
        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = m_solUVals.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            transformGradients(mapData, k, trialFunData[1], trialPhysGrad);

            // aGradTest = u^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = m_solUVals.col(k).transpose() * testPhysGrad;

            // localMat[d] += tau * aGradTest^T * (d psi_j / d x_d)
            for (size_t d = 0; d < localMat.size(); ++d)
                localMat[d].noalias() += (weight * tau) * aGradTest.transpose() * trialPhysGrad.row(d);
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGPressureGrad)
};

// ===================================================================================================================

/// @brief SUPG time discretization term for (non-ALE) incompressible Navier-Stokes:
///        (1/dt) * tau_SUPG * (u . grad phi_i) * phi_j
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGTimeDiscr : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    gsINSTerm_SUPGTimeDiscr(T viscosity, T timeStep, T scaleFactor = T(1))
    : Base(), m_viscosity(viscosity), m_timeStep(timeStep), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_VALUE;
    }

    /// @brief Assemble SUPG time discretization (single matrix, UU block)
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;

        this->computeCoeffSolU(mapData);
        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = m_solUVals.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);

            // aGradTest = u^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = m_solUVals.col(k).transpose() * testPhysGrad;

            // localMat += (1/dt) * tau * aGradTest^T * phi_j^T
            localMat.noalias() += (weight * tau / m_timeStep)
                * aGradTest.transpose() * trialFunData[0].col(k).transpose();
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGTimeDiscr)
};

// ===================================================================================================================

/// @brief SUPG viscous term for (non-ALE) incompressible Navier-Stokes:
///        -tau_SUPG * nu * (u . grad phi_i) * laplacian(phi_j)
///        Non-zero for C^(p-1) IGA splines.
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGViscous : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    gsINSTerm_SUPGViscous(T viscosity, T timeStep, T scaleFactor = T(1))
    : Base(), m_viscosity(viscosity), m_timeStep(timeStep), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV | NEED_DERIV2;
    }

    /// @brief Assemble the SUPG viscous term
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;
        if (trialFunData.size() < 3 || trialFunData[2].size() == 0)
            return;

        this->computeCoeffSolU(mapData);
        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);

        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysLaplace;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = m_solUVals.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            transformLaplaceHgrad(mapData, k, trialFunData[1], trialFunData[2], trialPhysLaplace);

            // aGradTest = u^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = m_solUVals.col(k).transpose() * testPhysGrad;

            // -tau * nu * (u . grad phi_i) * laplacian(phi_j)
            localMat.noalias() -= (weight * tau * m_viscosity) * aGradTest.transpose() * trialPhysLaplace;
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGViscous)
};

// ===================================================================================================================

/// @brief SUPG viscous term for ALE incompressible Navier-Stokes:
///        -tau_SUPG * nu * ((u-w) . grad phi_i) * laplacian(phi_j)
///        Non-zero for C^(p-1) IGA splines.
/// @tparam T coefficient type
template <class T>
class gsINSTerm_SUPGViscousALE : public gsFlowTermNonlin<T>
{
public:
    typedef gsFlowTermNonlin<T> Base;

protected:
    using Base::m_solUVals;

    index_t m_tarDim;
    gsMatrix<T> m_meshVelVals;
    const gsField<T>* m_meshVelField;
    T m_viscosity;
    T m_timeStep;
    T m_scaleFactor;

public:
    /// @brief Constructor
    gsINSTerm_SUPGViscousALE(const gsField<T>* meshVelField,
                              T viscosity,
                              T timeStep,
                              index_t tarDim = 2,
                              T scaleFactor = T(1))
    : Base(), m_meshVelField(meshVelField), m_viscosity(viscosity),
      m_timeStep(timeStep), m_tarDim(tarDim), m_scaleFactor(scaleFactor)
    {
        this->m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
        this->m_testFunFlags = NEED_DERIV;
        this->m_trialFunFlags = NEED_DERIV | NEED_DERIV2;
    }

    void setMeshVelocityField(const gsField<T>* meshVelField)
    { m_meshVelField = meshVelField; }

    void computeMeshVelocity(const gsMapData<T>& mapData)
    {
        if (m_meshVelField != nullptr)
            m_meshVelVals = m_meshVelField->value(mapData.points, mapData.patchId);
        else
            m_meshVelVals.setZero(m_tarDim, mapData.points.cols());
    }

    /// @brief Assemble the ALE SUPG viscous term
    virtual void assemble(const gsMapData<T>& mapData,
                         const gsVector<T>& quWeights,
                         const std::vector< gsMatrix<T> >& testFunData,
                         const std::vector< gsMatrix<T> >& trialFunData,
                         gsMatrix<T>& localMat) override
    {
        if (mapData.points.cols() == 0)
            return;

        if (testFunData.size() < 2 || testFunData[1].size() == 0)
            return;
        if (trialFunData.size() < 3 || trialFunData[2].size() == 0)
            return;

        computeMeshVelocity(mapData);
        this->computeCoeffSolU(mapData);

        gsMatrix<T> relativeVel = m_solUVals - m_meshVelVals;

        gsVector<T> coeffMeasure = this->getCoeffGeoMapProduct(mapData);
        const index_t nQuPoints = quWeights.rows();
        const index_t dim = mapData.dim.first;

        // Element size from physical element volume
        T elemVol = T(0);
        for (index_t k = 0; k < nQuPoints; ++k)
            elemVol += quWeights(k) * mapData.measure(k);
        const T h = math::pow(elemVol, T(1) / T(dim));

        gsMatrix<T> testPhysGrad, trialPhysLaplace;

        for (index_t k = 0; k < nQuPoints; ++k)
        {
            const T weight = quWeights(k) * coeffMeasure(k);

            const T aNorm = relativeVel.col(k).norm();

            // Stabilization parameter (Tezduyar 1992)
            const T inv_time = T(2) / m_timeStep;
            const T inv_conv = T(2) * aNorm / h;
            const T inv_diff = T(4) * m_viscosity / (h * h);
            const T tau = m_scaleFactor / math::sqrt(inv_time*inv_time + inv_conv*inv_conv + inv_diff*inv_diff);

            transformGradients(mapData, k, testFunData[1], testPhysGrad);
            transformLaplaceHgrad(mapData, k, trialFunData[1], trialFunData[2], trialPhysLaplace);

            // aGradTest = (u-w)^T * grad(phi_i)  -> (1 x nTest)
            gsMatrix<T> aGradTest = relativeVel.col(k).transpose() * testPhysGrad;

            // -tau * nu * ((u-w) . grad phi_i) * laplacian(phi_j)
            localMat.noalias() -= (weight * tau * m_viscosity) * aGradTest.transpose() * trialPhysLaplace;
        }
    }

    GISMO_CLONE_FUNCTION(gsINSTerm_SUPGViscousALE)
};

} // namespace gismo
