/** @file gsINSTerms.h
    
    @brief 
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova (Hornikova)
 */

#pragma once

namespace gismo
{

// ===================================================================================================================

/// @brief      A class computing individual terms of the weak formulation appearing in incompressible flow problems.
/// @tparam T   coefficient type
template <class T>
class gsINSTerm
{

protected: // *** Class members ***

    unsigned m_geoFlags, m_testFunFlags, m_shapeFunFlags; // evaluation flags
    gsVector<T> m_coeff; // coefficient of the term (size = number of quadrature points)


public: // *** Constructor/destructor ***

    gsINSTerm()
    {
        m_geoFlags = 0;
        m_testFunFlags = 0;
        m_shapeFunFlags = 0;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
    { GISMO_NO_IMPLEMENTATION }

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
    { GISMO_NO_IMPLEMENTATION }

    void updateEvalFlags(unsigned& geoFlags, unsigned& testFunFlags, unsigned& shapeFunFlags)
    { 
        geoFlags |= m_geoFlags;
        testFunFlags |= m_testFunFlags;
        shapeFunFlags |= m_shapeFunFlags;
    }


protected: // *** Member functions ***

    virtual void computeCoeff(const gsMapData<T>& mapData, real_t constValue = 1.0)
    { 
        m_coeff.resize(mapData.points.cols());
        m_coeff.setConstant(constValue);
    }

};

// ===================================================================================================================

/// @brief      
/// @tparam T   coefficient type
template <class T>
class gsINSTermNonlin : public gsINSTerm<T>
{

public:
    typedef gsINSTerm<T> Base;

protected: // *** Class members ***

    gsField<T> m_currentSolU;
    bool m_isCurrentSolSet;
    gsMatrix<T> m_solUVals;


public: // *** Constructor/destructor ***

    gsINSTermNonlin()
    { }


public: // *** Member functions ***

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_currentSolU = solutions.front();
        m_isCurrentSolSet = true;
    }


protected: // *** Member functions ***

    virtual void computeCoeffSolU(const gsMapData<T>& mapData)
    { 
        GISMO_ASSERT(m_isCurrentSolSet, "No velocity solution set in the gsINSTermNonlin visitor.");

        m_solUVals.resize(mapData.dim, mapData.points.cols());
        m_solUVals = m_currentSolU.value(mapData.points, mapData.patchID);
    }

};

// ===================================================================================================================

/// @brief 
/// @tparam T 
template <class T>
class gsINSTermValVal : public gsINSTerm<T>
{

public:
    typedef gsINSTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermValVal()
    {
        Base::m_geoFlags = NEED_MEASURE;
        Base::m_testFunFlags = NEED_VALUE;
        Base::m_shapeFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
    { 
        this->computeCoeff(mapData);

        const gsMatrix<T>& testFunVals = testFunData[0];
        const gsMatrix<T>& shapeFunVals = shapeFunData[0];

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = m_coeff(k) * quWeights(k) * mapData.measure(k);
            localMat += weight * (testFunVals.col(k) * shapeFunVals.col(k).transpose());
        }
    }

};

// ===================================================================================================================

/// @brief 
/// @tparam T 
template <class T>
class gsINSTermGradGrad : public gsINSTerm<T>
{

public:
    typedef gsINSTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermGradGrad()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_DERIV;
        Base::m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, gsMatrix<T>& localMat)
    { 
        this->computeCoeff(mapData);

        const gsMatrix<T>& testFunGrads = testFunData[1];
        const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

        const index_t nQuPoints = quWeights.rows();
        gsMatrix<T> testFunPhysGrad, shapeFunPhysGrad;

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = m_coeff(k) * quWeights(k) * mapData.measure(k);

            transformGradients(mapData, k, testFunGrads, testFunPhysGrad);
            transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

            localMat += weight * (testFunPhysGrad.transpose() * shapeFunPhysGrad);
        }
    }

};

// ===================================================================================================================

/// @brief 
/// @tparam T 
template <class T>
class gsINSTermPvalUdiv : public gsINSTerm<T> // order: shape, test
{

public:
    typedef gsINSTerm<T> Base;

protected: // *** Base class members ***

    using Base::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermPvalUdiv()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_DERIV;
        Base::m_shapeFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
    { 
        this->computeCoeff(mapData);

        const gsMatrix<T>& testFunGrads = testFunData[1];
        const gsMatrix<T>& shapeFunVals = shapeFunData[0];

        gsMatrix<T> testFunPhysGrad;

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = m_coeff(k) * quWeights(k) * mapData.measure(k);

            transformGradients(mapData, k, testFunGrads, testFunPhysGrad);

            for (size_t i = 0; i != localMat.size(); ++i)
                localMat[i].noalias() += weight * (shapeFunVals.col(k) * testFunPhysGrad.row(i));
        }
    }

};

// ===================================================================================================================

/// @brief 
/// @tparam T 
template <class T>
class gsINSTermUsolGradVal : public gsINSTermNonlin<T> // order: shape, test
{

public:
    typedef gsINSTermNonlin<T> Base;

protected: // *** Base class members ***

    using Base::m_solUVals;
    using gsINSTerm<T>::m_coeff;

public: // *** Constructor/destructor ***

    gsINSTermUsolGradVal()
    {
        Base::m_geoFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        Base::m_testFunFlags = NEED_VALUE;
        Base::m_shapeFunFlags = NEED_DERIV;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, const std::vector< gsMatrix<T> >& shapeFunData, std::vector< gsMatrix<T> >& localMat)
    { 
        this->computeCoeff(mapData);
        this->computeCoeffSolU(mapData);

        const gsMatrix<T>& testFunVals = testFunData[0];
        const gsMatrix<T>& shapeFunGrads = shapeFunData[1];

        gsMatrix<T> shapeFunPhysGrad;

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = m_coeff(k) * m_quWeights(k) * mapData.measure(k);

            transformGradients(mapData, k, shapeFunGrads, shapeFunPhysGrad);

            localMat.noalias() += weight * (testFunVals.col(k) * (m_solUVals.col(k).transpose() * shapeFunPhysGrad));
        }
    }

};

// ===================================================================================================================
// ===================================================================================================================

/// @brief      
/// @tparam T   coefficient type
template <class T>
class gsINSTermRhs
{

protected: // *** Class members ***

    unsigned m_geoFlags, m_testFunFlags; // evaluation flags
    const gsFunction<T>* m_pRhsFun;
    gsMatrix<T> m_rhsVals;


public: // *** Constructor/destructor ***

    gsINSTermRhs()
    {
        m_geoFlags = NEED_VALUE | NEED_MEASURE;
        m_testFunFlags = NEED_VALUE;
    }


public: // *** Member functions ***

    virtual void assemble(const gsMapData<T>& mapData, const gsVector<T>& quWeights, const std::vector< gsMatrix<T> >& testFunData, gsMatrix<T>& localRhs)
    { 
        m_pRhsFun->eval_into(mapData.values[0], m_rhsVals);

        const index_t nQuPoints = quWeights.rows();

        for (index_t k = 0; k < nQuPoints; k++)
        {
            const T weight = quWeights(k) * mapData.measure(k);

            localRhs.noalias() += weight * (testFunData[0].col(k) *  m_rhsVals.col(k).transpose());
        }
    }

    void updateEvalFlags(unsigned& geoFlags, unsigned& testFunFlags, unsigned& shapeFunFlags)
    { 
        geoFlags |= m_geoFlags;
        testFunFlags |= m_testFunFlags;
    }

};

} // namespace gismo