/** @file gsFlowBndEvaluators.h
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: H. Honnerova
 */

#pragma once
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>


namespace gismo
{

/**
 * @brief A base class for boundary evaluators.
 * 
 * Used for evaluating some integral quantity defined on a part of a patch boundary, e.g. flow rate.
 * 
 * @tparam T real number type
 */
template<class T>
class gsFlowBndEvaluator
{

protected: // *** Class members ***
    
    typename gsFlowSolverParams<T>::Ptr m_paramsPtr;
    std::vector< std::pair<int, boxSide> > m_bndPart;
    gsField<T> m_velocityField;
    gsField<T> m_pressureField;
    gsMapData<T> m_mapData;
    index_t m_unkID; // the index of variable, from which the quantitiy is computed (0 - velocity, 1 - pressure)
    T m_quantValue; // the resulting integral value


public: // *** Constructor/destructor ***

    /// @brief Constructor.
    /// @param[in] params the parameters of the given problem
    gsFlowBndEvaluator(const gsFlowSolverParams<T>& params):
    gsFlowBndEvaluator(memory::make_shared_not_owned(&params))
    { }

    /// @brief Constructor.
    /// @param[in] paramsPtr smart poiter to the parameters of the given problem
    gsFlowBndEvaluator(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    m_paramsPtr(paramsPtr)
    {
        initMembers();
    }

    /// @brief Constructor.
    /// @param[in] params  the parameters of the given problem
    /// @param[in] bndPart container of pairs (patch, side) defining the boundary part, over which the quantity will be integrated
    gsFlowBndEvaluator(const gsFlowSolverParams<T>& params, const std::vector< std::pair<int, boxSide> >& bndPart):
    gsFlowBndEvaluator(memory::make_shared_not_owned(&params), bndPart)
    { }

    /// @brief Constructor.
    /// @param[in] paramsPtr smart poiter to the parameters of the given problem
    /// @param[in] bndPart   container of pairs (patch, side) defining the boundary part, over which the quantity will be integrated
    gsFlowBndEvaluator(typename gsFlowSolverParams<T>::Ptr paramsPtr, const std::vector< std::pair<int, boxSide> >& bndPart):
    m_paramsPtr(paramsPtr), m_bndPart(bndPart)
    {
        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers()
    {
        m_unkID = 0;
        m_mapData.flags = 0;
        m_quantValue = 0.0;
    }

    /// @brief Evaluate the quantity in one element of the boundary part.
    /// @param[in] patchID      patch index
    /// @param[in] side         patch side
    /// @param[in] quNodes      quadrature nodes for the given element
    /// @param[in] quWeights    quadrature weights for the given element
    virtual void evalOnElement(index_t patchID, boxSide side, const gsMatrix<T>& quNodes, const gsVector<T>& quWeights)
    { GISMO_NO_IMPLEMENTATION }

    /// @brief Evaluate the quantity over one patch side.
    /// @param[in] patchID  patch index
    /// @param[in] side     patch side
    void evalOnPatchSide(index_t patchID, boxSide side);


public: // *** Member functions ***

    /// @brief Evaluate the quantity, i.e., perform the integration over the given boundary part.
    void evaluate()
    {
        m_quantValue = 0.0;

        for(size_t i = 0; i < m_bndPart.size(); i++)
            evalOnPatchSide(m_bndPart[i].first, m_bndPart[i].second);
    }


public: // *** Getters/setters ***

    /// @brief Set the velocity field for evaluation.
    /// @param[in] velocity the new velocity field
    void setVelocityField(const gsField<T>& velocity)
    { m_velocityField = velocity; }

    /// @brief Set the pressure field for evaluation.
    /// @param[in] pressure the new pressure field
    void setPressureField(const gsField<T>& pressure)
    { m_pressureField = pressure; }

    /// @brief Set the velocity and pressure fields for evaluation.
    /// @param[in] velocity the new velocity field
    /// @param[in] pressure the new pressure field
    void setSolutionFields(const gsField<T>& velocity, const gsField<T>& pressure)
    {
        m_velocityField = velocity;
        m_pressureField = pressure;
    }

    /// @brief Set the boundary part, over which the quantity will be integrated.
    /// @param[in] bndPart container of pairs (patch, side) defining the boundary part
    void setBndPart(const std::vector< std::pair<int, boxSide> >& bndPart)
    { m_bndPart = bndPart; }

    /// @brief Get the computed value.
    T getValue()
    { return m_quantValue; }


}; // gsFlowBndEvaluator


// ===================================================================================================================


/// @brief Flow rate evaluator.
/// @tparam T real number type
template<class T>
class gsFlowBndEvaluator_flowRate : public gsFlowBndEvaluator<T>
{

public:
    typedef gsFlowBndEvaluator<T> Base;


protected: // *** Base class members ***
    
    using Base::m_unkID;
    using Base::m_mapData;
    using Base::m_quantValue;
    using Base::m_velocityField;


public: // *** Constructor/destructor ***

    gsFlowBndEvaluator_flowRate(const gsFlowSolverParams<T>& params):
    gsFlowBndEvaluator_flowRate(memory::make_shared_not_owned(&params))
    { }

    gsFlowBndEvaluator_flowRate(typename gsFlowSolverParams<T>::Ptr paramsPtr):
    Base(paramsPtr)
    {
        initMembers();
    }

    gsFlowBndEvaluator_flowRate(const gsFlowSolverParams<T>& params, const std::vector< std::pair<int, boxSide> >& bndPart):
    gsFlowBndEvaluator_flowRate(memory::make_shared_not_owned(&params), bndPart)
    { }

    gsFlowBndEvaluator_flowRate(typename gsFlowSolverParams<T>::Ptr paramsPtr, const std::vector< std::pair<int, boxSide> >& bndPart):
    Base(paramsPtr, bndPart)
    {
        initMembers();
    }


protected: // *** Member functions ***

    /// @brief Initialize the class members.
    void initMembers()
    {
        m_unkID = 0;
        m_mapData.flags = NEED_VALUE | NEED_OUTER_NORMAL;
        m_quantValue = 0.0;
    }

    /// @brief Evaluate the flow rate in one element of the boundary part.
    /// @param[in] patchID      patch index
    /// @param[in] side         patch side
    /// @param[in] quNodes      quadrature nodes for the given element
    /// @param[in] quWeights    quadrature weights for the given element
    virtual void evalOnElement(index_t patchID, boxSide side, const gsMatrix<T>& quNodes, const gsVector<T>& quWeights);

}; // gsFlowBndEvaluator_flowRate


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowBndEvaluators.hpp)
#endif