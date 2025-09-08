/** @file gsINSVisitorsALE.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author: J. Li
*/

#pragma once

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsINSTermsALE.h>
#include <gsIncompressibleFlow/src/gsINSTerms.h>

namespace gismo
{

/// @brief ALE visitor for nonlinear velocity-velocity terms in incompressible Navier-Stokes
/// @tparam T coefficient type
/// @tparam MatOrder sparse matrix storage order
template <class T, int MatOrder = RowMajor>
class gsINSVisitorUUnonlinALE : public gsINSVisitorUU<T, MatOrder>
{
public:
    typedef gsINSVisitorUU<T, MatOrder> Base;
    
protected:
    using Base::m_dofMappers;
    using Base::m_terms;
    
    // Target dimension
    index_t m_tarDim;
    
    // Mesh velocity field
    const gsField<T>* m_meshVelField;
    
    // ALE convection term
    gsINSTerm_ALEConvection<T>* m_aleConvectionTerm;
    
public:
    /// @brief Constructor with paramsPtr (recommended)
    gsINSVisitorUUnonlinALE(typename gsFlowSolverParams<T>::Ptr paramsPtr,
                           const std::vector<gsDofMapper>& dofMappers,
                           index_t targetDim = 2,
                           const gsField<T>* meshVelField = nullptr)
    : Base(paramsPtr), 
      m_tarDim(targetDim),
      m_meshVelField(meshVelField)
    { 
        m_dofMappers = dofMappers;
        
        // Create ALE convection term
        m_aleConvectionTerm = new gsINSTerm_ALEConvection<T>(meshVelField, targetDim);
        
        // Add the term to the visitor
        m_terms.push_back(m_aleConvectionTerm);
    }
    
    /// @brief Constructor without paramsPtr (deprecated, for backward compatibility)
    gsINSVisitorUUnonlinALE(const std::vector<gsDofMapper>& dofMappers,
                           index_t targetDim = 2,
                           const gsField<T>* meshVelField = nullptr)
    : Base(), 
      m_tarDim(targetDim),
      m_meshVelField(meshVelField)
    { 
        m_dofMappers = dofMappers;
        
        // Create ALE convection term
        m_aleConvectionTerm = new gsINSTerm_ALEConvection<T>(meshVelField, targetDim);
        
        // Add the term to the visitor
        m_terms.push_back(m_aleConvectionTerm);
    }
    
    /// @brief Destructor
    ~gsINSVisitorUUnonlinALE()
    {
        // Terms are deleted by base class
    }
    
    /// @brief Set mesh velocity field
    void setMeshVelocityField(const gsField<T>* meshVelField)
    {
        m_meshVelField = meshVelField;
        if (m_aleConvectionTerm)
            m_aleConvectionTerm->setMeshVelocityField(meshVelField);
    }
    
    /// @brief Set current solution
    void setCurrentSolution(const gsField<T>& solution)
    {
        for (size_t i = 0; i < m_terms.size(); ++i)
        {
            gsFlowTermNonlin<T>* nonlinTerm = dynamic_cast<gsFlowTermNonlin<T>*>(m_terms[i]);
            if (nonlinTerm)
                nonlinTerm->setCurrentSolution(const_cast<gsField<T>&>(solution));
        }
    }
};

} // namespace gismo