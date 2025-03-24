/** @file gsINSAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>

namespace gismo
{

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();
    updateSizes();

    //m_TMsolverPtr = TMsolver;

    m_visitorRANSsymgrad = gsRANSVisitorUUSymmetricGradient<T, MatOrder>(m_paramsPtr);
    m_visitorRANSsymgrad.initialize();

    //m_visitorRANSsymgradoffdiag = gsRANSVisitorUUSymmetricGradientOffdiag<T, MatOrder>(m_paramsPtr);
    //m_visitorRANSsymgradoffdiag.initialize();
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::updateSizes()
{
    Base::updateSizes();

    m_matRANSsymgrad.resize(m_pshift, m_pshift);
    //m_matRANSsymgradoffdiag.resize(m_pshift, m_pshift);
    m_rhsRANS.setZero(m_pshift, 1);

    m_currentFieldU = Base::constructSolution(m_solution, 0);
    m_paramsPtr->setVelocitySolution(m_currentFieldU);
    
    m_oldTimeFieldU = m_currentFieldU;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::updateDofMappers()
{
    Base::updateDofMappers();

    m_visitorRANSsymgrad.updateDofMappers(m_dofMappers);
    //m_visitorRANSsymgradoffdiag.updateDofMappers(m_dofMappers);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    //m_visitorRANSsymgrad.setTurbulenceSolver(m_TMsolverPtr);
    //m_visitorRANSsymgrad.setRANSsolution(m_solution);
    
    // matrix cleaning
    //m_matRANSsymgrad.resize(m_pshift, m_pshift);
    //m_matRANSsymgrad.reserve(gsVector<index_t>::Constant(m_matRANSsymgrad.outerSize(), m_tarDim * m_nnzPerRowU));

    //m_rhsRANS.setZero(m_pshift, 1);

    //this->assembleBlock(m_visitorRANSsymgrad, 0, m_matRANSsymgrad, m_rhsRANS);

    // matrix cleaning
    //m_matRANSsymgradoffdiag.resize(m_pshift, m_pshift);
    //m_matRANSsymgradoffdiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgradoffdiag.outerSize(), m_tarDim * m_nnzPerRowU));

    //gsMatrix<T> dummyRhs;
    //dummyRhs.setZero(m_pshift, 1);

    //this->assembleBlock(m_visitorRANSsymgradoffdiag, 0, m_matRANSsymgradoffdiag, dummyRhs);

    //m_rhsRANSsymgradoffdiag = m_matRANSsymgradoffdiag * m_solution.topRows(m_pshift);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    //real_t sum = 0.0;
    //for (index_t outer = 0; outer < m_matRANSsymgrad.outerSize(); outer++)
    //    for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(m_matRANSsymgrad, outer); it; ++it)
    //        sum += it.value();
    //gsInfo << "m_matRANSsymgrad sum = " << m_matRANSsymgrad.sum() << std::endl;
    //gsInfo << "m_rhsRANS sum = " << m_rhsRANS.sum() << std::endl;

    this->fillGlobalMat_UU(m_matrix, m_matRANSsymgrad);
    m_rhs.topRows(m_pshift) += m_rhsRANS;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    //Base::update(solVector, updateSol);

    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
    {
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

        //gsInfo << "Updating symmetric gradient terms ... ";
        m_visitorRANSsymgrad.setTurbulenceSolver(m_TMsolverPtr);
        m_visitorRANSsymgrad.setRANSsolution(solVector);
        m_matRANSsymgrad.resize(m_pshift, m_pshift);
        m_matRANSsymgrad.reserve(gsVector<index_t>::Constant(m_matRANSsymgrad.outerSize(), m_tarDim * m_nnzPerRowU));
        m_rhsRANS.setZero(m_pshift, 1);
        this->assembleBlock(m_visitorRANSsymgrad, 0, m_matRANSsymgrad, m_rhsRANS);

        //gsMatrix<T> dummyRhs;
        //dummyRhs.setZero(m_pshift, 1);
        //this->assembleBlock(m_visitorRANSsymgradoffdiag, 0, m_matRANSsymgradoffdiag, dummyRhs);
        //m_rhsRANSsymgradoffdiag = m_matRANSsymgradoffdiag * m_solution.topRows(m_pshift);
        //gsInfo << "Done " << m_solution.rows() << ", " << m_solution.cols() << std::endl;
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
    {
        fillSystem();
    }
}

}