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
void gsRANSAssemblerUnsteady<T, MatOrder>::initMembers(/*gsTMSolverSST<T, MatOrder>* TMsolver*/)
{
    Base::initMembers();
    updateSizes();

    //m_TMsolver = TMsolver;

    m_visitorRANSsymgraddiag = gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder>(m_paramsPtr/*, m_TMsolver*/);
    m_visitorRANSsymgraddiag.initialize();

    m_visitorRANSsymgradoffdiag = gsRANSVisitorUUSymmetricGradientOffdiag<T, MatOrder>(m_paramsPtr/*, m_TMsolver*/);
    m_visitorRANSsymgradoffdiag.initialize();
    
    // neco s turbulentni viskozitou nebo modelem?
    //gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder>* m_visitorRANSsymgraddiag = dynamic_cast<gsRANSVisitorUUSymmetricGradientDiag<T, MatOrder>*>(m_visitorRANSsymgraddiag);
    m_visitorRANSsymgraddiag.setTurbulenceSolver(m_TMsolver);
    m_visitorRANSsymgradoffdiag.setTurbulenceSolver(m_TMsolver);

    //m_visitorRANSsymgraddiag.setRANSassembler(this);
    //m_visitorRANSsymgradoffdiag.setRANSassembler(this);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::updateSizes()
{
    Base::updateSizes();

    m_matRANSsymgraddiag.resize(m_pshift, m_pshift);
    m_matRANSsymgradoffdiag.resize(m_pshift, m_pshift);
    m_rhsRANS.setZero(m_pshift, 1);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::updateDofMappers()
{
    Base::updateDofMappers();

    m_visitorRANSsymgraddiag.updateDofMappers(m_dofMappers);
    m_visitorRANSsymgradoffdiag.updateDofMappers(m_dofMappers);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaning
    m_matRANSsymgraddiag.resize(m_pshift, m_pshift);
    m_matRANSsymgraddiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgraddiag.outerSize(), m_tarDim * m_nnzPerRowU));

    m_rhsRANS.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorRANSsymgraddiag, 0, m_matRANSsymgraddiag, m_rhsRANS);

    // matrix cleaning
    m_matRANSsymgradoffdiag.resize(m_pshift, m_pshift);
    m_matRANSsymgradoffdiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgradoffdiag.outerSize(), m_tarDim * m_nnzPerRowU));

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorRANSsymgradoffdiag, 0, m_matRANSsymgradoffdiag, dummyRhs);

    m_rhsRANSsymgradoffdiag = m_matRANSsymgradoffdiag * m_solution.topRows(m_pshift);
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_matRANSsymgraddiag);
    m_rhs.topRows(m_pshift) += m_rhsRANS;
    m_rhs.topRows(m_pshift) += m_rhsRANSsymgradoffdiag;
}

template<class T, int MatOrder>
void gsRANSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if(updateSol)
    {
        gsInfo << "Updatine symmetric gradient terms ... ";
        m_matRANSsymgraddiag.resize(m_pshift, m_pshift);
        m_matRANSsymgraddiag.reserve(gsVector<index_t>::Constant(m_matRANSsymgraddiag.outerSize(), m_tarDim * m_nnzPerRowU));
        m_rhsRANS.setZero(m_pshift, 1);
        this->assembleBlock(m_visitorRANSsymgraddiag, 0, m_matRANSsymgraddiag, m_rhsRANS);

        gsMatrix<T> dummyRhs;
        dummyRhs.setZero(m_pshift, 1);
        this->assembleBlock(m_visitorRANSsymgradoffdiag, 0, m_matRANSsymgradoffdiag, dummyRhs);
        m_rhsRANSsymgradoffdiag = m_matRANSsymgradoffdiag * m_solution.topRows(m_pshift);
        gsInfo << "Done" << std::endl;
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();
}

}