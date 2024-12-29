/** @file gsFlowPeriodicHelper.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowPeriodicHelper.h>

namespace gismo
{

template <class T>
void gsFlowPeriodicHelper<T>::initialize()
{
    m_numDofsFull = m_mapper.freeSize();

    m_globalPerDofsFree.clear();
    m_globalPerDofsElim.clear();
    m_isEliminated.setZero(m_numDofsFull);
    
    for (typename gsBoundaryConditions<T>::const_ppiterator it = m_bcInfo.periodicBegin(); it != m_bcInfo.periodicEnd(); ++it)
    {
        index_t patch1 = it->first().patch;
        index_t patch2 = it->second().patch;

        // local indices of periodic dofs on both sides
        gsMatrix<index_t> locPerDofs1, locPerDofs2;
        m_basis.basis(patch1).matchWith(*it, m_basis.basis(patch2), locPerDofs1, locPerDofs2);
    
        // global indices of the periodic dofs
        gsMatrix<index_t> globPerDofs1, globPerDofs2;
        m_mapper.localToGlobal(locPerDofs1, patch1, globPerDofs1);
        m_mapper.localToGlobal(locPerDofs2, patch2, globPerDofs2);

        // append the global indices to vectors of periodic dofs
        for (index_t i = 0; i < globPerDofs1.rows(); i++)
        {
            if ((m_mapper.is_free_index(globPerDofs1(i))) && (!m_isEliminated(globPerDofs2(i))) && (!m_isEliminated(globPerDofs1(i))))
            {
                m_globalPerDofsFree.push_back(globPerDofs1(i));
                m_globalPerDofsElim.push_back(globPerDofs2(i));
                m_isEliminated(globPerDofs2(i)) = true;
            }
        }

        m_numDofs = m_numDofsFull - m_globalPerDofsElim.size();

        m_map = gsVector<index_t>::LinSpaced(m_numDofsFull, 0, m_numDofsFull - 1);
        m_invMap = m_map.transpose();
        for (std::vector<index_t>::iterator it = m_globalPerDofsElim.begin(); it != m_globalPerDofsElim.end(); ++it)
        {
            m_invMap.removeCol(m_map(*it));
            m_map.bottomRows(m_numDofsFull - (*it) - 1) -= gsVector<index_t>::Ones(m_numDofsFull - (*it) - 1);
        }
        m_invMap.transposeInPlace();

        for (std::vector<index_t>::iterator it = m_globalPerDofsElim.begin(); it != m_globalPerDofsElim.end(); ++it)
            m_map(*it) = m_map(m_globalPerDofsFree[it - m_globalPerDofsElim.begin()]);
    }
}

} // namespace gismo