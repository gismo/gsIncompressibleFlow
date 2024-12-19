/** @file gsFlowPeriodicHelper.h
 
    @brief A helper class for periodic conditions in radially symmetric domains.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once

namespace gismo
{

/// @brief A helper class for periodic conditions in radially symmetric domains.
/// @tparam T real number type
template<class T>
class gsFlowPeriodicHelper
{

public: // *** Smart pointers ***

    typedef memory::shared_ptr<gsFlowPeriodicHelper> Ptr;
    typedef memory::unique_ptr<gsFlowPeriodicHelper> uPtr;

protected: // *** Class members ***

    const gsMultiBasis<T>& m_basis;
    const gsDofMapper& m_mapper;
    const gsBoundaryConditions<T>& m_bcInfo;

    index_t m_numDofsFull, m_numDofs;
    std::vector<index_t> m_globalPerDofsFree, m_globalPerDofsElim;
    gsMatrix<index_t> m_map; // index in full basis -> index in reduced basis (with periodic dofs eliminated)
    gsMatrix<index_t> m_invMap; // index in reduced basis -> index in full basis
    gsVector<bool> m_isEliminated;

public: // *** Constructor/destructor ***

    gsFlowPeriodicHelper() {}
    

    gsFlowPeriodicHelper(const gsMultiBasis<T>& basis, const gsDofMapper& mapper, const gsBoundaryConditions<T>& bcInfo):
    m_basis(basis), m_mapper(mapper), m_bcInfo(bcInfo)
    {
       initialize();
    }

    ~gsFlowPeriodicHelper()
    { }


public: // *** Static functions ***

    /// @brief Returns a shared pointer to a newly created instance.
    static Ptr make(const gsMultiBasis<T>& basis, const gsDofMapper& mapper, const gsBoundaryConditions<T>& bcInfo)
    { return memory::make_shared(new gsFlowPeriodicHelper<T>(basis, mapper, bcInfo)); }


protected: // *** Member functions ***

    void initialize();


public: // *** Member functions ***

    inline bool isEliminated(int i) const
    { return m_isEliminated(i); }

    inline index_t map(int i) const
    { return m_map(i); }

    inline index_t invMap(int i) const
    { return m_invMap(i); }


public: // *** Getters/setters ***

    index_t numFreeDofs()
    { return m_numDofs;}


}; // class gsFlowPeriodicHelper

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFlowPeriodicHelper.hpp)
#endif