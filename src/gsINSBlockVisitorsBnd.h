/** @file gsINSBlockVisitorsBnd.h

    @brief Visitors for the boundary terms in the Navier-Stokes equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s) : H.Hornikova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSBlockVisitors.h>

#include <gsCore/gsFuncData.h>

namespace gismo {

// === PARENT === //

/// @brief A parent class for incompressible Navier-Stokes boundary visitors.
/// @tparam T coefficient type
template <class T>
class gsINSBlockVisitorBnd : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsVector<T> m_unormal;
    boxSide m_side;

protected: // *** Base class members ***

    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockVisitorBnd(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

public: // *** Member functions ***

    virtual void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags,
        boxSide side);

    virtual void assemble(gsDomainIterator<T>& element,
    const gsMapData<T>& mapData,
    const gsVector<T>& quWeights)
    { GISMO_NO_IMPLEMENTATION }

};


// === PCD Robin === //

/// @brief The element visitor for the block arising from Robin boundary conditions in PCD preconditioner.
/// @tparam T coefficient type
template <class T>
class gsINSBlockVisitorRobinPCD : public gsINSBlockVisitorBnd<T>
{

public:
    typedef gsINSBlockVisitorBnd<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_actives;
    gsMatrix<T> m_bVals, m_solUVals, m_localMat;

protected: // *** Base class members ***

    using Base::m_side;
    using Base::m_unormal;
    using gsINSBlockVisitor<T>::m_patchIndex;
    using gsINSBlockVisitor<T>::m_viscosity;
    using gsINSBlockVisitor<T>::m_solU;
    using gsINSBlockVisitor<T>::m_Pmap;

public: // *** Constructor ***

    gsINSBlockVisitorRobinPCD(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    void evaluate(const gsBasisRefs<T>& basisRefs,
        const gsMapData<T>& mapData,
        gsMatrix<T>& quNodes);

    void assemble(gsDomainIterator<T>& element,
        const gsMapData<T>& mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

};


// === NITSCHE PARENT === //

/// @brief A parent class for incompressible Navier-Stokes Nitsche visitors.
/// @tparam T coefficient type
template <class T>
class gsINSBlockVisitorNitsche : public gsINSBlockVisitorBnd<T>
{

public:
    typedef gsINSBlockVisitorBnd<T> Base;

protected: // *** Class members ***

    T m_penalty;
    const gsFunctionSet<T>* m_pDirFcn;

public: // *** Constructor ***

    gsINSBlockVisitorNitsche(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

public: // *** Member functions ***

    void initialize(gsBasisRefs<T> const& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags,
        const boundary_condition<T> & dirCond);

};


// === NITSCHE BLOCK A === //

/// @brief The Nitsche element visitor for the viscous term of the incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template <class T>
class gsINSBlockAVisitorNitsche : public gsINSBlockVisitorNitsche<T>
{

public:
    typedef gsINSBlockVisitorNitsche<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_actives;
    std::vector<gsMatrix<T> > m_basisData;
    gsMatrix<T> m_physGrad;
    gsMatrix<T> m_dirData;

    gsMatrix<T> m_localMat;
    gsMatrix<T> m_localRhs;

protected: // *** Base class members ***

    using Base::m_penalty;
    using Base::m_pDirFcn;
    using gsINSBlockVisitorBnd<T>::m_side;
    using gsINSBlockVisitorBnd<T>::m_unormal;
    using gsINSBlockVisitor<T>::m_patchIndex;
    using gsINSBlockVisitor<T>::m_viscosity;
    using gsINSBlockVisitor<T>::m_Umap;
    using gsINSBlockVisitor<T>::m_dim;

public: // *** Constructor ***

    gsINSBlockAVisitorNitsche(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    void evaluate(const gsBasisRefs<T>& basisRefs,
        const gsMapData<T>& mapData,
        gsMatrix<T>& quNodes);

    void assemble(gsDomainIterator<T>& element,
        const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal( gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

};


// === NITSCHE BLOCK B === //

/// @brief The Nitsche element visitor for the pressure term of the incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template <class T>
class gsINSBlockBVisitorNitsche : public gsINSBlockVisitorNitsche<T>
{

public:
    typedef gsINSBlockVisitorNitsche<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU, m_activesP;
    gsMatrix<T> m_bValsU, m_bValsP;
    gsMatrix<T> m_dirData;
    std::vector<gsMatrix<T> > m_localMatVec;
    gsMatrix<T> m_localRhs;

protected: // *** Base class members ***

    using Base::m_pDirFcn;
    using gsINSBlockVisitorBnd<T>::m_side;
    using gsINSBlockVisitorBnd<T>::m_unormal;
    using gsINSBlockVisitor<T>::m_patchIndex;
    using gsINSBlockVisitor<T>::m_viscosity;
    using gsINSBlockVisitor<T>::m_Umap;
    using gsINSBlockVisitor<T>::m_Pmap;
    using gsINSBlockVisitor<T>::m_dim;

public: // *** Constructor ***

    gsINSBlockBVisitorNitsche(std::vector<gsDofMapper> & mappers, const T viscosity) :
        Base(mappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        // Set Geometry evaluation flags
        evFlags = NEED_VALUE | NEED_JACOBIAN;
    }

    inline void evaluate(const gsBasisRefs<T>& basisRefs,
        const gsMapData<T>& mapData,
        gsMatrix<T>& quNodes);

    inline void assemble(gsDomainIterator<T>& element,
        const gsMapData<T>& mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

}; 

} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSBlockVisitorsBnd.hpp)
#endif

