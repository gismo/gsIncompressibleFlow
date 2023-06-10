/** @file gsINSBlockVisitors.h

    @brief Visitors for individual terms in the incompressible Navier-Stokes equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Sourek, H. Hornikova
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsField.h>

namespace gismo
{

// === SUPER CLASS ==== //

/// @brief A base class for incompressible Navier-Stokes element visitors.
/// @tparam T coefficient type
template <class T>
class gsINSVisitorBase
{

protected: // *** Class members ***

    bool m_bSolutionSet; // current solution is already set

public: // *** Constructor/destructor ***

    gsINSVisitorBase()
    {
        m_bSolutionSet = false;
    }

public: // *** Member functions ***

    virtual void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { GISMO_NO_IMPLEMENTATION }

};


// === PARENT === //

/// @brief A parent class for incompressible Navier-Stokes element visitors.
/// @tparam T coefficient type
template <class T>
class gsINSBlockVisitor : public gsINSVisitorBase<T>
{

protected: // *** Class members ***

    const T m_viscosity;
    const gsDofMapper & m_Umap;
    const gsDofMapper & m_Pmap;
    gsField<T> m_solU;
    index_t m_dim;
    int m_patchIndex;
    gsMatrix<T> m_diffusionCoeff;

public: // *** Constructor ***

    gsINSBlockVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        m_viscosity(viscosity),
        m_Umap(dofMappers.front()),
        m_Pmap(dofMappers.back())
    { }

public: // *** Member functions ***

    virtual void initialize(const gsBasisRefs<T>& basisRefs,
        const index_t patchIndex,
        const gsAssemblerOptions& options,
        gsQuadRule<T>& rule,
        unsigned& evFlags);

    virtual void initialize(const gsBasisRefs<T>& basisRefs,
        const index_t patchIndex,
        gsQuadRule<T>& rule,
        unsigned& evFlags);

    virtual void evaluate(const gsBasisRefs<T>& basisRefs,
        gsMatrix<T>            & quNodes)
    { GISMO_NO_IMPLEMENTATION }

    virtual void evaluate(const gsBasisRefs<T>& basisRefs,
        const gsMapData<T>& mapData,
        gsMatrix<T>& quNodes)
    { GISMO_NO_IMPLEMENTATION }

    virtual void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights)
    { GISMO_NO_IMPLEMENTATION }

    virtual void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_solU = solutions.front();
        this->m_bSolutionSet = true;
    }

protected:

    /// Set geometry evaluation flags
    virtual void initializeSpecific(unsigned & evFlags)
    { GISMO_NO_IMPLEMENTATION }
};


// === BLOCK A === //

/// @brief The element visitor for the viscous term of the incompressible N.-S. equations.
/// @tparam T coefficient type
template <class T>
class gsINSBlockAVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_basisGradsU;

    gsMatrix<T> m_physGradU;
    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_viscosity;
    using Base::m_Umap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockAVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCKS B === //

/// @brief The element visitor for the pressure term of the incompressible N.-S. equations.
/// @tparam T coefficient type
template <class T>
class gsINSBlocksBVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU;
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsU;
    gsMatrix<T> m_basisValsP;

    gsMatrix<T> m_physGradU;
    std::vector<gsMatrix<T> > m_localMat;

protected: // *** Base class members ***

    using Base::m_Umap;
    using Base::m_Pmap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlocksBVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }


    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCKS C === //

/// @brief 
/// @tparam T coefficient type
template <class T>
class gsINSBlocksCVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU;
    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsP;
    gsMatrix<T> m_basisValsU;

    gsMatrix<T> m_physGradP;
    std::vector<gsMatrix<T> > m_localMat;

protected: // *** Base class members ***

    using Base::m_Umap;
    using Base::m_Pmap;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;


public: // *** Constructor ***

    gsINSBlocksCVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCK N === //

/// @brief The element visitor for the convective term of the incompressible N.-S. equations.
/// @tparam T coefficient type
template <class T>
class gsINSBlockNVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;
    gsMatrix<T> solActUCoeffs;
    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_physGradU;
    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_Umap;
    using Base::m_solU;
    using Base::m_dim; // Velocity vector dimension
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockNVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCK M === //

/// @brief The element visitor for the velocity mass matrix.
/// @tparam T coefficient type
template <class T>
class gsINSBlockMVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesU;
    std::vector<gsMatrix<T> > m_basisDataU;

    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_Umap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockMVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE;
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCK Ap === //

/// @brief The element visitor for the Poisson operator on pressure space.
/// @tparam T coefficient type
template <class T>
class gsINSBlockApVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesP;
    gsMatrix<T> m_basisGradsP;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_Pmap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockApVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCK Np === //

/// @brief The element visitor for the convective operator on pressure space.
/// @tparam T coefficient type
template <class T>
class gsINSBlockNpVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesP;
    std::vector<gsMatrix<T> > m_basisDataP;
    gsMatrix<T> m_solUVals;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_Pmap;
    using Base::m_solU;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockNpVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE | NEED_GRAD_TRANSFORM;

        GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === BLOCK Mp === //

/// @brief The element visitor for the pressure mass matrix.
/// @tparam T coefficient type
template <class T>
class gsINSBlockMpVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    gsMatrix<index_t> m_activesP;
    gsMatrix<T> bVals_p;

    gsMatrix<T> m_physGradP;
    gsMatrix<T> m_localMat;

protected: // *** Base class members ***

    using Base::m_Pmap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSBlockMpVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity)
    { }

public: // *** Member functions ***

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_MEASURE;
    }

    void evaluate(gsBasisRefs<T> const   & basisRefs,
        gsMatrix<T>            & quNodes);

    void assemble(const gsMapData<T> & mapData,
        const gsVector<T>& quWeights);

    void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & sysBlock,
        gsMatrix<T>           & rhs);
};


// === RHS === //

/// @brief The element visitor for the right-hand side of the incompressible Navier-Stokes equations.
/// @tparam T coefficient type
template <class T>
class gsINSRhsVisitor : public gsINSBlockVisitor<T>
{

public:
    typedef gsINSBlockVisitor<T> Base;

protected: // *** Class members ***

    const gsFunction<T>* m_pRhsFcn;

    gsMatrix<index_t> m_activesU;
    gsMatrix<T> m_basisValsU;
    gsMatrix<T> m_rhsVals;

    gsMatrix<T> m_localRhsU;

protected: // *** Base class members ***

    using Base::m_Umap;
    using Base::m_dim;
    using Base::m_patchIndex;

public: // *** Constructor ***

    gsINSRhsVisitor(std::vector< gsDofMapper > & dofMappers, const T viscosity) :
        Base(dofMappers, viscosity), m_pRhsFcn(NULL)
    { }

public: // *** Member functions ***

    void setRhsFunction(const gsFunction<T>* rhsFcn)
    { m_pRhsFcn = rhsFcn; }

    void initializeSpecific(unsigned & evFlags)
    {
        evFlags = NEED_VALUE | NEED_MEASURE;

        GISMO_ASSERT(m_pRhsFcn != NULL, "No rhs function set in gsINSRhsVisitor!");
    }

    void evaluate(const gsBasisRefs<T>& basisRefs,
        const gsMapData<T>& mapData,
        gsMatrix<T>& quNodes);

    void assemble(const gsMapData<T>& mapData,
        const gsVector<T>& quWeights);

    void localToGlobal( gsMatrix<T> & rhs );
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsINSBlockVisitors.hpp)
#endif