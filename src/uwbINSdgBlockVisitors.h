/** @file uwbINSdgBlockVisitors.h

Author(s): H. Hornikova
*/

#pragma once
#include "uwbINSBlockVisitors.h"

namespace gismo {

template <class T>
class uwbINSdgBlockVisitor : public uwbVisitorBase<T>
{   

protected:
    index_t m_dim;
    const T m_penalty;
    const T m_viscosity;
    const gsDofMapper & m_Umap;
    const gsDofMapper & m_Pmap;
    gsField<T> m_solU;
    boxSide m_side1;
    boxSide m_side2;
    int m_patch1;
    int m_patch2;
    gsVector<T> m_unormal;
    
public:

    uwbINSdgBlockVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        m_dim(0),
        m_penalty(penalty),
        m_viscosity(viscosity),
        m_Umap(mappers.front()),
        m_Pmap(mappers.back())
    { }

    void initialize(gsBasisRefs<T> const& basisRefs1,
        boxSide s1,
        boxSide s2,
        const int patchIndex1,
        const int patchIndex2,
        gsQuadRule<T> & rule,
        unsigned & evFlags)
    {   
        const gsBasis<T>& basis1 = basisRefs1.front();

        m_side1 = s1;
        m_side2 = s2;
        m_patch1 = patchIndex1;
        m_patch2 = patchIndex2;

        m_dim = basis1.dim();
        const int dir = m_side1.direction();
        gsVector<int> numQuadNodes(m_dim);
        for (int i = 0; i < basis1.dim(); ++i)
            numQuadNodes[i] = basis1.degree(i) + 1;
        numQuadNodes[dir] = 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        initializeSpecific(evFlags, numQuadNodes.prod());
    };

    virtual inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2)
    { GISMO_NO_IMPLEMENTATION }

    // assemble on element
    virtual inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights)
    { GISMO_NO_IMPLEMENTATION }

    virtual inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs)
    { GISMO_NO_IMPLEMENTATION }

    void setCurrentSolution(std::vector<gsField<T> >& solutions)
    { 
        m_solU = solutions.front();
        this->m_bSolutionSet = true;
    }

protected:
    virtual inline void activesPatch2()
    { GISMO_NO_IMPLEMENTATION }

    virtual inline void initializeSpecific(unsigned & evFlags, int numQuadNodes)
    { GISMO_NO_IMPLEMENTATION }
        
};


template <class T>
class uwbINSdgBlockAVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveU1, m_numActiveU2, m_numActPerNodeU2;
    std::vector<gsMatrix<T> > m_basisDataU1, m_basisDataU2;
    gsMatrix<T> m_physGrad1, m_physGrad2;
    gsMatrix<index_t> m_activesU1, m_activesU2, m_elementsU2;
    gsMatrix<T> K11, K12, K21, K22, M11, M12, M21, M22, N1, N2;

protected:

    using Base::m_dim;
    using Base::m_penalty;
    using Base::m_viscosity;
    using Base::m_Umap;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;
    
public:

    uwbINSdgBlockAVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }

    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);
    
    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);
   
    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);
    
protected:

    inline void activesPatch2();
    
    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSdgBlockBVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveU1, m_numActiveU2, m_numActiveP1, m_numActiveP2, m_numActPerNodeU2, m_numActPerNodeP2;
    gsMatrix<T> m_basisValsU1, m_basisValsU2, m_basisValsP1, m_basisValsP2;
    gsMatrix<index_t> m_activesU1, m_activesP1, m_activesU2, m_activesP2, m_elementsU2, m_elementsP2;
    std::vector<gsMatrix<T> > L11, L12, L21, L22;

protected:

    using Base::m_dim;
    using Base::m_Umap;
    using Base::m_Pmap;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;

public:

    uwbINSdgBlockBVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }


    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

protected:

    inline void activesPatch2();

    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSdgBlockNVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveU1, m_numActiveU2, m_numActPerNodeU2;
    gsMatrix<T> m_basisValsU1, m_basisValsU2;
    gsMatrix<index_t> m_activesU1, m_activesU2, m_elementsU2;
    gsMatrix<T> P11, P12, P21, P22, Q11, Q12, Q21, Q22;
    gsMatrix<T> m_solActUCoeffs1, m_solActUCoeffs2;
    gsMatrix<T> m_solUVals1, m_solUVals2;

protected:

    using Base::m_dim;
    using Base::m_Umap;
    using Base::m_solU;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;

public:

    uwbINSdgBlockNVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }

    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

protected:

    inline void activesPatch2();

    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSdgBlockPVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveP1, m_numActiveP2, m_numActPerNodeP2;
    gsMatrix<T> m_basisValsP1, m_basisValsP2;
    gsMatrix<index_t> m_activesP1, m_activesP2, m_elementsP2;
    gsMatrix<T> P11, P12, P21, P22;

protected:

    using Base::m_dim;
    using Base::m_penalty;
    using Base::m_Pmap;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;

public:

    uwbINSdgBlockPVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }

    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

protected:

    inline void activesPatch2();

    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSdgBlockApVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveP1, m_numActiveP2, m_numActPerNodeP2;
    std::vector<gsMatrix<T> > basisData1_p, basisData2_p;
    gsMatrix<T> m_physGrad1, m_physGrad2;
    gsMatrix<index_t> m_activesP1, m_activesP2, m_elementsP2;
    gsMatrix<T> K11, K12, K21, K22, M11, M12, M21, M22, N1, N2;

protected:

    using Base::m_dim;
    using Base::m_penalty;
    using Base::m_viscosity;
    using Base::m_Pmap;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;

public:

    uwbINSdgBlockApVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }

    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

protected:

    inline void activesPatch2();

    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSdgBlockNpVisitor : public uwbINSdgBlockVisitor<T>
{
public:
    typedef uwbINSdgBlockVisitor<T> Base;

protected:

    index_t m_numActiveU1, m_numActiveU2, m_numActPerNodeU2;
    index_t m_numActiveP1, m_numActiveP2, m_numActPerNodeP2;
    gsMatrix<T> m_basisValsU1, m_basisValsU2, m_basisValsP1, m_basisValsP2;
    gsMatrix<index_t> m_activesU1, m_activesU2, m_elementsU2;
    gsMatrix<index_t> m_activesP1, m_activesP2, m_elementsP2;
    gsMatrix<T> P11, P12, P21, P22, Q11, Q12, Q21, Q22;
    gsMatrix<T> m_solActUCoeffs1, m_solActUCoeffs2;
    gsMatrix<T> m_solUVals1, m_solUVals2;

protected:

    using Base::m_dim;
    using Base::m_Pmap;
    using Base::m_solU;
    using Base::m_side1;
    using Base::m_side2;
    using Base::m_patch1;
    using Base::m_patch2;
    using Base::m_unormal;

public:

    uwbINSdgBlockNpVisitor(const T penalty, const T viscosity, const std::vector<gsDofMapper> & mappers) :
        Base(penalty, viscosity, mappers)
    { }

    inline void evaluate(gsBasisRefs<T> const       & bases1,
        gsBasisRefs<T> const       & bases2,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsMatrix<T>            & quNodes1,
        gsMatrix<T>            & quNodes2);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element1,
        gsGeometryEvaluator<T> & geoEval1,
        gsGeometryEvaluator<T> & geoEval2,
        gsVector<T>            & quWeights);

    inline void localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
        gsSparseMatrix<T>     & blockMatrix,
        gsMatrix<T>           & blockRhs);

protected:

    inline void activesPatch2();

    inline void initializeSpecific(unsigned & evFlags, int numQuadNodes);

};

template <class T>
class uwbINSiFaceSizeVisitor
{

protected:
    index_t m_dim;
    boxSide m_side;
    int m_patch;
    gsVector<T> m_unormal;
    T m_size;

public:

    uwbINSiFaceSizeVisitor() { }

    void initialize(const gsBasis<T> & basis,
        boxSide s,
        const int patchIndex,
        gsQuadRule<T> & rule,
        unsigned & evFlags);

    // assemble on element
    inline void assemble(gsDomainIterator<T>    & element,
        gsGeometryEvaluator<T> & geoEval,
        gsMatrix<T>            & quNodes,
        gsVector<T>            & quWeights);

    T getSize() { return m_size; }

};


}

#include "uwbINSdgBlockVisitors.hpp"

// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(uwbINSdgBlockVisitors.hpp)
// #endif
