/** @file gsINSVisitors.h
    
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

template <class T>
class gsINSVisitor
{

protected: // *** Class members ***

    // zvenku
    index_t m_patchID;
    gsINSSolverParams<T> m_params; // pde, bases, assemblerOptions, options, precOptions
    // pde members: viscosity, rhs (f,g), gsBoundaryConditions, patches, unknownDim    // ulozit pointer/referenci?
 
    // definuje se tady, pak nemenne
    gsQuadRule<T> m_quRule;
    std::vector< gsINSTerm<T> > m_terms;
    unsigned m_geoFlags, m_testFunFlags, m_shapeFunFlags;
    std::vector< gsDofMapper > m_dofMappers;

    // aktualizuje se
    index_t m_currentTestFunID; // updated in evaluate()
    gsMatrix<T> m_localMat;
    gsVector<T> m_quWeights;
    gsMapData<T> m_mapData; // members: points, dim, patchID
    std::vector< gsMatrix<index_t> > m_shapeFunActives;
    std::vector< gsMatrix<T> > m_testFunData; // 0 - value, 1 - deriv, 2 - deriv2
    std::vector< gsMatrix<T> > m_shapeFunData; 
    

public: // *** Constructor/destructor ***

    gsINSVisitor(index_t patchID, const gsINSSolverParams<T>& params) :
    m_patchID(patchID), m_params(params)
    {
        this->defineTerms();
        gatherEvalFlags();
        setupQuadRule();
        m_params.createDofMappers(m_dofMappers);
    }


protected: // *** Member functions ***

// === defineTerms ===
// virtualni, podle params se rozhodne, ktere cleny vytvorit

// === gatherEvalFlags ===
// neni virtualni
// loop over m_terms, m_evFlags |= getEvalFlags
// jak zajistit, aby se to pridalo jen v pripade, ze to tam jeste neni: flag = flag|NEW_FLAG;

// === setupQuadRule ===
// neni virtualni
// gsVector<index_t> numQuadNodes(dim); 
// numQuadNodes[i] = basis.maxDegree() + 1; // take quadrature from highest degree
// m_quRule = gsGaussRule<T>(numQuadNodes);

public: // *** Member functions ***

// === evaluate(testFcnID) ===
// nekde driv spocitat pocet kvad. bodu na element
// support fce, elementy v supportu
// resize m_quNodes, m_quWeights (pocet elem * pocet bodu na elem)
// loop pres elementy, na kazdem spocitat kvad body a vahy, pridat do m_quNodes, m_quWeights
// actives
// ve vsech bodech vyhodnotit, co je potreba (fci s testFcnID nevyhodnocovat dvakrat, pokud je i mezi target fcn) - hodnoty fci a mapData (patch.computeMap(mapData))


// === assemble ===
// loop pres terms, na kazdy assemble(&locMatrix)
// do terms: basisVals, basisGrads

// === localToGlobal ===
// param: const &eliminatedDofs, &globalBlock, &globalRhs
// radek matice zname - testFcnID, mapper pro sloupce

};

} // namespace gismo