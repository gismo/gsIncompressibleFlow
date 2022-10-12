/** @file uwbINSdgBlockVisitors.hpp

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSdgBlockVisitors.h"

namespace gismo {

template <class T>
inline void uwbINSdgBlockAVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM;

    m_elementsU2.setZero(1, numQuadNodes);
}

template <class T>
inline void uwbINSdgBlockBVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN;

    L11.resize(m_dim);
    L12.resize(m_dim);
    L21.resize(m_dim);
    L22.resize(m_dim);

    m_elementsU2.setZero(1, numQuadNodes);
    m_elementsP2.setZero(1, numQuadNodes);
}

template <class T>
inline void uwbINSdgBlockNVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN;

    m_elementsU2.setZero(1, numQuadNodes);

    GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
}

template <class T>
inline void uwbINSdgBlockPVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN;

    m_elementsP2.setZero(1, numQuadNodes);
}

template <class T>
inline void uwbINSdgBlockApVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM;

    m_elementsP2.setZero(1, numQuadNodes);
}

template <class T>
inline void uwbINSdgBlockNpVisitor<T>::initializeSpecific(unsigned & evFlags, int numQuadNodes)
{
    // Set Geometry evaluation flags
    evFlags = NEED_VALUE | NEED_JACOBIAN;

    m_elementsU2.setZero(1, numQuadNodes);
    m_elementsP2.setZero(1, numQuadNodes);

    GISMO_ASSERT(this->m_bSolutionSet, "No velocity solution set in the visitor.");
}

template <class T>
inline void uwbINSdgBlockAVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    // Compute the active basis functions
    bases1.front().active_into(quNodes1.col(0), m_activesU1); // the active functions are the same for all quNodes on patch1
    bases2.front().active_into(quNodes2, m_activesU2); // the active functions can differ on patch2
        
    m_numActPerNodeU2 = m_activesU2.rows();

    // Prepare actives2 vector
    activesPatch2();

    m_numActiveU1 = m_activesU1.rows();
    m_numActiveU2 = m_activesU2.rows();

    // Evaluate basis functions and their first derivatives
    bases1.front().evalAllDers_into(quNodes1, 1, m_basisDataU1);
    bases2.front().evalAllDers_into(quNodes2, 1, m_basisDataU2);

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Initialize local matrices
    K11.setZero(m_numActiveU1, m_numActiveU1); K12.setZero(m_numActiveU1, m_numActiveU2);
    K22.setZero(m_numActiveU2, m_numActiveU2); K21.setZero(m_numActiveU2, m_numActiveU1);
    M11.setZero(m_numActiveU1, m_numActiveU1); M12.setZero(m_numActiveU1, m_numActiveU2);
    M22.setZero(m_numActiveU2, m_numActiveU2); M21.setZero(m_numActiveU2, m_numActiveU1);

}

template <class T>
inline void uwbINSdgBlockBVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    // Reducing the dimension of quadrature nodes to the interface dimension
    gsMatrix<T> bQuNodes1 = quNodes1.transpose();
    gsMatrix<T> bQuNodes2 = quNodes2.transpose();
    bQuNodes1.removeCol(m_side1.direction());
    bQuNodes2.removeCol(m_side2.direction());
    bQuNodes1.transposeInPlace();
    bQuNodes2.transposeInPlace();

    // Compute the active basis functions
    gsMatrix<index_t> boundaryBasisIDs1_u = bases1.front().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_u = bases2.front().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> boundaryBasisIDs1_p = bases1.back().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_p = bases2.back().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> tmp_actives1_u, tmp_actives2_u, tmp_actives1_p, tmp_actives2_p; // local indices of active functions in boundaryBasis
    bases1.front().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_u); // the active functions are the same for all quNodes on patch1
    bases2.front().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_u); // the active functions can differ on patch2
    bases1.back().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_p); // the active functions are the same for all quNodes on patch1
    bases2.back().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_p); // the active functions can differ on patch2

    std::vector<unsigned> tmpAct1U_vec, tmpAct1P_vec;

    for (index_t i = 0; i < tmp_actives1_u.rows(); i++)
        tmpAct1U_vec.push_back(tmp_actives1_u(i,0));

    for (index_t i = 0; i < tmp_actives1_p.rows(); i++)
        tmpAct1P_vec.push_back(tmp_actives1_p(i,0));

    // Global indices of active functions
    boundaryBasisIDs1_u.submatrixRows(tmpAct1U_vec, m_activesU1);
    boundaryBasisIDs1_p.submatrixRows(tmpAct1P_vec, m_activesP1);
    m_activesU2.setZero(tmp_actives2_u.rows(), tmp_actives2_u.cols());
    m_activesP2.setZero(tmp_actives2_p.rows(), tmp_actives2_p.cols());

    for (index_t j = 0; j < tmp_actives2_u.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_u.rows(); i++)
            m_activesU2(i, j) = boundaryBasisIDs2_u(tmp_actives2_u(i, j));

    for (index_t j = 0; j < tmp_actives2_p.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_p.rows(); i++)
            m_activesP2(i, j) = boundaryBasisIDs2_p(tmp_actives2_p(i, j));

    // Prepare actives2 vector and some info about the elements on patch2
    activesPatch2();

    m_numActPerNodeU2 = tmp_actives2_u.rows();
    m_numActPerNodeP2 = tmp_actives2_p.rows();
    m_numActiveU1 = m_activesU1.rows();
    m_numActiveU2 = m_activesU2.rows();
    m_numActiveP1 = m_activesP1.rows();
    m_numActiveP2 = m_activesP2.rows();

    // Evaluate basis functions
    m_basisValsU1.setZero(m_numActiveU1, quNodes1.cols());
    m_basisValsU2.setZero(m_numActPerNodeU2, quNodes2.cols());
    m_basisValsP1.setZero(m_numActiveP1, quNodes1.cols());
    m_basisValsP2.setZero(m_numActPerNodeP2, quNodes2.cols());

    for (index_t i = 0; i < m_numActiveU1; i++)
        m_basisValsU1.row(i) = bases1.front().evalSingle(m_activesU1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeU2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsU2.row(i).col(j) = bases2.front().evalSingle(m_activesU2(m_elementsU2(j) + i), quNodes2.col(j));

    for (index_t i = 0; i < m_numActiveP1; i++)
        m_basisValsP1.row(i) = bases1.back().evalSingle(m_activesP1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeP2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsP2.row(i).col(j) = bases2.back().evalSingle(m_activesP2(m_elementsP2(j) + i), quNodes2.col(j));

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Initialize local matrices
    for (index_t i = 0; i != m_dim; i++)
    {
        L11[i].setZero(m_numActiveP1, m_numActiveU1);
        L12[i].setZero(m_numActiveP2, m_numActiveU1);
        L21[i].setZero(m_numActiveP1, m_numActiveU2);
        L22[i].setZero(m_numActiveP2, m_numActiveU2);
    }
}


template <class T>
inline void uwbINSdgBlockNVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    // Reducing the dimension of quadrature nodes to the interface dimension
    gsMatrix<T> bQuNodes1 = quNodes1.transpose();
    gsMatrix<T> bQuNodes2 = quNodes2.transpose();
    bQuNodes1.removeCol(m_side1.direction());
    bQuNodes2.removeCol(m_side2.direction());
    bQuNodes1.transposeInPlace();
    bQuNodes2.transposeInPlace();

    // Compute the active basis functions
    gsMatrix<index_t> boundaryBasisIDs1_u = bases1.front().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_u = bases2.front().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> tmp_actives1_u, tmp_actives2_u; // local indices of active functions in boundaryBasis
    bases1.front().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_u); // the active functions are the same for all quNodes on patch1
    bases2.front().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_u); // the active functions can differ on patch2

    std::vector<unsigned> tmpAct1U_vec;
    for (index_t i = 0; i < tmp_actives1_u.rows(); i++)
        tmpAct1U_vec.push_back(tmp_actives1_u(i,0));

    // Global indices of active functions
    boundaryBasisIDs1_u.submatrixRows(tmpAct1U_vec, m_activesU1);
    m_activesU2.setZero(tmp_actives2_u.rows(), tmp_actives2_u.cols());

    for (index_t j = 0; j < tmp_actives2_u.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_u.rows(); i++)
            m_activesU2(i, j) = boundaryBasisIDs2_u(tmp_actives2_u(i, j));

    // Prepare actives2 vector and some info about the elements on patch2
    activesPatch2();

    m_numActPerNodeU2 = tmp_actives2_u.rows();
    m_numActiveU1 = m_activesU1.rows();
    m_numActiveU2 = m_activesU2.rows();

    // Evaluate basis functions
    m_basisValsU1.setZero(m_numActiveU1, quNodes1.cols());
    m_basisValsU2.setZero(m_numActPerNodeU2, quNodes2.cols());

    for (index_t i = 0; i < m_numActiveU1; i++)
        m_basisValsU1.row(i) = bases1.front().evalSingle(m_activesU1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeU2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsU2.row(i).col(j) = bases2.front().evalSingle(m_activesU2(m_elementsU2(j) + i), quNodes2.col(j));

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Initialize local matrices
    P11.setZero(m_numActiveU1, m_numActiveU1); P12.setZero(m_numActiveU1, m_numActiveU2);
    P22.setZero(m_numActiveU2, m_numActiveU2); P21.setZero(m_numActiveU2, m_numActiveU1);
    Q11.setZero(m_numActiveU1, m_numActiveU1); Q12.setZero(m_numActiveU1, m_numActiveU2);
    Q22.setZero(m_numActiveU2, m_numActiveU2); Q21.setZero(m_numActiveU2, m_numActiveU1);

    // Evaluate solution on element nodes
    m_solActUCoeffs1.setZero(m_dim, m_numActiveU1);
    m_solActUCoeffs2.setZero(m_dim, m_numActiveU2);
    m_solUVals2.setZero(m_dim, m_basisValsU2.cols());

    for (int j = 0; j < m_numActiveU1; j++)
    {
        m_solActUCoeffs1.col(j) = m_solU.coefficientVector(m_patch1).row(m_activesU1(j)).transpose();
    }
    for (int j = 0; j < m_numActiveU2; j++)
    {
        m_solActUCoeffs2.col(j) = m_solU.coefficientVector(m_patch2).row(m_activesU2(j)).transpose();
    }

    m_solUVals1.noalias() = m_solActUCoeffs1 * m_basisValsU1;
    for (int k = 0; k < m_basisValsU2.cols(); k++)
    {
        m_solUVals2.col(k) = m_solActUCoeffs2.block(0, m_elementsU2(k), m_dim, m_numActPerNodeU2) * m_basisValsU2.col(k);
    }
}

template <class T>
inline void uwbINSdgBlockPVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    // Reducing the dimension of quadrature nodes to the interface dimension
    gsMatrix<T> bQuNodes1 = quNodes1.transpose();
    gsMatrix<T> bQuNodes2 = quNodes2.transpose();
    bQuNodes1.removeCol(m_side1.direction());
    bQuNodes2.removeCol(m_side2.direction());
    bQuNodes1.transposeInPlace();
    bQuNodes2.transposeInPlace();

    // Compute the active basis functions
    gsMatrix<index_t> boundaryBasisIDs1_p = bases1.back().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_p = bases2.back().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> tmp_actives1_p, tmp_actives2_p; // local indices of active functions in boundaryBasis
    bases1.back().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_p); // the active functions are the same for all quNodes on patch1
    bases2.back().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_p); // the active functions can differ on patch2

    std::vector<unsigned> tmpAct1P_vec;
    for (index_t i = 0; i < tmp_actives1_p.rows(); i++)
        tmpAct1P_vec.push_back(tmp_actives1_p(i,0));

    // Global indices of active functions
    boundaryBasisIDs1_p.submatrixRows(tmpAct1P_vec, m_activesP1);
    m_activesP2.setZero(tmp_actives2_p.rows(), tmp_actives2_p.cols());

    for (index_t j = 0; j < tmp_actives2_p.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_p.rows(); i++)
            m_activesP2(i, j) = boundaryBasisIDs2_p(tmp_actives2_p(i, j));

    // Prepare actives2 vector and some info about the elements on patch2
    activesPatch2();

    m_numActPerNodeP2 = tmp_actives2_p.rows();
    m_numActiveP1 = m_activesP1.rows();
    m_numActiveP2 = m_activesP2.rows();

    // Evaluate basis functions
    m_basisValsP1.setZero(m_numActiveP1, quNodes1.cols());
    m_basisValsP2.setZero(m_numActPerNodeP2, quNodes2.cols());

    for (index_t i = 0; i < m_numActiveP1; i++)
        m_basisValsP1.row(i) = bases1.back().evalSingle(m_activesP1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeP2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsP2.row(i).col(j) = bases2.back().evalSingle(m_activesP2(m_elementsP2(j) + i), quNodes2.col(j));

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Pressure penalty
    P11.setZero(m_numActiveP1, m_numActiveP1); P12.setZero(m_numActiveP1, m_numActiveP2);
    P22.setZero(m_numActiveP2, m_numActiveP2); P21.setZero(m_numActiveP2, m_numActiveP1);
}

template <class T>
inline void uwbINSdgBlockApVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    
    bases1.back().active_into(quNodes1.col(0), m_activesP1); // the active functions are the same for all quNodes on patch1
    bases2.back().active_into(quNodes2, m_activesP2); // the active functions can differ on patch2

    m_numActPerNodeP2 = m_activesP2.rows();

    // Prepare actives2 vector and some info about the elements on patch2
    activesPatch2();

    m_numActiveP1 = m_activesP1.rows();
    m_numActiveP2 = m_activesP2.rows();

    bases1.back().evalAllDers_into(quNodes1, 1, basisData1_p);
    bases2.back().evalAllDers_into(quNodes2, 1, basisData2_p);

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Initialize local matrices
    K11.setZero(m_numActiveP1, m_numActiveP1); K12.setZero(m_numActiveP1, m_numActiveP2);
    K22.setZero(m_numActiveP2, m_numActiveP2); K21.setZero(m_numActiveP2, m_numActiveP1);
    M11.setZero(m_numActiveP1, m_numActiveP1); M12.setZero(m_numActiveP1, m_numActiveP2);
    M22.setZero(m_numActiveP2, m_numActiveP2); M21.setZero(m_numActiveP2, m_numActiveP1);
}

template <class T>
inline void uwbINSdgBlockNpVisitor<T>::evaluate(gsBasisRefs<T> const       & bases1,
    gsBasisRefs<T> const       & bases2,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsMatrix<T>            & quNodes1,
    gsMatrix<T>            & quNodes2)
{
    // Reducing the dimension of quadrature nodes to the interface dimension
    gsMatrix<T> bQuNodes1 = quNodes1.transpose();
    gsMatrix<T> bQuNodes2 = quNodes2.transpose();
    bQuNodes1.removeCol(m_side1.direction());
    bQuNodes2.removeCol(m_side2.direction());
    bQuNodes1.transposeInPlace();
    bQuNodes2.transposeInPlace();

    // Compute the active basis functions
    gsMatrix<index_t> boundaryBasisIDs1_u = bases1.front().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_u = bases2.front().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> boundaryBasisIDs1_p = bases1.back().boundary(m_side1); // indices of basis functions on m_side1
    gsMatrix<index_t> boundaryBasisIDs2_p = bases2.back().boundary(m_side2); // indices of basis functions on m_side2
    gsMatrix<index_t> tmp_actives1_u, tmp_actives2_u, tmp_actives1_p, tmp_actives2_p; // local indices of active functions in boundaryBasis
    bases1.front().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_u); // the active functions are the same for all quNodes on patch1
    bases2.front().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_u); // the active functions can differ on patch2
    bases1.back().boundaryBasis(m_side1)->active_into(bQuNodes1.col(0), tmp_actives1_p); // the active functions are the same for all quNodes on patch1
    bases2.back().boundaryBasis(m_side2)->active_into(bQuNodes2, tmp_actives2_p); // the active functions can differ on patch2

    std::vector<unsigned> tmpAct1U_vec, tmpAct1P_vec;

    for (index_t i = 0; i < tmp_actives1_u.rows(); i++)
        tmpAct1U_vec.push_back(tmp_actives1_u(i,0));

    for (index_t i = 0; i < tmp_actives1_p.rows(); i++)
        tmpAct1P_vec.push_back(tmp_actives1_p(i,0));

    // Global indices of active functions
    boundaryBasisIDs1_u.submatrixRows(tmpAct1U_vec, m_activesU1);
    m_activesU2.setZero(tmp_actives2_u.rows(), tmp_actives2_u.cols());
    boundaryBasisIDs1_p.submatrixRows(tmpAct1P_vec, m_activesP1);
    m_activesP2.setZero(tmp_actives2_p.rows(), tmp_actives2_p.cols());

    for (index_t j = 0; j < tmp_actives2_u.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_u.rows(); i++)
            m_activesU2(i, j) = boundaryBasisIDs2_u(tmp_actives2_u(i, j));

    for (index_t j = 0; j < tmp_actives2_p.cols(); j++)
        for (index_t i = 0; i < tmp_actives2_p.rows(); i++)
            m_activesP2(i, j) = boundaryBasisIDs2_p(tmp_actives2_p(i, j));

    // Prepare actives2 vector and some info about the elements on patch2
    activesPatch2();

    m_numActPerNodeU2 = tmp_actives2_u.rows();
    m_numActiveU1 = m_activesU1.rows();
    m_numActiveU2 = m_activesU2.rows();
    m_numActPerNodeP2 = tmp_actives2_p.rows();
    m_numActiveP1 = m_activesP1.rows();
    m_numActiveP2 = m_activesP2.rows();

    // Evaluate basis functions
    m_basisValsU1.setZero(m_numActiveU1, quNodes1.cols());
    m_basisValsU2.setZero(m_numActPerNodeU2, quNodes2.cols());
    m_basisValsP1.setZero(m_numActiveP1, quNodes1.cols());
    m_basisValsP2.setZero(m_numActPerNodeP2, quNodes2.cols());

    for (index_t i = 0; i < m_numActiveU1; i++)
        m_basisValsU1.row(i) = bases1.front().evalSingle(m_activesU1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeU2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsU2.row(i).col(j) = bases2.front().evalSingle(m_activesU2(m_elementsU2(j) + i), quNodes2.col(j));

    for (index_t i = 0; i < m_numActiveP1; i++)
        m_basisValsP1.row(i) = bases1.back().evalSingle(m_activesP1(i), quNodes1);

    for (index_t i = 0; i < m_numActPerNodeP2; i++)
        for (index_t j = 0; j < quNodes2.cols(); j++)
            m_basisValsP2.row(i).col(j) = bases2.back().evalSingle(m_activesP2(m_elementsP2(j) + i), quNodes2.col(j));

    // Compute image of Gauss nodes under geometry mapping as well as Jacobians
    geoEval1.evaluateAt(quNodes1);
    geoEval2.evaluateAt(quNodes2);

    // Initialize local matrices
    P11.setZero(m_numActiveP1, m_numActiveP1); P12.setZero(m_numActiveP1, m_numActiveP2);
    P22.setZero(m_numActiveP2, m_numActiveP2); P21.setZero(m_numActiveP2, m_numActiveP1);
    Q11.setZero(m_numActiveP1, m_numActiveP1); Q12.setZero(m_numActiveP1, m_numActiveP2);
    Q22.setZero(m_numActiveP2, m_numActiveP2); Q21.setZero(m_numActiveP2, m_numActiveP1);

    // Evaluate solution on element nodes
    m_solActUCoeffs1.setZero(m_dim, m_numActiveU1);
    m_solActUCoeffs2.setZero(m_dim, m_numActiveU2);
    m_solUVals2.setZero(m_dim, m_basisValsU2.cols());

    for (int j = 0; j < m_numActiveU1; j++)
    {
        m_solActUCoeffs1.col(j) = m_solU.coefficientVector(m_patch1).row(m_activesU1(j)).transpose();
    }
    for (int j = 0; j < m_numActiveU2; j++)
    {
        m_solActUCoeffs2.col(j) = m_solU.coefficientVector(m_patch2).row(m_activesU2(j)).transpose();
    }

    m_solUVals1.noalias() = m_solActUCoeffs1 * m_basisValsU1;
    for (int k = 0; k < m_basisValsU2.cols(); k++)
    {
        m_solUVals2.col(k) = m_solActUCoeffs2.block(0, m_elementsU2(k), m_dim, m_numActPerNodeU2) * m_basisValsU2.col(k);
    }
}


template <class T>
inline void uwbINSdgBlockAVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        const gsMatrix<T> & basisValsU1 = m_basisDataU1[0].col(k);
        const gsMatrix<T> & bGrads1_u = m_basisDataU1[1];
        const gsMatrix<T> & basisValsU2 = m_basisDataU2[0].col(k);
        const gsMatrix<T> & bGrads2_u = m_basisDataU2[1];

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides
        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Transform the basis gradients
        geoEval1.transformGradients(k, bGrads1_u, m_physGrad1);
        geoEval2.transformGradients(k, bGrads2_u, m_physGrad2);

        // Compute element matrices

         // matrix K from the penalty term mu*[u][v]
        const T c1 = weight * m_penalty * (1 / element1.getCellSize());
        K11.noalias() += c1 * (basisValsU1 * basisValsU1.transpose());
        K22.block(m_elementsU2(k), m_elementsU2(k), m_numActPerNodeU2, m_numActPerNodeU2).noalias()
            += c1 * (basisValsU2 * basisValsU2.transpose());
        K12.block(0, m_elementsU2(k), m_numActiveU1, m_numActPerNodeU2).noalias()
            -= c1 * (basisValsU1 * basisValsU2.transpose());
        K21.block(m_elementsU2(k), 0, m_numActPerNodeU2, m_numActiveU1).noalias()
            -= c1 * (basisValsU2 * basisValsU1.transpose());

        // matrix Mx + My (+Mz) from the term {grad u}[n x v]
        const T c2 = weight * T(0.5) * m_viscosity;

        N1.noalias() = m_unormal.transpose() * m_physGrad1;
        N2.noalias() = m_unormal.transpose() * m_physGrad2;

        M11.noalias() += c2 * (basisValsU1 * N1);
        M22.block(m_elementsU2(k), m_elementsU2(k), m_numActPerNodeU2, m_numActPerNodeU2).noalias()
            -= c2 * (basisValsU2 * N2);
        M12.block(0, m_elementsU2(k), m_numActiveU1, m_numActPerNodeU2).noalias()
            += c2 * (basisValsU1 * N2);
        M21.block(m_elementsU2(k), 0, m_numActPerNodeU2, m_numActiveU1).noalias()
            -= c2 * (basisValsU2 * N1);

    }
}

template <class T>
inline void uwbINSdgBlockBVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides
        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Compute element matrices

        // matrices Lx, Ly, (Lz) from the term {p}[n . v]
        for (index_t i = 0; i != m_dim; i++)
        {
            const T c = weight * m_unormal(i) * T(0.5);

            L11[i].noalias() -= c * (m_basisValsP1.col(k)*m_basisValsU1.col(k).transpose());
            L22[i].block(m_elementsP2(k), m_elementsU2(k), m_numActPerNodeP2, m_numActPerNodeU2).noalias()
                += c * (m_basisValsP2.col(k)*m_basisValsU2.col(k).transpose());
            L12[i].block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveU1).noalias()
                -= c * (m_basisValsP2.col(k)*m_basisValsU1.col(k).transpose());
            L21[i].block(0, m_elementsU2(k), m_numActiveP1, m_numActPerNodeU2).noalias()
                += c * (m_basisValsP1.col(k)*m_basisValsU2.col(k).transpose());
            
        }
    }
}

template <class T>
inline void uwbINSdgBlockNVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides
        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Compute element matrices 

        // values of old u.n for patch 1 and 2
        T Un1 = m_unormal.transpose() * m_solUVals1.col(k);
        T Un2 = m_unormal.transpose() * m_solUVals2.col(k);
        T Un = T(0.5) * (Un1 + Un2);
  
        const T c1 = weight * Un * T(0.5);

        P11.noalias() += c1 * (m_basisValsU1.col(k) * m_basisValsU1.col(k).transpose());
        P22.block(m_elementsU2(k), m_elementsU2(k), m_numActPerNodeU2, m_numActPerNodeU2).noalias()
            += c1 * (m_basisValsU2.col(k) * m_basisValsU2.col(k).transpose());
        P12.block(0, m_elementsU2(k), m_numActiveU1, m_numActPerNodeU2).noalias()
            -= c1 * (m_basisValsU1.col(k) * m_basisValsU2.col(k).transpose());
        P21.block(m_elementsU2(k), 0, m_numActPerNodeU2, m_numActiveU1).noalias()
            -= c1 * (m_basisValsU2.col(k) * m_basisValsU1.col(k).transpose());

        const T c2 = weight * math::abs(Un) * T(0.5);

        Q11.noalias() += c2 * (m_basisValsU1.col(k) * m_basisValsU1.col(k).transpose());
        Q22.block(m_elementsU2(k), m_elementsU2(k), m_numActPerNodeU2, m_numActPerNodeU2).noalias()
            += c2 * (m_basisValsU2.col(k) * m_basisValsU2.col(k).transpose());
        Q12.block(0, m_elementsU2(k), m_numActiveU1, m_numActPerNodeU2).noalias()
            -= c2 * (m_basisValsU1.col(k) * m_basisValsU2.col(k).transpose());
        Q21.block(m_elementsU2(k), 0, m_numActPerNodeU2, m_numActiveU1).noalias()
            -= c2 * (m_basisValsU2.col(k) * m_basisValsU1.col(k).transpose());
        // what is correct??
    }
}

template <class T>
inline void uwbINSdgBlockPVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides
        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Compute element matrices

        // matrix K from the penalty term mu*[u][v]
        const T c = weight * m_penalty * (1 / element1.getCellSize());
        P11.noalias() += c * (m_basisValsP1.col(k) * m_basisValsP1.col(k).transpose());
        P22.block(m_elementsP2(k), m_elementsP2(k), m_numActPerNodeP2, m_numActPerNodeP2).noalias()
            += c * (m_basisValsP2.col(k) * m_basisValsP2.col(k).transpose());
        P12.block(0, m_elementsP2(k), m_numActiveP1, m_numActPerNodeP2).noalias()
            -= c * (m_basisValsP1.col(k) * m_basisValsP2.col(k).transpose());
        P21.block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveP1).noalias()
            -= c * (m_basisValsP2.col(k) * m_basisValsP1.col(k).transpose());

    }
}

template <class T>
inline void uwbINSdgBlockApVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        const gsMatrix<T> & basisValsP1 = basisData1_p[0].col(k);
        const gsMatrix<T> & bGrads1_p = basisData1_p[1];
        const gsMatrix<T> & basisValsP2 = basisData2_p[0].col(k);
        const gsMatrix<T> & bGrads2_p = basisData2_p[1];

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides

        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Transform the basis gradients
        geoEval1.transformGradients(k, bGrads1_p, m_physGrad1);
        geoEval2.transformGradients(k, bGrads2_p, m_physGrad2);

        // Compute element matrices

        // matrix K from the penalty term mu*[p][q]
        const T c1 = weight * m_penalty * (1 / element1.getCellSize());
        K11.noalias() += c1 * (basisValsP1 * basisValsP1.transpose());
        K22.block(m_elementsP2(k), m_elementsP2(k), m_numActPerNodeP2, m_numActPerNodeP2).noalias()
            += c1 * (basisValsP2 * basisValsP2.transpose());
        K12.block(0, m_elementsP2(k), m_numActiveP1, m_numActPerNodeP2).noalias()
            -= c1 * (basisValsP1 * basisValsP2.transpose());
        K21.block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveP1).noalias()
            -= c1 * (basisValsP2 * basisValsP1.transpose());

        // matrix Mx + My (+Mz) from the term {grad p}[n x q]
        const T c2 = weight * T(0.5) * m_viscosity;

        N1.noalias() = m_unormal.transpose() * m_physGrad1;
        N2.noalias() = m_unormal.transpose() * m_physGrad2;

        M11.noalias() += c2 * (basisValsP1 * N1);
        M22.block(m_elementsP2(k), m_elementsP2(k), m_numActPerNodeP2, m_numActPerNodeP2).noalias()
            -= c2 * (basisValsP2 * N2);
        M12.block(0, m_elementsP2(k), m_numActiveP1, m_numActPerNodeP2).noalias()
            += c2 * (basisValsP1 * N2);
        M21.block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveP1).noalias()
            -= c2 * (basisValsP2 * N1);

    }
}

template <class T>
inline void uwbINSdgBlockNpVisitor<T>::assemble(gsDomainIterator<T>    & element1,
    gsGeometryEvaluator<T> & geoEval1,
    gsGeometryEvaluator<T> & geoEval2,
    gsVector<T>            & quWeights)
{
    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        // Compute the outer normal vector from patch1
        geoEval1.outerNormal(k, m_side1, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides
        const T weight = quWeights[k] * m_unormal.norm();
        m_unormal.normalize();

        // Compute element matrices 

        // values of old u.n for patch 1 and 2
        T Un1 = m_unormal.transpose() * m_solUVals1.col(k);
        T Un2 = m_unormal.transpose() * m_solUVals2.col(k);
        T Un = T(0.5) * (Un1 + Un2);

        const T c1 = weight * Un * T(0.5);

        P11.noalias() += c1 * (m_basisValsP1.col(k) * m_basisValsP1.col(k).transpose());
        P22.block(m_elementsP2(k), m_elementsP2(k), m_numActPerNodeP2, m_numActPerNodeP2).noalias()
            += c1 * (m_basisValsP2.col(k) * m_basisValsP2.col(k).transpose());
        P12.block(0, m_elementsP2(k), m_numActiveP1, m_numActPerNodeP2).noalias()
            -= c1 * (m_basisValsP1.col(k) * m_basisValsP2.col(k).transpose());
        P21.block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveP1).noalias()
            -= c1 * (m_basisValsP2.col(k) * m_basisValsP1.col(k).transpose());

        const T c2 = weight * math::abs(Un) * T(0.5);

        Q11.noalias() += c2 * (m_basisValsP1.col(k) * m_basisValsP1.col(k).transpose());
        Q22.block(m_elementsP2(k), m_elementsP2(k), m_numActPerNodeP2, m_numActPerNodeP2).noalias()
            += c2 * (m_basisValsP2.col(k) * m_basisValsP2.col(k).transpose());
        Q12.block(0, m_elementsP2(k), m_numActiveP1, m_numActPerNodeP2).noalias()
            -= c2 * (m_basisValsP1.col(k) * m_basisValsP2.col(k).transpose());
        Q21.block(m_elementsP2(k), 0, m_numActPerNodeP2, m_numActiveP1).noalias()
            -= c2 * (m_basisValsP2.col(k) * m_basisValsP1.col(k).transpose());
        // what is correct??
    }
}

template <class T>
inline void uwbINSdgBlockAVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    
    const index_t usz = m_Umap.freeSize();
    // Local DoFs to global DoFs
    m_Umap.localToGlobal(m_activesU1, m_patch1, m_activesU1);
    m_Umap.localToGlobal(m_activesU2, m_patch2, m_activesU2);
    
    // Push element contributions 1-2 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU1; ++i)
    {
        const index_t ii1 = m_activesU1(i);  //N1_i
        if (m_Umap.is_free_index(ii1))
        {
            for (index_t j = 0; j != m_numActiveU1; ++j)
            {
                const index_t  jj1 = m_activesU1(j); //N1_j
                const T tmp = K11(i, j) - M11(i, j) - M11(j, i);

                if (tmp != 0)
                {
                    if (m_Umap.is_free_index(jj1))
                    {
                        blockMatrix.coeffRef(ii1, jj1) += tmp;
                    }
                    else //jj1 is boundary index
                    {
                        const int bb = m_Umap.global_to_bindex(jj1);
                        for (index_t s = 0; s != m_dim; s++)
                        {
                            blockRhs(ii1 + s*usz, 0) -= tmp * eliminatedDofs[0](bb, s);
                        }
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveU2; ++j)
            {
                const index_t  jj2 = m_activesU2(j); //N2_j
                const T tmp = K12(i, j) - M12(i, j) - M21(j, i);

                if (tmp != 0)
                {
                    if (m_Umap.is_free_index(jj2))
                    {
                        blockMatrix.coeffRef(ii1, jj2) += tmp;

                    }
                    else //jj2 is boundary index
                    {
                        const int bb = m_Umap.global_to_bindex(jj2);
                        for (index_t s = 0; s != m_dim; s++)
                        {
                            blockRhs(ii1 + s*usz, 0) -= tmp * eliminatedDofs[0](bb, s);
                        }
                    }
                }
            }
        }
    }

    // Push element contributions 2-1 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU2; ++i)
    {
        const index_t ii2 = m_activesU2(i);
        if (m_Umap.is_free_index(ii2))
        {
            for (index_t j = 0; j != m_numActiveU2; ++j)
            {
                const index_t  jj1 = m_activesU2(j);
                const T tmp = K22(i, j) - M22(i, j) - M22(j, i);

                if (tmp != 0)
                {
                    if (m_Umap.is_free_index(jj1))
                    {
                        blockMatrix.coeffRef(ii2, jj1) += tmp;
                    }
                    else //jj1 is boundary index
                    {
                        const int bb = m_Umap.global_to_bindex(jj1);
                        for (index_t s = 0; s != m_dim; s++)
                        {
                            blockRhs(ii2 + s*usz, 0) -= tmp * eliminatedDofs[0](bb, s);
                        }
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveU1; ++j)
            {
                const index_t  jj2 = m_activesU1(j);
                const T tmp = K21(i, j) - M21(i, j) - M12(j, i);

                if (tmp != 0)
                {
                    if (m_Umap.is_free_index(jj2))
                    {
                        blockMatrix.coeffRef(ii2, jj2) += tmp;
                    }
                    else //jj2 is boundary index
                    {
                        const int bb = m_Umap.global_to_bindex(jj2);
                        for (index_t s = 0; s != m_dim; s++)
                        {
                            blockRhs(ii2 + s*usz, 0) -= tmp * eliminatedDofs[0](bb, s);
                        }
                    }
                }
            }
        }
    }

}

template <class T>
inline void uwbINSdgBlockBVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    
    const index_t usz = m_Umap.freeSize();
    const index_t ps = m_dim*usz;

    // Local DoFs to global DoFs
    m_Umap.localToGlobal(m_activesU1, m_patch1, m_activesU1);
    m_Umap.localToGlobal(m_activesU2, m_patch2, m_activesU2);
    m_Pmap.localToGlobal(m_activesP1, m_patch1, m_activesP1);
    m_Pmap.localToGlobal(m_activesP2, m_patch2, m_activesP2);

    // Push element contributions 1-2 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU1; ++i)
    {
        const index_t ii1 = m_activesU1(i);  //N1_i
        if (m_Umap.is_free_index(ii1))
        {
            
            for (index_t j = 0; j < m_numActiveP1; ++j)
            {
                const int jj1 = m_activesP1(j); //Q1_j

                if (m_Pmap.is_free_index(jj1))
                {
                    for (index_t s = 0; s != m_dim; s++)
                        blockMatrix.coeffRef(jj1, ii1 + s*usz) += L11[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj1)
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    for (index_t s = 0; s < m_dim; ++s)
                        blockRhs(ii1 + s * usz, 0) += L11[s](i, j) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j < m_numActiveP2; ++j)
            {
                const int jj2 = m_activesP2(j); //Q2_j

                if (m_Pmap.is_free_index(jj2))
                {
                    for (index_t s = 0; s != m_dim; s++)
                        blockMatrix.coeffRef(jj2, ii1 + s*usz) += L12[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj2)
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    for (index_t s = 0; s < m_dim; ++s)
                        blockRhs(ii1 + s * usz, 0) += L12[s](i, j) * eliminatedDofs[1](bb, 0);
                }
            }
        }
        else //ii1 is boundary index
        {
            const int bb = m_Umap.global_to_bindex(ii1);
            for (index_t k = 0; k < m_numActiveP1; ++k)
            {
                const int kk = m_activesP1(k);
                if (m_Pmap.is_free_index(kk))
                {
                    T tmp = L11[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += L11[s](k, i) * eliminatedDofs[0](bb, s);
                    blockRhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
            for (index_t k = 0; k < m_numActiveP2; ++k)
            {
                const int kk = m_activesP2(k);
                if (m_Pmap.is_free_index(kk))
                {
                    T tmp = L12[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += L12[s](k, i) * eliminatedDofs[0](bb, s);
                    blockRhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
        }
    }

    // Push element contributions 2-1 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU2; ++i)
    {
        const index_t ii2 = m_activesU2(i);
        if (m_Umap.is_free_index(ii2))
        {
            
            for (index_t j = 0; j < m_numActiveP2; ++j)
            {
                const int jj1 = m_activesP2(j);

                if (m_Pmap.is_free_index(jj1))
                {
                    for (index_t s = 0; s != m_dim; s++)
                        blockMatrix.coeffRef(jj1, ii2 + s*usz) += L22[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj1)
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    for (index_t s = 0; s < m_dim; ++s)
                        blockRhs(ii2 + s * usz, 0) += L22[s](i, j) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j < m_numActiveP1; ++j)
            {
                const int jj2 = m_activesP1(j);

                if (m_Pmap.is_free_index(jj2))
                {
                    for (index_t s = 0; s != m_dim; s++)
                        blockMatrix.coeffRef(jj2, ii2 + s*usz) += L21[s](j, i);
                }
                else //m_Pmap.is_boundary_index(jj2)
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    for (index_t s = 0; s < m_dim; ++s)
                        blockRhs(ii2 + s * usz, 0) += L21[s](i, j) * eliminatedDofs[1](bb, 0);
                }
            }
        }
        else //ii2 is boundary index
        {
            const int bb = m_Umap.global_to_bindex(ii2);
            for (index_t k = 0; k < m_numActiveP2; ++k)
            {
                const int kk = m_activesP2(k);
                if (m_Pmap.is_free_index(kk))
                {
                    T tmp = L22[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += L22[s](k, i) * eliminatedDofs[0](bb, s);
                    blockRhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
            for (index_t k = 0; k < m_numActiveP1; ++k)
            {
                const int kk = m_activesP1(k);
                if (m_Pmap.is_free_index(kk))
                {
                    T tmp = L21[0](k, i)*eliminatedDofs[0](bb, 0);
                    for (index_t s = 1; s != m_dim; ++s)
                        tmp += L21[s](k, i) * eliminatedDofs[0](bb, s);
                    blockRhs(ps + kk, 0) -= tmp;// assuming single rhs
                }
            }
        }
    }
}

template <class T>
inline void uwbINSdgBlockNVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    const index_t usz = m_Umap.freeSize();

    // Local DoFs to global DoFs
    m_Umap.localToGlobal(m_activesU1, m_patch1, m_activesU1);
    m_Umap.localToGlobal(m_activesU2, m_patch2, m_activesU2);
    
    // Push element contributions 1-2 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU1; ++i)
    {
        const index_t ii1 = m_activesU1(i);  //N1_i
        if (m_Umap.is_free_index(ii1))
        {
            for (index_t j = 0; j != m_numActiveU1; ++j)
            {
                const index_t  jj1 = m_activesU1(j); //N1_j
                if (m_Umap.is_free_index(jj1))
                {
                    blockMatrix.coeffRef(ii1, jj1) += P11(i, j) - Q11(i, j);
                }
                else //jj1 is boundary index
                {
                    const int bb = m_Umap.global_to_bindex(jj1);
                    for (index_t s = 0; s != m_dim; s++)
                    {
                        blockRhs(ii1 + s*usz, 0) -= (P11(i, j) - Q11(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveU2; ++j)
            {
                const index_t  jj2 = m_activesU2(j); //N2_j
                if (m_Umap.is_free_index(jj2))
                {
                    blockMatrix.coeffRef(ii1, jj2) += P12(i, j) - Q12(i, j);
                }
                else //jj2 is boundary index
                {
                    const int bb = m_Umap.global_to_bindex(jj2);
                    for (index_t s = 0; s != m_dim; s++)
                    {
                        blockRhs(ii1 + s*usz, 0) -= (P12(i, j) - Q12(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

    // Push element contributions 2-1 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveU2; ++i)
    {
        const index_t ii2 = m_activesU2(i);
        if (m_Umap.is_free_index(ii2))
        {
            for (index_t j = 0; j != m_numActiveU2; ++j)
            {
                const index_t  jj1 = m_activesU2(j);
                if (m_Umap.is_free_index(jj1))
                {
                    blockMatrix.coeffRef(ii2, jj1) += P22(i, j) - Q22(i, j);
                }
                else //jj1 is boundary index
                {
                    const int bb = m_Umap.global_to_bindex(jj1);
                    for (index_t s = 0; s != m_dim; s++)
                    {
                        blockRhs(ii2 + s*usz, 0) -= (P22(i, j) - Q22(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveU1; ++j)
            {
                const index_t  jj2 = m_activesU1(j);
                if (m_Umap.is_free_index(jj2))
                {
                    blockMatrix.coeffRef(ii2, jj2) += P21(i, j) - Q21(i, j);
                }
                else //jj2 is boundary index
                {
                    const int bb = m_Umap.global_to_bindex(jj2);
                    for (index_t s = 0; s != m_dim; s++)
                    {
                        blockRhs(ii2 + s*usz, 0) -= (P21(i, j) - Q21(i, j)) * eliminatedDofs[0](bb, s);
                    }
                }
            }
        }
    }

}

template <class T>
inline void uwbINSdgBlockPVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    
    // Local DoFs to global DoFs
    m_Pmap.localToGlobal(m_activesP1, m_patch1, m_activesP1);
    m_Pmap.localToGlobal(m_activesP2, m_patch2, m_activesP2);

    for (index_t i = 0; i != m_numActiveP1; ++i)
    {
        const index_t ii1 = m_activesP1(i);
        if (m_Pmap.is_free_index(ii1))
        {
            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t  jj1 = m_activesP1(j);
                if (m_Pmap.is_free_index(jj1))
                    blockMatrix.coeffRef(ii1, jj1) += P11(i, j);
                else // m_Pmap.is_boundary_index(jj1)
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    blockRhs(ii1, 0) -= P11(i, j) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t  jj2 = m_activesP2(j);
                if (m_Pmap.is_free_index(jj2))
                    blockMatrix.coeffRef(ii1, jj2) += P12(i, j);
                else // m_Pmap.is_boundary_index(jj2)
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    blockRhs(ii1, 0) -= P12(i, j) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }
    for (index_t i = 0; i != m_numActiveP2; ++i)
    {
        const index_t ii2 = m_activesP2(i);
        if (m_Pmap.is_free_index(ii2))
        {
            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t jj1 = m_activesP1(j);
                if (m_Pmap.is_free_index(jj1))
                    blockMatrix.coeffRef(ii2, jj1) += P21(i, j);
                else // m_Pmap.is_boundary_index(jj1)
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    blockRhs(ii2, 0) -= P21(i, j) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t jj2 = m_activesP2(j);
                if (m_Pmap.is_free_index(jj2))
                    blockMatrix.coeffRef(ii2, jj2) += P22(i, j);
                else // m_Pmap.is_boundary_index(jj2)
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    blockRhs(ii2, 0) -= P22(i, j) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }

}

template <class T>
inline void uwbINSdgBlockApVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs, 
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    // Local DoFs to global DoFs
    m_Pmap.localToGlobal(m_activesP1, m_patch1, m_activesP1);
    m_Pmap.localToGlobal(m_activesP2, m_patch2, m_activesP2);

    // Push element contributions 1-2 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveP1; ++i)
    {
        const index_t ii1 = m_activesP1(i);  //N1_i
        if (m_Pmap.is_free_index(ii1))
        {
            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t  jj1 = m_activesP1(j); //N1_j
                const T tmp = K11(i, j) - M11(i, j) - M11(j, i);

                if (tmp != 0)
                {
                    if (m_Pmap.is_free_index(jj1))
                    {
                        blockMatrix.coeffRef(ii1, jj1) += tmp;
                    }
                    else //jj1 is boundary index
                    {
                        const int bb = m_Pmap.global_to_bindex(jj1);
                        blockRhs(ii1, 0) -= tmp * eliminatedDofs[1](bb, 0);
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t  jj2 = m_activesP2(j); //N2_j
                const T tmp = K12(i, j) - M12(i, j) - M21(j, i);

                if (tmp != 0)
                {
                    if (m_Pmap.is_free_index(jj2))
                    {
                        blockMatrix.coeffRef(ii1, jj2) += tmp;
                    }
                    else //jj2 is boundary index
                    {
                        const int bb = m_Pmap.global_to_bindex(jj2);
                        blockRhs(ii1, 0) -= tmp * eliminatedDofs[1](bb, 0);
                    }
                }
            }
        }
    }

    // Push element contributions 2-1 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveP2; ++i)
    {
        const index_t ii2 = m_activesP2(i);
        if (m_Pmap.is_free_index(ii2))
        {
            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t  jj1 = m_activesP2(j);
                const T tmp = K22(i, j) - M22(i, j) - M22(j, i);

                if (tmp != 0)
                {
                    if (m_Pmap.is_free_index(jj1))
                    {
                        blockMatrix.coeffRef(ii2, jj1) += tmp;
                    }
                    else //jj1 is boundary index
                    {
                        const int bb = m_Pmap.global_to_bindex(jj1);
                        blockRhs(ii2, 0) -= tmp * eliminatedDofs[1](bb, 0);
                    }
                }
            }

            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t  jj2 = m_activesP1(j);
                const T tmp = K21(i, j) - M21(i, j) - M12(j, i);

                if (tmp != 0)
                {
                    if (m_Pmap.is_free_index(jj2))
                    {
                        blockMatrix.coeffRef(ii2, jj2) += tmp;
                    }
                    else //jj2 is boundary index
                    {
                        const int bb = m_Pmap.global_to_bindex(jj2);
                        blockRhs(ii2, 0) -= tmp * eliminatedDofs[1](bb, 0);
                    }
                }
            }
        }
    }
}

template <class T>
inline void uwbINSdgBlockNpVisitor<T>::localToGlobal(const std::vector<gsMatrix<T> >& eliminatedDofs,
    gsSparseMatrix<T>     & blockMatrix,
    gsMatrix<T>           & blockRhs)
{
    // Local DoFs to global DoFs
    m_Pmap.localToGlobal(m_activesP1, m_patch1, m_activesP1);
    m_Pmap.localToGlobal(m_activesP2, m_patch2, m_activesP2);

    // Push element contributions 1-2 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveP1; ++i)
    {
        const index_t ii1 = m_activesP1(i);  //N1_i
        if (m_Pmap.is_free_index(ii1))
        {
            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t  jj1 = m_activesP1(j); //N1_j
                if (m_Pmap.is_free_index(jj1))
                {
                    blockMatrix.coeffRef(ii1, jj1) += P11(i, j) - Q11(i, j);
                }
                else //jj1 is boundary index
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    blockRhs(ii1, 0) -= (P11(i, j) - Q11(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t  jj2 = m_activesP2(j); //N2_j
                if (m_Pmap.is_free_index(jj2))
                {
                    blockMatrix.coeffRef(ii1, jj2) += P12(i, j) - Q12(i, j);
                }
                else //jj2 is boundary index
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    blockRhs(ii1, 0) -= (P12(i, j) - Q12(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }

    // Push element contributions 2-1 to the global matrix and load vector
    for (index_t i = 0; i != m_numActiveP2; ++i)
    {
        const index_t ii2 = m_activesP2(i);
        if (m_Pmap.is_free_index(ii2))
        {
            for (index_t j = 0; j != m_numActiveP2; ++j)
            {
                const index_t  jj1 = m_activesP2(j);
                if (m_Pmap.is_free_index(jj1))
                {
                    blockMatrix.coeffRef(ii2, jj1) += P22(i, j) - Q22(i, j);
                }
                else //jj1 is boundary index
                {
                    const int bb = m_Pmap.global_to_bindex(jj1);
                    blockRhs(ii2, 0) -= (P22(i, j) - Q22(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }

            for (index_t j = 0; j != m_numActiveP1; ++j)
            {
                const index_t  jj2 = m_activesP1(j);
                if (m_Pmap.is_free_index(jj2))
                {
                    blockMatrix.coeffRef(ii2, jj2) += P21(i, j) - Q21(i, j);
                }
                else //jj2 is boundary index
                {
                    const int bb = m_Pmap.global_to_bindex(jj2);
                    blockRhs(ii2, 0) -= (P21(i, j) - Q21(i, j)) * eliminatedDofs[1](bb, 0);
                }
            }
        }
    }
}

template <class T>
inline void uwbINSdgBlockAVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_u = m_activesU2.row(0);

    // creating active function vectors and removing duplicate indices
    m_activesU2.resize(1, m_activesU2.cols() * m_activesU2.rows());
    for (index_t k = 1; k < m_activesU2.cols(); k++)
    {
        if (m_activesU2(k) <= m_activesU2(k - 1))
        {
            m_activesU2.removeCol(k--);
        }
    }
    m_activesU2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesU2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_u.cols(); i++)
        {
            if (m_activesU2(k) == basisFuncIndices_u(i))
                m_elementsU2(i) = k;
        }
    }

}

template <class T>
inline void uwbINSdgBlockBVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_u = m_activesU2.row(0);
    gsMatrix<index_t> basisFuncIndices_p = m_activesP2.row(0);

    // creating active function vectors and removing duplicate indices
    m_activesU2.resize(1, m_activesU2.cols() * m_activesU2.rows());
    for (index_t k = 1; k < m_activesU2.cols(); k++)
    {
        if (m_activesU2(k) <= m_activesU2(k - 1))
        {
            m_activesU2.removeCol(k--);
        }
    }
    m_activesU2.transposeInPlace(); // back to the column vector

    m_activesP2.resize(1, m_activesP2.cols() * m_activesP2.rows());
    for (index_t k = 1; k < m_activesP2.cols(); k++)
    {
        if (m_activesP2(k) <= m_activesP2(k - 1))
        {
            m_activesP2.removeCol(k--);
        }
    }
    m_activesP2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesU2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_u.cols(); i++)
        {
            if (m_activesU2(k) == basisFuncIndices_u(i))
                m_elementsU2(i) = k;
        }
    }

    for (index_t k = 0; k < m_activesP2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_p.cols(); i++)
        {
            if (m_activesP2(k) == basisFuncIndices_p(i))
                m_elementsP2(i) = k;
        }
    }
}

template <class T>
inline void uwbINSdgBlockNVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_u = m_activesU2.row(0);
    
    // creating active function vectors and removing duplicate indices
    m_activesU2.resize(1, m_activesU2.cols() * m_activesU2.rows());
    for (index_t k = 1; k < m_activesU2.cols(); k++)
    {
        if (m_activesU2(k) <= m_activesU2(k - 1))
        {
            m_activesU2.removeCol(k--);
        }
    }
    m_activesU2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesU2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_u.cols(); i++)
        {
            if (m_activesU2(k) == basisFuncIndices_u(i))
                m_elementsU2(i) = k;
        }
    }
}

template <class T>
inline void uwbINSdgBlockPVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_p = m_activesP2.row(0);

    // creating active function vectors and removing duplicate indices
    m_activesP2.resize(1, m_activesP2.cols() * m_activesP2.rows());
    for (index_t k = 1; k < m_activesP2.cols(); k++)
    {
        if (m_activesP2(k) <= m_activesP2(k - 1))
        {
            m_activesP2.removeCol(k--);
        }
    }
    m_activesP2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesP2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_p.cols(); i++)
        {
            if (m_activesP2(k) == basisFuncIndices_p(i))
                m_elementsP2(i) = k;
        }
    }
}

template <class T>
inline void uwbINSdgBlockApVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_p = m_activesP2.row(0);

    // creating active function vectors and removing duplicate indices
    m_activesP2.resize(1, m_activesP2.cols() * m_activesP2.rows());
    for (index_t k = 1; k < m_activesP2.cols(); k++)
    {
        if (m_activesP2(k) <= m_activesP2(k - 1))
        {
            m_activesP2.removeCol(k--);
        }
    }
    m_activesP2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesP2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_p.cols(); i++)
        {
            if (m_activesP2(k) == basisFuncIndices_p(i))
                m_elementsP2(i) = k;
        }
    }
}

template <class T>
inline void uwbINSdgBlockNpVisitor<T>::activesPatch2()
{
    // vectors of first basis function index in each quadrature node
    gsMatrix<index_t> basisFuncIndices_u = m_activesU2.row(0);
    gsMatrix<index_t> basisFuncIndices_p = m_activesP2.row(0);

    // creating active function vectors and removing duplicate indices
    m_activesU2.resize(1, m_activesU2.cols() * m_activesU2.rows());
    for (index_t k = 1; k < m_activesU2.cols(); k++)
    {
        if (m_activesU2(k) <= m_activesU2(k - 1))
        {
            m_activesU2.removeCol(k--);
        }
    }
    m_activesU2.transposeInPlace(); // back to the column vector

    m_activesP2.resize(1, m_activesP2.cols() * m_activesP2.rows());
    for (index_t k = 1; k < m_activesP2.cols(); k++)
    {
        if (m_activesP2(k) <= m_activesP2(k - 1))
        {
            m_activesP2.removeCol(k--);
        }
    }
    m_activesP2.transposeInPlace(); // back to the column vector

    // saving indices (in actives2 vectors) of the first basis functions 
    for (index_t k = 0; k < m_activesU2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_u.cols(); i++)
        {
            if (m_activesU2(k) == basisFuncIndices_u(i))
                m_elementsU2(i) = k;
        }
    }

    for (index_t k = 0; k < m_activesP2.rows(); k++)
    {
        for (index_t i = 0; i < basisFuncIndices_p.cols(); i++)
        {
            if (m_activesP2(k) == basisFuncIndices_p(i))
                m_elementsP2(i) = k;
        }
    }
}

template <class T>
void uwbINSiFaceSizeVisitor<T>::initialize(const gsBasis<T> & basis,
    boxSide s,
    const int patchIndex,
    gsQuadRule<T> & rule,
    unsigned & evFlags)
{
    m_side = s;
    m_patch = patchIndex;
    m_size = 0;

    m_dim = basis.dim();
    const int dir = m_side.direction();
    gsVector<int> numQuadNodes(m_dim);
    for (int i = 0; i < basis.dim(); ++i)
        numQuadNodes[i] = basis.degree(i) + 1;
    numQuadNodes[dir] = 1;

    // Setup Quadrature
    rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

    evFlags = NEED_VALUE | NEED_JACOBIAN;
}

template <class T>
void uwbINSiFaceSizeVisitor<T>::assemble(gsDomainIterator<T>    & element,
    gsGeometryEvaluator<T> & geoEval,
    gsMatrix<T>            & quNodes,
    gsVector<T>            & quWeights)
{
    geoEval.evaluateAt(quNodes);

    for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
    {

        // Compute the outer normal vector from patch1
        geoEval.outerNormal(k, m_side, m_unormal);

        // Integral transformation and quadrature weight (patch1)
        // assumed the same on both sides

        const T weight = quWeights[k] * m_unormal.norm();

        m_size += weight;

    }
}

}
