/** @file gsINSAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSAssembler.h>

namespace gismo
{

template<class T>
void gsINSAssembler<T>::initMembers()
{ 
    Base::initMembers();

    m_viscosity = m_params.getPde().viscosity();

    m_dofMappers.resize(2);
    m_ddof.resize(2);
    getBases().front().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.front(), 0);
    getBases().back().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.back(), 1);

    updateSizes();

    m_nnzPerRowU = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowU *= 2 * getBases().front().maxDegree(i) + 1;

    m_nnzPerRowP = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowP *= 2 * getBases().back().maxDegree(i) + 1;

    m_visitorUUlin = gsINSVisitorUUlin<T>(m_params);
    m_visitorUUlin.initialize();

    m_visitorUUnonlin = gsINSVisitorUUnonlin<T>(m_params);
    m_visitorUUnonlin.initialize();
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);

    m_visitorUP = gsINSVisitorPU_withUPrhs<T>(m_params);
    m_visitorUP.initialize();

    m_visitorF = gsINSVisitorRhsU<T>(m_params);
    m_visitorF.initialize();

    m_visitorG = gsINSVisitorRhsP<T>(m_params);
    m_visitorG.initialize();
}


template<class T>
void gsINSAssembler<T>::updateSizes()
{
    m_udofs = m_dofMappers.front().freeSize();
    m_pdofs = m_dofMappers.back().freeSize();
    m_pshift = m_tarDim * m_udofs;
    m_dofs = m_pshift + m_pdofs;

    m_ddof[0].setZero(m_dofMappers.front().boundarySize(), m_tarDim);
    m_ddof[1].setZero(m_dofMappers.back().boundarySize(), 1);

    if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
    {
        this->computeDirichletDofs(0, 0, m_ddof[0]);
        this->computeDirichletDofs(1, 1, m_ddof[1]);
    }

    m_solution.setZero(m_dofs, 1);

    m_currentVelField = constructSolution(m_solution, 0);

    m_blockUUlin.resize(m_pshift, m_pshift);
    m_blockUUnonlin.resize(m_pshift, m_pshift);
    m_blockUP.resize(m_pshift, m_pdofs);

    // memory allocation
    m_blockUUlin.reserve(gsVector<index_t>::Constant(m_blockUUlin.rows(), m_nnzPerRowU));
    m_blockUUnonlin.reserve(gsVector<index_t>::Constant(m_blockUUnonlin.rows(), m_nnzPerRowU));
    m_blockUP.reserve(gsVector<index_t>::Constant(m_blockUP.rows(), m_nnzPerRowU));

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_rhsUlin.setZero(m_pshift, 1);
    m_rhsUnonlin.setZero(m_pshift, 1);
    m_rhsBtB.setZero(m_dofs, 1);
    m_rhsFG.setZero(m_dofs, 1);
    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);
}


template<class T>
void gsINSAssembler<T>::updateDofMappers()
{
    m_visitorUUlin.updateDofMappers(m_dofMappers);
    m_visitorUUnonlin.updateDofMappers(m_dofMappers);
    m_visitorUP.updateDofMappers(m_dofMappers);
    m_visitorF.updateDofMappers(m_dofMappers);
    m_visitorG.updateDofMappers(m_dofMappers);
}


template<class T>
void gsINSAssembler<T>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
        m_solution = solVector;

    m_currentVelField = constructSolution(solVector, 0);
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);
}


template<class T>
void gsINSAssembler<T>::assembleLinearPart()
{
    // matrix and rhs cleaning
    m_blockUUlin.setZero();
    m_blockUP.setZero();
    m_rhsUlin.setZero();
    m_rhsBtB.setZero();
    m_rhsFG.setZero();

    this->assembleBlock(m_visitorUUlin, 0, m_blockUUlin, m_rhsUlin);
    this->assembleBlock(m_visitorUP, 0, m_blockUP, m_rhsBtB);
    this->assembleRhs(m_visitorF, 0, m_rhsFG);

    if(m_params.getPde().source()) // if the continuity eqn rhs is given
        this->assembleRhs(m_visitorG, 1, m_rhsFG);
}


template<class T>
void gsINSAssembler<T>::assembleNonlinearPart()
{
    // matrix and rhs cleaning
    m_blockUUnonlin.setZero();
    m_rhsUnonlin.setZero();

    this->assembleBlock(m_visitorUUnonlin, 0, m_blockUUnonlin, m_rhsUnonlin);
}


template<class T>
void gsINSAssembler<T>::fillBaseSystem() 
{
    gsVector<index_t> nnzPerRowVector;
    nnzPerRowVector.setZero(m_dofs);

    gsSparseMatrix<T, RowMajor> blockPU = gsSparseMatrix<T, RowMajor>(-m_blockUP.transpose());

    for (index_t i = 0; i < m_pshift; i++)
        nnzPerRowVector(i) = m_blockUUlin.row(i).nonZeros() + m_blockUP.row(i).nonZeros();

    for (index_t i = 0; i < m_pdofs; i++)
        nnzPerRowVector(m_pshift + i) = blockPU.row(i).nonZeros();

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_baseMatrix.reserve(nnzPerRowVector);

    for (index_t row = 0; row < m_pshift; row++)
    {
        for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_blockUUlin, row); it; ++it)
            m_baseMatrix.insert(row, it.col()) = it.value();

        for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_blockUP, row); it; ++it)
            m_baseMatrix.insert(row, m_pshift + it.col()) = it.value();
    }

    for (index_t row = 0; row < m_pdofs; row++)
        for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(blockPU, row); it; ++it)
            m_baseMatrix.insert(m_pshift + row, it.col()) = it.value();


    m_baseRhs.noalias() = m_rhsFG + m_rhsBtB;
    m_baseRhs.topRows(m_pshift) += m_rhsUlin;

    m_isBaseReady = true;
    m_isSystemReady = false;
}


template<class T>
void gsINSAssembler<T>::fillSystem()
{
    m_matrix = m_baseMatrix;

    for (index_t row = 0; row < m_pshift; row++)
        for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_blockUUnonlin, row); it; ++it)
            m_matrix.coeffRef(row, it.col()) += it.value();

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_rhs = m_baseRhs;
    m_rhs.topRows(m_pshift) += m_rhsUnonlin;

    m_isSystemReady = true;
}


template<class T>
void gsINSAssembler<T>::initialize()
{
    Base::initialize();

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
}


template<class T>
void gsINSAssembler<T>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


template<class T>
void gsINSAssembler<T>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
{
    m_dofMappers[unk] = gsDofMapper(getBases().at(unk), getBCs(), unk);

    if (getAssemblerOptions().intStrategy == iFace::conforming)
        for (gsBoxTopology::const_iiterator it = getPatches().iBegin(); it != getPatches().iEnd(); ++it)
        {
            getBases().at(unk).matchInterface(*it, m_dofMappers[unk]);
        }

    for (size_t i = 0; i < boundaryDofs.size(); i++)
        m_dofMappers[unk].markBoundary(i, boundaryDofs[i]);

    m_dofMappers[unk].finalize();

    this->updateDofMappers();
    this->updateSizes();
}


template<class T>
gsField<T> gsINSAssembler<T>::constructSolution(const gsMatrix<T>& solVector, index_t unk) const
{
    GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

    gsMultiPatch<T>* result = new gsMultiPatch<T>;

    const gsDofMapper& mapper = m_dofMappers[unk];

    const index_t dim = (unk == 0 ? m_tarDim : 1);
    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    gsAsConstMatrix<T> solV = (unk == 0 ?
        gsAsConstMatrix<T>(solVector.data(), m_udofs, dim)
        :
        gsAsConstMatrix<T>(solVector.data() + m_pshift, m_pdofs, 1)
        );

    for (size_t p = 0; p < this->getPatches().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const index_t sz = this->getBases().at(unk).piece(p).size();
        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, p)) // DoF value is in the solVector
            {
                coeffs.row(i) = solV.row(mapper.index(i, p));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof[unk].row(mapper.bindex(i, p));
            }
        }

        result->addPatch(this->getBases().at(unk).piece(p).makeGeometry(coeffs));
    }

    return gsField<T>(this->getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
}


template<class T>
T gsINSAssembler<T>::computeFlowRate(index_t patch, boxSide side, gsMatrix<T> solution) const
{
    T flowRate = 0;

    gsField<T> solutionField = constructSolution(solution, 0); // velocity field

    const gsGeometry<T>& geo = this->getPatches().patch(patch);
    const gsBasis<T>& basis = this->getBases().at(0).basis(patch);

    gsVector<index_t> numQuadNodes(m_tarDim);
    const index_t dir = side.direction();
    for (short_t i = 0; i < m_tarDim; ++i)
        numQuadNodes[i] = (2 * basis.degree(i) + 1);
    numQuadNodes[dir] = 1;

    // Setup Quadrature
    gsGaussRule<T> QuRule(numQuadNodes);

    gsMatrix<T> quNodes; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights
   
    // Initialize geometry evaluator
    gsMapData<T> mapData;
    mapData.flags = NEED_VALUE | NEED_OUTER_NORMAL;

    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(side);
    for (; domIt->good(); domIt->next())
    {
        // Compute the quadrature rule on patch1
        QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

        // Compute image of Gauss nodes under geometry mapping as well as Jacobians
        mapData.points = quNodes;
        geo.computeMap(mapData);

        // Evaluate solution on element nodes
        gsMatrix<T> solUVals = solutionField.value(quNodes, patch);

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector from patch1
            gsVector<T> normal;
            outerNormal(mapData, k, side, normal);

            // the normal norm is equal to integral measure
            flowRate += quWeights[k] * normal.dot(solUVals.col(k));
        }
    }
    return flowRate;
}

template<class T>
void gsINSAssembler<T>::fillStokesSystem(gsSparseMatrix<T, RowMajor>& stokesMat, gsMatrix<T>& stokesRhs)
{
    if (!m_isBaseReady)
        this->fillBaseSystem();

    stokesMat = m_baseMatrix;
    stokesRhs = m_baseRhs;
}


template<class T>
index_t gsINSAssembler<T>::numDofsUnk(index_t i)
{
    if (i == 0)
        return m_udofs;
    else if (i == 1)
        return m_pdofs;
    else
        GISMO_ERROR("numDofsUnk(i): i must be 0 or 1.");
}


template<class T>
gsSparseMatrix<T, RowMajor> gsINSAssembler<T>::getBlockUU() const
{
    if (m_isSystemReady)
        return m_matrix.topLeftCorner(m_pshift, m_pshift);
    else
        return m_blockUUlin + m_blockUUnonlin;
}


template<class T>
gsMatrix<T> gsINSAssembler<T>::getRhsU() const
{ 
    if (m_isSystemReady)
        return m_rhs.topRows(m_pshift);
    else
    {
        gsMatrix<T> rhsUpart = (m_rhsFG + m_rhsBtB).topRows(m_pshift);
        return (rhsUpart + m_rhsUlin + m_rhsUnonlin);
    }
}


template<class T>
gsMatrix<T> gsINSAssembler<T>::getRhsP() const
{
    if (m_isSystemReady)
        return m_rhs.bottomRows(m_pdofs);
    else
        return (m_rhsFG + m_rhsBtB).bottomRows(m_pdofs);
}

// =============================================================================


template<class T>
void gsINSAssemblerUnsteady<T>::initMembers()
{
    Base::initMembers();
    updateSizes();

    m_visitorTimeDiscr = gsINSVisitorUUtimeDiscr<T>(m_params);
    m_visitorTimeDiscr.initialize();
}


template<class T>
void gsINSAssemblerUnsteady<T>::updateSizes()
{
    Base::updateSizes();

    m_oldTimeVelField = m_currentVelField;

    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_rhsTimeDiscr.setZero(m_pshift, 1);

    // memory allocation
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.rows(), m_nnzPerRowU));
}


template<class T>
void gsINSAssemblerUnsteady<T>::updateDofMappers()
{
    Base::updateDofMappers();

    m_visitorTimeDiscr.updateDofMappers(m_dofMappers);
}


template<class T>
void gsINSAssemblerUnsteady<T>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::updateCurrentSolField(solVector, updateSol);

    if (updateSol)
        m_oldTimeVelField = this->constructSolution(solVector, 0);
}


template<class T>
void gsINSAssemblerUnsteady<T>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaning
    m_blockTimeDiscr.setZero();

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorTimeDiscr, 0, m_blockTimeDiscr, dummyRhs);

    m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);
}


template<class T>
void gsINSAssemblerUnsteady<T>::fillBaseSystem() 
{
    Base::fillBaseSystem();

    for (index_t row = 0; row < m_pshift; row++)
        for (typename gsSparseMatrix<T, RowMajor>::InnerIterator it(m_blockTimeDiscr, row); it; ++it)
            m_baseMatrix.coeffRef(row, it.col()) += it.value();

}


template<class T>
void gsINSAssemblerUnsteady<T>::fillSystem()
{
    Base::fillSystem();

    m_rhs.topRows(m_pshift) += m_rhsTimeDiscr;
}


template<class T>
void gsINSAssemblerUnsteady<T>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    gsFlowAssemblerBase<T>::update(solVector, updateSol);

    if(updateSol)
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


} // namespace gismo