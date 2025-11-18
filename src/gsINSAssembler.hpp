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

template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::initMembers()
{ 
    Base::initMembers();

    m_viscosity = m_paramsPtr->getPde().viscosity();

    m_ddof.resize(2);

    m_nnzPerOuterU = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerOuterU *= 2 * getBasis(0).maxDegree(i) + 1;

    m_nnzPerOuterP = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerOuterP *= 2 * getBasis(1).maxDegree(i) + 1;

    if (MatOrder == RowMajor)
        m_nnzPerOuterUP = m_nnzPerOuterP;
    else 
    {
        m_nnzPerOuterUP = 1;
        for (short_t i = 0; i < m_tarDim; i++)
            m_nnzPerOuterUP *= 3 * getBasis(0).maxDegree(i) - 1;

        m_nnzPerOuterUP *= m_tarDim;
    }

    m_visitorUUlin = gsINSVisitorUUlin<T, MatOrder>(m_paramsPtr);
    m_visitorUUlin.initialize();

    m_visitorUUnonlin = gsINSVisitorUUnonlin<T, MatOrder>(m_paramsPtr);
    m_visitorUUnonlin.initialize();
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);

    m_visitorUP = gsINSVisitorPU<T, MatOrder>(m_paramsPtr);
    m_visitorUP.initialize();

    m_visitorPU = gsINSVisitorUP<T, MatOrder>(m_paramsPtr);
    m_visitorPU.initialize();

    m_visitorF = gsINSVisitorRhsU<T, MatOrder>(m_paramsPtr);
    m_visitorF.initialize();

    m_visitorG = gsINSVisitorRhsP<T, MatOrder>(m_paramsPtr);
    m_visitorG.initialize();

    m_isMassMatReady = false;
    m_massMatBlocks.resize(2); // [velocity, pressure]
    m_massMatRhs.resize(2); // [velocity, pressure]

    updateSizes();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::updateSizes()
{
    if (m_paramsPtr->hasPeriodicBC())
        m_udofs = m_paramsPtr->getPerHelperPtr()->numFreeDofs();
    else
        m_udofs = this->getMapper(0).freeSize();

    m_pdofs = this->getMapper(1).freeSize();
    m_pshift = m_tarDim * m_udofs;
    m_dofs = m_pshift + m_pdofs;

    m_ddof[0].setZero(this->getMapper(0).boundarySize(), m_tarDim);
    m_ddof[1].setZero(this->getMapper(1).boundarySize(), 1);

    if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
    {
        this->computeDirichletDofs(0, 0, m_ddof[0]);
        this->computeDirichletDofs(1, 1, m_ddof[1]);
    }

    m_solution.setZero(m_dofs, 1);
    m_currentVelField = constructSolution(m_solution, 0);

    m_blockUUlin_comp.resize(m_udofs, m_udofs);
    m_blockUUnonlin_comp.resize(m_udofs, m_udofs);
    m_blockUUlin_whole.resize(m_pshift, m_pshift);
    m_blockUUnonlin_whole.resize(m_pshift, m_pshift);
    m_blockUP.resize(m_pshift, m_pdofs);
    m_blockPU.resize(m_pdofs, m_pshift);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_rhsUlin.setZero(m_pshift, 1);
    m_rhsUnonlin.setZero(m_pshift, 1);
    m_rhsPlin.setZero(m_pdofs, 1);
    m_rhsF.setZero(m_pshift, 1);
    m_rhsG.setZero(m_pdofs, 1);
    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);

    this->initParallel(); 
}

template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::initParallel()
{
    std::pair<index_t, index_t> locInfoU, locInfoP;

    petsc_computeMatLayout(m_udofs, locInfoU, m_paramsPtr->getMpiComm());
    petsc_computeMatLayout(m_pdofs, locInfoP, m_paramsPtr->getMpiComm());

    m_globalStartEnd.resize(2, 2);
    m_globalStartEnd(0, 0) = locInfoU.second; // first local velocity dof
    m_globalStartEnd(0, 1) = locInfoU.second + locInfoU.first; // last local velocity dof + 1
    m_globalStartEnd(1, 0) = locInfoP.second; // first local pressure dof
    m_globalStartEnd(1, 1) = locInfoP.second + locInfoP.first; // last local pressure dof + 1

    for (index_t unk = 0; unk <= 1; unk++)
    {
        for (index_t i = m_globalStartEnd(unk, 0); i < m_globalStartEnd(unk, 1); i++)
        {
            std::vector< std::pair<index_t, index_t> > localIDs;
            this->getMapper(unk).preImage(i, localIDs); // returns (patch, dof)

            for (size_t k = 0; k < localIDs.size(); k++)
                this->getMapper(unk).markTagged(localIDs[k].second, localIDs[k].first); // (dof, patch)
        }
    }

    m_ownedLocalDofs.resize(m_paramsPtr->getPde().patches().nPatches());
    for (size_t p = 0; p < m_ownedLocalDofs.size(); p++)
    {
        m_ownedLocalDofs[p].resize(2);

        for (index_t unk = 0; unk <= 1; unk++)
        {
            gsVector<index_t> ownedGlobalDofs = this->getMapper(unk).findTagged(p); // returns global indices
            gsVector<index_t> bndDofs = this->getMapper(unk).findBoundary(p); // returns local indices

            index_t numOwned = ownedGlobalDofs.size();
            index_t numBnd = bndDofs.size();
            m_ownedLocalDofs[p][unk].resize(numOwned + numBnd);

            for (index_t i = 0; i < numOwned; i++)
                m_ownedLocalDofs[p][unk](i) = this->globalToLocalOnPatch(p, unk, ownedGlobalDofs(i));

            for (index_t i = 0; i < numBnd; i++)
                m_ownedLocalDofs[p][unk](i + numOwned) = bndDofs(i);
        }
    }

    m_isParallelInitialized = true;

}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
        m_solution = solVector;

    m_currentVelField = constructSolution(solVector, 0, m_paramsPtr->isRotation()); // construct relative velocity if isRotation() = true
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::assembleLinearPart()
{
    // matrix and rhs cleaning
    m_blockUUlin_comp.resize(m_udofs, m_udofs);
    m_blockUUlin_whole.resize(m_pshift, m_pshift);
    m_blockUP.resize(m_pshift, m_pdofs);
    m_blockPU.resize(m_pdofs, m_pshift);
    m_rhsUlin.setZero();
    m_rhsPlin.setZero();
    m_rhsF.setZero();
    m_rhsG.setZero();

    bool hasPeriodicBC = m_paramsPtr->hasPeriodicBC();
    bool isRotation = m_paramsPtr->isRotation();

    if (hasPeriodicBC || isRotation)
    {
        // nnz multiplication:
        // rotation: velocity mass matrix is added into two off-diagonal blocks
        // periodic: a row/column of a diag. block is added into (some) rows/columns of (some) zero off-diagonal blocks
        // assuming periodicity/rotation wrt one of the coordinate axes

        gsVector<index_t> nnzVectorU = gsVector<index_t>::Constant(m_blockUUlin_whole.outerSize(), m_nnzPerOuterU);

        // if (isRotation)
        //     nnzVectorU *= 2;
        // else
        // {
        //     // this is not correct, TODO fix:
        //     std::vector<index_t> freePerDofsU = m_paramsPtr->getPerHelperPtr(0)->getFreePeriodicDofs();
        //     for (size_t i = 0; i < freePerDofsU.size(); i++)
        //         for (short_t d = 0; d < m_tarDim; d++)
        //             nnzVectorU(m_paramsPtr->getPerHelperPtr(0)->map(freePerDofsU[i]) + d*m_udofs) *=2;
        // }

        // TODO: improve
        nnzVectorU *= 2; // this is a lot more than needed for most of the rows/columns

        m_blockUUlin_whole.reserve(nnzVectorU);
        bool compressMat = isRotation ? false : true;

        this->assembleBlock(m_visitorUUlin, 0, m_blockUUlin_whole, m_rhsUlin, compressMat);
    
        if (isRotation)
        {
            gsINSVisitorUUrotation<T, MatOrder> visitorUUrot(m_paramsPtr);
            visitorUUrot.initialize();
            this->assembleBlock(visitorUUrot, 0, m_blockUUlin_whole, m_rhsUlin);
        }
    }
    else
    {
        m_blockUUlin_comp.reserve(gsVector<index_t>::Constant(m_blockUUlin_comp.outerSize(), m_nnzPerOuterU));
        this->assembleBlock(m_visitorUUlin, 0, m_blockUUlin_comp, m_rhsUlin);
    }

    m_blockUP.reserve(gsVector<index_t>::Constant(m_blockUP.outerSize(), m_nnzPerOuterUP));
    m_blockPU.reserve(gsVector<index_t>::Constant(m_blockPU.outerSize(), m_tarDim * m_nnzPerOuterU));
    this->assembleBlock(m_visitorUP, 0, m_blockUP, m_rhsUlin);
    this->assembleBlock(m_visitorPU, 1, m_blockPU, m_rhsPlin);
    this->assembleRhs(m_visitorF, 0, m_rhsF);

    if (m_paramsPtr->getPde().source()) // if the continuity eqn rhs is given
        this->assembleRhs(m_visitorG, 1, m_rhsG);

    // velocity and pressure mass matrix (needed for preconditioners)
    if ( m_paramsPtr->options().getString("lin.solver") == "iter" )
    {
        m_massMatRhs[0].setZero(m_pshift, 1);
        m_massMatRhs[1].setZero(m_pdofs, 1);

        gsINSVisitorUUmass<T, MatOrder> visitorUUmass(m_paramsPtr);
        gsINSVisitorPPmass<T, MatOrder> visitorPPmass(m_paramsPtr);
        visitorUUmass.initialize();
        visitorPPmass.initialize();
        
        index_t nnzU = m_nnzPerOuterU;
        index_t nnzP = m_nnzPerOuterP;

        if (hasPeriodicBC)
        {
            m_massMatBlocks[0].resize(m_pshift, m_pshift);

            // this is a lot more than needed for most of the rows / columns:
            nnzU *= 2; // TODO: improve
            nnzP *= 2; // TODO: improve
        }
        else
            m_massMatBlocks[0].resize(m_udofs, m_udofs);

        m_massMatBlocks[0].reserve(gsVector<index_t>::Constant(m_massMatBlocks[0].outerSize(), nnzU));
        this->assembleBlock(visitorUUmass, 0, m_massMatBlocks[0], m_massMatRhs[0]);
        
        m_massMatBlocks[1].resize(m_pdofs, m_pdofs);
        m_massMatBlocks[1].reserve(gsVector<index_t>::Constant(m_pdofs, nnzP));
        this->assembleBlock(visitorPPmass, 1, m_massMatBlocks[1], m_massMatRhs[1]);

        m_isMassMatReady = true;
    }
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::assembleNonlinearPart()
{
    // matrix and rhs cleaning
    m_blockUUnonlin_comp.resize(m_udofs, m_udofs);
    m_blockUUnonlin_whole.resize(m_pshift, m_pshift);
    m_rhsUnonlin.setZero();

    if (!m_paramsPtr->hasPeriodicBC())
    {
        m_blockUUnonlin_comp.reserve(gsVector<index_t>::Constant(m_blockUUnonlin_comp.outerSize(), m_nnzPerOuterU));
        this->assembleBlock(m_visitorUUnonlin, 0, m_blockUUnonlin_comp, m_rhsUnonlin);
    }
    else
    {
        // TODO: improve nnz
        m_blockUUnonlin_whole.reserve(gsVector<index_t>::Constant(m_blockUUnonlin_whole.outerSize(), 2 * m_nnzPerOuterU));
        this->assembleBlock(m_visitorUUnonlin, 0, m_blockUUnonlin_whole, m_rhsUnonlin);
    }
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::makeBlockUU(gsSparseMatrix<T, MatOrder>& result, bool linPartOnly)
{
    if (m_blockUUlin_comp.nonZeros() > 0)
    {
        gsSparseMatrix<T, MatOrder> blockUUcomp = m_blockUUlin_comp;

        if(!linPartOnly)
            blockUUcomp += m_blockUUnonlin_comp;

        gsVector<index_t> nonZerosVector(m_pshift);
        gsVector<index_t> nonZerosVector_UUcomp = getNnzVectorPerOuter(blockUUcomp);

        for (short_t d = 0; d < m_tarDim; d++)
            nonZerosVector.middleRows(d * m_udofs, m_udofs) = nonZerosVector_UUcomp;

        result.reserve(nonZerosVector);
        fillGlobalMat_UU(result, blockUUcomp);
    }
    
    result += m_blockUUlin_whole;

    if(!linPartOnly)
        result += m_blockUUnonlin_whole;

    result.makeCompressed();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::makeRhsU(gsMatrix<T>& result, bool linPartOnly)
{
    result = m_rhsF + m_rhsUlin;

    if(!linPartOnly)
        result += m_rhsUnonlin;
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillGlobalMat_UU(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    if (sourceMat.rows() == m_udofs) // sourceMat is a block for one velocity component
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                for (short_t d = 0; d < m_tarDim; d++)
                    globalMat.coeffRef(it.row() + d*m_udofs, it.col() + d*m_udofs) += it.value();
    }
    else // sourceMat is the whole UU block
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row(), it.col()) += it.value();
    }
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillGlobalMat_UP(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row(), it.col() + m_pshift) += it.value();

}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillGlobalMat_PU(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
            globalMat.coeffRef(it.row() + m_pshift, it.col()) += it.value();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillGlobalMat_PP(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
            globalMat.coeffRef(it.row() + m_pshift, it.col() + m_pshift) += it.value();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillBaseSystem() 
{
    gsVector<index_t> nonZerosVector(m_dofs);
    gsVector<index_t> nonZerosVector_UUcomp = getNnzVectorPerOuter(m_blockUUlin_comp);

    for (short_t d = 0; d < m_tarDim; d++)
        nonZerosVector.middleRows(d * m_udofs, m_udofs) = nonZerosVector_UUcomp;

    nonZerosVector.topRows(m_pshift) += getNnzVectorPerOuter(m_blockUUlin_whole);

    gsSparseMatrix<T, MatOrder> blockPU = gsSparseMatrix<T, MatOrder>(-m_blockUP.transpose());

    if (MatOrder == ColMajor)
    {
        nonZerosVector.topRows(m_pshift) += getNnzVectorPerOuter(blockPU);
        nonZerosVector.bottomRows(m_pdofs) = getNnzVectorPerOuter(m_blockUP);
    }
    else
    {
        nonZerosVector.topRows(m_pshift) += getNnzVectorPerOuter(m_blockUP);
        nonZerosVector.bottomRows(m_pdofs) = getNnzVectorPerOuter(blockPU);
    }

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_baseMatrix.reserve(nonZerosVector);

    fillGlobalMat_UU(m_baseMatrix, m_blockUUlin_comp);
    fillGlobalMat_UU(m_baseMatrix, m_blockUUlin_whole);
    fillGlobalMat_UP(m_baseMatrix, m_blockUP);
    fillGlobalMat_PU(m_baseMatrix, blockPU);

    m_baseRhs.topRows(m_pshift) += m_rhsF + m_rhsUlin;
    m_baseRhs.bottomRows(m_pdofs) += m_rhsG + m_rhsPlin;

    m_isBaseReady = true;
    m_isSystemReady = false;
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillSystem()
{
    m_matrix = m_baseMatrix;

    fillGlobalMat_UU(m_matrix, m_blockUUnonlin_comp);
    fillGlobalMat_UU(m_matrix, m_blockUUnonlin_whole);

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_rhs = m_baseRhs;
    m_rhs.topRows(m_pshift) += m_rhsUnonlin;

    m_isSystemReady = true;
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::initialize()
{
    Base::initialize();

    if (m_paramsPtr->isRotation())
        computeOmegaXrCoeffs();

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillStokesSystem(gsSparseMatrix<T, MatOrder>& stokesMat, gsMatrix<T>& stokesRhs)
{
    if (!m_isBaseReady)
        this->fillBaseSystem();

    stokesMat = m_baseMatrix;
    stokesRhs = m_baseRhs;

    if (!stokesMat.isCompressed())
        stokesMat.makeCompressed();
}


template<class T, int MatOrder>
gsField<T> gsINSAssembler<T, MatOrder>::constructSolution(const gsMatrix<T>& solVector, index_t unk, bool customSwitch) const
{
    // if there is no rotation, absolute and relative velocity are equal, m_omegaXrCoeffs not computed
    if (!m_paramsPtr->isRotation())
        customSwitch = false;

    index_t fullUdofs = this->getMapper(0).freeSize();
    index_t fullPdofs = this->getMapper(1).freeSize();
    index_t fullPshift = m_tarDim * fullUdofs;
    index_t fullDofs = fullPshift + fullPdofs;

    gsMatrix<T> solution;

    if (m_paramsPtr->hasPeriodicBC() && solVector.rows() == m_dofs) // the passed solVector is periodic
        per2nonper_into(solVector, solution);
    else
    {
        GISMO_ASSERT(solVector.rows() == fullDofs, "Something went wrong, is solution vector valid?");
        solution = solVector;
    }

    typename gsMultiPatch<T>::Ptr result(new gsMultiPatch<T>);

    const gsDofMapper& mapper = this->getMapper(unk);
    gsDofMapper rotMapper; // mapper for omega x r coeffs

    if (unk == 0 && customSwitch)
        this->getBasis(0).getMapper(getAssemblerOptions().intStrategy, rotMapper);

    const index_t dim = (unk == 0 ? m_tarDim : 1);
    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    gsAsConstMatrix<T> solV = (unk == 0 ?
        gsAsConstMatrix<T>(solution.data(), fullUdofs, dim)
        :
        gsAsConstMatrix<T>(solution.data() + fullPshift, fullPdofs, 1)
        );

    for (size_t p = 0; p < this->getPatches().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const index_t sz = this->getBasis(unk).piece(p).size();
        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, p)) // DoF value is in the solution vector
                coeffs.row(i) = solV.row(mapper.index(i, p));
            else // eliminated DoF: fill with Dirichlet data
                coeffs.row(i) = m_ddof[unk].row(mapper.bindex(i, p));

            if (unk == 0 && customSwitch)
                coeffs.row(i) -= m_omegaXrCoeffs.row(rotMapper.index(i, p));
        }

        result->addPatch(this->getBasis(unk).piece(p).makeGeometry(coeffs));
    }

    return gsField<T>(this->getPatches(), result, true);
}


template<class T, int MatOrder>
T gsINSAssembler<T, MatOrder>::computeFlowRate(index_t patch, boxSide side, gsMatrix<T> solution) const
{
    T flowRate = 0;

    gsField<T> solutionField = constructSolution(solution, 0); // velocity field

    const gsGeometry<T>& geo = this->getPatches().patch(patch);
    const gsBasis<T>& basis = this->getBasis(0).basis(patch);

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
    mapData.side = side;

    //typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(side);
    typename gsBasis<T>::domainIter domIt = basis.domain()->beginBdr(side);
    typename gsBasis<T>::domainIter domItEnd = basis.domain()->endBdr(side);
    for (; domIt!=domItEnd; ++domIt)
    {
        // Compute the quadrature rule on patch1
        QuRule.mapTo(domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights);

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


template<class T, int MatOrder>
index_t gsINSAssembler<T, MatOrder>::numDofsUnk(index_t i) const
{
    if (i == 0)
        return m_udofs;
    else if (i == 1)
        return m_pdofs;
    else
        GISMO_ERROR("numDofsUnk(i): i must be 0 or 1.");
}


template<class T, int MatOrder>
gsSparseMatrix<T, MatOrder> gsINSAssembler<T, MatOrder>::getBlockUU(bool linPartOnly)
{
    if (m_isSystemReady && !linPartOnly)
    {
        return m_matrix.topLeftCorner(m_pshift, m_pshift);
    }
    else if (m_isBaseReady && linPartOnly)
    {
        return m_baseMatrix.topLeftCorner(m_pshift, m_pshift);
    }
    else
    {
        gsSparseMatrix<T, MatOrder> blockUU(m_pshift, m_pshift);
        makeBlockUU(blockUU, linPartOnly);
        return blockUU;
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsINSAssembler<T, MatOrder>::getRhsU(bool linPartOnly)
{ 
    if (m_isSystemReady && !linPartOnly)
        return m_rhs.topRows(m_pshift);
    else if (m_isBaseReady && linPartOnly)
        return m_baseRhs.topRows(m_pshift);
    else
    {
        gsMatrix<T> rhsU;
        makeRhsU(rhsU, linPartOnly);
        return rhsU;
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsINSAssembler<T, MatOrder>::getRhsP() const
{
    if (m_isSystemReady)
        return m_rhs.bottomRows(m_pdofs);
    else
        return m_rhsG + m_rhsPlin;
}


// --- periodic BC functions ---

template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::nonper2per_into(const gsMatrix<T>& fullVector, gsMatrix<T>& perVector) const
{
    perVector.setZero(m_dofs, 1);

    index_t fullUdofs = this->getMapper(0).freeSize();
    index_t fullPshift = m_tarDim * fullUdofs;

    // mapping velocity DOFs
    for (index_t row = 0; row < fullPshift; row++)
    {
        if (!isEliminatedPeriodic(row))
            perVector(mapPeriodic(row), 0) += fullVector(row, 0);
        else
        {
            for (int d = 0; d < m_tarDim; d++)
                perVector( m_paramsPtr->getPerHelperPtr()->map(row % fullUdofs) + d * m_udofs, 0) += m_paramsPtr->getPde().bc().getTransformMatrix()(d, row / fullUdofs) * fullVector(row, 0);
        }
    }

    // pressure DOFS
    perVector.bottomRows(m_pdofs) = fullVector.bottomRows(m_pdofs);
}

template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::per2nonper_into(const gsMatrix<T>& perVector, gsMatrix<T>& fullVector) const
{
    index_t fullUdofs = this->getMapper(0).freeSize();
    index_t fullPshift = m_tarDim * fullUdofs;
    index_t fullDofs = fullPshift + m_pdofs;
    
    fullVector.setZero(fullDofs, 1);

    for (index_t row = 0; row < m_dofs; row++)
        fullVector(invMapPeriodic(row), 0) += perVector(row, 0);

    std::vector<index_t> elimDofsU = m_paramsPtr->getPerHelperPtr()->getElimPeriodicDofs();

    for (size_t i = 0; i < elimDofsU.size(); i++)
    {
        gsMatrix<T> u(m_tarDim, 1);
        for (int d = 0; d < m_tarDim; d++)
            u(d, 0) = perVector(m_paramsPtr->getPerHelperPtr()->map(elimDofsU[i]) + d * m_udofs, 0);

        u = m_paramsPtr->getPde().bc().getTransformMatrix().transpose() * u; // multiplication by an inverse tranform matrix

        for (int d = 0; d < m_tarDim; d++)
            fullVector(elimDofsU[i] + d * fullUdofs, 0) = u(d, 0);
    }
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::computeOmegaXrCoeffs()
{
    real_t omega = m_paramsPtr->options().getReal("omega");

    gsFunctionExpr<T> omegaXr;
    std::ostringstream s1, s2;

    s1 << -omega << " * y";
    s2 << omega << " * x";

    switch (m_tarDim)
    {
    case 2:
        omegaXr = gsFunctionExpr<T>(s1.str(), s2.str(), 2);
        break;

    case 3:
        omegaXr = gsFunctionExpr<T>(s1.str(), s2.str(), "0", 3);
        break;

    default:
        GISMO_ERROR("Wrong space dimension.");
        break;
    }

    gsDofMapper mapperOxR;
    getBasis(0).getMapper(getAssemblerOptions().intStrategy, mapperOxR);
    int ndofsOxR = mapperOxR.freeSize();

    gsBoundaryConditions<T> dummyBC;
    std::vector<gsMatrix<T> > dummyDirichletDofs;
    gsMatrix<T> dummyRhs;

    gsNavStokesPde<real_t> pde(this->getPatches(), dummyBC, &omegaXr, m_viscosity);
    gsFlowSolverParams<real_t> params(pde, m_paramsPtr->getBases());

    gsSparseMatrix<T, MatOrder> projMatrix(ndofsOxR, ndofsOxR);
    projMatrix.reserve(gsVector<int>::Constant(ndofsOxR, m_nnzPerOuterU));
    gsMatrix<T> projRhs(ndofsOxR * m_tarDim, 1);
    projRhs.setZero();

    gsINSVisitorUUmass<T, MatOrder> visitorMass(memory::make_shared_not_owned(&params));
    gsINSVisitorRhsU<T, MatOrder> visitorRhs(memory::make_shared_not_owned(&params));
    visitorMass.initialize();
    visitorRhs.initialize();

    // TODO: parallelize
    for(size_t p = 0; p < this->getPatches().nPatches() ; p++)
    {
        visitorMass.initOnPatch(p);
        visitorRhs.initOnPatch(p);

        typename gsBasis<T>::domainIter domIt = m_paramsPtr->getBasis(0).piece(p).domain()->beginAll();
        typename gsBasis<T>::domainIter domItEnd = m_paramsPtr->getBasis(0).piece(p).domain()->endAll();

        while (domIt!=domItEnd)
        {
            visitorMass.evaluate(domIt.get());
            visitorMass.assemble();
            visitorMass.localToGlobal(dummyDirichletDofs, projMatrix, dummyRhs);

            visitorRhs.evaluate(domIt.get());
            visitorRhs.assemble();
            visitorRhs.localToGlobal(projRhs);

            ++domIt;
        }
    }

    projMatrix.makeCompressed();

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT linSolver;
    #else
    typename gsSparseSolver<T>::SimplicialLDLT linSolver;
    #endif

    linSolver.analyzePattern(projMatrix);
    linSolver.factorize(projMatrix);

    m_omegaXrCoeffs.resize(ndofsOxR, m_tarDim);

    for (short_t i = 0; i < m_tarDim; i++)
        m_omegaXrCoeffs.col(i) = linSolver.solve(projRhs.middleRows(i*ndofsOxR, ndofsOxR));

}


// =============================================================================


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();

    m_visitorTimeDiscr = gsINSVisitorUUtimeDiscr<T, MatOrder>(m_paramsPtr);
    m_visitorTimeDiscr.initialize();
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::updateSizes()
{
    Base::updateSizes();

    m_oldTimeVelField = m_currentVelField;

    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_rhsTimeDiscr.setZero(m_pshift, 1);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::updateCurrentSolField(solVector, updateSol);

    if (updateSol)
        m_oldTimeVelField = this->constructSolution(solVector, 0, m_paramsPtr->isRotation());
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaningm_nnzPerOuterU
    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.outerSize(), m_nnzPerOuterU));

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorTimeDiscr, 0, m_blockTimeDiscr, dummyRhs);

    m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    this->fillGlobalMat_UU(m_matrix, m_blockTimeDiscr);
    m_rhs.topRows(m_pshift) += m_rhsTimeDiscr;
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


} // namespace gismo
