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

    m_viscosity = m_params.getPde().viscosity();

    m_dofMappers.resize(2);
    m_ddof.resize(2);
    getBases().front().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.front(), 0);
    getBases().back().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.back(), 1);

    m_nnzPerRowU = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowU *= 2 * getBases().front().maxDegree(i) + 1;

    m_nnzPerRowP = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowP *= 2 * getBases().back().maxDegree(i) + 1;

    updateSizes();

    m_visitorUUlin = gsINSVisitorUUlin<T, MatOrder>(m_params);
    m_visitorUUlin.initialize();

    m_visitorUUnonlin = gsINSVisitorUUnonlin<T, MatOrder>(m_params);
    m_visitorUUnonlin.initialize();
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);

    m_visitorUP = gsINSVisitorPU_withUPrhs<T, MatOrder>(m_params);
    m_visitorUP.initialize();

    m_visitorF = gsINSVisitorRhsU<T, MatOrder>(m_params);
    m_visitorF.initialize();

    m_visitorG = gsINSVisitorRhsP<T, MatOrder>(m_params);
    m_visitorG.initialize();

    m_isMassMatReady = false;
    m_massMatBlocks.resize(2); // [velocity, pressure]
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::updateSizes()
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

    m_blockUUlin_comp.resize(m_udofs, m_udofs);
    m_blockUUnonlin_comp.resize(m_udofs, m_udofs);
    m_blockUUlin_whole.resize(m_pshift, m_pshift);
    m_blockUUnonlin_whole.resize(m_pshift, m_pshift);
    m_blockUP.resize(m_pshift, m_pdofs);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_rhsUlin.setZero(m_pshift, 1);
    m_rhsUnonlin.setZero(m_pshift, 1);
    m_rhsBtB.setZero(m_dofs, 1);
    m_rhsFG.setZero(m_dofs, 1);
    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::updateDofMappers()
{
    m_visitorUUlin.updateDofMappers(m_dofMappers);
    m_visitorUUnonlin.updateDofMappers(m_dofMappers);
    m_visitorUP.updateDofMappers(m_dofMappers);
    m_visitorF.updateDofMappers(m_dofMappers);
    m_visitorG.updateDofMappers(m_dofMappers);
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
        m_solution = solVector;

    m_currentVelField = constructSolution(solVector, 0);
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::assembleLinearPart()
{
    // matrix and rhs cleaning
    m_blockUUlin_comp.resize(m_udofs, m_udofs);
    m_blockUP.resize(m_pshift, m_pdofs);
    m_rhsUlin.setZero();
    m_rhsBtB.setZero();
    m_rhsFG.setZero();

    int nnzPerOuterUP = 0;

    if (MatOrder == RowMajor)
        nnzPerOuterUP = m_nnzPerRowP;
    else 
        nnzPerOuterUP = m_tarDim * m_nnzPerRowU;

    // memory allocation
    m_blockUUlin_comp.reserve(gsVector<index_t>::Constant(m_blockUUlin_comp.outerSize(), m_nnzPerRowU));
    m_blockUP.reserve(gsVector<index_t>::Constant(m_blockUP.outerSize(), nnzPerOuterUP));

    this->assembleBlock(m_visitorUUlin, 0, m_blockUUlin_comp, m_rhsUlin);
    this->assembleBlock(m_visitorUP, 0, m_blockUP, m_rhsBtB);
    this->assembleRhs(m_visitorF, 0, m_rhsFG);

    if(m_params.getPde().source()) // if the continuity eqn rhs is given
        this->assembleRhs(m_visitorG, 1, m_rhsFG);

    // mass matrices for velocity and pressure (needed for preconditioners)
    if ( m_params.options().getString("lin.solver") == "iter" )
    {
        gsMatrix<T> dummyRhsU(m_pshift, 1);
        gsMatrix<T> dummyRhsP(m_pdofs, 1);

        gsINSVisitorUUmass<T, MatOrder> visitorUUmass(m_params);
        gsINSVisitorPPmass<T, MatOrder> visitorPPmass(m_params);
        
        visitorUUmass.initialize();
        visitorPPmass.initialize();

        m_massMatBlocks[0].resize(m_udofs, m_udofs);
        m_massMatBlocks[1].resize(m_pdofs, m_pdofs);
        m_massMatBlocks[0].reserve(gsVector<index_t>::Constant(m_udofs, m_nnzPerRowU));
        m_massMatBlocks[1].reserve(gsVector<index_t>::Constant(m_pdofs, m_nnzPerRowP));

        this->assembleBlock(visitorUUmass, 0, m_massMatBlocks[0], dummyRhsU);
        this->assembleBlock(visitorPPmass, 0, m_massMatBlocks[1], dummyRhsP);

        m_isMassMatReady = true;
        

        // linear operators needed for PCD preconditioner
        // if ( m_paramsRef.options().getString("lin.precType").substr(0, 3) == "PCD" )
        // {
        //     // pressure Laplace operator

        //     gsINSVisitorPPlaplace<T, MatOrder> visitorPPlaplace(m_params);
        //     visitorPPlaplace.initialize();

        //     m_pcdBlocks[0].resize(m_pdofs, m_pdofs);
        //     m_pcdBlocks[0].reserve(gsVector<index_t>::Constant(m_pdofs, m_nnzPerRowP));

        //     this->assembleBlock(visitorPPlaplace, 0, m_pcdBlocks[0], dummyRhsP);

        //     // prepare for BCs
        //     findPressureBoundaryIDs();
        // }
    }
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::assembleNonlinearPart()
{
    // matrix and rhs cleaning
    m_blockUUnonlin_comp.resize(m_udofs, m_udofs);
    m_rhsUnonlin.setZero();

    // memory allocation
    m_blockUUnonlin_comp.reserve(gsVector<index_t>::Constant(m_blockUUnonlin_comp.outerSize(), m_nnzPerRowU));

    this->assembleBlock(m_visitorUUnonlin, 0, m_blockUUnonlin_comp, m_rhsUnonlin);

    // linear operators needed for PCD preconditioner
    // if ( m_params.options().getString("lin.solver") == "iter" &&  m_paramsRef.options().getString("lin.precType").substr(0, 3) == "PCD" )
    // {
    //     gsMatrix<T> dummyRhsP(m_pdofs, 1);

    //     // Robin BC

    //     // TODO

    //     // pressure convection operator

    //     gsINSVisitorPPconvection<T, MatOrder> visitorPPconv(m_params);
    //     visitorPPconv.initialize();
    //     visitorPPconv.setCurrentSolution(m_currentVelField);

    //     m_pcdBlocks[2].resize(m_pdofs, m_pdofs);
    //     m_pcdBlocks[2].reserve(gsVector<index_t>::Constant(m_pdofs, m_nnzPerRowP));

    //     this->assembleBlock(visitorPPconv, 0, m_pcdBlocks[2], dummyRhsP);
    // }
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

    m_baseRhs.noalias() = m_rhsFG + m_rhsBtB;
    m_baseRhs.topRows(m_pshift) += m_rhsUlin;

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

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
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


template<class T, int MatOrder>
void gsINSAssembler<T, MatOrder>::fillStokesSystem(gsSparseMatrix<T, MatOrder>& stokesMat, gsMatrix<T>& stokesRhs)
{
    if (!m_isBaseReady)
        this->fillBaseSystem();

    stokesMat = m_baseMatrix;
    stokesRhs = m_baseRhs;
}


template<class T, int MatOrder>
gsField<T> gsINSAssembler<T, MatOrder>::constructSolution(const gsMatrix<T>& solVector, index_t unk) const
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


template<class T, int MatOrder>
T gsINSAssembler<T, MatOrder>::computeFlowRate(index_t patch, boxSide side, gsMatrix<T> solution) const
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


template<class T, int MatOrder>
index_t gsINSAssembler<T, MatOrder>::numDofsUnk(index_t i)
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
        gsSparseMatrix<T, MatOrder> blockUUcomp = m_blockUUlin_comp;

        if(!linPartOnly)
            blockUUcomp += m_blockUUnonlin_comp;

        gsSparseMatrix<T, MatOrder> blockUU(m_pshift, m_pshift);

        gsVector<index_t> nonZerosVector(m_pshift);
        gsVector<index_t> nonZerosVector_UUcomp = getNnzVectorPerOuter(blockUUcomp);

        for (short_t d = 0; d < m_tarDim; d++)
            nonZerosVector.middleRows(d * m_udofs, m_udofs) = nonZerosVector_UUcomp;

        blockUU.reserve(nonZerosVector);
        fillGlobalMat_UU(blockUU, blockUUcomp);

        blockUU += m_blockUUlin_whole;

        if(!linPartOnly)
            blockUU += m_blockUUnonlin_whole;

        blockUU.makeCompressed();

        return blockUU;
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsINSAssembler<T, MatOrder>::getRhsU() const
{ 
    if (m_isSystemReady)
        return m_rhs.topRows(m_pshift);
    else
    {
        gsMatrix<T> rhsUpart = (m_rhsFG + m_rhsBtB).topRows(m_pshift);
        return (rhsUpart + m_rhsUlin + m_rhsUnonlin);
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsINSAssembler<T, MatOrder>::getRhsP() const
{
    if (m_isSystemReady)
        return m_rhs.bottomRows(m_pdofs);
    else
        return (m_rhsFG + m_rhsBtB).bottomRows(m_pdofs);
}


// --- PCD functions ---

// template<class T, int MatOrder>
// void gsINSAssembler<T, MatOrder>::findPressureBoundaryIDs()
// {
//     findPressureBoundaryPartIDs(m_params.getBndIn(), m_presInIDs);
//     findPressureBoundaryPartIDs(m_params.getBndOut(), m_presOutIDs);
//     findPressureBoundaryPartIDs(m_params.getBndWall(), m_presWallIDs);
// }

// template<class T, int MatOrder>
// void gsINSAssembler<T, MatOrder>::findPressureBoundaryPartIDs(std::vector<std::pair<int, boxSide> > bndPart, std::vector<index_t>& idVector)
//     {
//         const gsMultiBasis<T>& pBasis = getBases()[1];

//         gsMatrix<index_t> bnd;
//         idVector.clear();

//         for (std::vector<std::pair<int, boxSide> >::iterator it = bndPart.begin(); it != bndPart.end(); it++)
//         {
//             bnd = pBasis.piece(it->first).boundary(it->second);

//             for (int i = 0; i < bnd.rows(); i++)
//                 idVector.push_back(this->getMappers()[1].index(bnd(i, 0), it->first));
//         }
//     }

// =============================================================================


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::initMembers()
{
    Base::initMembers();
    updateSizes();

    m_visitorTimeDiscr = gsINSVisitorUUtimeDiscr<T, MatOrder>(m_params);
    m_visitorTimeDiscr.initialize();
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::updateSizes()
{
    Base::updateSizes();

    m_oldTimeVelField = m_currentVelField;

    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_rhsTimeDiscr.setZero(m_pshift, 1);

    // memory allocation
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.outerSize(), m_nnzPerRowU));
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::updateDofMappers()
{
    Base::updateDofMappers();

    m_visitorTimeDiscr.updateDofMappers(m_dofMappers);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::updateCurrentSolField(solVector, updateSol);

    if (updateSol)
        m_oldTimeVelField = this->constructSolution(solVector, 0);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix cleaning
    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.outerSize(), m_nnzPerRowU));

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorTimeDiscr, 0, m_blockTimeDiscr, dummyRhs);

    m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::fillBaseSystem() 
{
    Base::fillBaseSystem();
    this->fillGlobalMat_UU(m_baseMatrix, m_blockTimeDiscr);
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::fillSystem()
{
    Base::fillSystem();

    m_rhs.topRows(m_pshift) += m_rhsTimeDiscr;
}


template<class T, int MatOrder>
void gsINSAssemblerUnsteady<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


} // namespace gismo