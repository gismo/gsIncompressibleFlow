/** @file gsINSBlockAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Hornikova, J. Sourek
*/

#pragma once
#include <gsIncompressibleFlow/src/gsINSBlockAssembler.h>
#include <gsIncompressibleFlow/src/gsINSUtils.h>

#include <gsAssembler/gsAssembler.h>
#include <gsCore/gsFuncData.h>
#include <gsAssembler/gsGaussRule.h>

namespace gismo
{

template<class T>
void gsINSBlockAssembler<T>::initMembers()
{
    m_udofs = m_dofMappers.front().freeSize();
    m_pdofs = m_dofMappers.back().freeSize();
    m_pshift = m_tarDim * m_udofs;
    m_dofs = m_tarDim * m_udofs + m_pdofs;

    m_nonZerosPerColU = 1;
    for (int i = 0; i < m_tarDim; i++)
        m_nonZerosPerColU *= 2 * getBases().front().maxDegree(i) + 1;

    m_nonZerosPerColP = 1;
    for (int i = 0; i < m_tarDim; i++)
        m_nonZerosPerColP *= 2 * getBases().back().maxDegree(i) + 1;

    m_ddof[0].setZero(m_dofMappers.front().boundarySize(), m_tarDim);
    m_ddof[1].setZero(m_dofMappers.back().boundarySize(), 1);

    m_blockA.resize(m_udofs, m_udofs);

    m_blockB.resize(m_tarDim);
    m_blockMinusBT.resize(m_tarDim);
    for (int s = 0; s < m_tarDim; ++s)
    {
        m_blockB[s].resize(m_pdofs, m_udofs);
        m_blockMinusBT[s].resize(m_udofs, m_pdofs);
    }

    m_blockCT.resize(m_tarDim);
    for (int s = 0; s < m_tarDim; ++s)
    {
        m_blockCT[s].resize(m_udofs, m_pdofs);
    }

    m_blockN.resize(m_udofs, m_udofs);
    m_blockNpattern.resize(m_udofs, m_udofs);
    m_blockM.resize(m_udofs, m_udofs);
    m_blockAp.resize(m_pdofs, m_pdofs);
    m_blockMp.resize(m_pdofs, m_pdofs);

    m_rhsA.setZero(m_dofs, 1);
    m_rhsB.setZero(m_dofs, 1);
    m_rhsC.setZero(m_dofs, 1);
    m_rhsN.setZero(m_dofs, 1);
    m_rhsM.setZero(m_dofs, 1);
    m_rhsAp.setZero(m_dofs, 1);
    m_rhsMp.setZero(m_dofs, 1);
    m_rhsF.setZero(m_dofs, 1);

    m_solution.setZero(m_dofs, 1);

    // if (getAssemblerOptions().intStrategy == iFace::dg)
    // {
    //     m_blockPdg.resize(m_pdofs, m_pdofs);
    //     m_rhsPdg.setZero(m_dofs, 1);
    // }

    if (getAssemblerOptions().dirStrategy == dirichlet::elimination)
    {
        computeDirichletDofs(0, 0, m_ddof[0]);
        computeDirichletDofs(1, 1, m_ddof[1]);
    }

    m_currentVelField = constructSolution(m_solution, 0);
    m_oldTimeVelField = m_currentVelField;

    m_bPCDbndPrepared = false;
}

template<class T>
void gsINSBlockAssembler<T>::constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result, int unk) const
{
    GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

    const gsDofMapper & mapper = m_dofMappers[unk];

    result.clear(); // result is cleared first

    const index_t dim = (unk == 0 ? m_tarDim : 1);
    gsMatrix<T> coeffs;

    // Point to the correct entries of the solution vector
    gsAsConstMatrix<T> solV = (unk == 0 ?
        gsAsConstMatrix<T>(solVector.data(), m_udofs, dim)
        :
        gsAsConstMatrix<T>(solVector.data() + m_pshift, m_pdofs, 1)
        );

    for (unsigned int p = 0; p < getPatches().nPatches(); ++p)
    {
        // Reconstruct solution coefficients on patch p
        const int sz = getBases().at(unk).piece(p).size();
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

        result.addPatch(getBases().at(unk).piece(p).makeGeometry(coeffs));
    }
}


template<class T>
gsField<T> gsINSBlockAssembler<T>::constructSolution(const gsMatrix<T>& solVector, int unk) const
{
    gsMultiPatch<T> * result = new gsMultiPatch<T>;
    constructSolution(solVector, *result, unk);
    return gsField<T>(getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
}

template<class T>
void gsINSBlockAssembler<T>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
    {
        m_solution = solVector;

        if(isUnsteady())
            m_oldTimeVelField = constructSolution(solVector, 0);
    }

    m_currentVelField = constructSolution(solVector, 0);
}


template<class T>
T gsINSBlockAssembler<T>::computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const
{
    T flowRate = 0;

    gsField<T> solutionField = constructSolution(solution, 0); // velocity field

    const gsGeometry<T>& geo = getPatches().patch(patch);
    const gsBasis<T>& basis = getBases().at(0).basis(patch);

    gsVector<int> numQuadNodes(m_tarDim);
    const int dir = side.direction();
    for (int i = 0; i < m_tarDim; ++i)
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
void gsINSBlockAssembler<T>::assembleLinearStokesPart()
{
    // matrix and rhs cleaning
    m_rhsA.setZero();
    m_blockA.resize(m_udofs, m_udofs);
    gsSparseMatrix<T> blockAsym(m_udofs, m_udofs);

    m_rhsB.setZero();

    gsSparseMatrix<T> blocksB(m_pdofs, m_pshift);

    // if (getAssemblerOptions().intStrategy == iFace::dg)
    // {
    //     m_rhsPdg.setZero();
    //     m_blockPdg.resize(m_pdofs, m_pdofs);
    // }

    m_rhsF.setZero();

    // blocks assembly
    blockAsym.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    Base::template assembleBlock< gsINSBlockAVisitor<T> >(blockAsym, m_rhsA, m_nonZerosPerColU);
    m_blockA = blockAsym.template selfadjointView<gsEigen::Upper>();
    // if (getAssemblerOptions().intStrategy == iFace::dg) 
    // {
    //     gsSparseMatrix<T> blockAdg(m_udofs, m_udofs);
    //     blockAdg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    //     Base::template assembleBlockDg< gsINSBlockAVisitorDG<T> >(blockAdg, m_rhsA, m_nonZerosPerColU);
    //     m_blockA += blockAdg;
    // }

    m_blockA.makeCompressed();

    blocksB.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
    Base::template assembleBlock< gsINSBlocksBVisitor<T> >(blocksB, m_rhsB, m_nonZerosPerColP);
    // if (getAssemblerOptions().intStrategy == iFace::dg) 
    // {
    //     gsSparseMatrix<T> blocksBdg(m_pdofs, m_pshift);
    //     blocksBdg.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
    //     Base::template assembleBlockDg< gsINSBlockBVisitorDG<T> >(blocksBdg, m_rhsB, m_nonZerosPerColP);
    //     blocksB += blocksBdg;
    // }
    blocksB.makeCompressed();

    for (int s = 0; s < m_tarDim; ++s) {
        m_blockB[s] = gsSparseMatrix<T>(blocksB.middleCols(s * m_udofs, m_udofs));
        m_blockB[s].makeCompressed();
    }

    for (int s = 0; s < m_tarDim; ++s) {
        m_blockMinusBT[s] = gsSparseMatrix<T, gsEigen::ColMajor>(-m_blockB[s].transpose());
        m_blockMinusBT[s].makeCompressed();
    }

    // if (getAssemblerOptions().intStrategy == iFace::dg)
    // {
    //     m_blockPdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    //     Base::template assembleBlockDg< gsINSBlockPVisitorDG<T> >(m_blockPdg, m_rhsPdg, m_nonZerosPerColP);
    //     m_blockPdg.makeCompressed();
    // }

    Base::template assembleRhs< gsINSRhsVisitor<T> >(m_rhsF);

    if (isUnsteady())
    {
        assembleMassMatrix();
    }
}


template<class T>
void gsINSBlockAssembler<T>::assembleNonlinearPart()
{
    // matrix and rhs cleaning
    m_blockN.resize(m_udofs, m_udofs);
    m_rhsN.setZero();

    // blocks assembly
    m_blockN.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    Base::template assembleBlock< gsINSBlockNVisitor<T> >(m_blockN, m_rhsN, m_nonZerosPerColU);
    // if (getAssemblerOptions().intStrategy == iFace::dg) {
    //     gsSparseMatrix<T> m_blockNdg(m_udofs, m_udofs);
    //     m_blockNdg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    //     Base::template assembleBlockDg< gsINSBlockNVisitorDG<T> >(m_blockNdg, m_rhsN, m_nonZerosPerColU);
    //     m_blockN += m_blockNdg;
    // }
    // m_blockN.makeCompressed();
}


template<class T>
void gsINSBlockAssembler<T>::assembleMassMatrix()
{
    // matrix and rhs cleaning
    m_rhsM.setZero();
    m_blockM.resize(m_udofs, m_udofs);
    gsSparseMatrix<T> blockMsym(m_udofs, m_udofs);

    // blocks assembly
    blockMsym.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    Base::template assembleBlock< gsINSBlockMVisitor<T> >(blockMsym, m_rhsM, m_nonZerosPerColU);
    m_blockM = blockMsym.template selfadjointView<gsEigen::Upper>();
    m_blockM.makeCompressed();
}


template<class T>
void gsINSBlockAssembler<T>::assemblePressureMassMatrix()
{
    // matrix and rhs cleaning
    m_rhsMp.setZero();
    m_blockMp.resize(m_pdofs, m_pdofs);
    gsSparseMatrix<T> blockMpsym(m_pdofs, m_pdofs);

    // blocks assembly
    blockMpsym.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    Base::template assembleBlock< gsINSBlockMpVisitor<T> >(blockMpsym, m_rhsMp, m_nonZerosPerColP);
    m_blockMp = blockMpsym.template selfadjointView<gsEigen::Upper>();
    m_blockMp.makeCompressed();
}


template<class T>
void gsINSBlockAssembler<T>::assembleBlockNpattern()
{
    gsMatrix<T> current_solution = m_solution;

    m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

    // matrix and rhs cleaning
    m_blockNpattern.resize(m_udofs, m_udofs);

    gsMatrix<T> dummy_rhsN(m_dofs, 1);
    dummy_rhsN.setZero();

    // blocks assembly
    m_blockNpattern.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    Base::template assembleBlock< gsINSBlockNVisitor<T> >(m_blockNpattern, dummy_rhsN, m_nonZerosPerColU);
    // if (getAssemblerOptions().intStrategy == iFace::dg) {
    //     gsSparseMatrix<T> blockNpatterndg(m_udofs, m_udofs);
    //     blockNpatterndg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
    //     Base::template assembleBlockDg< gsINSBlockNVisitorDG<T> >(blockNpatterndg, dummy_rhsN, m_nonZerosPerColU);
    //     m_blockNpattern += blockNpatterndg;
    // }
    m_blockNpattern.makeCompressed();
    m_blockNpattern = 0. * m_blockNpattern;

    m_solution = current_solution;
}


template<class T>
void gsINSBlockAssembler<T>::fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
{
    gsVector<int> nonZerosPerColumnVector;
    nonZerosPerColumnVector.setZero(m_dofs);
    for (int s = 0; s < m_tarDim; ++s)
        for (int i = 0; i < m_udofs; i++)
        {
            nonZerosPerColumnVector(i + s * m_udofs) = m_blockA.col(i).nonZeros() + m_blockB[s].col(i).nonZeros();
        }

    for (int i = 0; i < m_pdofs; i++)
        for (int s = 0; s < m_tarDim; ++s)
            nonZerosPerColumnVector(i + m_pshift) += m_blockMinusBT[s].col(i).nonZeros();

    // if (getAssemblerOptions().intStrategy == iFace::dg)
    // {
    //     for (int i = 0; i < m_pdofs; i++)
    //         nonZerosPerColumnVector(i + m_pshift) += m_blockPdg.col(i).nonZeros();
    // }

    stokesMatrix.resize(m_dofs, m_dofs); // Clean matrix
    stokesMatrix.reserve(nonZerosPerColumnVector);

#pragma omp parallel for num_threads(m_numThreads)
    for (index_t col = 0; col < m_udofs; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_blockA, col); it; ++it)
            for (index_t s = 0; s < m_tarDim; ++s)
                stokesMatrix.insert(it.row() + s*m_udofs, it.col() + s*m_udofs) = it.value();

#pragma omp parallel for num_threads(m_numThreads)
    for (index_t col = 0; col < m_pshift; ++col)
        for (typename gsSparseMatrix<T>::InnerIterator it(m_blockB[col / m_udofs], col % m_udofs); it; ++it)
            stokesMatrix.insert(it.row() + m_pshift, col) = it.value();

#pragma omp parallel for num_threads(m_numThreads)
    for (index_t col = 0; col < m_pdofs; ++col)
        for (index_t s = 0; s < m_tarDim; ++s)
            for (typename gsSparseMatrix<T>::InnerIterator it(m_blockMinusBT[s], col); it; ++it)
                stokesMatrix.insert(it.row() + s*m_udofs, it.col() + m_pshift) = it.value();

//     if (getAssemblerOptions().intStrategy == iFace::dg)
//     {
// #pragma omp parallel for num_threads(m_numThreads)
//         for (index_t col = 0; col < m_pdofs; ++col)
//             for (typename gsSparseMatrix<T>::InnerIterator it(m_blockPdg, col); it; ++it)
//                 stokesMatrix.insert(it.row() + m_pshift, it.col() + m_pshift) = it.value();
//     }

    stokesMatrix.makeCompressed();

    stokesRhs = m_rhsA + m_rhsB + m_rhsF;

    // if (getAssemblerOptions().intStrategy == iFace::dg)
    //     stokesRhs += m_rhsPdg;
}


template<class T>
void gsINSBlockAssembler<T>::assemblePressurePoisson()
{
    // matrix and rhs cleaning
    m_rhsAp.setZero();
    m_blockAp.resize(m_pdofs, m_pdofs);
    gsSparseMatrix<T> blockApsym(m_pdofs, m_pdofs);

    // block assembly
    blockApsym.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    Base::template assembleBlock< gsINSBlockApVisitor<T> >(blockApsym, m_rhsAp, m_nonZerosPerColP);
    m_blockAp = blockApsym.template selfadjointView<gsEigen::Upper>();
    // if (getAssemblerOptions().intStrategy == iFace::dg) {
    //     gsSparseMatrix<T> blocksApdg(m_pdofs, m_pdofs);
    //     blocksApdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    //     Base::template assembleBlockDg< gsINSBlockApVisitorDG<T> >(blocksApdg, m_rhsAp, m_nonZerosPerColP);
    //     m_blockAp += blocksApdg;
    // }
    m_blockAp.makeCompressed();
}


template<class T>
void gsINSBlockAssembler<T>::assemblePressureConvection()
{
    // matrix and rhs cleaning
    m_blockNp.resize(m_pdofs, m_pdofs);
    gsMatrix<T> dummyRhsNp(m_dofs, 1);
    dummyRhsNp.setZero();

    // blocks assembly
    m_blockNp.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    Base::template assembleBlock< gsINSBlockNpVisitor<T> >(m_blockNp, dummyRhsNp, m_nonZerosPerColP);
    // if (getAssemblerOptions().intStrategy == iFace::dg) {
    //     gsSparseMatrix<T> m_blockNpdg(m_pdofs, m_pdofs);
    //     m_blockNpdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
    //     Base::template assembleBlockDg< gsINSBlockNpVisitorDG<T> >(m_blockNpdg, dummyRhsNp, m_nonZerosPerColP);
    //     m_blockNp += m_blockNpdg;
    // }
    m_blockNp.makeCompressed();
}


template<class T>
void gsINSBlockAssembler<T>::assembleApRobinBlock(std::vector<std::pair<int, boxSide> > bndPart)
{
    // matrix and rhs cleaning
    m_blockRobin.resize(m_pdofs, m_pdofs);
    gsMatrix<T> dummyRhs(m_dofs, 1);
    dummyRhs.setZero();

    m_blockRobin.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));

    gsINSBlockVisitorRobinPCD<T> visitor(m_dofMappers, m_viscosity);
    this->blockSpecificSettings(&visitor);

    gsMatrix<T> quNodes, physNodes;
    gsVector<T> quWeights; 
    unsigned evFlags(0);
    gsQuadRule<T> QuRule;
    gsMapData<T> mapData;
    typename gsBasis<T>::domainIter domIt;
    memory::shared_ptr< gsBasisRefs<T> > basesPtr;

    for (size_t sideIndex = 0; sideIndex < bndPart.size(); sideIndex++) // loop over inflow patch sides
    {
        int patch = bndPart[sideIndex].first;
        boxSide side = bndPart[sideIndex].second;

        basesPtr.reset(new gsBasisRefs<T>(getBases(), patch));

        visitor.initialize(*basesPtr, patch, QuRule, evFlags, side);

        mapData.flags = evFlags;

        domIt = (*basesPtr)[1].makeDomainIterator(side);

        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            mapData.points = quNodes;
            getPatches().patch(patch).computeMap(mapData);

            visitor.evaluate(*basesPtr, mapData, quNodes);
            visitor.assemble(*domIt, mapData, quWeights);
            visitor.localToGlobal(m_blockRobin, dummyRhs);
        }
    }
}


template<class T>
void gsINSBlockAssembler<T>::preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType)
{
    findPressureBoundaryIDs(bndIn, bndOut, bndWall);

    // Robin condition at inflow boundary
    if (bcType == 1 || bcType == 2 || bcType == 3 || bcType == 5)
        assembleApRobinBlock(bndIn);

    m_bPCDbndPrepared = true;
}


template<class T>
gsSparseMatrix<T> gsINSBlockAssembler<T>::getPressurePoissonMatrix(bool assemb, bool lumping)
{
    if (assemb)
    {
        if (!m_blockAp.nonZeros())
            assemblePressurePoisson();

        return m_blockAp;
    }
    else
    {
        gsSparseMatrix<T> velMinv(m_udofs, m_udofs);
        diagInvMatrix_into(m_blockM, velMinv, 1, lumping);

        gsSparseMatrix<T> Ap = m_blockB[0] * velMinv * m_blockB[0].transpose();
        for (size_t i = 1; i < m_blockB.size(); i++)
            Ap += m_blockB[i] * (velMinv * m_blockB[i].transpose());

        return Ap;
    }
}


template<class T>
void gsINSBlockAssembler<T>::applyPCDboundaryConditions(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType)
{
    switch (bcType)
    {
    case 0:
    default:
    {
        applyPCDboundaryDirichletPart(Ap, m_presInIDs);
        applyPCDboundaryDirichletPart(Fp, m_presInIDs);
        break;
    }
    case 1:
    {
        applyPCDboundaryRobinPart(Fp);
        break;
    }
    case 2:
    {
        applyPCDboundaryRobinPart(Fp);
        applyPCDboundaryDirichletPart(Fp, m_presOutIDs);
        break;
    }
    case 3:
    {
        applyPCDboundaryRobinPart(Fp);
        applyPCDboundaryDirichletPart(Ap, m_presOutIDs);
        applyPCDboundaryDirichletPart(Fp, m_presOutIDs);
        break;
    }
    case 4:
        break;
    }
}


template<class T>
void gsINSBlockAssembler<T>::findPressureBoundaryPartIDs(std::vector<std::pair<int, boxSide> > bndPart, std::vector<index_t>& idVector)
{
    const gsMultiBasis<T>& pBasis = getBases()[1];

    gsMatrix<index_t> bnd;
    idVector.clear();

    for (std::vector<std::pair<int, boxSide> >::iterator it = bndPart.begin(); it != bndPart.end(); it++)
    {
        bnd = pBasis.piece(it->first).boundary(it->second);

        for (int i = 0; i < bnd.rows(); i++)
            idVector.push_back(this->getMappers()[1].index(bnd(i, 0), it->first));
    }
}


template<class T>
void gsINSBlockAssembler<T>::applyPCDboundaryDirichletPart(gsSparseMatrix<T>& block, const std::vector<index_t>& bndPartIDs, bool scaleDiag, bool exclInflowPart)
{
    GISMO_ASSERT(m_bPCDbndPrepared, "Pressure boundary DOF IDs not prepared!");

    gsSparseMatrix<T> mat(block.rows(), block.cols());
    mat.setIdentity();

    for (size_t i = 0; i < bndPartIDs.size(); i++)
        mat.coeffRef(bndPartIDs[i], bndPartIDs[i]) = 0;

    mat.prune(0,0);

    block = mat * block * mat;

    gsMatrix<T> diag = block.diagonal();
    
    if (exclInflowPart)
    {
        for (size_t i = 0; i < m_presInIDs.size(); i++)
            diag(m_presInIDs[i]) = 0;
    }

    real_t alpha = 1.0;

    if(scaleDiag)
        alpha = diag.sum() / diag.count(); // average of the diagonal entries not corresponding to bndPart (+ inflow)        

    block.reserve(gsVector<int>::Constant(block.cols(), 1));

    for (size_t i = 0; i < bndPartIDs.size(); i++)
        block.coeffRef(bndPartIDs[i], bndPartIDs[i]) = alpha;

    block.makeCompressed();
}


} // namespace gismo
