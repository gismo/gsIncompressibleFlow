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
void gsINSAssemblerBase<T>::initMembers()
{ 
    m_dofs = 0;
    m_tarDim = getPatches().dim();
    m_viscosity = m_params.getPde().viscosity();
    m_isInitialized = false;
    m_isSystemReady = false;

    m_dofMappers.resize(2);
    m_ddof.resize(2);

    getBases().front().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.front(), 0);
    getBases().back().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.back(), 1);
}

template<class T>
void gsINSAssemblerBase<T>::computeDirichletDofs(const index_t unk, const index_t basisID, gsMatrix<T>& ddofVector)
{
    GISMO_ASSERT(ddofVector.rows() == m_dofMappers[basisID].boundarySize(), "Dirichlet DOF vector has wrong size.");

    const gsDofMapper & mapper = m_dofMappers[basisID];
    const gsMultiBasis<T> & mbasis = getBases().at(basisID);

    switch (getAssemblerOptions().dirValues)
    {
    case dirichlet::homogeneous:
        ddofVector.setZero();
        break;
    case dirichlet::interpolation:
        computeDirichletDofsIntpl(unk, mapper, mbasis, ddofVector);
        break;
    case dirichlet::l2Projection:
        computeDirichletDofsL2Proj(unk, mapper, mbasis, ddofVector);
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values.");
    }

    // Corner values
    for (typename gsBoundaryConditions<T>::const_citerator
        it = getBCs().cornerBegin();
        it != getBCs().cornerEnd(); ++it)
    {
        if (it->unknown == unk)
        {
            const index_t i = mbasis[it->patch].functionAtCorner(it->corner);
            const index_t ii = mapper.bindex(i, it->patch);
            ddofVector.row(ii).setConstant(it->value);
        }
        else
            continue;
    }
}


template<class T>
void gsINSAssemblerBase<T>::computeDirichletDofsIntpl(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
{
    for (typename gsBoundaryConditions<T>::const_iterator
        it = getBCs().dirichletBegin();
        it != getBCs().dirichletEnd(); ++it)
    {
        if (it->unknown() != unk)
            continue;

        const index_t patchID = it->patch();
        const gsBasis<T> & basis = mbasis.piece(patchID);

        // Get dofs on this boundary
        gsMatrix<index_t> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if (it->isHomogeneous())
        {
            for (index_t i = 0; i != boundary.size(); ++i)
            {
                const index_t ii = mapper.bindex((boundary)(i), it->patch());
                ddofVector.row(ii).setZero();
            }
            continue;
        }

        // Get the side information
        short_t dir = it->side().direction();
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve(getPatches().parDim());

        for (short_t i = 0; i < getPatches().parDim(); ++i)
        {
            if (i == dir)
            {
                gsVector<T> b(1);
                b[0] = (basis.component(i).support()) (0, param);
                rr.push_back(b);
            }
            else
            {
                rr.push_back(basis.component(i).anchors().transpose());
            }
        }

        GISMO_ASSERT(it->function()->targetDim() == ddofVector.cols(),
            "Given Dirichlet boundary function does not match problem dimension."
            << it->function()->targetDim() << " != " << ddofVector.cols() << "\n");

        // Compute dirichlet values
        gsMatrix<T> pts = getPatches().patch(it->patch()).eval(gsPointGrid<T>(rr));
        gsMatrix<T> fpts;
        if (!it->parametric())
        {
            fpts = it->function()->eval(pts);
        }
//            else if (it->parametric() && m_params.settings().isIgaDirichletGeometry())
//            {
//                gsMatrix<T> parPts;
//                getIgaBCGeom().invertPoints(pts, parPts);
//                fpts = it->function()->eval(parPts);
//            }
        else
            fpts = it->function()->eval(gsPointGrid<T>(rr));

        // Interpolate dirichlet boundary 
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals = geo->coefs();

        // Save corresponding boundary dofs
        for (index_t i = 0; i != boundary.size(); ++i)
        {
            const index_t ii = mapper.bindex(boundary.at(i), it->patch());
            ddofVector.row(ii) = dVals.row(i);
        }
    }
}


template<class T>
void gsINSAssemblerBase<T>::computeDirichletDofsL2Proj(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
{
    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero(ddofVector.rows(), ddofVector.cols());

    // Temporaries
    gsVector<T> quWeights;
    gsMatrix<T> rhsVals;
    gsMatrix<index_t> globIdxAct;
    gsMatrix<T> basisVals;
    gsMapData<T> mapData(NEED_MEASURE);

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for (typename gsBoundaryConditions<T>::const_iterator
        it = getBCs().dirichletBegin();
        it != getBCs().dirichletEnd(); ++it)
    {
        if (it->unknown() != unk)
            continue;

        if (it->isHomogeneous())
            continue;

        GISMO_ASSERT(it->function()->targetDim() == ddofVector.cols(),
            "Given Dirichlet boundary function does not match problem dimension."
            << it->function()->targetDim() << " != " << ddofVector.cols() << "\n");

        const index_t patchID = it->patch();
        const gsBasis<T> & basis = mbasis.piece(patchID);
        const gsGeometry<T> & patch = getPatches()[patchID];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on it->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, it->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(it->side());

        for (; bdryIter->good(); bdryIter->next())
        {
            bdQuRule.mapTo(bdryIter->lowerCorner(), bdryIter->upperCorner(),
                mapData.points, quWeights);

            patch.computeMap(mapData);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = it->function()->eval(getPatches()[patchID].eval(mapData.points));

            basis.eval_into(mapData.points, basisVals);

            // active basis (first line) functions/DOFs:
            basis.active_into(mapData.points.col(0), globIdxAct);
            mapper.localToGlobal(globIdxAct, patchID, globIdxAct);

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for (index_t i = 0; i < globIdxAct.rows(); i++)
                if (mapper.is_boundary_index(globIdxAct(i, 0)))
                    eltBdryFcts.push_back(i);

            // Do the actual assembly:
            for (index_t k = 0; k < mapData.points.cols(); k++)
            {
                const T weight_k = quWeights[k] * mapData.measure(k);

                // Only run through the active boundary functions on the element:
                for (size_t i0 = 0; i0 < eltBdryFcts.size(); i0++)
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const index_t i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const index_t ii = mapper.global_to_bindex(globIdxAct(i));

                    for (size_t j0 = 0; j0 < eltBdryFcts.size(); j0++)
                    {
                        const index_t j = eltBdryFcts[j0];
                        const index_t jj = mapper.global_to_bindex(globIdxAct(j));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i, k) * basisVals(j, k));
                    } // for j

                    globProjRhs.row(ii) += weight_k * basisVals(i, k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat(mapper.boundarySize(), mapper.boundarySize());
    globProjMat.setFrom(projMatEntries);
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    ddofVector = solver.compute(globProjMat).solve(globProjRhs);
}

template<class T>
void gsINSAssemblerBase<T>::assembleBlock(gsINSVisitor<T>& visitor, index_t testBasisID, gsSparseMatrix<T, RowMajor>& block, gsMatrix<T>& blockRhs)
{
    for(size_t p = 0; p < getPatches().nPatches(); p++)
    {
        index_t nBases = m_params.getBases()[testBasisID].piece(p).size();

        visitor.initOnPatch(p);

        for(index_t i = 0; i < nBases; i++)
        {
            visitor.evaluate(i);
            visitor.assemble();
            visitor.localToGlobal(m_ddof, block, blockRhs);
        }
    }

    block.makeCompressed();
}

template<class T>
void gsINSAssemblerBase<T>::assembleRhs(gsINSVisitor<T>& visitor, index_t testBasisID, gsMatrix<T>& rhs)
{
    for(size_t p = 0; p < getPatches().nPatches(); p++)
    {
        index_t nBases = m_params.getBases()[testBasisID].piece(p).size();

        visitor.initOnPatch(p);

        for(index_t i = 0; i < nBases; i++)
        {
            visitor.evaluate(i);
            visitor.assemble();
            visitor.localToGlobal(rhs);
        }
    }
}


template<class T>
void gsINSAssemblerBase<T>::initialize()
{
    assembleLinearPart();
    m_isInitialized = true;

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
}


template<class T>
void gsINSAssemblerBase<T>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    GISMO_ASSERT(m_isInitialized, "Assembler must be initialized first, call initialize()");

    updateCurrentSolField(solVector, updateSol);

    updateAssembly();

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


template<class T>
void gsINSAssemblerBase<T>::updateAssembly()
{
    assembleNonlinearPart();
}


template<class T>
void gsINSAssemblerBase<T>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
{
    m_dofMappers[unk] = gsDofMapper(getBases().at(unk), getBCs(), unk);

    if (getAssemblerOptions().intStrategy == iFace::conforming)
        for (gsBoxTopology::const_iiterator it = getPatches().iBegin(); it != getPatches().iEnd(); ++it)
        {
            getBases().at(unk).matchInterface(*it, m_dofMappers[unk]);
        }

    for (size_t i = 0; i < boundaryDofs.size(); i++)
    {
        m_dofMappers[unk].markBoundary(i, boundaryDofs[i]);
    }

    m_dofMappers[unk].finalize();

    initMembers();
}


// =============================================================================


template<class T>
void gsINSAssembler<T>::initMembers()
{
    m_udofs = m_dofMappers.front().freeSize();
    m_pdofs = m_dofMappers.back().freeSize();
    m_pshift = m_tarDim * m_udofs;
    m_dofs = m_pshift + m_pdofs;

    m_nnzPerRowU = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowU *= 2 * this->getBases().front().maxDegree(i) + 1;

    m_nnzPerRowP = 1;
    for (short_t i = 0; i < m_tarDim; i++)
        m_nnzPerRowP *= 2 * this->getBases().back().maxDegree(i) + 1;

    m_ddof[0].setZero(m_dofMappers.front().boundarySize(), m_tarDim);
    m_ddof[1].setZero(m_dofMappers.back().boundarySize(), 1);

    if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
    {
        this->computeDirichletDofs(0, 0, m_ddof[0]);
        this->computeDirichletDofs(1, 1, m_ddof[1]);
    }

    m_solution.setZero(m_dofs, 1);

    m_currentVelField = constructSolution(m_solution, 0);
    m_oldTimeVelField = m_currentVelField;

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

    m_blockUUlin.resize(m_pshift, m_pshift);
    m_blockUUnonlin.resize(m_pshift, m_pshift);
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

template<class T>
void gsINSAssembler<T>::assembleLinearPart()
{
    // matrix and rhs cleaning
    m_blockUUlin.resize(m_pshift, m_pshift);
    m_blockUP.resize(m_pshift, m_pdofs);
    m_rhsUlin.setZero(m_pshift, 1);
    m_rhsBtB.setZero(m_dofs, 1);
    m_rhsFG.setZero(m_dofs, 1);

    // memory allocation
    m_blockUUlin.reserve(gsVector<index_t>::Constant(m_blockUUlin.rows(), m_nnzPerRowU));
    m_blockUP.reserve(gsVector<index_t>::Constant(m_blockUP.rows(), m_nnzPerRowU));

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
    m_blockUUnonlin.resize(m_pshift, m_pshift);
    m_rhsUnonlin.setZero(m_pshift, 1);

    // memory allocation
    m_blockUUnonlin.reserve(gsVector<index_t>::Constant(m_blockUUnonlin.rows(), m_nnzPerRowU));

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
gsField<T> gsINSAssembler<T>::constructSolution(const gsMatrix<T>& solVector, index_t unk) const
{
    GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

    gsMultiPatch<T> * result = new gsMultiPatch<T>;

    const gsDofMapper & mapper = m_dofMappers[unk];

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
void gsINSAssembler<T>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
    {
        m_solution = solVector;

        if(this->isUnsteady())
            m_oldTimeVelField = constructSolution(solVector, 0);
    }

    m_currentVelField = constructSolution(solVector, 0);
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);
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

// =============================================================================


template<class T>
void gsINSAssemblerUnsteady<T>::initMembers()
{
    Base::initMembers();

    m_visitorTimeDiscr = gsINSVisitorUUtimeDiscr<T>(m_params);
    m_visitorTimeDiscr.initialize();

    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_rhsTimeDiscr.setZero(m_pshift, 1);
}

template<class T>
void gsINSAssemblerUnsteady<T>::assembleLinearPart()
{
    Base::assembleLinearPart();

    // matrix and rhs cleaning
    m_blockTimeDiscr.resize(m_pshift, m_pshift);

    // memory allocation
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.rows(), m_nnzPerRowU));

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
    GISMO_ASSERT(m_isInitialized, "Assembler must be initialized first, call initialize()");

    this->updateCurrentSolField(solVector, updateSol);

    this->updateAssembly();

    if(updateSol)
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

    if (m_params.options().getSwitch("fillGlobalSyst"))
        fillSystem();
}


} // namespace gismo