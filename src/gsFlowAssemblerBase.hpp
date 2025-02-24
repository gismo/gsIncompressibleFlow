/** @file gsFlowAssemblerBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova)
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>

namespace gismo
{

template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::initMembers()
{ 
    m_dofs = 0;
    m_tarDim = getPatches().dim();
    m_isInitialized = false;
    m_isBaseReady = false;
    m_isSystemReady = false;
}


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::computeDirichletDofs(const index_t unk, const index_t basisID, gsMatrix<T>& ddofVector)
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


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::computeDirichletDofsIntpl(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
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
//            else if (it->parametric() && m_paramsPtr->settings().isIgaDirichletGeometry())
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


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::computeDirichletDofsL2Proj(const index_t unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
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


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::assembleBlock(gsFlowVisitor<T, MatOrder>& visitor, index_t testBasisID, gsSparseMatrix<T, MatOrder>& block, gsMatrix<T>& blockRhs, bool compressMat)
{
    for(size_t p = 0; p < getPatches().nPatches(); p++)
    {
        visitor.initOnPatch(p);

        if (m_paramsPtr->options().getString("assemb.loop") == "RbR")
        {
            index_t nBases = m_paramsPtr->getBases()[testBasisID].piece(p).size();

            for(index_t i = 0; i < nBases; i++)
            {
                visitor.evaluate(i);
                visitor.assemble();
                visitor.localToGlobal(m_ddof, block, blockRhs);
            }
        }
        else
        {
            if (m_paramsPtr->options().getString("assemb.loop") != "EbE")
                gsWarn << "Unknown matrix formation method, using EbE (element by element)!\n";

            typename gsBasis<T>::domainIter domIt = m_paramsPtr->getBases().front().piece(p).makeDomainIterator(boundary::none);

            while (domIt->good())
            {
                visitor.evaluate(domIt.get());
                visitor.assemble();
                visitor.localToGlobal(m_ddof, block, blockRhs);

                domIt->next();
            }
        }    
    }

    if (compressMat)
        block.makeCompressed();
}


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::assembleRhs(gsFlowVisitor<T, MatOrder>& visitor, index_t testBasisID, gsMatrix<T>& rhs)
{
    for(size_t p = 0; p < getPatches().nPatches(); p++)
    {
        visitor.initOnPatch(p);

        if (m_paramsPtr->options().getString("assemb.loop") == "RbR")
        {
            index_t nBases = m_paramsPtr->getBases()[testBasisID].piece(p).size();

            for(index_t i = 0; i < nBases; i++)
            {
                visitor.evaluate(i);
                visitor.assemble();
                visitor.localToGlobal(rhs);
            }
        }
        else
        {
            typename gsBasis<T>::domainIter domIt = m_paramsPtr->getBases().front().piece(p).makeDomainIterator(boundary::none);

            while (domIt->good())
            {
                visitor.evaluate(domIt.get());
                visitor.assemble();
                visitor.localToGlobal(rhs);

                domIt->next();
            }
        }
    }
}


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::updateAssembly()
{
    assembleNonlinearPart();
}


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::initialize()
{
    assembleLinearPart();
    m_isInitialized = true;
}


template<class T, int MatOrder>
void gsFlowAssemblerBase<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    GISMO_ASSERT(m_isInitialized, "Assembler must be initialized first, call initialize()");

    updateCurrentSolField(solVector, updateSol);
    updateAssembly();
}


} // namespace gismo