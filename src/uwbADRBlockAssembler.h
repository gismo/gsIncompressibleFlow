/** @file uwbADRBlockAssembler.h

Author(s): E. Turnerova
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include <gsCore/gsGeometryEvaluator.h>
#include "uwbADRSolverParams.h"
#include "uwbADRBlockVisitors.h"
#include "uwbADRSUPGBlockVisitors.h"
#include "uwbADRCrosswindVisitors.h"
#include "uwbADRArtificialDiffusionVisitors.h"

namespace gismo
{

template<class T>
class uwbADRBlockAssembler
{

public:

    uwbADRBlockAssembler(const uwbADRSolverParams<T>& params, int numVar) : m_params(params), m_numVar(numVar)
    {
        initMembers();
    }

    virtual ~uwbADRBlockAssembler()
    {}

protected:

    virtual void initMembers()
    {
        m_tarDim = getPatches().dim();
        setNumThreads(m_params.getNumThreads());

        getBases().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMapper, 0);

        m_varDofs = getMapper().freeSize();
        m_dofs = m_numVar * m_varDofs;

        m_nonZerosPerCol = 1;
        for (int i = 0; i < m_tarDim; i++)
            m_nonZerosPerCol *= 2 * getBases().maxDegree(i) + 1;

        m_ddof.resize(m_numVar);

        initElementList();

        if (getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            for (size_t i = 0; i < m_ddof.size(); i++)
            {
                m_ddof[i].setZero(getMapper().boundarySize(), 1);
                computeDirichletDofs(i, m_ddof[i]);
            }
        }
            
        m_solution.setZero(m_varDofs, m_numVar);
        m_oldTimeSol.setZero(m_varDofs, m_numVar);

        m_mMass.resize(m_varDofs, m_varDofs);
        m_mADR.resize(m_varDofs, m_dofs);
        m_mADRpattern.resize(m_varDofs, m_dofs);

        m_rhsMass.setZero(m_varDofs, m_numVar);
        m_rhsADR.setZero(m_varDofs, m_numVar);

        m_mMass_SUPG.resize(m_varDofs, m_varDofs);
        m_mADR_SUPG.resize(m_varDofs, m_dofs);
        m_mADRPattern_SUPG.resize(m_varDofs, m_dofs);

        m_rhsMass_SUPG.setZero(m_varDofs, m_numVar);
        m_rhsADR_SUPG.setZero(m_varDofs, m_numVar);

        m_mCrosswindPattern.resize(m_varDofs, m_dofs);
        m_mCrosswind.resize(m_varDofs, m_dofs);

        m_rhsCrosswind.setZero(m_varDofs, m_numVar);

        m_mIsoArtificialDiffPattern.resize(m_varDofs, m_dofs);
        m_mIsoArtificialDiff.resize(m_varDofs, m_dofs);
        m_mArtificialDiffPattern.resize(m_varDofs, m_dofs);
        m_mArtificialDiff.resize(m_varDofs, m_dofs);

        m_rhsIsoArtificialDiff.setZero(m_varDofs, m_numVar);
        m_rhsArtificialDiff.setZero(m_varDofs, m_numVar);

        m_bLinearADRCoeffsSet = false;
        m_bNonlinADcoeffsSet = false;
        m_bNonlinADRcoeffsSet = false;
    }

public:
    // assumes that solVector is either solution vector of one variable
    // or each column of solVector contains solution vector of one variable
    void constructSolution(const gsMatrix<T>& solVector, gsMultiPatch<T>& result) const
    {
        GISMO_ASSERT(m_varDofs == solVector.rows(), "Something went wrong, is solution vector valid?");

        result.clear(); // result is cleared first

        gsMatrix<T> ddofMat(m_ddof[0].rows(), m_ddof.size());

        if (ddofMat.rows()) // if number of eliminated dofs is not zero
        {
            for (int j = 0; j < ddofMat.cols(); j++)
                ddofMat.col(j) = m_ddof[j];
        }

        const index_t dim = solVector.cols();
        gsMatrix<T> coeffs;

        for (unsigned int p = 0; p < getPatches().nPatches(); ++p)
        {
            // Reconstruct solution coefficients on patch p
            const int sz = getBases().piece(p).size();
            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (getMapper().is_free(i, p)) // DoF value is in the solVector
                    coeffs.row(i) = solVector.row(getMapper().index(i, p));
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs.row(i) = ddofMat.row(getMapper().bindex(i, p));
                }
            }

            result.addPatch(getBases().piece(p).makeGeometry(coeffs));
        }
    }

    gsField<T> constructSolution(const gsMatrix<T>& solVector) const
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructSolution(solVector, *result);
        return gsField<T>(getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
    }

    void setSolution(const gsMatrix<T> & solVector)
    { m_solution = solVector; }

    void setLinearADRCoefficients(T diffusionCoeff, gsVector<T> advectionCoeff, T reactionCoeff)
    {
        m_diffusionCoeff = diffusionCoeff;
        m_advectionCoeff = advectionCoeff;
        m_reactionCoeff = reactionCoeff;
        m_bLinearADRCoeffsSet = true;
    }

    void setAdvectionDiffusionCoeffFields(gsField<T>& advectionField, gsField<T>& diffusionField)
    {
        m_advectionField = advectionField;
        m_diffusionField = diffusionField;
        m_bNonlinADcoeffsSet = true;
    }
    void setADRcoeffFields(gsField<T>& advectionField, gsField<T>& diffusionField, gsField<T>& reactionField)
    {
        m_advectionField = advectionField;
        m_diffusionField = diffusionField;
        m_reactionField = reactionField;
        m_bNonlinADRcoeffsSet = true;
    }

    void setCoeffGeometry(const gsMultiPatch<T>& patchesCoeffs, const gsMultiBasis<T>& basesCoeffs)
    { m_patchesCoeffs = patchesCoeffs; m_basesCoeffs = basesCoeffs; }

protected:
    void initElementList()
    {
        m_elementList.setZero(getBases().totalElements(), 2);

        int nelementglobal = 0;
        for (size_t npatch = 0; npatch < getPatches().nPatches(); ++npatch) {
            for (size_t nelement = 0; nelement < getBases().basis(npatch).numElements(); ++nelement) {
                m_elementList(nelementglobal, 0) = npatch;
                m_elementList(nelementglobal, 1) = nelement;
                ++nelementglobal;
            }
        }
    }

    void computeDirichletDofs(const int unk, gsMatrix<T>& ddofVector)
    {
        GISMO_ASSERT(ddofVector.rows() == m_dofMapper.boundarySize(), "Dirichlet DOF vector has wrong size.");

        const gsDofMapper & mapper = m_dofMapper;
        const gsMultiBasis<T> & mbasis = getBases();

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
                const int i = mbasis[it->patch].functionAtCorner(it->corner);
                const int ii = mapper.bindex(i, it->patch);
                ddofVector.row(ii).setConstant(it->value);
            }
            else
                continue;
        }
    }

    // Compute Dirichlet dofs using interpolation
    void computeDirichletDofsIntpl(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
    {
        for (typename gsBoundaryConditions<T>::const_iterator
            it = getBCs().dirichletBegin();
            it != getBCs().dirichletEnd(); ++it)
        {
            if (it->unknown() != unk)
                continue;

            const int patchIdx = it->patch();
            const gsBasis<T> & basis = mbasis.piece(patchIdx);

            // Get dofs on this boundary
            gsMatrix<index_t> boundary = basis.boundary(it->side());

            // If the condition is homogeneous then fill with zeros
            if (it->isHomogeneous())
            {
                for (index_t i = 0; i != boundary.size(); ++i)
                {
                    const int ii = mapper.bindex((boundary)(i), it->patch());
                    ddofVector.row(ii).setZero();
                }
                continue;
            }

            // Get the side information
            int dir = it->side().direction();
            index_t param = (it->side().parameter() ? 1 : 0);

            // Compute grid of points on the face ("face anchors")
            std::vector< gsVector<T> > rr;
            rr.reserve(getPatches().parDim());

            for (int i = 0; i < getPatches().parDim(); ++i)
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
            if (it->parametric())
            {
                gsMatrix<T> parPts;
                getIgaBCGeom().invertPoints(pts, parPts);
                fpts = it->function()->eval(parPts);
            }
            else
                fpts = it->function()->eval(pts);

            // Interpolate dirichlet boundary
            typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
            typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
            const gsMatrix<T> & dVals = geo->coefs();

            // Save corresponding boundary dofs
            for (index_t i = 0; i != boundary.size(); ++i)
            {
                const int ii = mapper.bindex(boundary.at(i), it->patch());
                ddofVector.row(ii) = dVals.row(i);
            }
        }
    }

    void computeDirichletDofsL2Proj(const int unk, const gsDofMapper & mapper, const gsMultiBasis<T> & mbasis, gsMatrix<T>& ddofVector)
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
        gsMapData<T> md(NEED_MEASURE);

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

            const int patchIdx = it->patch();
            const gsBasis<T> & basis = mbasis.piece(patchIdx);
            const gsGeometry<T> & patch = getPatches()[patchIdx];

            // Set up quadrature to degree+1 Gauss points per direction,
            // all lying on it->side() except from the direction which
            // is NOT along the element

            gsGaussRule<T> bdQuRule(basis, 1.0, 1, it->side().direction());

            // Create the iterator along the given part boundary.
            typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(it->side());

            for (; bdryIter->good(); bdryIter->next())
            {
                bdQuRule.mapTo(bdryIter->lowerCorner(), bdryIter->upperCorner(),
                    md.points, quWeights);

                //geoEval->evaluateAt( md.points );
                patch.computeMap(md);

                // the values of the boundary condition are stored
                // to rhsVals. Here, "rhs" refers to the right-hand-side
                // of the L2-projection, not of the PDE.
                rhsVals = it->function()->eval(getPatches()[patchIdx].eval(md.points));

                basis.eval_into(md.points, basisVals);

                // active basis (first line) functions/DOFs:
                basis.active_into(md.points.col(0), globIdxAct);
                mapper.localToGlobal(globIdxAct, patchIdx, globIdxAct);

                // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
                // something like a "element-wise index"
                std::vector<index_t> eltBdryFcts;
                eltBdryFcts.reserve(mapper.boundarySize());
                for (index_t i = 0; i < globIdxAct.rows(); i++)
                    if (mapper.is_boundary_index(globIdxAct(i, 0)))
                        eltBdryFcts.push_back(i);

                // Do the actual assembly:
                for (index_t k = 0; k < md.points.cols(); k++)
                {
                    const T weight_k = quWeights[k] * md.measure(k);

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

    template<class ElementBlockVisitor>
    void assembleBlock(gsSparseMatrix<T> & matrixBlock,
        gsMatrix<T> & rhs,
        const int nonZerosPerCol)
    {
        std::vector< gsSparseMatrix<T> *> threadBlockPtr;
        threadBlockPtr.resize(m_numThreads);
        threadBlockPtr[0] = &matrixBlock;
        for (int t = 1; t < m_numThreads; ++t) {
            threadBlockPtr[t] = new gsSparseMatrix<T>(matrixBlock.rows(), matrixBlock.cols());
            threadBlockPtr[t]->reserve(gsVector<int>::Constant(matrixBlock.cols(), nonZerosPerCol));
        }

        std::vector< gsMatrix<T> * > threadRhsPtr;
        threadRhsPtr.resize(m_numThreads);
        threadRhsPtr[0] = &rhs;
        for (int t = 1; t < m_numThreads; ++t) {
            threadRhsPtr[t] = new gsMatrix<T>(rhs.rows(), rhs.cols());
            threadRhsPtr[t]->setZero();
        }

        #pragma omp parallel num_threads(m_numThreads)
        {
            #ifdef _OPENMP
            const int threadNumber(omp_get_thread_num());
            #else
            const int threadNumber = 0;
            #endif

            // deep copies
            const gsMultiPatch<T> patches(getPatches());
            //gsMultiBasis<T> bases(getBases());
            //std::vector< gsMultiBasis<T> > bases(getBases());
            std::vector< gsMultiBasis<T> > bases;
            bases.push_back(getBases());
            gsDofMapper dofMapper(m_dofMapper);
            const std::vector<gsMatrix<T> > ddof(m_ddof);
            const int totalElementsNumber(m_elementList.rows());
            const gsMatrix<index_t> elementList(m_elementList);

            gsQuadRule<T> QuRule;
            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            unsigned evFlags(0);
            typename gsGeometryEvaluator<T>::uPtr geoEval;
            memory::shared_ptr< gsBasisRefs<T> > basesPtr;
            typename gsBasis<T>::domainIter domIt;

            gsSparseMatrix<T> * blockPtr = threadBlockPtr[threadNumber];
            gsMatrix<T> * rhsPtr = threadRhsPtr[threadNumber];

            ElementBlockVisitor blockVisitor(dofMapper, getADREvaluator());

            blockSpecificSettings(&blockVisitor);

            index_t patchIndex = -1;
            index_t elementIndex = -1;

            // Parallel assembly procedure.
            #pragma omp for
            for (int nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal) {

                if (patchIndex != elementList(nElementGlobal, 0)) {
                    patchIndex = elementList(nElementGlobal, 0);

                    //basesPtr.reset(&bases.basis(patchIndex));
                    basesPtr.reset(new gsBasisRefs<T>(bases, patchIndex));
                    domIt = basesPtr->front().makeDomainIterator(boundary::none);
                    //domIt = basesPtr->makeDomainIterator(boundary::none);
                    elementIndex = 0;

                    blockVisitor.initialize(*basesPtr, patchIndex, QuRule, evFlags);

                    geoEval = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patchIndex]));
                }

                if (elementIndex != elementList(nElementGlobal, 1)) {
                    while ((elementIndex < elementList(nElementGlobal, 1)) && (domIt->good())) {
                        domIt->next();
                        ++elementIndex;
                    }
                }

                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                if (checkStabMethod())
                {
                    T elemDiam = computeElementDiameter(domIt, patchIndex);
                    //gsInfo << "elemDiam = " << elemDiam << "\n";
                    blockVisitor.setElementLength(elemDiam);
                }

                blockVisitor.evaluate(*basesPtr, *geoEval, quNodes);

                blockVisitor.assemble(*geoEval, quWeights);

                blockVisitor.localToGlobal(ddof, *blockPtr, *rhsPtr);

                domIt->next();
                ++elementIndex;
            }
        }

        // Summing up block matrices and rhs from all threads
        // into threadBlockPtr[0], threadRhsPtr[0] which are
        // pointers to the input variables & matrixBlock, & rhs
        for (int levelSkip = 1; levelSkip < m_numThreads; levelSkip *= 2) {
            #pragma omp parallel for num_threads(m_numThreads)
            for (int i = 0; i < m_numThreads; i += 2 * levelSkip) {
                if (i + levelSkip < m_numThreads) {
                    *threadBlockPtr[i] += *threadBlockPtr[i + levelSkip]; //matrix loses allocated space!!
                    *threadRhsPtr[i] += *threadRhsPtr[i + levelSkip];
                }
            }
        }

        for (int t = 1; t < m_numThreads; ++t) {
            delete threadBlockPtr[t];
            delete threadRhsPtr[t];
        }
    }

    T computeElementDiameter(typename gsBasis<T>::domainIter& domIt, int patchIndex)
    {
        gsMatrix<T> nodes;
        gsVector<T> nodeA = domIt->lowerCorner();
        gsVector<T> nodeC = domIt->upperCorner();
        nodes.setZero(getPatches().dim(), 4);
        nodes.middleCols(0, 1) = nodeA; //point A of the element in parametric domain
        nodes.middleCols(2, 1) = nodeC; //point C of the element in parametric domain
        nodes(0, 1) = nodeC(0); //x-coordinate of point B of the element in parametric domain
        nodes(1, 1) = nodeA(1); //y-coordinate of point B of the element in parametric domain
        nodes(0, 3) = nodeA(0); //x-coordinate of point D of the element in parametric domain
        nodes(1, 3) = nodeC(1); //y-coordinate of point D of the element in parametric domain
        gsMatrix<T> physNodes;
        getPatches().patch(patchIndex).eval_into(nodes, physNodes);
        T lengthE = (physNodes.col(2) - physNodes.col(0)).norm();
        T lengthF = (physNodes.col(3) - physNodes.col(1)).norm();

        T elemDiam = math::max(lengthE, lengthF);
        //gsInfo << "elemDiam = " << elemDiam << "\n";
        return elemDiam;
    }

    virtual void blockSpecificSettings(uwbADRBlockVisitor<T>* blockVisitor)
    {
        gsField<T> solution = constructSolution(m_solution);
        blockVisitor->setCurrentSolution(solution);

        if (getADREvaluator() == "linConstCoeffs")
        {
            GISMO_ASSERT(m_bLinearADRCoeffsSet, "Coefficients of the linear advection-diffusion-reaction equation not set in assembler.");
            blockVisitor->setADRCoefficients(m_diffusionCoeff, m_advectionCoeff, m_reactionCoeff);
        }
        else if (getADREvaluator() == "nonlinCoeffsField" && m_bNonlinADcoeffsSet)
        {
            blockVisitor->setAdvectionDiffusionCoeffFields(m_advectionField, m_diffusionField);
            //blockVisitor->setCoeffGeometry(m_patchesCoeffs, m_basesCoeffs);
        }
        else if (getADREvaluator() == "nonlinCoeffsField" && m_bNonlinADRcoeffsSet)
        {
            blockVisitor->setADRcoeffFields(m_advectionField, m_diffusionField, m_reactionField);
            //blockVisitor->setCoeffGeometry(m_patchesCoeffs, m_basesCoeffs);
        }

        if (isSUPG())
        {
            uwbADRSUPGBlockVisitor<T>* pVisitorSUPG = dynamic_cast<uwbADRSUPGBlockVisitor<T>*>(blockVisitor);
            if (pVisitorSUPG != NULL)
                pVisitorSUPG->setSUPG(getTauStabType(), getTimeStep());
        }

        if (isCrosswind())
        {
            uwbADRCrosswindVisitor<T>* pVisitorCrosswind = dynamic_cast<uwbADRCrosswindVisitor<T>*>(blockVisitor);
            if (pVisitorCrosswind != NULL)
            {
                pVisitorCrosswind->setCrosswind(getCrosswindType(), getTauStabType(), isUnsteady(), getTimeStep());
                if (isUnsteady())
                {
                    gsField<T> oldSolField = constructSolution(m_oldTimeSol);
                    pVisitorCrosswind->setOldSolutionField(oldSolField);
                }
            }
        }

        if (isIsoArtificialDiff())
        {
            uwbADRIsoArtificialDiffusionVisitor<T>* pVisitorIsoAD = dynamic_cast<uwbADRIsoArtificialDiffusionVisitor<T>*>(blockVisitor);
            if (pVisitorIsoAD != NULL)
            {
                pVisitorIsoAD->setIsoArtificialDiffusion(isUnsteady(), getTimeStep());
                if (isUnsteady())
                {
                    gsField<T> oldSolField = constructSolution(m_oldTimeSol);
                    pVisitorIsoAD->setOldSolutionField(oldSolField);
                }
            }
        }

        if (isArtificialDiff())
        {
            uwbADRArtificialDiffusionVisitor<T>* pVisitorAD = dynamic_cast<uwbADRArtificialDiffusionVisitor<T>*>(blockVisitor);
            if (pVisitorAD != NULL)
                pVisitorAD->setArtificialDiffusion(getTauStabType(), getTimeStep());
        }
    }

public:
    void assembleMassMatrix()
    {
        // mass matrix
        // matrix and rhs cleaning
        m_mMass.resize(m_varDofs, m_varDofs);
        m_rhsMass.setZero();
        
        m_mMass.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
        if (isFCTlowOrder())
            assembleBlock< uwbMassMatrixFCTLowOrderVisitor<T> >(m_mMass, m_rhsMass, m_nonZerosPerCol);
        else
            assembleBlock< uwbMassMatrixVisitor<T> >(m_mMass, m_rhsMass, m_nonZerosPerCol);
        m_mMass.makeCompressed();

        if (isSUPG())
            assembleMassMatrixSUPG();
    }

    void assemblePatternBlocks()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        // advection-diffusion-reaction pattern matrix
        // matrix cleaning
        m_mADRpattern.resize(m_varDofs, m_dofs);

        m_mADRpattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        if (isFCTlowOrder())
            assembleBlock< uwbADRFCTLowOrderVisitor<T> >(m_mADRpattern, dummyRhs, m_nonZerosPerCol);
        else
            assembleBlock< uwbADRVisitor<T> >(m_mADRpattern, dummyRhs, m_nonZerosPerCol);
        m_mADRpattern.makeCompressed();
        m_mADRpattern = 0. * m_mADRpattern;

        if (isSUPG())
            assembleSUPGpatternPart();
        if (isCrosswind())
            assembleCrosswindPatternPart();
        if (isIsoArtificialDiff())
            assembleIsoArtificialDiffPatternPart();
        if (isArtificialDiff())
            assembleArtificialDiffPatternPart();
    }

    void assembleNonlinearBlocks(const gsMatrix<T> & solVector, bool updateSol = true)
    {
        if (updateSol)
            m_oldTimeSol = solVector;

        gsMatrix<T> current_solution = m_solution;
        m_solution = solVector;

        assembleADR();

        if (!updateSol)
            m_solution = current_solution;
    }

    void assembleADR()
    {
        // advection-diffusion-reaction matrix
        // matrix cleaning
        m_mADR.resize(m_varDofs, m_dofs);
        m_rhsADR.setZero();

        m_mADR.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        if (isFCTlowOrder())
            assembleBlock< uwbADRFCTLowOrderVisitor<T> >(m_mADR, m_rhsADR, m_nonZerosPerCol);
        else
            assembleBlock< uwbADRVisitor<T> >(m_mADR, m_rhsADR, m_nonZerosPerCol);
        m_mADR.makeCompressed();

        if (isSUPG())
            assembleNonlinearSUPGPart();
        if (isCrosswind())
            assembleCrosswindPart();
        if (isIsoArtificialDiff())
            assembleIsoArtificialDiffPart();
        if (isArtificialDiff())
            assembleArtificialDiffPart();
    }

    //======================== assemble SUPG pattern part =============================
    void assembleSUPGpatternPart()
    {
        gsMatrix<T> dummyRhs_SUPG(m_varDofs, m_numVar);
        dummyRhs_SUPG.setZero();

        // SUPG block Mass matrix pattern
        // matrix cleaning
        m_mMassPattern_SUPG.resize(m_varDofs, m_varDofs);
        m_mMassPattern_SUPG.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
        assembleBlock< uwbMassMatrixSUPGVisitor<T> >(m_mMassPattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
        m_mMassPattern_SUPG.makeCompressed();
        m_mMassPattern_SUPG = 0. * m_mMassPattern_SUPG;
            
        // SUPG advection-diffusion-reaction pattern
        // matrix cleaning
        m_mADRPattern_SUPG.resize(m_varDofs, m_dofs);
        m_mADRPattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRSUPGVisitor<T> >(m_mADRPattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
        m_mADRPattern_SUPG.makeCompressed();
        m_mADRPattern_SUPG = 0. * m_mADRPattern_SUPG;
    }

    //======================== assemble CROSSWIND pattern part =============================
    void assembleCrosswindPatternPart()
    {
        gsMatrix<T> dummyRhs_CW(m_varDofs, m_numVar);
        dummyRhs_CW.setZero();

        // CROSSWIND matrix pattern
        // matrix cleaning
        m_mCrosswindPattern.resize(m_varDofs, m_dofs);
        m_mCrosswindPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRCrosswindVisitor<T> >(m_mCrosswindPattern, dummyRhs_CW, m_nonZerosPerCol);
        m_mCrosswindPattern.makeCompressed();
        m_mCrosswindPattern = 0. * m_mCrosswindPattern;
    }

    //=============== assemble isotropic artificial diffusion pattern part ================
    void assembleIsoArtificialDiffPatternPart()
    {
        gsMatrix<T> dummyRhs_CW(m_varDofs, m_numVar);
        dummyRhs_CW.setZero();

        // isoAD matrix pattern
        // matrix cleaning
        m_mIsoArtificialDiffPattern.resize(m_varDofs, m_dofs);
        m_mIsoArtificialDiffPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRIsoArtificialDiffusionVisitor<T> >(m_mIsoArtificialDiffPattern, dummyRhs_CW, m_nonZerosPerCol);
        m_mIsoArtificialDiffPattern.makeCompressed();
        m_mIsoArtificialDiffPattern = 0. * m_mIsoArtificialDiffPattern;
    }

    //=============== assemble artificial diffusion pattern part ================
    void assembleArtificialDiffPatternPart()
    {
        gsMatrix<T> dummyRhs_CW(m_varDofs, m_numVar);
        dummyRhs_CW.setZero();

        // ADstab matrix pattern
        // matrix cleaning
        m_mArtificialDiffPattern.resize(m_varDofs, m_dofs);
        m_mArtificialDiffPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRArtificialDiffusionVisitor<T> >(m_mArtificialDiffPattern, dummyRhs_CW, m_nonZerosPerCol);
        m_mArtificialDiffPattern.makeCompressed();
        m_mArtificialDiffPattern = 0. * m_mArtificialDiffPattern;
    }

    //======================== assemble SUPG terms =============================
    void assembleNonlinearSUPGPart()
    {
        // SUPG advection-diffusion-reaction terms
        // matrix cleaning
        m_mADR_SUPG.resize(m_varDofs, m_dofs);
        m_rhsADR_SUPG.setZero();

        m_mADR_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRSUPGVisitor<T> >(m_mADR_SUPG, m_rhsADR_SUPG, m_nonZerosPerCol);
        m_mADR_SUPG.makeCompressed();
    }

    //======================== assemble CROSSWIND term =========================
    void assembleCrosswindPart()
    {
        // matrix cleaning
        m_mCrosswind.resize(m_varDofs, m_dofs);
        m_rhsCrosswind.setZero();

        m_mCrosswind.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRCrosswindVisitor<T> >(m_mCrosswind, m_rhsCrosswind, m_nonZerosPerCol);
        m_mCrosswind.makeCompressed();
    }

    //================= assemble isotropi artificial diffusion term ==============
    void assembleIsoArtificialDiffPart()
    {
        // matrix cleaning
        m_mIsoArtificialDiff.resize(m_varDofs, m_dofs);
        m_rhsIsoArtificialDiff.setZero();

        m_mIsoArtificialDiff.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRIsoArtificialDiffusionVisitor<T> >(m_mIsoArtificialDiff, m_rhsIsoArtificialDiff, m_nonZerosPerCol);
        m_mIsoArtificialDiff.makeCompressed();
    }

    //================= assemble isotropi artificial diffusion term ==============
    void assembleArtificialDiffPart()
    {
        // matrix cleaning
        m_mArtificialDiff.resize(m_varDofs, m_dofs);
        m_rhsArtificialDiff.setZero();

        m_mArtificialDiff.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        assembleBlock< uwbADRArtificialDiffusionVisitor<T> >(m_mArtificialDiff, m_rhsArtificialDiff, m_nonZerosPerCol);
        m_mArtificialDiff.makeCompressed();
    }

    void assembleMassMatrixSUPG()
    {
        // SUPG block M
        // matrix and rhs cleaning
        m_mMass_SUPG.resize(m_varDofs, m_varDofs);
        m_rhsMass_SUPG.setZero();

        m_mMass_SUPG.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
        assembleBlock< uwbMassMatrixSUPGVisitor<T> >(m_mMass_SUPG, m_rhsMass_SUPG, m_nonZerosPerCol);
        m_mMass_SUPG.makeCompressed();
    }


public:
    int getNumVar() const { return m_numVar; }
    int numVarDofs() const { return m_varDofs; }

    std::string getADREvaluator() {return m_params.settings().getADREvaluator(); }

    const gsSparseMatrix<T>& getMassMatrix() const { return m_mMass; }
    const gsSparseMatrix<T>& getADRMatrix() const { return m_mADR; }
    const gsSparseMatrix<T>& getADRMatrixPattern() const { return m_mADRpattern; }
    const gsMatrix<T>& getRhsMass() const { return m_rhsMass; }
    const gsMatrix<T>& getRhsADR() const { return m_rhsADR; }

    //--------- SUPG ------------------
    bool    isSUPG() const { return m_params.settings().get(constantsADR::SUPG); }
    int     getTauStabType() const { return m_params.settings().get(constantsADR::tauStabType); }
    T       getTimeStep() const { return m_params.settings().get(constantsADR::timeStep); }

    //--------- CROSSWIND --------------
    bool    isCrosswind() const { return m_params.settings().get(constantsADR::CROSSWIND); }
    int     getCrosswindType() const { return m_params.settings().get(constantsADR::crosswindType); }

    //--------- artificial diffusion --------------
    bool    isIsoArtificialDiff() const { return m_params.settings().get(constantsADR::isoArtificialDiffusion); }
    bool    isArtificialDiff() const { return m_params.settings().get(constantsADR::artificialDiffusion); }

    //--------- FCT ------------------
    bool    isFCTlowOrder() const { return m_params.settings().get(constantsADR::FCT_lowOrder); }

    //SUPG getters
    const gsSparseMatrix<T>& getMassMatrix_SUPG() const { return m_mMass_SUPG; }
    const gsSparseMatrix<T>& getMassMatrixPattern_SUPG() const { return m_mMassPattern_SUPG; }
    const gsSparseMatrix<T>& getADRMatrix_SUPG() const { return m_mADR_SUPG; }
    const gsSparseMatrix<T>& getADRMatrixPattern_SUPG() const { return m_mADRPattern_SUPG; }

    const gsMatrix<T>& getRhsMass_SUPG() const { return m_rhsMass_SUPG; }
    const gsMatrix<T>& getRhsADR_SUPG() const { return m_rhsADR_SUPG; }

    //CROSSWIND getters
    const gsSparseMatrix<T>& getCrosswindMatrix() const { return m_mCrosswind; }
    const gsSparseMatrix<T>& getCrosswindPatternMatrix() const { return m_mCrosswindPattern; }

    const gsMatrix<T>& getRhsCrosswind() const { return m_rhsCrosswind; }

    //artificial diffusion getters
    const gsSparseMatrix<T>& getIsoArtificialDiffMatrix() const { return m_mIsoArtificialDiff; }
    const gsSparseMatrix<T>& getIsoArtificialDiffPatternMatrix() const { return m_mIsoArtificialDiffPattern; }
    const gsSparseMatrix<T>& getArtificialDiffMatrix() const { return m_mArtificialDiff; }
    const gsSparseMatrix<T>& getArtificialDiffPatternMatrix() const { return m_mArtificialDiffPattern; }

    const gsMatrix<T>& getRhsIsoArtificialDiff() const { return m_rhsIsoArtificialDiff; }
    const gsMatrix<T>& getRhsArtificialDiff() const { return m_rhsArtificialDiff; }

    //SUPG members
    gsSparseMatrix<T> m_mMass_SUPG;
    gsSparseMatrix<T> m_mMassPattern_SUPG;
    gsSparseMatrix<T> m_mADR_SUPG;
    gsSparseMatrix<T> m_mADRPattern_SUPG;

    gsMatrix<T> m_rhsMass_SUPG;
    gsMatrix<T> m_rhsADR_SUPG;

    //CROSSWIND members
    gsSparseMatrix<T> m_mCrosswind;
    gsSparseMatrix<T> m_mCrosswindPattern;

    gsMatrix<T> m_rhsCrosswind;

    //artificial diffusion members
    gsSparseMatrix<T> m_mIsoArtificialDiff;
    gsSparseMatrix<T> m_mIsoArtificialDiffPattern;
    gsSparseMatrix<T> m_mArtificialDiff;
    gsSparseMatrix<T> m_mArtificialDiffPattern;

    gsMatrix<T> m_rhsIsoArtificialDiff;
    gsMatrix<T> m_rhsArtificialDiff;

    // member functions from uwbBlockAssemblerBase
public:
    int numDofs() const
    {
        GISMO_ASSERT(m_dofs > 0, "Something went wrong, number of DOFs is zero!");
        return m_dofs;
    }

    int getTarDim() const { return m_tarDim; }
    int getNumThreads() const { return m_numThreads; }

    void setNumThreads(const int numThreads)
    {
        #ifdef _OPENMP
        const int maxThreads = omp_get_max_threads();
        if ((numThreads > 0) && (numThreads <= maxThreads))
            m_numThreads = numThreads;
        else
        {
            gsWarn << "The maximum number of threads ( " << maxThreads << " ) will be used.\n";
            m_numThreads = maxThreads;
        }
        #else
        //gsWarn << "_OPENMP not defined.\n";
        m_numThreads = 1;
        #endif
    }

    bool isUnsteady() const { return m_params.settings().get(constantsADR::unsteady); }

    const gsDofMapper &    getMapper() const { return m_dofMapper; }
    const std::vector<gsMatrix<T> > &   getDirichletDofs() const { return m_ddof; }
    const gsMatrix<T> &                 getSolution() const { return m_solution; }
    const gsMultiPatch<T>&              getPatches() const { return m_params.getPde().patches(); }
    const gsMultiBasis<T>& getBases() const { return m_params.getBases(); }
    gsMultiBasis<T>& getBases() { return m_params.getBases(); }
    const gsBoundaryConditions<T>&        getBCs() const { return m_params.getBCs(); }
    //const gsFunction<T>&      getRhsFcn() const { return m_params.getPde().getRhs(); }
    gsAssemblerOptions  getAssemblerOptions() const { return m_params.getAssemblerOptions(); }
    gsGeometry<T>& getIgaBCGeom() const { return m_params.settings().getIgaDirichletGeometry(); }
    int     getParam(constantsADR::intConst name) const { return m_params.settings().get(name); }
    T       getParam(constantsADR::realConst name) const { return m_params.settings().get(name); }
    bool    getParam(constantsADR::boolConst name) const { return m_params.settings().get(name); }

    bool checkStabMethod()
    {
        return(isSUPG() || isCrosswind() || isIsoArtificialDiff() || isArtificialDiff());
    }



protected:
    uwbADRSolverParams<T> m_params;

    int m_numVar;
    int m_varDofs;
    int m_nonZerosPerCol;

    int m_dofs;
    int m_numThreads;
    int m_tarDim;

    gsDofMapper m_dofMapper;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution, m_oldTimeSol;

    gsMatrix<index_t> m_elementList;
    std::vector<boundaryInterface> m_iFaceList;
    gsMatrix<index_t> m_iFaceElementList;

    T m_diffusionCoeff;
    gsVector<T> m_advectionCoeff;
    T m_reactionCoeff;
    bool m_bLinearADRCoeffsSet;
    bool m_bNonlinADcoeffsSet;
    bool m_bNonlinADRcoeffsSet;

    gsField<T> m_advectionField;
    gsField<T> m_diffusionField;
    gsField<T> m_reactionField;
    gsMultiPatch<T> m_patchesCoeffs;
    gsMultiBasis<T> m_basesCoeffs;

    gsSparseMatrix<T> m_mMass;
    gsSparseMatrix<T> m_mADRpattern;
    gsSparseMatrix<T> m_mADR;
    gsMatrix<T> m_rhsMass;
    gsMatrix<T> m_rhsADR;
};

} // namespace gismo
