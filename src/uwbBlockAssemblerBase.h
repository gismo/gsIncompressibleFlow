/** @file uwbBlockAssemblerBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include <gsCore/gsGeometryEvaluator.h>
#include "uwbINSSolverParams.h"
#include "uwbINSBlockVisitors.h"
#include "uwbINSdgBlockVisitors.h"

#include "uwbResiduumEvaluators.h"

namespace gismo
{

template<class T>
class uwbBlockAssemblerBase
{
public:
    uwbBlockAssemblerBase(const uwbINSSolverParams<T>& params) : m_params(params)
    {
        initMembers();
    }

    virtual ~uwbBlockAssemblerBase()
    {}

protected:

    virtual void initMembers()
    { 
        m_dofs = 0;
        m_tarDim = getPatches().dim();
        m_viscosity = m_params.getPde().getViscosity();
        setNumThreads(m_params.getNumThreads());
        m_bUnsteady = false;

        m_dofMappers.resize(2);
        m_ddof.resize(2);

        getBases().front().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.front(), 0);
        getBases().back().getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers.back(), 1);
       
        initElementList();

        if (getAssemblerOptions().intStrategy == iFace::dg)
            initiFaceElementList();

        /*if ((getBases().front().maxCwiseDegree() == getBases().back().maxCwiseDegree()) && isPSPG() == false)
        {
            gsInfo << "Equal order of the basis functions for velicity and pressure is chosen. PSPG stabilization is set.";
            m_params.settings().set(constantsINS::PSPG, true);
        }*/
    }

    // Create list of all elements for purpose of parallelization. 
    // First column contains patch indices and 
    // second column contains element indices on the given patch.
    void initElementList()
    {
        for (size_t npatch = 0; npatch < getPatches().nPatches(); ++npatch)
            for (size_t nelement = 0; nelement < getBases().front().basis(npatch).numElements(); ++nelement)
                m_elementList.push_back(std::make_pair(npatch, nelement));
    }

    //void initElementVolumeList()
    //{
    //    m_elementVolumeList.setZero(getBases().front().totalElements(), 1);
    //
    //    gsMatrix<T> quNodes;            // Mapped nodes
    //    gsVector<T> quWeights;         // Mapped weights
    //    unsigned evFlags(0);           // Evaluation flags for the Geometry map
    //    gsQuadRule<T> QuRule;
    //    uwbINSelementVolumeVisitor<T> visitor;
    //
    //    int nelementglobal = 0;
    //    for (unsigned npatch = 0; npatch < getPatches().nPatches(); ++npatch)
    //    {
    //        const gsBasis<T> * basis(&getBases().front().basis(npatch));
    //        visitor.initialize(*basis, npatch, QuRule, evFlags);
    //        // Initialize geometry evaluator
    //        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, getPatches().patch(npatch)));
    //        // Initialize domain element iterator
    //        typename gsBasis<T>::domainIter domIt = basis->makeDomainIterator();
    //
    //        for (; domIt->good(); domIt->next()) {
    //                
    //            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
    //            visitor.assemble(*domIt, *geoEval, quNodes, quWeights);
    //
    //            m_elementVolumeList(nelementglobal) = visitor.getVolume();
    //
    //            visitor.newElement();
    //                
    //            ++nelementglobal;
    //        }
    //    }
    //}

    // Create list of all interface elements for purpose of parallelization. 
    // First column contains patch indices and 
    // second column contains element indices on the given patch.
    void initiFaceElementList()
    {
        m_iFaceList.clear();
        std::vector<int> iFaceElementNumberList;
        int totaliFaceElementNumber = 0;

        for (typename gsMultiPatch<T>::const_iiterator it = getPatches().iBegin(); it != getPatches().iEnd(); ++it)
        {

            T size1 = 0;
            T size2 = 0;

            const int patch1 = it->first().patch;
            const int patch2 = it->second().patch;
            const gsBasis<T> * basis1(&getPatches().basis(patch1));
            const gsBasis<T> * basis2(&getPatches().basis(patch2));

            gsMatrix<T> quNodes;// Mapped nodes
            gsVector<T> quWeights;         // Mapped weights
            unsigned evFlags(0);           // Evaluation flags for the Geometry map
            gsQuadRule<T> QuRule;

            uwbINSiFaceSizeVisitor<T> visitor;

            // compute size1 of patch1
            visitor.initialize(*basis1, it->first().side(), patch1, QuRule, evFlags);

            // Initialize geometry evaluator
            typename gsGeometryEvaluator<T>::uPtr geoEval1(getEvaluator(evFlags, getPatches().patch(patch1)));

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt1 = basis1->makeDomainIterator(it->first().side());
            for (; domIt1->good(); domIt1->next())
            {
                // Compute the quadrature rule on patch1
                QuRule.mapTo(domIt1->lowerCorner(), domIt1->upperCorner(), quNodes, quWeights);

                // Assemble on element
                visitor.assemble(*domIt1, *geoEval1, quNodes, quWeights);

            }
            size1 = visitor.getSize();

            // compute size2 of patch2
            visitor.initialize(*basis2, it->second().side(), patch2, QuRule, evFlags);

            // Initialize geometry evaluators
            typename gsGeometryEvaluator<T>::uPtr geoEval2(getEvaluator(evFlags, getPatches().patch(patch2)));

            // Initialize domain element iterators
            typename gsBasis<T>::domainIter domIt2 = basis2->makeDomainIterator(it->second().side());

            for (; domIt2->good(); domIt2->next())
            {
                // Compute the quadrature rule on patch1
                QuRule.mapTo(domIt2->lowerCorner(), domIt2->upperCorner(), quNodes, quWeights);

                // Assemble on element
                visitor.assemble(*domIt2, *geoEval2, quNodes, quWeights);
            }
            size2 = visitor.getSize();

            boundaryInterface iFace;

            if (math::almostEqual<10, T>(size1, size2)) // patches are the same size, first will be the patch with more elements
            {
                iFace = (getBases().at(0).piece(it->first().patch).numElements(it->first().side()) <
                    getBases().at(0).piece(it->second().patch).numElements(it->second().side()) ?
                    it->getInverse() : *it);
            }
            else if (size1 < size2)
            {
                iFace = *it;
            }
            else
            {
                iFace = it->getInverse();
            }

            //gsInfo << "Interface " << iFace.first().patch << " - " << iFace.second().patch << "\n";

            m_iFaceList.push_back(iFace);
            int iFaceElementNumber = getBases().at(0).piece(iFace.first().patch).numElements(iFace.first().side());
            iFaceElementNumberList.push_back(iFaceElementNumber);
            totaliFaceElementNumber += iFaceElementNumber;

        }

        for (size_t iFaceIndex = 0; iFaceIndex < m_iFaceList.size(); iFaceIndex++)
            for (int elementIndex = 0; elementIndex < iFaceElementNumberList[iFaceIndex]; elementIndex++)
                m_iFaceElementList.push_back(std::make_pair(iFaceIndex, elementIndex));

    }

    // Compute Dirichlet dofs for unknown unk which corresponds to the basis basisID
    // and save the coefficients in ddofVector.
    // Assuming unk and basisID to be 0 or 1.
    void computeDirichletDofs(const int unk, const int basisID, gsMatrix<T>& ddofVector)
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
            if (!it->parametric())
            {
                fpts = it->function()->eval(pts);
            }
            else if (it->parametric() && m_params.settings().isIgaDirichletGeometry())
            {
                gsMatrix<T> parPts;
                getIgaBCGeom().invertPoints(pts, parPts);
                fpts = it->function()->eval(parPts);
            }
            else
                fpts = it->function()->eval(gsPointGrid<T>(rr));

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
            std::vector< gsMultiBasis<T> > bases(getBases());
            std::vector< gsDofMapper > dofMappers(m_dofMappers);
            const std::vector<gsMatrix<T> > ddof(m_ddof);
            const index_t totalElementsNumber(m_elementList.size());
            const std::vector< std::pair<index_t, index_t> > elementList(m_elementList);

            gsQuadRule<T> QuRule;
            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            unsigned evFlags(0);
            typename gsGeometryEvaluator<T>::uPtr geoEval;
            memory::shared_ptr< gsBasisRefs<T> > basesPtr;
            typename gsBasis<T>::domainIter domIt;

            gsSparseMatrix<T> * blockPtr = threadBlockPtr[threadNumber];
            gsMatrix<T> * rhsPtr = threadRhsPtr[threadNumber];

            index_t patchIndex = -1;
            index_t elementIndex = -1;

            // Parallel assembly procedure.
            #pragma omp for
            for (index_t nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal) {

                ElementBlockVisitor blockVisitor(dofMappers, m_viscosity);

                blockSpecificSettings(&blockVisitor);

                if (patchIndex != elementList[nElementGlobal].first)
                {
                    patchIndex = elementList[nElementGlobal].first;

                    basesPtr.reset(new gsBasisRefs<T>(bases, patchIndex));
                    domIt = basesPtr->front().makeDomainIterator(boundary::none);
                    elementIndex = 0;
                }

                blockVisitor.initialize(*basesPtr, patchIndex, QuRule, evFlags);

                geoEval = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patchIndex]));

                if (elementIndex != elementList[nElementGlobal].second)
                {
                    while ((elementIndex < elementList[nElementGlobal].second) && (domIt->good())) {
                        domIt->next();
                        ++elementIndex;
                    }
                }

                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                if (isStabilization())
                {
                    T elemDiam = computeElementDiameter(domIt, patchIndex);
                    blockVisitor.setElementLength(elemDiam);
                }

                blockVisitor.evaluate(*basesPtr, *geoEval, quNodes);//, *domIt);

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

    virtual void blockSpecificSettings(uwbVisitorBase<T>* blockVisitor)
    { GISMO_NO_IMPLEMENTATION }

    template<class RhsVisitor>
    void assembleRhs(gsMatrix<T> & rhs)
    {
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
            std::vector< gsMultiBasis<T> > bases(getBases());
            std::vector< gsDofMapper > dofMappers(m_dofMappers);
            const index_t totalElementsNumber(m_elementList.size());
            const std::vector< std::pair<index_t, index_t> > elementList(m_elementList);

            gsQuadRule<T> QuRule;
            gsMatrix<T> quNodes;
            gsVector<T> quWeights;
            unsigned evFlags(0);
            typename gsGeometryEvaluator<T>::uPtr geoEval;
            memory::shared_ptr< gsBasisRefs<T> > basesPtr;
            typename gsBasis<T>::domainIter domIt;

            gsMatrix<T> * rhsPtr = threadRhsPtr[threadNumber];

            index_t patchIndex = -1;
            index_t elementIndex = -1;

            // Parallel assembly procedure.
            #pragma omp for
            for (index_t nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal) {

                RhsVisitor rhsVisitor(dofMappers, m_viscosity);

                rhsSpecificSettings(&rhsVisitor);

                if (patchIndex != elementList[nElementGlobal].first) {
                    patchIndex = elementList[nElementGlobal].first;

                    basesPtr.reset(new gsBasisRefs<T>(bases, patchIndex));
                    domIt = basesPtr->front().makeDomainIterator(boundary::none);
                    elementIndex = 0;
                }

                rhsVisitor.initialize(*basesPtr, patchIndex, QuRule, evFlags);

                geoEval = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patchIndex]));

                if (elementIndex != elementList[nElementGlobal].second) {
                    while ((elementIndex < elementList[nElementGlobal].second) && (domIt->good())) {
                        domIt->next();
                        ++elementIndex;
                    }
                }

                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

                if (isStabilization())
                {
                    T elemDiam = computeElementDiameter(domIt, patchIndex);
                    rhsVisitor.setElementLength(elemDiam);
                }

                rhsVisitor.evaluate(*basesPtr, *geoEval, quNodes);

                rhsVisitor.assemble(*geoEval, quWeights);

                rhsVisitor.localToGlobal(*rhsPtr);

                domIt->next();
                ++elementIndex;
            }

        }

        // Summing up rhs from all threads 
        // into threadRhsPtr[0] which is a 
        // pointer to the input variable & rhs
        for (int t = 1; t < m_numThreads; ++t)
            *threadRhsPtr[0] += *threadRhsPtr[t];

        for (int t = 1; t < m_numThreads; ++t) {
            delete threadRhsPtr[t];
        }
    }

    virtual void rhsSpecificSettings(uwbVisitorBase<T>* blockVisitor)
    { GISMO_NO_IMPLEMENTATION }

    template<class InterfaceBlockVisitor>
    void assembleBlockDg(gsSparseMatrix<T> & matrixBlock,
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
            threadRhsPtr[t] = new gsMatrix<T>(m_dofs, 1);
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
            std::vector< gsMultiBasis<T> > bases(getBases());
            std::vector< gsDofMapper > dofMappers(m_dofMappers);
            const std::vector<gsMatrix<T> > ddof(m_ddof);
            const index_t totalElementsNumber(m_iFaceElementList.size());
            const std::vector< std::pair<index_t, index_t> > elementList(m_iFaceElementList);
            std::vector<boundaryInterface> iFaceList(m_iFaceList);

            gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
            gsMatrix<T> physNodes;         // Physical coordinates of quadrature nodes
            gsVector<T> quWeights;         // Mapped weights
            unsigned evFlags(0);           // Evaluation flags for the Geometry map
            gsQuadRule<T> QuRule;
            typename gsGeometryEvaluator<T>::uPtr geoEval1;
            typename gsGeometryEvaluator<T>::uPtr geoEval2;
            memory::shared_ptr< gsBasisRefs<T> > bases1Ptr;
            memory::shared_ptr< gsBasisRefs<T> > bases2Ptr;
            typename gsBasis<T>::domainIter domIt1;
            int patch1 = -1;
            int patch2 = -1;

            gsSparseMatrix<T> * blockPtr = threadBlockPtr[threadNumber];
            gsMatrix<T> * rhsPtr = threadRhsPtr[threadNumber];

            const T pp = 1000.0;

            InterfaceBlockVisitor blockDgVisitor(pp, m_viscosity, dofMappers);
            blockSpecificSettings(&blockDgVisitor);

            index_t iFaceIndex = -1;
            index_t elementIndex = -1;

            #pragma omp for 
            for (index_t nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal) {

                // visitor iFace initialization
                if (iFaceIndex != elementList[nElementGlobal].first) {
                    iFaceIndex = elementList[nElementGlobal].first;

                    const boundaryInterface & iFace = iFaceList[iFaceIndex];

                    patch1 = iFace.first().patch;
                    patch2 = iFace.second().patch;

                    bases1Ptr.reset(new gsBasisRefs<T>(bases, patch1));
                    bases2Ptr.reset(new gsBasisRefs<T>(bases, patch2));

                    domIt1 = (*bases1Ptr)[0].makeDomainIterator(iFace.first().side());
                    elementIndex = 0;

                    blockDgVisitor.initialize(*bases1Ptr, iFace.first().side(), iFace.second().side(), patch1, patch2, QuRule, evFlags);

                    geoEval1 = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patch1]));
                    geoEval2 = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patch2]));
                }

                // visitor proceeds to correct element
                if (elementIndex != elementList[nElementGlobal].second) {
                    while ((elementIndex < elementList[nElementGlobal].second) && (domIt1->good())) {
                        domIt1->next();
                        ++elementIndex;
                    }
                }

                // Compute the quadrature rule on patch1
                QuRule.mapTo(domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);

                // Compute quadrature nodes on patch2
                patches.patch(patch1).eval_into(quNodes1, physNodes); // physical coordinates of quNodes1
                patches.patch(patch2).invertPoints(physNodes, quNodes2); // parameters for physNodes on patch2

                // Perform required evaluations on the quadrature nodes
                blockDgVisitor.evaluate(*bases1Ptr, *bases2Ptr, *geoEval1, *geoEval2, quNodes1, quNodes2);

                // Assemble on element
                blockDgVisitor.assemble(*domIt1, *geoEval1, *geoEval2, quWeights);

                // Push to global patch matrix (m_rhs is filled in place)
                blockDgVisitor.localToGlobal(ddof, *blockPtr, *rhsPtr);

                domIt1->next();
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

    template<class ElementResiduumVisitor>
    void assembleProblemResiduum(std::vector<T> & residuum)
    {
        // deep copies
        const gsMultiPatch<T> patches(getPatches());
        std::vector< gsMultiBasis<T> > bases(getBases());
        std::vector< gsDofMapper > dofMappers(m_dofMappers);
        const size_t totalElementsNumber(m_elementList.size());
        const std::vector< std::pair<index_t, index_t> > elementList(m_elementList);

        gsQuadRule<T> QuRule;
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        unsigned evFlags(0);
        typename gsGeometryEvaluator<T>::uPtr geoEval;
        memory::shared_ptr< gsBasisRefs<T> > basesPtr;
        typename gsBasis<T>::domainIter domIt;

        index_t patchIndex = -1;
        index_t elementIndex = -1;

        // Parallel assembly procedure.
        #pragma omp for
        for (index_t nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal) {

            ElementResiduumVisitor residuumVisitor(dofMappers, m_viscosity);

            residuumSpecificSettings(&residuumVisitor);

            if (patchIndex != elementList[nElementGlobal].first)
            {
                patchIndex = elementList[nElementGlobal].first;

                basesPtr.reset(new gsBasisRefs<T>(bases, patchIndex));
                domIt = basesPtr->front().makeDomainIterator(boundary::none);
                elementIndex = 0;
            }

            residuumVisitor.initialize(*basesPtr, patchIndex, QuRule, evFlags);

            geoEval = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, patches[patchIndex]));

            if (elementIndex != elementList[nElementGlobal].second)
            {
                while ((elementIndex < elementList[nElementGlobal].second) && (domIt->good())) {
                    domIt->next();
                    ++elementIndex;
                }
            }

            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            T elemDiam = computeElementDiameter(domIt, patchIndex);
            residuumVisitor.setElementLength(elemDiam);

            residuumVisitor.evaluate(*basesPtr, *geoEval, quNodes);//, *domIt);

            residuumVisitor.assemble(*geoEval, quWeights);

            residuum.push_back(residuumVisitor.getElementResiduum());

            domIt->next();
            ++elementIndex;
        }

    }

    virtual void residuumSpecificSettings(uwbEvaluatorBase<T>* residuumVisitor)
    { GISMO_NO_IMPLEMENTATION }

public:

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

    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
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

    //works only in 2D.To do: generalize!!!
    T computeAspectRatio(bool minAR = false)
    {
            // deep copies
            std::vector< gsMultiBasis<T> > bases(getBases());
            const index_t totalElementsNumber(m_elementList.size());
            const std::vector< std::pair<index_t, index_t> > elementList(m_elementList);

            //typename gsGeometry<T>::Evaluator geoEval;
            //unsigned evFlags = NEED_VALUE | NEED_MEASURE;
            gsMatrix<T> nodes;
            memory::shared_ptr< gsBasisRefs<T> > basesPtr;
            typename gsBasis<T>::domainIter domIt;

            gsVector<T> aspectRatios;
            aspectRatios.setZero(totalElementsNumber);
            //gsVector<T> aspectRatiosParam;
            //aspectRatiosParam.setZero(totalElementsNumber);

            index_t patchIndex = -1;
            index_t elementIndex = -1;

            for (index_t nElementGlobal = 0; nElementGlobal < totalElementsNumber; ++nElementGlobal)
            {
                if (patchIndex != elementList[nElementGlobal].first)
                {
                    patchIndex = elementList[nElementGlobal].first;

                    basesPtr.reset(new gsBasisRefs<T>(bases, patchIndex));
                    domIt = basesPtr->front().makeDomainIterator(boundary::none);
                    elementIndex = 0;
                }

                if (elementIndex != elementList[nElementGlobal].second)
                {
                    while ((elementIndex < elementList[nElementGlobal].second) && (domIt->good()))
                    {
                        domIt->next();
                        ++elementIndex;
                    }
                }

                gsVector<T> nodeA = domIt->lowerCorner();
                gsVector<T> nodeC = domIt->upperCorner();
                nodes.setZero(getPatches().dim(), 4);
                nodes.middleCols(0, 1) = nodeA; //point A of the element in parametric domain
                nodes.middleCols(2, 1) = nodeC; //point C of the element in parametric domain
                nodes(0, 1) = nodeC(0); //x-coordinate of point B of the element in parametric domain
                nodes(1, 1) = nodeA(1); //y-coordinate of point B of the element in parametric domain
                nodes(0, 3) = nodeA(0); //x-coordinate of point D of the element in parametric domain
                nodes(1, 3) = nodeC(1); //y-coordinate of point D of the element in parametric domain
                //gsMatrix<T> physNodes;
                //getPatches().patch(patchIndex).eval_into(nodes, physNodes); // physical coordinates of paramBoundaryNodes

                //---
                index_t numSubIntervals = 10;
                T stepABx = math::abs(nodes(0, 1) - nodes(0, 0)) / numSubIntervals;
                T stepADy = math::abs(nodes(1, 3) - nodes(1, 0)) / numSubIntervals;
                T lengthABphys = 0.;
                T lengthADphys = 0.;
                gsMatrix<T> integralPointsAB;
                gsMatrix<T> integralPointsAD;
                for (int k = 0; k < numSubIntervals; k++)
                {
                    integralPointsAB.setZero(2, 2);
                    integralPointsAD.setZero(2, 2);
                    integralPointsAB(1, 0) = nodeA(1);
                    integralPointsAB(1, 1) = nodeA(1);
                    integralPointsAD(0, 0) = nodeA(0);
                    integralPointsAD(0, 1) = nodeA(0);
                    integralPointsAB(0, 0) = nodeA(0) + k * stepABx;
                    integralPointsAB(0, 1) = nodeA(0) + (k + 1) * stepABx;
                    integralPointsAD(1, 0) = nodeA(1) + k * stepADy;
                    integralPointsAD(1, 1) = nodeA(1) + (k + 1) * stepADy;

                    gsMatrix<T> integralPointsABphys;
                    gsMatrix<T> integralPointsADphys;
                    getPatches().patch(patchIndex).eval_into(integralPointsAB, integralPointsABphys);
                    getPatches().patch(patchIndex).eval_into(integralPointsAD, integralPointsADphys);

                    lengthABphys += (integralPointsABphys.col(1) - integralPointsABphys.col(0)).norm();
                    lengthADphys += (integralPointsADphys.col(1) - integralPointsADphys.col(0)).norm();
                }
                if (lengthABphys < lengthADphys)
                    aspectRatios(elementIndex) = lengthADphys / lengthABphys;
                else
                    aspectRatios(elementIndex) = lengthABphys / lengthADphys;

                //aspectRatiosParam(elementIndex) = domIt->getMaxCellLength() / domIt->getMinCellLength();

                domIt->next();
                ++elementIndex;
            }
            //gsInfo << "aspectRatiosParam.maxCoeff() = " << aspectRatiosParam.maxCoeff() << "\n";
        if (minAR)
            return aspectRatios.minCoeff();
        else
            return aspectRatios.maxCoeff();
    }

    bool isUnsteady() const { return m_bUnsteady; }
    void setUnsteady(bool unsteady = true) { m_bUnsteady = unsteady; }

    int     getTauStabTypeTanhCSD() const { return m_params.settings().get(constantsINS::tauStabTypeTanhCSD); }
    int     getTanhCSDtype() const { return m_params.settings().get(constantsINS::tanhCSDtype); }
    int     getTanhCSDresidualType() const { return m_params.settings().get(constantsINS::tanhCSDresidualType); }
    T       getTanhCSDalpha() const { return m_params.settings().get(constantsINS::tanhCSDalpha); }
    
    T getViscosity() const { return m_viscosity; }

    const std::vector<gsDofMapper> &    getMappers() const { return m_dofMappers; }
    const std::vector<gsMatrix<T> > &   getDirichletDofs() const { return m_ddof; }
    const gsMatrix<T> &                 getSolution() const { return m_solution; }
    const gsMultiPatch<T>&              getPatches() const { return m_params.getPde().patches(); }
    const std::vector< gsMultiBasis<T> >& getBases() const { return m_params.getBases(); }
    std::vector< gsMultiBasis<T> >& getBases() { return m_params.getBases(); }
    const gsBoundaryConditions<T>&        getBCs() const { return m_params.getBCs(); }
    const gsFunction<T>&      getRhsFcn() const { return m_params.getPde().getRhs(); }
    gsAssemblerOptions  getAssemblerOptions() const { return m_params.getAssemblerOptions(); }
    gsGeometry<T>& getIgaBCGeom() const { return m_params.settings().getIgaDirichletGeometry(); }
    int     getParam(constantsINS::intConst name) const { return m_params.settings().get(name); }
    T       getParam(constantsINS::realConst name) const { return m_params.settings().get(name); }

    virtual bool isStabilization() { GISMO_NO_IMPLEMENTATION }

protected:    
    int m_dofs;
    int m_numThreads;
    int m_tarDim;
    T m_viscosity;
    bool m_bUnsteady;

    uwbINSSolverParams<T> m_params;
    std::vector<gsDofMapper> m_dofMappers;
    std::vector<gsMatrix<T> > m_ddof;
    gsMatrix<T> m_solution;

    std::vector< std::pair<index_t, index_t> > m_elementList, m_iFaceElementList;
    std::vector<boundaryInterface> m_iFaceList;

}; // class uwbBlockAssemblerBase

} // namespace gismo
