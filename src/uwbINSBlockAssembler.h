/** @file uwbINSBlockAssembler.h

Author(s): H. Hornikova, J. Sourek
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include <gsCore/gsGeometryEvaluator.h>
#include "uwbBlockAssemblerBase.h"
#include "uwbINSSUPGBlockVisitors.h"
#include "uwbINSPSPGBlockVisitors.h"
#include "uwbINSCrosswindVisitors.h"
#include "uwbRANSSUPGBlockVisitors.h"
#include "uwbRANSCrosswindVisitors.h"
#include "uwbRANSADBlockVisitor.h"
#include "uwbRANSCSDVisitors.h"
#include "uwbRANSisoADVisitor.h"
#include "uwbRANSBlockVisitors.h"
#include "uwbINSBlockVisitorsBnd.h"

namespace gismo
{

template<class T>
class uwbINSBlockAssembler : public uwbBlockAssemblerBase<T>
{
public:
    typedef uwbBlockAssemblerBase<T> Base;

public:

    uwbINSBlockAssembler(const uwbINSSolverParams<T>& params) :
        Base(params)
    {
        initMembers();
    }

    virtual ~uwbINSBlockAssembler()
    {}

protected:

    virtual void initMembers()
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

        if (getAssemblerOptions().intStrategy == iFace::dg)
        {
            m_blockPdg.resize(m_pdofs, m_pdofs);
            m_rhsPdg.setZero(m_dofs, 1);
        }

        if (getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            computeDirichletDofs(0, 0, m_ddof[0]);
            computeDirichletDofs(1, 1, m_ddof[1]);
        }

        m_currentVelField = constructSolution(m_solution, 0);
        m_oldTimeVelField = m_currentVelField;

        if (isCROSSWIND() || isRANScrosswind() || isPSPG() || isRANStanhCSD() || isRANSisoAD())
        {
            m_currentPresField = constructSolution(m_solution, 1);
            //gsInfo << "Currect solution set.\n";
        }

        m_bPCDbndPrepared = false;

        // rotation:
        if (isRotation())
            computeOmegaXrCoeffs();

        m_bRANSsupg = false;
        //SUPG
        if (isSUPG())
        {
            m_blockA_SUPG.resize(m_udofs, m_udofs);
            m_rhsA_SUPG.setZero(m_dofs, 1);
            m_blockB_SUPG.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockB_SUPG[s].resize(m_udofs, m_pdofs);
            m_rhsB_SUPG.setZero(m_dofs, 1);
            m_blockN_SUPG.resize(m_udofs, m_udofs);
            m_rhsN_SUPG.setZero(m_dofs, 1);
            m_blockM_SUPG.resize(m_udofs, m_udofs);
            m_rhsM_SUPG.setZero(m_dofs, 1);

            m_blockApattern_SUPG.resize(m_udofs, m_udofs);
            m_blockNpattern_SUPG.resize(m_udofs, m_udofs);
            m_blockMpattern_SUPG.resize(m_udofs, m_udofs);
            m_blockBpattern_SUPG.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockBpattern_SUPG[s].resize(m_udofs, m_pdofs);
        }
        if (isTCSD())
        {
            m_blockN_SUPG.resize(m_udofs, m_udofs);
            m_rhsN_SUPG.setZero(m_dofs, 1);
            m_blockM_SUPG.resize(m_udofs, m_udofs);
            m_rhsM_SUPG.setZero(m_dofs, 1);

            m_blockNpattern_SUPG.resize(m_udofs, m_udofs);
            m_blockMpattern_SUPG.resize(m_udofs, m_udofs);
        }

        if (isCROSSWIND())
        {
            m_blockCrosswind.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockCrosswind[s].resize(m_udofs, m_udofs);
            m_blockCrosswindPattern.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockCrosswindPattern[s].resize(m_udofs, m_udofs);
            m_rhsCrosswind.setZero(m_dofs, 1);
        }

        if (isPSPG())
            m_rhsPSPG.setZero(m_dofs, 1);

        m_bResiduumEval = false;
    }

    void initRANSMembers()
    {
        m_blockA_RANS.resize(m_udofs, m_udofs);
        m_blockEdiag_RANS.resize(m_tarDim);
        for (int s = 0; s < m_tarDim; ++s)
            m_blockEdiag_RANS[s].resize(m_udofs, m_udofs);
        m_rhsA_RANS.setZero(m_dofs, 1);
        m_rhsEdiag_RANS.setZero(m_dofs, 1);
        m_rhsEnondiag_RANS.setZero(m_dofs, 1);

        if (isSUPG())
        {
            m_blockA_RANS_SUPG.resize(m_udofs, m_udofs);
            m_rhsA_RANS_SUPG.setZero(m_dofs, 1);
        }

        if (isRANScrosswind())
        {
            m_blockCrosswind.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockCrosswind[s].resize(m_udofs, m_udofs);
            m_blockCrosswindPattern.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockCrosswindPattern[s].resize(m_udofs, m_udofs);
            m_rhsCrosswind.setZero(m_dofs, 1);
        }

        if (isRANStanhCSD())
        {
            m_blockTanhCSD.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockTanhCSD[s].resize(m_udofs, m_udofs);
            m_blockTanhCSDpattern.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockTanhCSDpattern[s].resize(m_udofs, m_udofs);
            m_rhsTanhCSD.setZero(m_dofs, 1);
        }

        if (isRANSisoAD())
        {
            m_blockIsoAD.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockIsoAD[s].resize(m_udofs, m_udofs);
            m_blockIsoADpattern.resize(m_tarDim);
            for (int s = 0; s < m_tarDim; ++s)
                m_blockIsoADpattern[s].resize(m_udofs, m_udofs);
            m_rhsIsoAD.setZero(m_dofs, 1);
        }

        if (isRANSad())
        {
            m_block_RANS_AD.resize(m_udofs, m_udofs);
            m_rhs_RANS_AD.setZero(m_dofs, 1);
        }

        /*if (isPSPG())
        {
            m_rhsPSPG_RANS.setZero(m_dofs, 1);
        }*/
    }

    // computation of (omega x r) coefficients with L2-projection  
    void computeOmegaXrCoeffs()
    {
        gsFunctionExpr<T> * omegaXr;
        std::ostringstream s1, s2;

        switch (m_tarDim)
        {
        case 2:
            s1 << -getOmega() << " * y";
            s2 << getOmega() << " * x";
            omegaXr = new gsFunctionExpr<T>(s1.str(), s2.str(), 2);
            break;

        case 3:
            s1 << -getOmega() << " * z";
            s2 << getOmega() << " * y";
            omegaXr = new gsFunctionExpr<T>("0", s1.str(), s2.str(), 3);
            break;

        default:
            GISMO_ERROR("Rotation is implemented only for 2D and 3D.");
            break;
        }

        gsDofMapper mapper;
        getBases().at(0).getMapper(getAssemblerOptions().intStrategy, mapper);
        int matSize = mapper.freeSize();
        m_omegaXrCoeffs.resize(matSize, m_tarDim);

        gsSparseMatrix<T> projMatrix(matSize, matSize);
        gsMatrix<T> projRhs(matSize, m_tarDim);
        projRhs.setZero();

        // nonZeroPerCols evaluation
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_tarDim; i++)
            nonZerosPerCol *= 2 * getBases().at(0).maxDegree(i) + 1;

        projMatrix.reserve(gsVector<int>::Constant(matSize, nonZerosPerCol));

        gsMatrix<T> quNodes;
        gsVector<T> quWeights;
        gsQuadRule<T> QuRule;

        gsMatrix<T> basisVals;
        gsMatrix<T> rhsVals;
        gsMatrix<index_t> globIdxAct;

        for (unsigned int npatch = 0; npatch < getPatches().nPatches(); ++npatch)
        {
            const gsBasis<T> & basis = (getBases().at(0)).piece(npatch);

            gsVector<index_t> numQuadNodes(basis.dim());
            for (int i = 0; i < basis.dim(); ++i)
                numQuadNodes[i] = basis.maxDegree() + 1;
            QuRule = gsGaussRule<T>(numQuadNodes);

            typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_MEASURE, getPatches().patch(npatch)));

            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

            for (; domIt->good(); domIt->next())
            {
                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
                geoEval->evaluateAt(quNodes);
                rhsVals = omegaXr->eval(getPatches().patch(npatch).eval(quNodes));
                basis.eval_into(quNodes, basisVals);
                basis.active_into(quNodes.col(0), globIdxAct);
                mapper.localToGlobal(globIdxAct, npatch, globIdxAct);

                for (index_t k = 0; k < quNodes.cols(); k++)
                {
                    const T weight_k = quWeights[k] * geoEval->measure(k);

                    for (int i = 0; i < globIdxAct.size(); i++)
                    {
                        const unsigned ii = globIdxAct(i);

                        for (int j = 0; j < globIdxAct.size(); j++)
                        {
                            const unsigned jj = globIdxAct(j);

                            projMatrix.coeffRef(ii, jj) += weight_k * basisVals(i, k) * basisVals(j, k);
                        }

                        projRhs.row(ii) += weight_k *  basisVals(i, k) * rhsVals.col(k).transpose();
                    }
                }
            }
        }

        projMatrix.makeCompressed();
        typename gsSparseSolver<T>::CGDiagonal solver;
        m_omegaXrCoeffs = solver.compute(projMatrix).solve(projRhs);
    }

public:

    void constructSolution(const gsMatrix<T>& solVector,
        gsMultiPatch<T>& result, int unk, bool relative = false) const
    {
        GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

        if (!isRotation())
            relative = false;

        const gsDofMapper & mapper = m_dofMappers[unk];

        gsDofMapper diffMapper;

        if (unk == 0 && relative)
        {
            getBases().at(0).getMapper(getAssemblerOptions().intStrategy, diffMapper);
        }

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

                    if (unk == 0 && relative)
                        coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(i, p));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    // Assuming Taylor-Hood 
                    if (!unk) //velocity
                    {
                        coeffs.row(i) = m_ddof[unk].row(mapper.bindex(i, p));

                        if (relative)
                            coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(i, p));
                    }
                    else //pressure
                        coeffs.row(i).setZero();
                }
            }

            result.addPatch(getBases().at(unk).piece(p).makeGeometry(coeffs));
        }
    }

    gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, bool relative = false) const
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructSolution(solVector, *result, unk, relative);
        return gsField<T>(getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
    }

    void constructSolution(const gsMatrix<T>& solVector,
        gsMultiPatch<T>& result, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, int unk, bool relative = false) const
    {
        GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

        if (!isRotation())
            relative = false;

        const gsDofMapper & mapper = m_dofMappers[unk];

        gsDofMapper diffMapper;

        if (unk == 0 && relative)
        {
            getBases().at(0).getMapper(getAssemblerOptions().intStrategy, diffMapper);
        }

        result.clear(); // result is cleared first

        const index_t dim = (unk == 0 ? m_tarDim : 1);
        gsMatrix<T> coeffs;
        gsMatrix<index_t> sliceCoefs;

        // Point to the correct entries of the solution vector
        gsAsConstMatrix<T> solV = (unk == 0 ?
            gsAsConstMatrix<T>(solVector.data(), m_udofs, dim)
            :
            gsAsConstMatrix<T>(solVector.data() + m_pshift, m_pdofs, 1)
            );

        for (int j = 0; j < patchNumbers.size(); ++j)
        {
            index_t p = patchNumbers(j);

            sliceCoefs = getBases().at(unk).piece(p).boundary(sides[j]);

            // Reconstruct solution coefficients on patch p
            const int sz = sliceCoefs.rows();
            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(sliceCoefs.at(i), p)) // DoF value is in the solVector
                {
                    coeffs.row(i) = solV.row(mapper.index(sliceCoefs.at(i), p));

                    if (unk == 0 && relative)
                        coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(sliceCoefs.at(i), p));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    // Assuming Taylor-Hood 
                    if (!unk) //velocity
                    {
                        coeffs.row(i) = m_ddof[unk].row(mapper.bindex(sliceCoefs.at(i), p));

                        if (relative)
                            coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(sliceCoefs.at(i), p));
                    }
                    else //pressure
                        coeffs.row(i).setZero();
                }
            }

            result.addPatch(getBases().at(unk).piece(p).boundaryBasis(sides[j])->makeGeometry(coeffs));
        }
    }

    gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructSolution(solVector, *result, patchNumbers, sides, unk, relative);

        gsMultiPatch<T> * slice = new gsMultiPatch<T>; // Memory leak! How to do it better?
        for (int p = 0; p < patchNumbers.size(); p++)
        {
            slice->addPatch(getPatches().patch(patchNumbers(p)).boundary(sides[p]));
        }
        return gsField<T>(*slice, typename gsFunctionSet<T>::Ptr(result), true);
    }

    void constructSolutionCombined(const gsMatrix<T>& solVector,
        gsMultiPatch<T>& result, gsVector<size_t> relPatches) const
    {
        GISMO_ASSERT(m_dofs == solVector.rows(), "Something went wrong, is solution vector valid?");

        const gsDofMapper & mapper = m_dofMappers[0];

        gsDofMapper diffMapper;
        getBases().at(0).getMapper(getAssemblerOptions().intStrategy, diffMapper);

        result.clear(); // result is cleared first

        const index_t dim = m_tarDim;
        gsMatrix<T> coeffs;

        // Point to the correct entries of the solution vector
        gsAsConstMatrix<T> solV = gsAsConstMatrix<T>(solVector.data(), m_udofs, dim);


        bool rel = false;

        for (unsigned int p = 0; p < getPatches().nPatches(); ++p)
        {

            if ((relPatches.array() == p).any() && isRotation())
                rel = true;
            else
                rel = false;

            // Reconstruct solution coefficients on patch p
            const int sz = getBases().at(0).piece(p).size();
            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, p)) // DoF value is in the solVector
                {
                    coeffs.row(i) = solV.row(mapper.index(i, p));

                    if (rel)
                        coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(i, p));
                }
                else // eliminated DoF: fill with Dirichlet data
                {

                    coeffs.row(i) = m_ddof[0].row(mapper.bindex(i, p));

                    if (rel)
                        coeffs.row(i) -= m_omegaXrCoeffs.row(diffMapper.index(i, p));

                }
            }

            result.addPatch(getBases().at(0).piece(p).makeGeometry(coeffs));
        }
    }

    gsField<T> constructSolutionCombined(const gsMatrix<T>& solVector, gsVector<size_t> relPatches) const
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructSolutionCombined(solVector, *result, relPatches);
        return gsField<T>(getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
    }

    void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
    {
        if (updateSol)
        {
            m_solution = solVector;

            if(isUnsteady())
                m_oldTimeVelField = constructSolution(solVector, 0);
        }

        m_currentVelField = constructSolution(solVector, 0, isRotation()); // use relative velocity for computation of block N (if rotational)

        if (isCROSSWIND() || isRANScrosswind() || isPSPG() || isRANStanhCSD() || isRANSisoAD() || isResiduumEval())
            m_currentPresField = constructSolution(solVector, 1);
    }

    // computes flow rate through a side of a patch
    T computeFlowRate(int patch, boxSide side, gsMatrix<T> solution) const
    {
        T flowRate = 0;

        gsField<T> solutionField = constructSolution(solution, 0); // velocity field

        const gsBasis<T>& basis = getBases().at(0).basis(patch);

        gsVector<int> numQuadNodes(m_tarDim);
        const int dir = side.direction();
        for (int i = 0; i < m_tarDim; ++i)
            numQuadNodes[i] = (2 * basis.degree(i) + 1);
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        gsGaussRule<T> QuRule(numQuadNodes);// harmless slicing occurs here

        gsMatrix<T> quNodes; // Mapped nodes
        gsVector<T> quWeights; // Mapped weights
        unsigned evFlags = NEED_VALUE | NEED_OUTER_NORMAL;   // Evaluation flags for the Geometry map

        // Initialize geometry evaluator
        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, getPatches().patch(patch)));

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(side);
        for (; domIt->good(); domIt->next())
        {
            // Compute the quadrature rule on patch1
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

            // Compute image of Gauss nodes under geometry mapping as well as Jacobians
            geoEval->evaluateAt(quNodes);

            // Evaluate solution on element nodes
            gsMatrix<T> solUVals = solutionField.value(quNodes, patch);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Compute the outer normal vector from patch1
                gsVector<T> normal;
                geoEval->outerNormal(k, side, normal); 
                
                // the normal norm is equal to integral measure
                flowRate += quWeights[k] * normal.dot(solUVals.col(k));
            }
        }
        return flowRate;
    }

protected:
    virtual void blockSpecificSettings(uwbVisitorBase<T>* blockVisitor)
    {
        std::vector<gsField<T> > sols;
        sols.push_back(m_currentVelField);
        blockVisitor->setCurrentSolution(sols);

        uwbRANSBlockVisitor<T>* pVisitor = dynamic_cast<uwbRANSBlockVisitor<T>*>(blockVisitor);

        if (pVisitor != NULL)
            pVisitor->setTurbulenceSolver(m_params.settings().getTurbulenceSolver());

        if (isSUPG() || isTCSD())
        {
            uwbINSSUPGBlockVisitor<T>* pVisitorSUPG = dynamic_cast<uwbINSSUPGBlockVisitor<T>*>(blockVisitor);
            if (pVisitorSUPG != NULL)
            {
                pVisitorSUPG->setTauStabType(getTauStabTypeSUPG(), isUnsteady(), getTimeStep());
                if (m_bRANSsupg)
                    pVisitorSUPG->setTurbulenceSolver(m_params.settings().getTurbulenceSolver());
            }

            uwbRANSSUPGBlockVisitor<T>* pVisitorRANSsupg = dynamic_cast<uwbRANSSUPGBlockVisitor<T>*>(blockVisitor);
            if (pVisitorRANSsupg != NULL)
                pVisitorRANSsupg->setTauStabType(getTauStabTypeSUPG(), getTimeStep(), getDiamPowerSUPG());
        }
        if (isCROSSWIND())
        {
            uwbINSCrosswindVisitor<T>* pVisitorCROSSWIND = dynamic_cast<uwbINSCrosswindVisitor<T>*>(blockVisitor);
            if (pVisitorCROSSWIND != NULL)
            {
                pVisitorCROSSWIND->setCrosswind(getCrosswindType(), getCWresidualType(), getTauStabTypeCW(), isUnsteady(), getTimeStep(), getDiamPowerCW());
                pVisitorCROSSWIND->setPressureSolutionField(m_currentPresField);
                if (isUnsteady())
                    pVisitorCROSSWIND->setOldSolutionField(m_oldTimeVelField);
            }
        }
        if (isPSPG())
        {
            uwbINSPSPGBlockVisitors<T>* pVisitorPSPG = dynamic_cast<uwbINSPSPGBlockVisitors<T>*>(blockVisitor);
            if (pVisitorPSPG != NULL)
            {
                pVisitorPSPG->setPSPG(getTauStabTypePSPG(), isUnsteady(), getTimeStep());
                pVisitorPSPG->setPressureSolutionField(m_currentPresField);
                if (isUnsteady())
                    pVisitorPSPG->setOldSolutionField(m_oldTimeVelField);
            }
        }
        if (isRANScrosswind())
        {
            uwbRANSCrosswindVisitor<T>* pVisitorCROSSWIND = dynamic_cast<uwbRANSCrosswindVisitor<T>*>(blockVisitor);
            if (pVisitorCROSSWIND != NULL)
            {
                pVisitorCROSSWIND->setCrosswind(getCrosswindType(), getCWresidualType(), getTauStabTypeCW(), isUnsteady(), getTimeStep(), getDiamPowerCW());
                pVisitorCROSSWIND->setPressureSolutionField(m_currentPresField);
                if (isUnsteady())
                    pVisitorCROSSWIND->setOldSolutionField(m_oldTimeVelField);
            }
        }
        if (isRANSad())
        {
            uwbRANSADBlockVisitor<T>* pVisitorAD = dynamic_cast<uwbRANSADBlockVisitor<T>*>(blockVisitor);
            if (pVisitorAD != NULL)
                pVisitorAD->setADType(getTauStabTypeAD(), getTimeStep(), getDiamPowerAD());
        }
        if (isRANStanhCSD())
        {
            uwbRANSCSDVisitor<T>* pVisitorCSD = dynamic_cast<uwbRANSCSDVisitor<T>*>(blockVisitor);
            if (pVisitorCSD != NULL)
            {
                pVisitorCSD->setCSD(getTanhCSDtype(), getTanhCSDresidualType(), getTauStabTypeTanhCSD(), isUnsteady(), getTimeStep(), getTanhCSDalpha());
                pVisitorCSD->setPressureSolutionField(m_currentPresField);
                if (isUnsteady())
                    pVisitorCSD->setOldSolutionField(m_oldTimeVelField);
            }
        }
        if (isRANSisoAD())
        {
            uwbRANSisoADVisitor<T>* pVisitorIsoAD = dynamic_cast<uwbRANSisoADVisitor<T>*>(blockVisitor);
            if (pVisitorIsoAD != NULL)
            {
                pVisitorIsoAD->setIsoAD(getIsoADtype(), getIsoADresidualType(), isUnsteady(), getTimeStep(), getIsoADalpha());
                pVisitorIsoAD->setPressureSolutionField(m_currentPresField);
                if (isUnsteady())
                    pVisitorIsoAD->setOldSolutionField(m_oldTimeVelField);
            }
        }
    }

    virtual void rhsSpecificSettings(uwbVisitorBase<T>* rhsVisitor)
    {
        uwbINSRhsVisitor<T>* pVisitor = dynamic_cast<uwbINSRhsVisitor<T>*>(rhsVisitor);

        if (pVisitor != NULL)
            pVisitor->setRhsFunction(getRhsFcn());
    }

public:

    void assembleLinearStokesPart()
    {
        // matrix and rhs cleaning
        m_rhsA.setZero();
        m_blockA.resize(m_udofs, m_udofs);
        gsSparseMatrix<T> blockAsym(m_udofs, m_udofs);

        m_rhsB.setZero();

        gsSparseMatrix<T> blocksB(m_pdofs, m_pshift);

        if (getAssemblerOptions().intStrategy == iFace::dg)
        {
            m_rhsPdg.setZero();
            m_blockPdg.resize(m_pdofs, m_pdofs);
        }

        m_rhsF.setZero();

        // blocks assembly
        blockAsym.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSBlockAsymVisitor<T> >(blockAsym, m_rhsA, m_nonZerosPerColU);
        m_blockA = blockAsym.template selfadjointView<Eigen::Upper>();
        if (getAssemblerOptions().intStrategy == iFace::dg) 
        {
            gsSparseMatrix<T> blockAdg(m_udofs, m_udofs);
            blockAdg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
            Base::template assembleBlockDg< uwbINSdgBlockAVisitor<T> >(blockAdg, m_rhsA, m_nonZerosPerColU);
            m_blockA += blockAdg;
        }

        m_blockA.makeCompressed();

        blocksB.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSBlocksBVisitor<T> >(blocksB, m_rhsB, m_nonZerosPerColP);
        if (getAssemblerOptions().intStrategy == iFace::dg) 
        {
            gsSparseMatrix<T> blocksBdg(m_pdofs, m_pshift);
            blocksBdg.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
            Base::template assembleBlockDg< uwbINSdgBlockBVisitor<T> >(blocksBdg, m_rhsB, m_nonZerosPerColP);
            blocksB += blocksBdg;
        }
        blocksB.makeCompressed();

#pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockB[s] = gsSparseMatrix<T>(blocksB.middleCols(s * m_udofs, m_udofs));
            m_blockB[s].makeCompressed();
        }

#pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockMinusBT[s] = gsSparseMatrix<T, Eigen::ColMajor>(-m_blockB[s].transpose());
            m_blockMinusBT[s].makeCompressed();
        }

        if (getAssemblerOptions().intStrategy == iFace::dg)
        {
            m_blockPdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
            Base::template assembleBlockDg< uwbINSdgBlockPVisitor<T> >(m_blockPdg, m_rhsPdg, m_nonZerosPerColP);
            m_blockPdg.makeCompressed();
        }

        Base::template assembleRhs< uwbINSRhsVisitor<T> >(m_rhsF);

        if ((isUnsteady() && isTimeDerTerm()) || isRotation())
        {
            assembleMassMatrix();
        }
    }

    //===================== Assemble Nonlinear Part =====================
    void assembleNonlinearPart()
    {
        // matrix and rhs cleaning
        m_blockN.resize(m_udofs, m_udofs);
        m_rhsN.setZero();

        // blocks assembly
        m_blockN.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSBlockNVisitor<T> >(m_blockN, m_rhsN, m_nonZerosPerColU);
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            gsSparseMatrix<T> m_blockNdg(m_udofs, m_udofs);
            m_blockNdg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
            Base::template assembleBlockDg< uwbINSdgBlockNVisitor<T> >(m_blockNdg, m_rhsN, m_nonZerosPerColU);
            m_blockN += m_blockNdg;
        }
        m_blockN.makeCompressed();
    }

    void assembleNonlinearSUPGPart()
    {
        assembleBlockA_SUPG();
        assembleBlockB_SUPG();
        assembleBlockN_SUPG();
        if ( (isUnsteady() && isTimeDerTerm()) || isRotation())
            assembleBlockM_SUPG();
    }

    void assembleNonlinearCSDPart()
    {
        assembleBlockN_SUPG();
        if ( (isUnsteady() && isTimeDerTerm()) || isRotation())
            assembleBlockM_SUPG();
    }

    //========================= A_SUPG ===================================
    void assembleBlockA_SUPG()
    {
        // matrix and rhs cleaning
        m_blockA_SUPG.resize(m_udofs, m_udofs);
        m_rhsA_SUPG.setZero();

        // blocks assembly
        m_blockA_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockAVisitor<T> >(m_blockA_SUPG, m_rhsA_SUPG, m_nonZerosPerColU);
        m_blockA_SUPG.makeCompressed();

    }
    //end A_SUPG

    //===================== B_SUPG ============================================
    void assembleBlockB_SUPG()
    {
        // matrix and rhs cleaning
        m_rhsB_SUPG.setZero();

        gsSparseMatrix<T> blocksB_SUPG(m_pdofs, m_pshift);

        blocksB_SUPG.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSSUPGBlockBVisitor<T> >(blocksB_SUPG, m_rhsB_SUPG, m_nonZerosPerColP);
        blocksB_SUPG.makeCompressed();

#pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockB_SUPG[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksB_SUPG.middleCols(s * m_udofs, m_udofs)).transpose());//(m_blockB_SUPG_pom[s].transpose());
            m_blockB_SUPG[s].makeCompressed();
        }

    }
    //end B_SUPG


    //=========================== N_SUPG ===============================================
    void assembleBlockN_SUPG()
    {
        // matrix and rhs cleaning
        m_blockN_SUPG.resize(m_udofs, m_udofs);
        m_rhsN_SUPG.setZero();

        // blocks assembly
        m_blockN_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockNVisitor<T> >(m_blockN_SUPG, m_rhsN_SUPG, m_nonZerosPerColU);
        m_blockN_SUPG.makeCompressed();

    }
    //end N_SUPG

    //=========================== M_SUPG ===============================================
    void assembleBlockM_SUPG()
    {
        // matrix and rhs cleaning
        m_blockM_SUPG.resize(m_udofs, m_udofs);
        m_rhsM_SUPG.setZero();

        // blocks assembly
        m_blockM_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockMVisitor<T> >(m_blockM_SUPG, m_rhsM_SUPG, m_nonZerosPerColU);
        m_blockM_SUPG.makeCompressed();

    }
    //end M_SUPG

    void assembleCrosswindPart()
    {
        // matrix and rhs cleaning
        m_rhsCrosswind.setZero();

        gsSparseMatrix<T> blocksCrosswind(m_udofs, m_pshift);

        blocksCrosswind.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSCrosswindVisitor<T> >(blocksCrosswind, m_rhsCrosswind, m_nonZerosPerColU);
        blocksCrosswind.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockCrosswind[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCrosswind.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockCrosswind[s].makeCompressed();
        }
    }

    void assembleBlockCROSSWINDpattern()
    {
        gsMatrix<T> current_solution = m_solution;
        m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

        gsMatrix<T> dummy_rhs(m_dofs, 1);
        dummy_rhs.setZero();

        gsSparseMatrix<T> blocksCROSSWINDpattern(m_udofs, m_pshift);

        blocksCROSSWINDpattern.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSCrosswindVisitor<T> >(blocksCROSSWINDpattern, dummy_rhs, m_nonZerosPerColU);
        blocksCROSSWINDpattern.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockCrosswindPattern[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCROSSWINDpattern.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockCrosswindPattern[s].makeCompressed();
            m_blockCrosswindPattern[s] = 0. * m_blockCrosswindPattern[s];
        }

        m_solution = current_solution;
    }

    void assemblePSPGPart()
    {
        // matrix and rhs cleaning
        m_rhsPSPG.setZero();

        gsSparseMatrix<T> blocksPSPG(m_pdofs, m_udofs);

        blocksPSPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSPSPGBlockVisitors<T> >(blocksPSPG, m_rhsPSPG, m_nonZerosPerColP);
        blocksPSPG.makeCompressed();
    }

    void assembleRANScrosswindPart()
    {
        // matrix and rhs cleaning
        m_rhsCrosswind.setZero();

        gsSparseMatrix<T> blocksCrosswind(m_udofs, m_pshift);

        blocksCrosswind.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSCrosswindVisitor<T> >(blocksCrosswind, m_rhsCrosswind, m_nonZerosPerColU);
        blocksCrosswind.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockCrosswind[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCrosswind.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockCrosswind[s].makeCompressed();
        }
    }

    void assembleBlockRANScrosswindPattern()
    {
        gsMatrix<T> current_solution = m_solution;
        m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

        gsMatrix<T> dummy_rhs(m_dofs, 1);
        dummy_rhs.setZero();

        gsSparseMatrix<T> blocksCROSSWINDpattern(m_udofs, m_pshift);

        blocksCROSSWINDpattern.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSCrosswindVisitor<T> >(blocksCROSSWINDpattern, dummy_rhs, m_nonZerosPerColU);
        blocksCROSSWINDpattern.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockCrosswindPattern[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCROSSWINDpattern.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockCrosswindPattern[s].makeCompressed();
            m_blockCrosswindPattern[s] = 0. * m_blockCrosswindPattern[s];
        }

        m_solution = current_solution;
    }

    void assembleRANSTanhCSDpart()
    {
        // matrix and rhs cleaning
        m_rhsTanhCSD.setZero();

        gsSparseMatrix<T> blocksCSD(m_udofs, m_pshift);

        blocksCSD.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSCSDVisitor<T> >(blocksCSD, m_rhsTanhCSD, m_nonZerosPerColU);
        blocksCSD.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockTanhCSD[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCSD.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockTanhCSD[s].makeCompressed();
        }
    }

    void assembleBlockRANStanhCSDpattern()
    {
        gsMatrix<T> current_solution = m_solution;
        m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

        gsMatrix<T> dummy_rhs(m_dofs, 1);
        dummy_rhs.setZero();

        gsSparseMatrix<T> blocksCSDpattern(m_udofs, m_pshift);

        blocksCSDpattern.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSCSDVisitor<T> >(blocksCSDpattern, dummy_rhs, m_nonZerosPerColU);
        blocksCSDpattern.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockTanhCSDpattern[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksCSDpattern.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockTanhCSDpattern[s].makeCompressed();
            m_blockTanhCSDpattern[s] = 0. * m_blockTanhCSDpattern[s];
        }

        m_solution = current_solution;
    }

    void assembleRANSisoADpart()
    {
        // matrix and rhs cleaning
        m_rhsIsoAD.setZero();

        gsSparseMatrix<T> blocksIsoAD(m_udofs, m_pshift);

        blocksIsoAD.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSisoADVisitor<T> >(blocksIsoAD, m_rhsIsoAD, m_nonZerosPerColU);
        blocksIsoAD.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockIsoAD[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksIsoAD.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockIsoAD[s].makeCompressed();
        }
    }

    void assembleBlockRANSisoADpattern()
    {
        gsMatrix<T> current_solution = m_solution;
        m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

        gsMatrix<T> dummy_rhs(m_dofs, 1);
        dummy_rhs.setZero();

        gsSparseMatrix<T> blocksIsoADpattern(m_udofs, m_pshift);

        blocksIsoADpattern.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSisoADVisitor<T> >(blocksIsoADpattern, dummy_rhs, m_nonZerosPerColU);
        blocksIsoADpattern.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockIsoADpattern[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksIsoADpattern.middleCols(s * m_udofs, m_udofs)).transpose());
            m_blockIsoADpattern[s].makeCompressed();
            m_blockIsoADpattern[s] = 0. * m_blockIsoADpattern[s];
        }

        m_solution = current_solution;
    }

    void assembleMassMatrix()
    {
        // matrix and rhs cleaning
        m_rhsM.setZero();
        m_blockM.resize(m_udofs, m_udofs);
        gsSparseMatrix<T> blockMsym(m_udofs, m_udofs);

        // blocks assembly
        blockMsym.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSBlockMVisitor<T> >(blockMsym, m_rhsM, m_nonZerosPerColU);
        m_blockM = blockMsym.template selfadjointView<Eigen::Upper>();
        m_blockM.makeCompressed();
    }

    void assemblePressureMassMatrix()
    {
        // matrix and rhs cleaning
        m_rhsMp.setZero();
        m_blockMp.resize(m_pdofs, m_pdofs);
        gsSparseMatrix<T> blockMpsym(m_pdofs, m_pdofs);

        // blocks assembly
        blockMpsym.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSBlockMpVisitor<T> >(blockMpsym, m_rhsMp, m_nonZerosPerColP);
        m_blockMp = blockMpsym.template selfadjointView<Eigen::Upper>();
        m_blockMp.makeCompressed();
    }

    void assembleBlocksC()
    {
        // matrix and rhs cleaning
        m_rhsC.setZero();

        gsSparseMatrix<T> blocksC(m_pdofs, m_pshift);

        blocksC.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSBlocksCVisitor<T> >(blocksC, m_rhsC, m_nonZerosPerColP);
        blocksC.makeCompressed();

#pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockCT[s] = gsSparseMatrix<T, Eigen::ColMajor>(blocksC.middleCols(s * m_udofs, m_udofs).transpose());
            m_blockCT[s].makeCompressed();
        }
    }

    //=================== assembleBlockRANS =========================
    void assembleRANSPart()
    {
        //------------ A_RANS ---------------
        m_rhsA_RANS.setZero();
        m_blockA_RANS.resize(m_udofs, m_udofs);
        gsSparseMatrix<T> blockAsym_RANS(m_udofs, m_udofs);

        blockAsym_RANS.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSBlockAsymVisitor<T> >(blockAsym_RANS, m_rhsA_RANS, m_nonZerosPerColU);
        m_blockA_RANS = blockAsym_RANS.template selfadjointView<Eigen::Upper>();
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            GISMO_ERROR("RANS is not implemented in DG case.");
        }
        m_blockA_RANS.makeCompressed();

        //------------ EdiagSym_RANS ---------------
        m_rhsEdiag_RANS.setZero();
        gsSparseMatrix<T> blocksEdiag_RANS(m_udofs, m_pshift);

        blocksEdiag_RANS.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSBlockEdiagsymVisitor<T> >(blocksEdiag_RANS, m_rhsEdiag_RANS, m_nonZerosPerColU);
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            GISMO_ERROR("RANS is not implemented in DG case.");
        }
        blocksEdiag_RANS.makeCompressed();

#pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockEdiag_RANS[s] = gsSparseMatrix<T, Eigen::ColMajor>(blocksEdiag_RANS.middleCols(s * m_udofs, m_udofs));
            m_blockEdiag_RANS[s].makeCompressed();
        }

        //------------ Erhs_RANS ---------------
        m_rhsEnondiag_RANS.setZero();
        int numBlocks;
        switch (m_tarDim)
        {
        case 2:
            numBlocks = m_tarDim - 1;
            break;
        case 3:
            numBlocks = m_tarDim;
            break;
        default:
            GISMO_ERROR("RANS are implemented only for 2D and 3D.");
            break;
        }

        gsSparseMatrix<T> blocksEnondiag_RANS(m_udofs, numBlocks * m_udofs);
        std::vector< gsSparseMatrix<T> > blockEnondiag_RANS;
        blockEnondiag_RANS.resize(numBlocks);
        for (int s = 0; s < numBlocks; ++s)
        {
            blockEnondiag_RANS[s].resize(m_udofs, m_udofs);
        }

        blocksEnondiag_RANS.reserve(gsVector<int>::Constant(numBlocks * m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSRhsEnondiagVisitor<T> >(blocksEnondiag_RANS, m_rhsEnondiag_RANS, m_nonZerosPerColU);
        blocksEnondiag_RANS.makeCompressed();

#pragma omp parallel for
        for (int s = 0; s < numBlocks; ++s) {
            blockEnondiag_RANS[s] = gsSparseMatrix<T, Eigen::ColMajor>(blocksEnondiag_RANS.middleCols(s * m_udofs, m_udofs));
            blockEnondiag_RANS[s].makeCompressed();
        }

        m_rhsEnondiag_RANS.middleRows(0, m_udofs).noalias() -= blockEnondiag_RANS[0] * m_solution.middleRows(m_udofs, m_udofs); //E12
        m_rhsEnondiag_RANS.middleRows(m_udofs, m_udofs).noalias() -= blockEnondiag_RANS[0].transpose() * m_solution.middleRows(0, m_udofs); //E21
        if (numBlocks == m_tarDim)
        {
            m_rhsEnondiag_RANS.middleRows(0, m_udofs).noalias() -= blockEnondiag_RANS[1] * m_solution.middleRows(2 * m_udofs, m_udofs); //E13
            m_rhsEnondiag_RANS.middleRows(2 * m_udofs, m_udofs).noalias() -= blockEnondiag_RANS[1].transpose() * m_solution.middleRows(0, m_udofs); //E31
            m_rhsEnondiag_RANS.middleRows(m_udofs, m_udofs).noalias() -= blockEnondiag_RANS[2] * m_solution.middleRows(2 * m_udofs, m_udofs); //E23
            m_rhsEnondiag_RANS.middleRows(2 * m_udofs, m_udofs).noalias() -= blockEnondiag_RANS[2].transpose() * m_solution.middleRows(m_udofs, m_udofs); //E32
        }

        if (isSUPG())
        {
            //------------ A_RANS_SUPG ---------------
            m_rhsA_RANS_SUPG.setZero();
            m_blockA_RANS_SUPG.resize(m_udofs, m_udofs);

            m_blockA_RANS_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
            Base::template assembleBlock< uwbRANSSUPGBlockAVisitor<T> >(m_blockA_RANS_SUPG, m_rhsA_RANS_SUPG, m_nonZerosPerColU);
            if (getAssemblerOptions().intStrategy == iFace::dg) {
                GISMO_ERROR("RANS is not implemented in DG case.");
            }
            m_blockA_RANS_SUPG.makeCompressed();
        }
    }

    //=========================== RANS AD ===============================================
    void assembleBlockRANSad()
    {
        // matrix and rhs cleaning
        m_block_RANS_AD.resize(m_udofs, m_udofs);
        m_rhs_RANS_AD.setZero();

        // blocks assembly
        m_block_RANS_AD.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbRANSADBlockVisitor<T> >(m_block_RANS_AD, m_rhs_RANS_AD, m_nonZerosPerColU);
        m_block_RANS_AD.makeCompressed();
    }

    //=================== assembleBlockNpattern =========================
    void assembleBlockNpattern()
    {
        gsMatrix<T> current_solution = m_solution;

        m_solution = gsMatrix<T>::Constant(m_dofs, 1, 1.);

        // matrix and rhs cleaning
        m_blockNpattern.resize(m_udofs, m_udofs);

        gsMatrix<T> dummy_rhsN(m_dofs, 1);
        dummy_rhsN.setZero();

        // blocks assembly
        m_blockNpattern.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSBlockNVisitor<T> >(m_blockNpattern, dummy_rhsN, m_nonZerosPerColU);
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            gsSparseMatrix<T> blockNpatterndg(m_udofs, m_udofs);
            blockNpatterndg.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
            Base::template assembleBlockDg< uwbINSdgBlockNVisitor<T> >(blockNpatterndg, dummy_rhsN, m_nonZerosPerColU);
            m_blockNpattern += blockNpatterndg;
        }
        m_blockNpattern.makeCompressed();
        m_blockNpattern = 0. * m_blockNpattern;

        m_solution = current_solution;
    }

    //=================== assembleBlockNonlinearSUPGpattern =========================
    void assembleBlockNonlinearSUPGpattern()
    {
        gsMatrix<T> current_solution = m_solution;

        gsSparseMatrix<T> blocksBpattern_SUPG(m_pdofs, m_pshift);

        gsMatrix<T> dummy_rhsA_SUPG(m_dofs, 1);
        gsMatrix<T> dummy_rhsN_SUPG(m_dofs, 1);
        gsMatrix<T> dummy_rhsB_SUPG(m_dofs, 1);

        m_blockApattern_SUPG.resize(m_udofs, m_udofs);
        dummy_rhsA_SUPG.setZero();
        dummy_rhsB_SUPG.setZero();

        m_blockApattern_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockAVisitor<T> >(m_blockApattern_SUPG, dummy_rhsA_SUPG, m_nonZerosPerColU);
        m_blockApattern_SUPG.makeCompressed();
        m_blockApattern_SUPG = 0. * m_blockApattern_SUPG;

        blocksBpattern_SUPG.reserve(gsVector<int>::Constant(m_pshift, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSSUPGBlockBVisitor<T> >(blocksBpattern_SUPG, dummy_rhsB_SUPG, m_nonZerosPerColP);
        blocksBpattern_SUPG.makeCompressed();

        #pragma omp parallel for
        for (int s = 0; s < m_tarDim; ++s) {
            m_blockBpattern_SUPG[s] = gsSparseMatrix<T, Eigen::ColMajor>((blocksBpattern_SUPG.middleCols(s * m_udofs, m_udofs)).transpose());//(m_blockB_SUPG_pom[s].transpose());
            m_blockBpattern_SUPG[s].makeCompressed();
            m_blockBpattern_SUPG[s] = 0. * m_blockBpattern_SUPG[s];
        }

        m_blockNpattern_SUPG.resize(m_udofs, m_udofs);
        dummy_rhsN_SUPG.setZero();

        m_blockNpattern_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockNVisitor<T> >(m_blockNpattern_SUPG, dummy_rhsN_SUPG, m_nonZerosPerColU);
        m_blockNpattern_SUPG.makeCompressed();
        m_blockNpattern_SUPG = 0. * m_blockNpattern_SUPG;

        if ((isUnsteady() && isTimeDerTerm()) || isRotation())
            assembleBlockMSUPGpattern();

        m_solution = current_solution;
    }

    //=================== assembleBlockNonlinearCSDpattern =========================
    void assembleBlockNonlinearCSDpattern()
    {
        gsMatrix<T> current_solution = m_solution;

        gsMatrix<T> dummy_rhsN_SUPG(m_dofs, 1);

        m_blockNpattern_SUPG.resize(m_udofs, m_udofs);
        dummy_rhsN_SUPG.setZero();

        m_blockNpattern_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockNVisitor<T> >(m_blockNpattern_SUPG, dummy_rhsN_SUPG, m_nonZerosPerColU);
        m_blockNpattern_SUPG.makeCompressed();
        m_blockNpattern_SUPG = 0. * m_blockNpattern_SUPG;

        if ((isUnsteady() && isTimeDerTerm()) || isRotation())
            assembleBlockMSUPGpattern();

        m_solution = current_solution;
    }

    //=================== assembleBlockMSUPGpattern =========================
    void assembleBlockMSUPGpattern()
    {
        m_blockMpattern_SUPG.resize(m_udofs, m_udofs);

        gsMatrix<T> dummy_rhsM_SUPG(m_dofs, 1);
        dummy_rhsM_SUPG.setZero();

        m_blockMpattern_SUPG.reserve(gsVector<int>::Constant(m_udofs, m_nonZerosPerColU));
        Base::template assembleBlock< uwbINSSUPGBlockMVisitor<T> >(m_blockMpattern_SUPG, dummy_rhsM_SUPG, m_nonZerosPerColU);
        m_blockMpattern_SUPG.makeCompressed();
        m_blockMpattern_SUPG = 0. * m_blockMpattern_SUPG;
    }

    void assembleAllLinearBlocks()
    {
        assembleLinearStokesPart();
        assemblePressurePoisson();
        assembleBlocksC();

        if (!isRotation())
            assembleMassMatrix();
    }

    void fillStokesSystem_into(gsSparseMatrix<T> & stokesMatrix, gsMatrix<T> & stokesRhs) const
    {
        gsVector<int> nonZerosPerColumnVector;
        nonZerosPerColumnVector.setZero(m_dofs);
        for (int s = 0; s < m_tarDim; ++s)
            for (int i = 0; i < m_udofs; i++)
            {
                nonZerosPerColumnVector(i + s * m_udofs) = m_blockA.col(i).nonZeros() + m_blockB[s].col(i).nonZeros();

                if (isRotation())
                    nonZerosPerColumnVector(i + s * m_udofs) += m_blockM.col(i).nonZeros();
            }

        for (int i = 0; i < m_pdofs; i++)
            for (int s = 0; s < m_tarDim; ++s)
                nonZerosPerColumnVector(i + m_pshift) += m_blockMinusBT[s].col(i).nonZeros();

        if (getAssemblerOptions().intStrategy == iFace::dg)
        {
            for (int i = 0; i < m_pdofs; i++)
                nonZerosPerColumnVector(i + m_pshift) += m_blockPdg.col(i).nonZeros();
        }

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

        if (getAssemblerOptions().intStrategy == iFace::dg)
        {
#pragma omp parallel for num_threads(m_numThreads)
            for (index_t col = 0; col < m_pdofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockPdg, col); it; ++it)
                    stokesMatrix.insert(it.row() + m_pshift, it.col() + m_pshift) = it.value();
        }

        if (isRotation())
        {
            GISMO_ENSURE(m_tarDim == 2 || m_tarDim == 3, "Rotation is implemented only for 2D and 3D.");

#pragma omp parallel for num_threads(m_numThreads)
            for (index_t col = 0; col < m_udofs; ++col)
                for (typename gsSparseMatrix<T>::InnerIterator it(m_blockM, col); it; ++it)
                {
                    stokesMatrix.insert(it.row() + ((m_tarDim - 2) * m_udofs), it.col() + ((m_tarDim - 1) * m_udofs)) = -getOmega() * it.value();
                    stokesMatrix.insert(it.row() + ((m_tarDim - 1) * m_udofs), it.col() + ((m_tarDim - 2) * m_udofs)) = getOmega() * it.value();
                }
        }

        stokesMatrix.makeCompressed();

        stokesRhs = m_rhsA + m_rhsB + m_rhsF;

        if (getAssemblerOptions().intStrategy == iFace::dg)
            stokesRhs += m_rhsPdg;

        if (isRotation())
        {
            stokesRhs.middleRows((m_tarDim - 2) * m_udofs, m_udofs) -= getOmega() * m_rhsM.middleRows((m_tarDim - 1) * m_udofs, m_udofs);
            stokesRhs.middleRows((m_tarDim - 1) * m_udofs, m_udofs) += getOmega() * m_rhsM.middleRows((m_tarDim - 2) * m_udofs, m_udofs);
        }
    }

    // ------------------------------ PCD blocks ------------------------------

public:
    void assemblePressurePoisson()
    {
        // matrix and rhs cleaning
        m_rhsAp.setZero();
        m_blockAp.resize(m_pdofs, m_pdofs);
        gsSparseMatrix<T> blockApsym(m_pdofs, m_pdofs);

        // block assembly
        blockApsym.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSBlockApVisitor<T> >(blockApsym, m_rhsAp, m_nonZerosPerColP);
        m_blockAp = blockApsym.template selfadjointView<Eigen::Upper>();
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            gsSparseMatrix<T> blocksApdg(m_pdofs, m_pdofs);
            blocksApdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
            Base::template assembleBlockDg< uwbINSdgBlockApVisitor<T> >(blocksApdg, m_rhsAp, m_nonZerosPerColP);
            m_blockAp += blocksApdg;
        }
        m_blockAp.makeCompressed();
    }

    void assemblePressureConvection()
    {
        // matrix and rhs cleaning
        m_blockNp.resize(m_pdofs, m_pdofs);
        gsMatrix<T> dummyRhsNp(m_dofs, 1);
        dummyRhsNp.setZero();

        // blocks assembly
        m_blockNp.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
        Base::template assembleBlock< uwbINSBlockNpVisitor<T> >(m_blockNp, dummyRhsNp, m_nonZerosPerColP);
        if (getAssemblerOptions().intStrategy == iFace::dg) {
            gsSparseMatrix<T> m_blockNpdg(m_pdofs, m_pdofs);
            m_blockNpdg.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));
            Base::template assembleBlockDg< uwbINSdgBlockNpVisitor<T> >(m_blockNpdg, dummyRhsNp, m_nonZerosPerColP);
            m_blockNp += m_blockNpdg;
        }
        m_blockNp.makeCompressed();

    }

    void assembleApRobinBlock(std::vector<std::pair<int, boxSide> > bndPart)
    {
        // matrix and rhs cleaning
        m_blockRobin.resize(m_pdofs, m_pdofs);
        gsMatrix<T> dummyRhs(m_dofs, 1);
        dummyRhs.setZero();

        m_blockRobin.reserve(gsVector<int>::Constant(m_pdofs, m_nonZerosPerColP));

        uwbINSBlockVisitorRobinPCD<T> visitor(m_dofMappers, m_viscosity);
        this->blockSpecificSettings(&visitor);

        gsMatrix<T> quNodes, physNodes;
        gsVector<T> quWeights; 
        unsigned evFlags(0);
        gsQuadRule<T> QuRule;
        typename gsGeometryEvaluator<T>::uPtr geoEval;
        typename gsBasis<T>::domainIter domIt;
        memory::shared_ptr< gsBasisRefs<T> > basesPtr;

        for (size_t sideIndex = 0; sideIndex < bndPart.size(); sideIndex++) // loop over inflow patch sides
        {
            int patch = bndPart[sideIndex].first;
            boxSide side = bndPart[sideIndex].second;

            basesPtr.reset(new gsBasisRefs<T>(getBases(), patch));

            visitor.initialize(*basesPtr, patch, QuRule, evFlags, side);

            geoEval = typename gsGeometryEvaluator<T>::uPtr(getEvaluator(evFlags, getPatches()[patch]));
            
            domIt = (*basesPtr)[1].makeDomainIterator(side);

            for (; domIt->good(); domIt->next())
            {
                QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
                visitor.evaluate(*basesPtr, *geoEval, quNodes);
                visitor.assemble(*domIt, *geoEval, quWeights);
                visitor.localToGlobal(m_blockRobin, dummyRhs);
            }
        }
    }

    void preparePCDboundary(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall, int bcType = 0)
    {
        findPressureBoundaryIDs(bndIn, bndOut, bndWall);

        // Robin condition at inflow boundary
        if (bcType == 1 || bcType == 2 || bcType == 3 || bcType == 5)
            assembleApRobinBlock(bndIn);

        m_bPCDbndPrepared = true;
    }

    gsSparseMatrix<T> getPressurePoissonMatrix(bool assemb = true, bool lumping = false)
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
            uwbINSPreconditioner<T>::diagInvMatrix_into(m_blockM, velMinv, 1, lumping);

            gsSparseMatrix<T> Ap = m_blockB[0] * velMinv * m_blockB[0].transpose();
            for (size_t i = 1; i < m_blockB.size(); i++)
                Ap += m_blockB[i] * (velMinv * m_blockB[i].transpose());

            return Ap;
        }
    }

    gsSparseMatrix<T> getPressureConvectionMatrix()
    {
        assemblePressureConvection();

        return m_blockNp;
    }

    gsSparseMatrix<T> getPressureMassMatrix()
    {
        if (!m_blockMp.nonZeros())
            assemblePressureMassMatrix();

        return m_blockMp;
    }

    // bcType:
    // 0 - original PCD: Dirichlet inflow + Neumann outflow and walls for both Ap, Fp
    // 1 - Elman, Tuminaro variant 1: Robin inflow + Neumann outflow and walls for Fp, nothing for Ap (Ap should not be assembled, i.e., Ap = B Mu^(-1) Bt)
    // 2 - Elman, Tuminaro variant 2: Robin inflow + Dirichlet outflow + Neumann walls for Fp, nothing for Ap (Ap should not be assembled, i.e., Ap = B Mu^(-1) Bt)
    // 3 - Blechta Y-variant: Robin inflow + Dirichlet outflow + Neumann walls for Fp, Dirichlet outflow + Neumann inflow and walls for Ap (Ap should be assembled)
    // 4 - no BCs set (for assembAp = assembFp = false)

    void applyPCDboundaryConditions(gsSparseMatrix<T>& Ap, gsSparseMatrix<T>& Fp, int bcType)
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

protected:
    void findPressureBoundaryIDs(std::vector<std::pair<int, boxSide> > bndIn, std::vector<std::pair<int, boxSide> > bndOut, std::vector<std::pair<int, boxSide> > bndWall)
    {
        findPressureBoundaryPartIDs(bndIn, m_presInIDs);
        findPressureBoundaryPartIDs(bndOut, m_presOutIDs);
        findPressureBoundaryPartIDs(bndWall, m_presWallIDs);
    }

    void findPressureBoundaryPartIDs(std::vector<std::pair<int, boxSide> > bndPart, std::vector<index_t>& idVector)
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

    // scaleDiag - scale the diagonal entries corresponding to the Dirichlet conditions (otherwise equal to 1)
    // exclInflowPart - exclude inflow bnd IDs from calculation of the average for alpha
    void applyPCDboundaryDirichletPart(gsSparseMatrix<T>& block, const std::vector<index_t>& bndPartIDs, bool scaleDiag = true, bool exclInflowPart = false)
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

    void applyPCDboundaryRobinPart(gsSparseMatrix<T>& block)
    {
        GISMO_ASSERT(m_bPCDbndPrepared, "Pressure boundary DOF IDs not prepared!");

        block += m_blockRobin;        
    }

    // ------------------------------ PCD blocks end ------------------------------

public:

    T computeDimensionlessWallDistance(gsMatrix<T> solution, gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        std::ofstream ofile;
        if (print)
        {
            time_t t = std::time(0);   // get time now
            struct tm * now = std::localtime(&t);
            std::stringstream filename;
            filename << "wallDistance_" << (now->tm_mon + 1) << "-" << now->tm_mday << "_"
                << now->tm_hour << "-" << now->tm_min << ".txt";

            ofile.open(filename.str().c_str());

            ofile << "\n             Dimensionless wall distance information                       \n";
            ofile << "=========================================================================================================\n";
            ofile << "patch index  |     side      |   number of unfit points (y+ > " << maxYplus << ") |   min y+   |   max y+   |   average y+   \n";
            ofile << "---------------------------------------------------------------------------------------------------------\n";
        }

        gsVector<T> minYPlusOverSides;
        minYPlusOverSides.setZero(distancePatches.rows());

        T skinFrCoeff = math::pow(2.0*math::log10(reynoldsNumber) - 0.65, -2.3);
        T tau_w = 0.5 * skinFrCoeff * math::pow(uFreeStream, 2);
        T u_tau = math::sqrt(tau_w);

        int unfitPointsNum;
        for (int numSide = 0; numSide < distancePatches.rows(); numSide++)
        {
            int patchIndex = distancePatches[numSide];
            boxSide side = distanceSides[numSide];

            typename gsGeometry<T>::uPtr geomBoundary = getPatches().patch(patchIndex).boundary(side);

            gsMatrix<T> ab = geomBoundary->support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);

            gsVector<unsigned> np = uniformSampleCount(a, b, npts);
            gsMatrix<T> paramBoundaryNodes = gsPointGrid(a, b, np);
            int numOfPoints = paramBoundaryNodes.cols();

            gsMatrix<T> paramBoundaryNodesFull;
            paramBoundaryNodesFull.setZero(getPatches().dim(), paramBoundaryNodes.cols());
            gsMatrix<T> paramFirstNodesFull = paramBoundaryNodesFull;

            int param = side.parameter();
            int index = 0;

            if (m_tarDim == 2)
            {
                const gsTensorBSplineBasis<2, T>*  basis = dynamic_cast<const gsTensorBSplineBasis<2, T>*>(&getBases().at(0).piece(patchIndex));

                for (int s = 0; s < getPatches().dim(); s++)
                {
                    if (s == side.direction())
                    {
                        if (param == 1) // east or north or back
                        {
                            paramBoundaryNodesFull.middleRows(s, 1).setOnes();

                            real_t nodeValue = basis->knot(s, basis->knots(s).size() - basis->degree(s) - 2);
                            //basis->knot(0, basis->knots(0).size() - basis->degree(0) - 2)//east
                            //basis->knot(1, basis->knots(1).size() - basis->degree(1) - 2)//north
                            //basis->knot(2, basis->knots(2).size() - basis->degree(2) - 2)//back

                            paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                        }
                        else // param = 0 ... i.e. west or south or front
                        {
                            real_t nodeValue = basis->knot(s, basis->degree(s) + 1);
                            //basis->knot(0, basis->degree(0) + 1)//west
                            //basis->knot(1, basis->degree(1) + 1)//south
                            //basis->knot(2, basis->degree(2) + 1)//front

                            paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                        }

                    }
                    else
                    {
                        paramBoundaryNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                        paramFirstNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                        index++;
                    }
                }

            }
            else {
                const gsTensorBSplineBasis<3, T>*  basis = dynamic_cast<const gsTensorBSplineBasis<3, T>*>(&getBases().at(0).piece(patchIndex));

                for (int s = 0; s < getPatches().dim(); s++)
                {
                    if (s == side.direction())
                    {
                        if (param == 1) // east or north or back
                        {
                            paramBoundaryNodesFull.middleRows(s, 1).setOnes();
                            real_t nodeValue = basis->knot(s, basis->knots(s).size() - basis->degree(s) - 2);
                            paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                        }
                        else // param = 0 ... i.e. west or south or front
                        {
                            real_t nodeValue = basis->knot(s, basis->degree(s) + 1);
                            paramFirstNodesFull.middleRows(s, 1).setConstant(nodeValue);
                        }

                    }
                    else
                    {
                        paramBoundaryNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                        paramFirstNodesFull.middleRows(s, 1) = paramBoundaryNodes.middleRows(index, 1);
                        index++;
                    }
                }
            }

            gsMatrix<T> physBoundaryNodesFull;
            gsMatrix<T> physFirstNodesFull;
            getPatches().patch(patchIndex).eval_into(paramBoundaryNodesFull, physBoundaryNodesFull); // physical coordinates of paramBoundaryNodes
            getPatches().patch(patchIndex).eval_into(paramFirstNodesFull, physFirstNodesFull); // physical coordinates of paramBoundaryNodes

            T averageValue = 0.;
            unfitPointsNum = 0;
            gsVector<T> yPlusAtSide;
            yPlusAtSide.setZero(numOfPoints);
            for (int nodeIndex = 0; nodeIndex < numOfPoints; nodeIndex++)
            {
                T distanceToWall = (physFirstNodesFull.col(nodeIndex) - physBoundaryNodesFull.col(nodeIndex)).norm();
                yPlusAtSide[nodeIndex] = distanceToWall * u_tau / viscosity;

                averageValue += yPlusAtSide[nodeIndex];
                if (yPlusAtSide[nodeIndex] > maxYplus)
                    unfitPointsNum++;
            }
            if (print)
                ofile << "      " << patchIndex << "      |   " << side << "    |                 " << unfitPointsNum << "                 |   " << yPlusAtSide.minCoeff() <<
                "  |   " << yPlusAtSide.maxCoeff() << "  |   " << averageValue / numOfPoints << "\n";

            minYPlusOverSides[numSide] = yPlusAtSide.minCoeff();
        }
        //ofile << "=========================================================================================================\n\n";
        if (print)
            ofile.close();

        return minYPlusOverSides.minCoeff() * viscosity / u_tau;
    }

protected:
    void blockCheck(const gsSparseMatrix<T> & block) const
    {
        GISMO_ASSERT(block.nonZeros(), "The required block is empty!");
    }

//==============================================================================================================================

public:
    void evaluateINSsteadyResiduumL2norm(const gsMatrix<T> & solVector, std::vector<T> & residuum)
    {
        m_bResiduumEval = true;
        updateCurrentSolField(solVector, false);
        Base::template assembleProblemResiduum< uwbINSsteadyResiduumEvaluatorL2norm<T> >(residuum);
    }

protected:
    virtual void residuumSpecificSettings(uwbEvaluatorBase<T>* residuumVisitor)
    {
        std::vector<gsField<T> > sols;
        sols.push_back(m_currentVelField);
        sols.push_back(m_currentPresField);
        residuumVisitor->setCurrentSolution(sols);

        //uwbINSsteadyResiduumEvaluatorL2norm<T>* pVisitor = dynamic_cast<uwbINSsteadyResiduumEvaluatorL2norm<T>*>(residuumVisitor);
        //if (pVisitor != NULL)
        //{
        //   if (isUnsteady())
        //    {
        //        gsField<T> oldSolField = constructSolution(m_oldTimeSol, 0); // TODO: improve, prepare the gsField into member
        //        pVisitor->setOldSolutionField(oldSolField);
        //    }
        //}


        //if (pVisitor != NULL)
        //    pVisitor->setTurbulenceSolver(m_params.settings().getTurbulenceSolver());
    }

//===================================================================================================================================

public:

    int getUdofs() const { return m_udofs; }
    int getPdofs() const { return m_pdofs; }
    int getPshift() const { return m_pshift; }

    const gsSparseMatrix<T> & getBlockA() const { return m_blockA; }
    const std::vector< gsSparseMatrix<T> > & getBlockB() const { return m_blockB; }
    const std::vector< gsSparseMatrix<T> > & getBlockMinusBT() const { return m_blockMinusBT; }
    const gsSparseMatrix<T> & getBlockB(const int i) const { return m_blockB[i]; }
    const gsSparseMatrix<T> & getBlockMinusBT(const int i) const { return m_blockMinusBT[i]; }
    const gsSparseMatrix<T> & getBlockCT(const int i) const { return m_blockCT[i]; }
    const gsSparseMatrix<T> & getBlockN() const { return m_blockN; }
    const gsSparseMatrix<T> & getBlockNpattern() const { return m_blockNpattern; }
    const gsSparseMatrix<T> & getBlockM() const { return m_blockM; }
    const gsSparseMatrix<T> & getBlockAp() const { blockCheck(m_blockAp); return m_blockAp; }
    const gsSparseMatrix<T> & getBlockNp() const { blockCheck(m_blockNp); return m_blockNp; }
    const gsSparseMatrix<T> & getBlockMp() const { blockCheck(m_blockMp); return m_blockMp; }

    const gsMatrix<T> & getRhsA() const { return m_rhsA; }
    const gsMatrix<T> & getRhsB() const { return m_rhsB; }
    const gsMatrix<T> & getRhsC() const { return m_rhsC; }
    const gsMatrix<T> & getRhsN() const { return m_rhsN; }
    const gsMatrix<T> & getRhsM() const { return m_rhsM; }
    const gsMatrix<T> & getRhsAp() const { return m_rhsAp; }
    const gsMatrix<T> & getRhsMp() const { return m_rhsMp; }
    const gsMatrix<T> & getRhsF() const { return m_rhsF; }

    //---------------- SUPG --------------------
    const std::vector< gsSparseMatrix<T> > & getBlockB_SUPG() const { return m_blockB_SUPG; }
    const gsSparseMatrix<T> & getBlockB_SUPG(const int i) const { return m_blockB_SUPG[i]; }
    const gsSparseMatrix<T> & getBlockA_SUPG() const { return m_blockA_SUPG; }
    const gsSparseMatrix<T> & getBlockN_SUPG() const { return m_blockN_SUPG; }
    const gsSparseMatrix<T> & getBlockM_SUPG() const { return m_blockM_SUPG; }
    const gsMatrix<T> & getRhsA_SUPG() const { return m_rhsA_SUPG; }
    const gsMatrix<T> & getRhsN_SUPG() const { return m_rhsN_SUPG; }
    const gsMatrix<T> & getRhsM_SUPG() const { return m_rhsM_SUPG; }
    const gsMatrix<T> & getRhsB_SUPG() const { return m_rhsB_SUPG; }
    const gsSparseMatrix<T> & getBlockApattern_SUPG() const { return m_blockApattern_SUPG; }
    const gsSparseMatrix<T> & getBlockNpattern_SUPG() const { return m_blockNpattern_SUPG; }
    const gsSparseMatrix<T> & getBlockMpattern_SUPG() const { return m_blockMpattern_SUPG; }
    const std::vector< gsSparseMatrix<T> > & getBlockBpattern_SUPG() const { return m_blockBpattern_SUPG; }
    const gsSparseMatrix<T> & getBlockBpattern_SUPG(const int i) const { return m_blockBpattern_SUPG[i]; }
    //----------------- CROSSWIND ------------------
    const std::vector< gsSparseMatrix<T> > & getBlockCrosswindPattern() const { return m_blockCrosswindPattern; }
    const gsSparseMatrix<T> & getBlockCrosswindPattern(const int i) const { return m_blockCrosswindPattern[i]; }
    const std::vector< gsSparseMatrix<T> > & getBlockCrosswind() const { return m_blockCrosswind; }
    const gsSparseMatrix<T> & getBlockCrosswind(const int i) const { return m_blockCrosswind[i]; }
    const gsMatrix<T> & getRhsCrosswind() const { return m_rhsCrosswind; }
    //----------------- tanh-CSD ------------------
    const std::vector< gsSparseMatrix<T> > & getBlockTanhCSDPattern() const { return m_blockTanhCSDpattern; }
    const gsSparseMatrix<T> & getBlockTanhCSDpattern(const int i) const { return m_blockTanhCSDpattern[i]; }
    const std::vector< gsSparseMatrix<T> > & getBlockTanhCSD() const { return m_blockTanhCSD; }
    const gsSparseMatrix<T> & getBlockTanhCSD(const int i) const { return m_blockTanhCSD[i]; }
    const gsMatrix<T> & getRhsTanhCSD() const { return m_rhsTanhCSD; }
    //----------------- iso-AD ------------------
    const std::vector< gsSparseMatrix<T> > & getBlockIsoADpattern() const { return m_blockIsoADpattern; }
    const gsSparseMatrix<T> & getBlockIsoADpattern(const int i) const { return m_blockIsoADpattern[i]; }
    const std::vector< gsSparseMatrix<T> > & getBlockIsoAD() const { return m_blockIsoAD; }
    const gsSparseMatrix<T> & getBlockIsoAD(const int i) const { return m_blockIsoAD[i]; }
    const gsMatrix<T> & getRhsIsoAD() const { return m_rhsIsoAD; }
    //----------------- PSPG -----------------------
    const gsMatrix<T> & getRhsPSPG() const { return m_rhsPSPG; }
    //----------------- RANS -----------------------
    const gsMatrix<T> & getRhsA_RANS() const { return m_rhsA_RANS; }
    const gsMatrix<T> & getRhsEdiag_RANS() const { return m_rhsEdiag_RANS; }
    const gsMatrix<T> & getRhsEnondiag_RANS() const { return m_rhsEnondiag_RANS; }
    const gsSparseMatrix<T> & getBlockA_RANS() const { return m_blockA_RANS; }
    const std::vector< gsSparseMatrix<T> > & getBlockEdiag_RANS() const { return m_blockEdiag_RANS; }
    const gsSparseMatrix<T> & getBlockEdiag_RANS(const int i) const { return m_blockEdiag_RANS[i]; }

    const gsSparseMatrix<T> & getBlockA_RANS_SUPG() const { return m_blockA_RANS_SUPG; }
    const gsMatrix<T> & getRhsA_RANS_SUPG() const { return m_rhsA_RANS_SUPG; }

    const gsSparseMatrix<T> & getBlock_RANS_AD() const { return m_block_RANS_AD; }
    const gsMatrix<T> & getRhs_RANS_AD() const { return m_rhs_RANS_AD; }

    bool    isRotation() const { return m_params.settings().isRotation(); }
    bool    isTimeDerTerm() const { return m_params.settings().get(constantsINS::timeDerTerm); }
    T       getOmega() const { return this->getParam(constantsINS::omega); }

    bool    isSUPG() const { return m_params.settings().get(constantsINS::SUPG); }
    bool    isTCSD() const { return m_params.settings().get(constantsINS::TCSD); }
    bool    isCROSSWIND() const { return m_params.settings().get(constantsINS::CROSSWIND); }
    bool    isPSPG() const { return m_params.settings().get(constantsINS::PSPG); }
    bool    isRANScrosswind() const { return m_params.settings().get(constantsINS::RANScrosswind); }
    bool    isRANSad() const { return m_params.settings().get(constantsINS::RANSad); }
    bool    isRANStanhCSD() const { return m_params.settings().get(constantsINS::tanhCSD); }
    bool    isRANSisoAD() const { return m_params.settings().get(constantsINS::RANSisoAD); }

    int     getTauStabTypeSUPG() const { return m_params.settings().get(constantsINS::tauStabTypeSUPG); }
    int     getTauStabTypeCW() const { return m_params.settings().get(constantsINS::tauStabTypeCW); }
    int     getTauStabTypeAD() const { return m_params.settings().get(constantsINS::tauStabTypeAD); }
    int     getTauStabTypePSPG() const { return m_params.settings().get(constantsINS::tauStabTypePSPG); }
    int     getDiamPowerSUPG() const { return m_params.settings().get(constantsINS::diamPowerSUPG); }
    int     getDiamPowerCW() const { return m_params.settings().get(constantsINS::diamPowerCW); }
    int     getDiamPowerAD() const { return m_params.settings().get(constantsINS::diamPowerAD); }
    int     getCrosswindType() const { return m_params.settings().get(constantsINS::crosswindType); }
    int     getCWresidualType() const { return m_params.settings().get(constantsINS::CWresidualType); }
    int     getIsoADtype() const { return m_params.settings().get(constantsINS::isoADtype); }
    int     getIsoADresidualType() const { return m_params.settings().get(constantsINS::isoADresidualType); }
    T       getIsoADalpha() const { return m_params.settings().get(constantsINS::isoADalpha); }
    T       getTimeStep() const { return m_params.settings().get(constantsINS::timeStep); }

    bool isStabilization()
    {
        return(isSUPG() || isTCSD() || isCROSSWIND() || isRANScrosswind() || isRANSad() || isRANStanhCSD() || isPSPG());
    }

    bool isResiduumEval() const { return m_bResiduumEval; }

    void setRANS()
    {
        initRANSMembers();
        if (isSUPG() || isTCSD())
            m_bRANSsupg = true;
        if (isPSPG())
            m_bRANSpspg = true;
    }

protected:

    int m_udofs;
    int m_pdofs;
    int m_pshift;
    int m_nonZerosPerColU, m_nonZerosPerColP;

    bool m_bRANSsupg, m_bRANSpspg;
    gsSparseMatrix<T> m_blockA;
    std::vector< gsSparseMatrix<T> > m_blockB;
    std::vector< gsSparseMatrix<T> > m_blockMinusBT;
    std::vector< gsSparseMatrix<T> > m_blockCT;
    gsSparseMatrix<T> m_blockN;
    gsSparseMatrix<T> m_blockNpattern;
    gsSparseMatrix<T> m_blockM;
    gsSparseMatrix<T> m_blockMinv;
    gsSparseMatrix<T> m_blockAp;
    gsSparseMatrix<T> m_blockNp;
    gsSparseMatrix<T> m_blockMp;
    gsSparseMatrix<T> m_blockPdg;
    gsSparseMatrix<T> m_blockRobin;
    gsMatrix<T> m_rhsA;
    gsMatrix<T> m_rhsB;
    gsMatrix<T> m_rhsC;
    gsMatrix<T> m_rhsN;
    gsMatrix<T> m_rhsM;
    gsMatrix<T> m_rhsAp;
    gsMatrix<T> m_rhsMp;
    gsMatrix<T> m_rhsPdg;
    gsMatrix<T> m_rhsF;

    gsField<T>  m_currentVelField, m_currentPresField, m_oldTimeVelField;

    // rotation members
    gsMatrix<T> m_omegaXrCoeffs;

    //SUPG members
    gsSparseMatrix<T> m_blockA_SUPG;
    std::vector< gsSparseMatrix<T> > m_blockB_SUPG;
    gsSparseMatrix<T> m_blockN_SUPG;
    gsSparseMatrix<T> m_blockM_SUPG;
    gsMatrix<T> m_rhsA_SUPG;
    gsMatrix<T> m_rhsB_SUPG;
    gsMatrix<T> m_rhsN_SUPG;
    gsMatrix<T> m_rhsM_SUPG;
    gsSparseMatrix<T> m_blockApattern_SUPG;
    gsSparseMatrix<T> m_blockNpattern_SUPG;
    gsSparseMatrix<T> m_blockMpattern_SUPG;
    std::vector< gsSparseMatrix<T> > m_blockBpattern_SUPG;

    //Crosswind members
    std::vector< gsSparseMatrix<T> > m_blockCrosswind;
    std::vector< gsSparseMatrix<T> > m_blockCrosswindPattern;
    gsMatrix<T> m_rhsCrosswind;

    //tanh-CSD members
    std::vector< gsSparseMatrix<T> > m_blockTanhCSD;
    std::vector< gsSparseMatrix<T> > m_blockTanhCSDpattern;
    gsMatrix<T> m_rhsTanhCSD;

    //iso-AD members
    std::vector< gsSparseMatrix<T> > m_blockIsoAD;
    std::vector< gsSparseMatrix<T> > m_blockIsoADpattern;
    gsMatrix<T> m_rhsIsoAD;

    //PSPG members
    gsMatrix<T> m_rhsPSPG;

    //RANS members
    gsMatrix<T> m_rhsA_RANS;
    gsMatrix<T> m_rhsEdiag_RANS;
    gsMatrix<T> m_rhsEnondiag_RANS;
    gsSparseMatrix<T> m_blockA_RANS;
    std::vector< gsSparseMatrix<T> > m_blockEdiag_RANS;
    gsSparseMatrix<T> m_blockApattern_RANS;
    std::vector< gsSparseMatrix<T> > m_blockEdiagpattern_RANS;

    //RANS SUPG members
    gsSparseMatrix<T> m_blockA_RANS_SUPG;
    gsMatrix<T> m_rhsA_RANS_SUPG;

    //RANS AD members
    gsSparseMatrix<T> m_block_RANS_AD;
    gsMatrix<T> m_rhs_RANS_AD;

    // PCD assembly
    std::vector<index_t> m_presInIDs, m_presOutIDs, m_presWallIDs;
    bool m_bPCDbndPrepared;

    bool m_bResiduumEval;

    // members from uwbBlockAssemblerBase
    using Base::m_dofs;
    using Base::m_numThreads;
    using Base::m_tarDim;
    using Base::m_bUnsteady;
    using Base::m_params;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;
    using Base::m_viscosity;

    // member functions from uwbBlockAssemblerBase
public:
    using Base::isUnsteady;
    using Base::getPatches;
    using Base::getBases;
    using Base::getRhsFcn;
    using Base::getAssemblerOptions;

    using Base::getTauStabTypeTanhCSD;
    using Base::getTanhCSDtype;
    using Base::getTanhCSDresidualType;
    using Base::getTanhCSDalpha;

protected:
    using Base::computeDirichletDofs;
};

} // namespace gismo
