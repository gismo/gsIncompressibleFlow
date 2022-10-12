/** @file uwbTMBlockAssembler.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#ifdef _OPENMP 
#include <omp.h>
#endif

#include "uwbBlockAssemblerBase.h"
#include "uwbTMBlockVisitorsKOmega.h"
#include "uwbTMSUPGBlockVisitorsKOmega.h"
#include "uwbTMCrosswindVisitors.h"
#include "uwbTMADVisitors.h"
#include "uwbTMFCTBlockVisitorsKOmega.h"
#include "uwbTMCSDVisitor.h"

namespace gismo
{

template<class T>
class uwbTMBlockAssembler : public uwbBlockAssemblerBase<T>
{
public:
    typedef uwbBlockAssemblerBase<T> Base;

public:

    uwbTMBlockAssembler(const uwbINSSolverParams<T>& params, int numVar) :
        Base(params), m_numVar(numVar)
    {
        initMembers();
    }

    virtual ~uwbTMBlockAssembler()
    {}

protected:

    virtual void initMembers()
    {
        m_varDofs = getMapper().freeSize();
        m_dofs = m_numVar * m_varDofs;

        m_nonZerosPerCol = 1;
        for (int i = 0; i < m_tarDim; i++)
            m_nonZerosPerCol *= 2 * getSolBasis().maxDegree(i) + 1;

        m_ddof.resize(m_numVar);

        if (getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            for (size_t i = 0; i < m_ddof.size(); i++)
            {
                m_ddof[i].setZero(getMapper().boundarySize(), 1);
                computeDirichletDofs(i, 1, m_ddof[i]);
            }
        }
            
        m_solution.setZero(m_varDofs, m_numVar);

        m_currentSolField = constructSolution(m_solution);
        m_oldTimeSolField = m_currentSolField;

        initMembers_kOmega();
    }

    void initMembers_kOmega()
    {
        m_blockM.resize(m_varDofs, m_varDofs);
        m_blockN.resize(m_varDofs, m_varDofs);
        m_blockNpattern.resize(m_varDofs, m_varDofs);
        m_blocksAnonlin.resize(m_varDofs, m_dofs);
        m_blocksAnonlinPattern.resize(m_varDofs, m_dofs);
        m_blocksReaction.resize(m_varDofs, m_dofs);
        m_blocksReactionPattern.resize(m_varDofs, m_dofs);
        m_blockBlend.resize(m_varDofs, m_varDofs);

        m_rhsM.setZero(m_varDofs, m_numVar);
        m_rhsN.setZero(m_varDofs, m_numVar);
        m_rhsAnonlin.setZero(m_varDofs, m_numVar);
        m_rhsReaction.setZero(m_varDofs, m_numVar);
        m_rhsBlend.setZero(m_varDofs, 1);
        m_rhsF.setZero(m_varDofs, m_numVar);

        m_blocks.resize(m_varDofs, m_dofs);
        m_blocksPattern.resize(m_varDofs, m_dofs);
        m_rhs.setZero(m_varDofs, m_numVar);

        if(isTMsupg())
            initSUPGmembers_kOmega();
        if(isTMcrosswind())
        {
            initCROSSWINDmembers_kOmega();
            if(isTMsupg())
                initSUPGmembers_kOmega();
        }
        if(isTMad())
            initADmembers_kOmega();
        if(isTMisoAD())
            initIsoADmembers_kOmega();
        if(isTMtanhCSD())
            initTanhCSDmembers_kOmega();
    }

    void initSUPGmembers_kOmega()
    {
        m_blockM_SUPG.resize(m_varDofs, m_dofs);
        m_blockN_SUPG.resize(m_varDofs, m_dofs);
        m_blockMpattern_SUPG.resize(m_varDofs, m_dofs);
        m_blockNpattern_SUPG.resize(m_varDofs, m_dofs);
        m_blocksAnonlin_SUPG.resize(m_varDofs, m_dofs);
        m_blocksAnonlinPattern_SUPG.resize(m_varDofs, m_dofs);
        m_blocksReaction_SUPG.resize(m_varDofs, m_dofs);
        m_blocksReactionPattern_SUPG.resize(m_varDofs, m_dofs);
        m_blockBlend_SUPG.resize(m_varDofs, m_varDofs);

        m_rhsM_SUPG.setZero(m_varDofs, m_numVar);
        m_rhsN_SUPG.setZero(m_varDofs, m_numVar);
        m_rhsAnonlin_SUPG.setZero(m_varDofs, m_numVar);
        m_rhsReaction_SUPG.setZero(m_varDofs, m_numVar);
        m_rhsBlend_SUPG.setZero(m_varDofs, 1);
        m_rhsF_SUPG.setZero(m_varDofs, m_numVar);

        m_blocks_SUPG.resize(m_varDofs, m_dofs);
        m_blocksPattern_SUPG.resize(m_varDofs, m_dofs);
        m_rhs_SUPG.setZero(m_varDofs, m_numVar);
    }

    void initCROSSWINDmembers_kOmega()
    {
        m_blockCrosswind.resize(m_varDofs, m_dofs);
        m_blockCrosswindPattern.resize(m_varDofs, m_dofs);
        m_rhsCrosswind.setZero(m_varDofs, m_numVar);
    }

    void initADmembers_kOmega()
    {
        m_blockAD.resize(m_varDofs, m_dofs);
        m_blockADpattern.resize(m_varDofs, m_dofs);
        m_rhsAD.setZero(m_varDofs, m_numVar);
    }

    void initIsoADmembers_kOmega()
    {
        m_blockIsoAD.resize(m_varDofs, m_dofs);
        m_blockIsoADpattern.resize(m_varDofs, m_dofs);
        m_rhsIsoAD.setZero(m_varDofs, m_numVar);
    }

    void initTanhCSDmembers_kOmega()
    {
        m_blockTanhCSD.resize(m_varDofs, m_dofs);
        m_blockTanhCSDpattern.resize(m_varDofs, m_dofs);
        m_rhsTanhCSD.setZero(m_varDofs, m_numVar);
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
            const int sz = getSolBasis().piece(p).size();
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

            result.addPatch(getSolBasis().piece(p).makeGeometry(coeffs));
        }
    }

    gsField<T> constructSolution(const gsMatrix<T>& solVector) const
    {
        gsMultiPatch<T> * result = new gsMultiPatch<T>;
        constructSolution(solVector, *result);
        return gsField<T>(getPatches(), typename gsFunctionSet<T>::Ptr(result), true);
    }

    void updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
    { 
        m_currentSolField = constructSolution(solVector);

        if (updateSol)
        {
            m_solution = solVector;

            if (isUnsteady())
                m_oldTimeSolField = m_currentSolField;
        }
    }

    void updateVelocitySolution(const gsField<T>& uSolField)
    { m_uSolField = uSolField; }

    void setSolution(const gsMatrix<T> & solVector)
    { m_solution = solVector; }

    void setPoissonSolution(const gsField<T>& solPoisson, const gsMultiPatch<T>& patches, const gsMultiBasis<T>& bases)
    { m_solPoisson = solPoisson; m_patchesPoisson = patches; m_basesPoisson = bases; }

    void setKOmegaDiffusionCoeffSolField(const std::vector<gsField<T> >& kOmegaDiffCoeffSolField)
    { m_kOmegaDiffCoeffSolField = kOmegaDiffCoeffSolField; }

protected:
    virtual void blockSpecificSettings(uwbVisitorBase<T>* blockVisitor)
    {
        std::vector<gsField<T> > sols;
        sols.push_back(m_uSolField);
        sols.push_back(m_currentSolField);
        blockVisitor->setCurrentSolution(sols);

        uwbTMBlockVisitor<T>* pVisitor = dynamic_cast<uwbTMBlockVisitor<T>*>(blockVisitor);

        if (pVisitor != NULL)
        {
            pVisitor->createTMEvaluator(getTMEvaluator());
            if (checkWallDistanceBasedTM())
                pVisitor->setPoissonSolution(m_solPoisson, m_patchesPoisson, m_basesPoisson);
            pVisitor->setTMProductionLimiter(isTMproductionLimited(), getStartXPointProduction());
        }

        if (isTMsupg())
        {
            uwbTMSUPGBlockVisitorKOmega<T>* pVisitorSUPG = dynamic_cast<uwbTMSUPGBlockVisitorKOmega<T>*>(blockVisitor);
            if (pVisitorSUPG != NULL)
            {
                //pVisitorSUPG->setKOmegaDiffusionCoeffSolField(m_kOmegaDiffCoeffSolField);
                pVisitorSUPG->setTauStabType(getTauStabTypeSUPG(), isUnsteady(), getTimeStep(), getDiamPowerSUPG());
            }
        }

        if (isTMcrosswind())
        {
            uwbTMCrosswindVisitorKOmega<T>* pVisitorCrosswind = dynamic_cast<uwbTMCrosswindVisitorKOmega<T>*>(blockVisitor);
            if (pVisitorCrosswind != NULL)
            {
                pVisitorCrosswind->setCrosswind(getCrosswindType(), getCWresidualType(), getTauStabTypeCW(), isUnsteady(), getTimeStep(), getDiamPowerCW());
                if (isUnsteady())
                    pVisitorCrosswind->setOldSolutionField(m_oldTimeSolField);
            }
        }

        if (isTMad())
        {
            uwbTMADVisitorKOmega<T>* pVisitorAD = dynamic_cast<uwbTMADVisitorKOmega<T>*>(blockVisitor);
            if (pVisitorAD != NULL)
                pVisitorAD->setArtificialDiffusion(getTauStabTypeAD(), getTimeStep(), getDiamPowerAD());
        }

        if (isTMisoAD())
        {
            uwbTMisoADVisitorKOmega<T>* pVisitorIsoAD = dynamic_cast<uwbTMisoADVisitorKOmega<T>*>(blockVisitor);
            if (pVisitorIsoAD != NULL)
            {
                pVisitorIsoAD->setIsoArtificialDiffusion(isUnsteady(), getTimeStep(), getTauStabTypeSUPG(), getIsoADtype());
                if (isUnsteady())
                    pVisitorIsoAD->setOldSolutionField(m_oldTimeSolField);
            }
        }

        if (isTMtanhCSD())
        {
            uwbTMCSDVisitorKOmega<T>* pVisitorCSD = dynamic_cast<uwbTMCSDVisitorKOmega<T>*>(blockVisitor);
            if (pVisitorCSD != NULL)
            {
                pVisitorCSD->setTanhCSD(isUnsteady(), getTimeStep(), getTauStabTypeTanhCSD(), getTanhCSDtype(), getTanhCSDresidualType(), getTanhCSDalpha());
                if (isUnsteady())
                    pVisitorCSD->setOldSolutionField(m_oldTimeSolField);
            }
        }
    }

    virtual void rhsSpecificSettings(uwbVisitorBase<T>* rhsVisitor)
    { 
        blockSpecificSettings(rhsVisitor);
    }

public:
    void assembleLinearPart_kOmega()
    {
        // block M
        // matrix and rhs cleaning
        m_blockM.resize(m_varDofs, m_varDofs);
        m_rhsM.setZero();
        
        m_blockM.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
        if (isFCTlowOrder())
            Base::template assembleBlock< uwbTMFCTBlockMVisitorKOmega<T> >(m_blockM, m_rhsM, m_nonZerosPerCol);
        else
            Base::template assembleBlock< uwbTMBlockMVisitorKOmega<T> >(m_blockM, m_rhsM, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockM.makeCompressed();
    }

    void assemblePatternBlocks_kOmega(bool assembleAll = false)
    {
        if (assembleAll)
            assembleAllPatternBlocks_kOmega();
        else
        {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        // block nonlinear A pattern
        m_blocksAnonlinPattern.resize(m_varDofs, m_dofs);
        m_blocksAnonlinPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        Base::template assembleBlock< uwbTMBlockANonlinearVisitorKOmega<T> >(m_blocksAnonlinPattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksAnonlinPattern.makeCompressed();
        m_blocksAnonlinPattern = 0. * m_blocksAnonlinPattern;

        if (isUnsteady())
        {
            // block N pattern 
            // matrix cleaning
            m_blockNpattern.resize(m_varDofs, m_varDofs);

            m_blockNpattern.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
            if (isFCTlowOrder())
                Base::template assembleBlock< uwbTMFCTBlockNVisitorKOmega<T> >(m_blockNpattern, dummyRhs, m_nonZerosPerCol);
            else
                Base::template assembleBlock< uwbTMBlocksNVisitorKOmega<T> >(m_blockNpattern, dummyRhs, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blockNpattern.makeCompressed();
            m_blockNpattern = 0. * m_blockNpattern;

            // reaction blocks pattern
            m_blocksReactionPattern.resize(m_varDofs, m_dofs);
            m_blocksReactionPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            Base::template assembleBlock< uwbTMBlocksReactionVisitorKOmega<T> >(m_blocksReactionPattern, dummyRhs, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blocksReactionPattern.makeCompressed();
            m_blocksReactionPattern = 0. * m_blocksReactionPattern;
        }
        }
    }

    void assembleAllPatternBlocks_kOmega()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        // block nonlinear A pattern + reaction pattern + blend pattern
        m_blocksPattern.resize(m_varDofs, m_dofs);
        m_blocksPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        Base::template assembleBlock< uwbTMBlockVisitorKOmega<T> >(m_blocksPattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksPattern.makeCompressed();
        m_blocksPattern = 0. * m_blocksPattern;
    }

    void assembleNonlinearPart_kOmega(bool assembleAll = false)
    {
        if (assembleAll)
            assembleAllNonlinearBlocks_kOmega();
        else
        {

        // block nonlinear A
        // matrix and rhs cleaning
        m_blocksAnonlin.resize(m_varDofs, m_dofs);
        m_rhsAnonlin.setZero();

        m_blocksAnonlin.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        Base::template assembleBlock< uwbTMBlockANonlinearVisitorKOmega<T> >(m_blocksAnonlin, m_rhsAnonlin, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksAnonlin.makeCompressed();

        if (isUnsteady())
        {
            // reaction blocks
            // matrix and rhs cleaning
            m_blocksReaction.resize(m_varDofs, m_dofs);
            m_rhsReaction.setZero();

            m_blocksReaction.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            Base::template assembleBlock< uwbTMBlocksReactionVisitorKOmega<T> >(m_blocksReaction, m_rhsReaction, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blocksReaction.makeCompressed();

            // blend block
            // matrix and rhs cleaning
            m_blockBlend.resize(m_varDofs, m_varDofs);
            m_rhsBlend.setZero();

            m_blockBlend.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
            Base::template assembleBlock< uwbTMBlockBlendVisitorKOmega<T> >(m_blockBlend, m_rhsBlend, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blockBlend.makeCompressed();
        }

        // rhs
        m_rhsF.setZero();
        Base::template assembleRhs< uwbTMRhsVisitorKOmega<T> >(m_rhsF);

        }
    }

    void assembleAllNonlinearBlocks_kOmega()
    {

        // block nonlinear A + reaction blocks + blend block + rhs
        // matrix and rhs cleaning
        m_blocks.resize(m_varDofs, m_dofs);
        m_rhs.setZero();

        m_blocks.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        Base::template assembleBlock< uwbTMBlockVisitorKOmega<T> >(m_blocks, m_rhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocks.makeCompressed();
    }

    void assembleNewTimestepPart_kOmega()
    {
        // block N 
        // matrix cleaning
        m_blockN.resize(m_varDofs, m_varDofs);
        m_rhsN.setZero();

        m_blockN.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
        if (isFCTlowOrder())
            Base::template assembleBlock< uwbTMFCTBlockNVisitorKOmega<T> >(m_blockN, m_rhsN, m_nonZerosPerCol);
        else
            Base::template assembleBlock< uwbTMBlocksNVisitorKOmega<T> >(m_blockN, m_rhsN, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockN.makeCompressed();
    }

    //======================== assemble SUPG pattern part =============================
    /*void assembleSUPGpatternBlocks_kOmega()
    {
        gsMatrix<T> dummyRhs_SUPG(m_varDofs, m_numVar);
        dummyRhs_SUPG.setZero();

        // SUPG block nonlinear A pattern
        m_blocksAnonlinPattern_SUPG.resize(m_varDofs, m_dofs);
        m_blocksAnonlinPattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockANonlinearVisitorKOmega<T> >(m_blocksAnonlinPattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksAnonlinPattern_SUPG.makeCompressed();
        m_blocksAnonlinPattern_SUPG = 0. * m_blocksAnonlinPattern_SUPG;

        if (isUnsteady())
        {
            if (isTimeDerTerm())
            {
                // SUPG block M pattern
                // matrix cleaning
                m_blockMpattern_SUPG.resize(m_varDofs, m_dofs);
                m_blockMpattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
                uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockMVisitorKOmega<T> >(m_blockMpattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
                //if (getAssemblerOptions().intStrategy == iFace::dg) { }
                m_blockMpattern_SUPG.makeCompressed();
                m_blockMpattern_SUPG = 0. * m_blockMpattern_SUPG;
            }

            // SUPG block N pattern (in steady case, N is not nonlinear, velocity does not change)
            // matrix cleaning
            m_blockNpattern_SUPG.resize(m_varDofs, m_dofs);
            m_blockNpattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksNVisitorKOmega<T> >(m_blockNpattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blockNpattern_SUPG.makeCompressed();
            m_blockNpattern_SUPG = 0. * m_blockNpattern_SUPG;

            // reaction SUPG blocks pattern
            m_blocksReactionPattern_SUPG.resize(m_varDofs, m_dofs);
            m_blocksReactionPattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksReactionVisitorKOmega<T> >(m_blocksReactionPattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blocksReactionPattern_SUPG.makeCompressed();
            m_blocksReactionPattern_SUPG = 0. * m_blocksReactionPattern_SUPG;
        }
    }*/

    void assembleSUPGpatternBlocks_kOmega()
    {
        gsMatrix<T> dummyRhs_SUPG(m_varDofs, m_numVar);
        dummyRhs_SUPG.setZero();

        // SUPG block nonlinear A pattern + reaction + blend + mass matrix
        m_blocksPattern_SUPG.resize(m_varDofs, m_dofs);
        m_blocksPattern_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksVisitorKOmega<T> >(m_blocksPattern_SUPG, dummyRhs_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksPattern_SUPG.makeCompressed();
        m_blocksPattern_SUPG = 0. * m_blocksPattern_SUPG;
    }

    //======================== assemble SUPG part =============================
    void assembleNonlinearPartSUPG_kOmega()
    {
        // SUPG block nonlinear A
        // matrix and rhs cleaning
        m_blocksAnonlin_SUPG.resize(m_varDofs, m_dofs);
        m_rhsAnonlin_SUPG.setZero();

        m_blocksAnonlin_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockANonlinearVisitorKOmega<T> >(m_blocksAnonlin_SUPG, m_rhsAnonlin_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocksAnonlin_SUPG.makeCompressed();

        if (isUnsteady())
        {
            // reaction blocks
            // matrix and rhs cleaning
            m_blocksReaction_SUPG.resize(m_varDofs, m_dofs);
            m_rhsReaction_SUPG.setZero();

            m_blocksReaction_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksReactionVisitorKOmega<T> >(m_blocksReaction_SUPG, m_rhsReaction_SUPG, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blocksReaction_SUPG.makeCompressed();

            // blend block
            // matrix and rhs cleaning
            m_blockBlend_SUPG.resize(m_varDofs, m_varDofs);
            m_rhsBlend_SUPG.setZero();

            m_blockBlend_SUPG.reserve(gsVector<int>::Constant(m_varDofs, m_nonZerosPerCol));
            uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockBlendVisitorKOmega<T> >(m_blockBlend_SUPG, m_rhsBlend_SUPG, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blockBlend_SUPG.makeCompressed();
        }

        // rhs
        m_rhsF_SUPG.setZero();
        uwbBlockAssemblerBase<T>::template assembleRhs< uwbTMSUPGRhsFVisitorKOmega<T> >(m_rhsF_SUPG);
    }

    void assembleSUPGblocks_kOmega()
    {
        // SUPG block M
        // matrix and rhs cleaning
        m_blockM_SUPG.resize(m_varDofs, m_dofs);
        m_rhsM_SUPG.setZero();

        m_blockM_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockMVisitorKOmega<T> >(m_blockM_SUPG, m_rhsM_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockM_SUPG.makeCompressed();


        // SUPG block nonlinear A + reaction + blend + mass matrix
        // matrix and rhs cleaning
        m_blocks_SUPG.resize(m_varDofs, m_dofs);
        m_rhs_SUPG.setZero();

        m_blocks_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksVisitorKOmega<T> >(m_blocks_SUPG, m_rhs_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blocks_SUPG.makeCompressed();
    }

    //======================== assemble SUPG new time step part =============================
    void assembleNewTimestepPartSUPG_kOmega()
    {
        // SUPG block N
        // matrix cleaning
        m_blockN_SUPG.resize(m_varDofs, m_dofs);
        m_rhsN_SUPG.setZero();

        m_blockN_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlocksNVisitorKOmega<T> >(m_blockN_SUPG, m_rhsN_SUPG, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockN_SUPG.makeCompressed();

        if (isUnsteady() && isTimeDerTerm())
        {
            // SUPG block M
            // matrix and rhs cleaning
            m_blockM_SUPG.resize(m_varDofs, m_dofs);
            m_rhsM_SUPG.setZero();

            m_blockM_SUPG.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
            uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMSUPGBlockMVisitorKOmega<T> >(m_blockM_SUPG, m_rhsM_SUPG, m_nonZerosPerCol);
            //if (getAssemblerOptions().intStrategy == iFace::dg) { }
            m_blockM_SUPG.makeCompressed();
        }
    }

    //======================== assemble CROSSWIND pattern part =============================
    void assembleCrosswindPatternBlock_kOmega()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        // SUPG block nonlinear A pattern
        m_blockCrosswindPattern.resize(m_varDofs, m_dofs);
        m_blockCrosswindPattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMCrosswindVisitorKOmega<T> >(m_blockCrosswindPattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockCrosswindPattern.makeCompressed();
        m_blockCrosswindPattern = 0. * m_blockCrosswindPattern;
    }

    //======================== assemble CROSSWIND part =============================
    void assembleCrosswindBlock_kOmega()
    {
        m_blockCrosswind.resize(m_varDofs, m_dofs);
        m_rhsCrosswind.setZero();

        m_blockCrosswind.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMCrosswindVisitorKOmega<T> >(m_blockCrosswind, m_rhsCrosswind, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockCrosswind.makeCompressed();
    }

    //======================== assemble AD pattern part =============================
    void assembleADpatternBlock_kOmega()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        // SUPG block nonlinear A pattern
        m_blockADpattern.resize(m_varDofs, m_dofs);
        m_blockADpattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMADVisitorKOmega<T> >(m_blockADpattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockADpattern.makeCompressed();
        m_blockADpattern = 0. * m_blockADpattern;
    }

    //======================== assemble AD part =============================
    void assembleADBlock_kOmega()
    {
        m_blockAD.resize(m_varDofs, m_dofs);
        m_rhsAD.setZero();

        m_blockAD.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMADVisitorKOmega<T> >(m_blockAD, m_rhsAD, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockAD.makeCompressed();
    }

    //======================== assemble isoAD pattern part =============================
    void assembleIsoADpatternBlock_kOmega()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        m_blockIsoADpattern.resize(m_varDofs, m_dofs);
        m_blockIsoADpattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMisoADVisitorKOmega<T> >(m_blockIsoADpattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockIsoADpattern.makeCompressed();
        m_blockIsoADpattern = 0. * m_blockIsoADpattern;
    }

    //======================== assemble AD part =============================
    void assembleIsoADBlock_kOmega()
    {
        m_blockIsoAD.resize(m_varDofs, m_dofs);
        m_rhsIsoAD.setZero();

        m_blockIsoAD.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMisoADVisitorKOmega<T> >(m_blockIsoAD, m_rhsIsoAD, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockIsoAD.makeCompressed();
    }

    //======================== assemble tanh-CSD pattern part =============================
    void assembleTanhCSDpatternBlock_kOmega()
    {
        gsMatrix<T> dummyRhs(m_varDofs, m_numVar);
        dummyRhs.setZero();

        m_blockTanhCSDpattern.resize(m_varDofs, m_dofs);
        m_blockTanhCSDpattern.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMCSDVisitorKOmega<T> >(m_blockTanhCSDpattern, dummyRhs, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockTanhCSDpattern.makeCompressed();
        m_blockTanhCSDpattern = 0. * m_blockTanhCSDpattern;
    }

    //======================== assemble tanh-CSD part =============================
    void assembleTanhCSDBlock_kOmega()
    {
        m_blockTanhCSD.resize(m_varDofs, m_dofs);
        m_rhsTanhCSD.setZero();

        m_blockTanhCSD.reserve(gsVector<int>::Constant(m_dofs, m_nonZerosPerCol));
        uwbBlockAssemblerBase<T>::template assembleBlock< uwbTMCSDVisitorKOmega<T> >(m_blockTanhCSD, m_rhsTanhCSD, m_nonZerosPerCol);
        //if (getAssemblerOptions().intStrategy == iFace::dg) { }
        m_blockTanhCSD.makeCompressed();
    }

public:
    bool checkWallDistanceBasedTM()
    {
        return (getTMEvaluator() == "koSST" || getTMEvaluator() == "koSSTMenter2009" || getTMEvaluator() == "koSAS"
                || getTMEvaluator() == "koSAS_SS" || getTMEvaluator() == "koSAS_SO" || getTMEvaluator() == "koSAS_OO");
    }

public:
    int getNumVar() const { return m_numVar; }
    int numVarDofs() const { return m_varDofs; }
    const gsDofMapper& getMapper() const { return m_dofMappers.back(); }

    const gsMultiBasis<T>&   getSolBasis() const { return getBases().back(); }
    gsMultiBasis<T>&        getSolBasis() { return getBases().back(); }

    std::string getTMEvaluator() {return m_params.settings().getTMEvaluator(); }
    const gsField<T>& getPoissonSolution() { return m_solPoisson; }
    const gsMultiPatch<T>& getPoissonPatches() { return m_patchesPoisson; }
    const gsMultiBasis<T>& getPoissonBasis() { return m_basesPoisson; }

    const std::vector<gsField<T> >& getKOmegaCoeffSolField() { return m_kOmegaDiffCoeffSolField; }

    int     getTauStabTypeSUPG() const { return m_params.settings().get(constantsINS::tauStabTypeSUPG); }
    int     getTauStabTypeCW() const { return m_params.settings().get(constantsINS::tauStabTypeCW); }
    int     getTauStabTypeAD() const { return m_params.settings().get(constantsINS::tauStabTypeAD); }
    int     getDiamPowerSUPG() const { return m_params.settings().get(constantsINS::diamPowerSUPG); }
    int     getDiamPowerCW() const { return m_params.settings().get(constantsINS::diamPowerCW); }
    int     getDiamPowerAD() const { return m_params.settings().get(constantsINS::diamPowerAD); }
    int     getCrosswindType() const { return m_params.settings().get(constantsINS::crosswindType); }
    int     getCWresidualType() const { return m_params.settings().get(constantsINS::CWresidualType); }
    int     getIsoADtype() const { return m_params.settings().get(constantsINS::isoADtype); }

    bool    isFCTlowOrder() const { return m_params.settings().get(constantsINS::TMfct_lowOrder); }

    T       getTimeStep() const { return m_params.settings().get(constantsINS::timeStep); }
    bool    isTimeDerTerm() const { return m_params.settings().get(constantsINS::timeDerTerm); }

    T       getStartXPointProduction() const { return m_params.settings().get(constantsINS::productionXPoint); }
    bool    isTMproductionLimited() const { return m_params.settings().get(constantsINS::limitTMProduction); }

    const gsSparseMatrix<T>& getBlockM() const { return m_blockM; }
    const gsSparseMatrix<T>& getBlockN() const { return m_blockN; }
    const gsSparseMatrix<T>& getBlockNpattern() const { return m_blockNpattern; }
    const gsSparseMatrix<T>& getBlocksAnonlin() const { return m_blocksAnonlin; }
    const gsSparseMatrix<T>& getBlocksReaction() const { return m_blocksReaction; }
    const gsSparseMatrix<T>& getBlockBlend() const { return m_blockBlend; }
    const gsSparseMatrix<T>& getBlocksAnonlinPattern() const { return m_blocksAnonlinPattern; }
    const gsSparseMatrix<T>& getBlocksReactionPattern() const { return m_blocksReactionPattern; }
    const gsMatrix<T>& getRhsM() const { return m_rhsM; }
    const gsMatrix<T>& getRhsN() const { return m_rhsN; }
    const gsMatrix<T>& getRhsAnonlin() const { return m_rhsAnonlin; }
    const gsMatrix<T>& getRhsReaction() const { return m_rhsReaction; }
    const gsMatrix<T>& getRhsBlend() const { return m_rhsBlend; }
    const gsMatrix<T>& getRhsF() const { return m_rhsF; }

    const gsSparseMatrix<T>& getBlocks() const { return m_blocks; }
    const gsSparseMatrix<T>& getBlocksPattern() const { return m_blocksPattern; }
    const gsMatrix<T>& getRhs() const { return m_rhs; }

    //SUPG getters
    const gsSparseMatrix<T>& getBlockM_SUPG() const { return m_blockM_SUPG; }
    const gsSparseMatrix<T>& getBlockMpattern_SUPG() const { return m_blockMpattern_SUPG; }
    const gsSparseMatrix<T>& getBlockN_SUPG() const { return m_blockN_SUPG; }
    const gsSparseMatrix<T>& getBlockNpattern_SUPG() const { return m_blockNpattern_SUPG; }

    const gsSparseMatrix<T>& getBlocksAnonlin_SUPG() const { return m_blocksAnonlin_SUPG; }
    const gsSparseMatrix<T>& getBlocksAnonlinPattern_SUPG() const { return m_blocksAnonlinPattern_SUPG; }
    const gsSparseMatrix<T>& getBlocksReaction_SUPG() const { return m_blocksReaction_SUPG; }
    const gsSparseMatrix<T>& getBlocksReactionPattern_SUPG() const { return m_blocksReactionPattern_SUPG; }
    const gsSparseMatrix<T>& getBlockBlend_SUPG() const { return m_blockBlend_SUPG; }

    const gsMatrix<T>& getRhsM_SUPG() const { return m_rhsM_SUPG; }
    const gsMatrix<T>& getRhsN_SUPG() const { return m_rhsN_SUPG; }
    const gsMatrix<T>& getRhsAnonlin_SUPG() const { return m_rhsAnonlin_SUPG; }
    const gsMatrix<T>& getRhsReaction_SUPG() const { return m_rhsReaction_SUPG; }
    const gsMatrix<T>& getRhsBlend_SUPG() const { return m_rhsBlend_SUPG; }
    const gsMatrix<T>& getRhsF_SUPG() const { return m_rhsF_SUPG; }

    const gsSparseMatrix<T>& getBlocks_SUPG() const { return m_blocks_SUPG; }
    const gsSparseMatrix<T>& getBlocksPattern_SUPG() const { return m_blocksPattern_SUPG; }
    const gsMatrix<T>& getRhs_SUPG() const { return m_rhs_SUPG; }

    //CROSSWIND getters
    const gsSparseMatrix<T>& getBlockCrosswind() const { return m_blockCrosswind; }
    const gsSparseMatrix<T>& getBlockCrosswindPattern() const { return m_blockCrosswindPattern; }
    const gsMatrix<T>& getRhsCrosswind() const { return m_rhsCrosswind; }

    //AD getters
    const gsSparseMatrix<T>& getBlockAD() const { return m_blockAD; }
    const gsSparseMatrix<T>& getBlockADpattern() const { return m_blockADpattern; }
    const gsMatrix<T>& getRhsAD() const { return m_rhsAD; }

    //isoAD getters
    const gsSparseMatrix<T>& getBlockIsoAD() const { return m_blockIsoAD; }
    const gsSparseMatrix<T>& getBlockIsoADpattern() const { return m_blockIsoADpattern; }
    const gsMatrix<T>& getRhsIsoAD() const { return m_rhsIsoAD; }

    //tanh-CSD getters
    const gsSparseMatrix<T>& getBlockTanhCSD() const { return m_blockTanhCSD; }
    const gsSparseMatrix<T>& getBlockTanhCSDpattern() const { return m_blockTanhCSDpattern; }
    const gsMatrix<T>& getRhsTanhCSD() const { return m_rhsTanhCSD; }

    bool    isTMsupg() const { return m_params.settings().get(constantsINS::TMsupg); }
    bool    isTMcrosswind() const { return m_params.settings().get(constantsINS::TMcrosswind); }
    bool    isTMisoAD() const { return m_params.settings().get(constantsINS::TMisoAD); }
    bool    isTMad() const { return m_params.settings().get(constantsINS::TMad); }
    bool    isTMtanhCSD() const { return m_params.settings().get(constantsINS::tanhCSD); }

    bool isStabilization()
    {
        return(isTMsupg() || isTMcrosswind() || isTMisoAD() || isTMad() || isTMtanhCSD());
    }

protected:

    int m_numVar;
    int m_varDofs;
    int m_nonZerosPerCol;
    gsField<T> m_uSolField; // u_n (velocity from n-th time step)
    gsField<T> m_solPoisson;
    gsMultiPatch<T> m_patchesPoisson;
    gsMultiBasis<T> m_basesPoisson;
    std::vector<gsField<T> > m_kOmegaDiffCoeffSolField;

    gsField<T>  m_currentSolField, m_oldTimeSolField;

    // members for k-omega model
    gsSparseMatrix<T> m_blockM;
    gsSparseMatrix<T> m_blockNpattern;
    gsSparseMatrix<T> m_blockN;
    gsSparseMatrix<T> m_blocksAnonlin;
    gsSparseMatrix<T> m_blocksReaction;
    gsSparseMatrix<T> m_blockBlend;
    gsSparseMatrix<T> m_blocksAnonlinPattern;
    gsSparseMatrix<T> m_blocksReactionPattern;
    gsMatrix<T> m_rhsM;
    gsMatrix<T> m_rhsN;
    gsMatrix<T> m_rhsAnonlin;
    gsMatrix<T> m_rhsReaction;
    gsMatrix<T> m_rhsBlend;
    gsMatrix<T> m_rhsF;

    gsSparseMatrix<T> m_blocks;
    gsSparseMatrix<T> m_blocksPattern;
    gsMatrix<T> m_rhs;

    //SUPG members
    gsSparseMatrix<T> m_blockM_SUPG;
    gsSparseMatrix<T> m_blockMpattern_SUPG;
    gsSparseMatrix<T> m_blockN_SUPG;
    gsSparseMatrix<T> m_blockNpattern_SUPG;

    gsSparseMatrix<T> m_blocksAnonlin_SUPG;
    gsSparseMatrix<T> m_blocksAnonlinPattern_SUPG;
    gsSparseMatrix<T> m_blocksReaction_SUPG;
    gsSparseMatrix<T> m_blocksReactionPattern_SUPG;
    gsSparseMatrix<T> m_blockBlend_SUPG;

    gsMatrix<T> m_rhsM_SUPG;
    gsMatrix<T> m_rhsN_SUPG;
    gsMatrix<T> m_rhsAnonlin_SUPG;
    gsMatrix<T> m_rhsReaction_SUPG;
    gsMatrix<T> m_rhsBlend_SUPG;
    gsMatrix<T> m_rhsF_SUPG;

    gsSparseMatrix<T> m_blocks_SUPG;
    gsSparseMatrix<T> m_blocksPattern_SUPG;
    gsMatrix<T> m_rhs_SUPG;

    //CROSSWIND members
    gsSparseMatrix<T> m_blockCrosswind;
    gsSparseMatrix<T> m_blockCrosswindPattern;
    gsMatrix<T> m_rhsCrosswind;

    //AD members
    gsSparseMatrix<T> m_blockAD;
    gsSparseMatrix<T> m_blockADpattern;
    gsMatrix<T> m_rhsAD;

    //isoAD members
    gsSparseMatrix<T> m_blockIsoAD;
    gsSparseMatrix<T> m_blockIsoADpattern;
    gsMatrix<T> m_rhsIsoAD;

    //tanh-CSD members
    gsSparseMatrix<T> m_blockTanhCSD;
    gsSparseMatrix<T> m_blockTanhCSDpattern;
    gsMatrix<T> m_rhsTanhCSD;

    // members from uwbBlockAssemblerBase
    using Base::m_dofs;
    using Base::m_numThreads;
    using Base::m_tarDim;
    using Base::m_bUnsteady;
    using Base::m_params;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_solution;

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
