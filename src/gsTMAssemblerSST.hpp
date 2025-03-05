/** @file gsINSAssembler.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova (Hornikova), B. Bastl
*/

#pragma once
#include <gsIncompressibleFlow/src/gsTMAssemblerSST.h>

namespace gismo
{

template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::initMembers()
{ 
    Base::initMembers();

    real_t uFreeStream = m_paramsPtr->options().getReal("TM.uFreeStream");
    real_t turbIntensity = m_paramsPtr->options().getReal("TM.turbIntensity");
    real_t viscRatio = m_paramsPtr->options().getReal("TM.viscosityRatio");

    // boundary conditions are automatically set based on boundary conditions for RANS
    real_t kin = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    real_t kwall = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    gsInfo << "kin = " << util::to_string(kin) << std::endl;
    gsInfo << "kwall = " << util::to_string(kwall) << std::endl;
    gsFunctionExpr<T> Kin(util::to_string(kin), m_tarDim);
    gsFunctionExpr<T> Kwall(util::to_string(kwall), m_tarDim);
    std::vector<std::pair<int, boxSide> > bndIn = m_paramsPtr->getBndIn();
    std::vector<std::pair<int, boxSide> > bndWall = m_paramsPtr->getBndWall();
    addBCs(m_bc, bndIn, bndWall, Kin, Kwall);

    real_t oin = kin / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at inlet
    real_t owall = kwall / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at the wall
    //real_t oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2)); // other way for prescribing this boundary condition
    gsInfo << "oin = " << util::to_string(oin) << std::endl;
    gsInfo << "owall = " << util::to_string(owall) << std::endl;
    gsFunctionExpr<T> Oin(util::to_string(oin), m_tarDim);
    gsFunctionExpr<T> Owall(util::to_string(owall), m_tarDim);
    addBCs(m_bc, bndIn, bndWall, Oin, Owall);

    m_dofMappers.resize(numTMvars);
    for (short_t i = 0; i < numTMvars; i++)
    {
        m_bases[i].getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers[i], 1);
    }

    updateSizes();

    
    m_visitorLinearSST = gsTMVisitorLinearSST<T, MatOrder>(m_paramsPtr);
    m_visitorLinearSST.initialize();

    m_visitorTimeIterationSST_K = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, 0);
    m_visitorTimeIterationSST_O = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, 1);
    m_visitorTimeIterationSST_K.initialize();
    m_visitorTimeIterationSST_O.initialize();

    real_t sigmaK1 = m_paramsPtr->getSSTModel().get_sigmaK1();
    real_t sigmaK2 = m_paramsPtr->getSSTModel().get_sigmaK2();
    real_t sigmaO1 = m_paramsPtr->getSSTModel().get_sigmaO1();
    real_t sigmaO2 = m_paramsPtr->getSSTModel().get_sigmaO2();
    real_t viscosity = m_paramsPtr->getPde().viscosity();
    m_visitorNonlinearSST_K = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, sigmaK1, sigmaK2, viscosity, 0);
    m_visitorNonlinearSST_O = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, sigmaO1, sigmaO2, viscosity, 1);
    m_visitorNonlinearSST_K.initialize();
    m_visitorNonlinearSST_O.initialize();
    m_visitorNonlinearSST_K.setCurrentSolution(m_solution);
    m_visitorNonlinearSST_O.setCurrentSolution(m_solution);

    //m_isMassMatReady = false;
    //m_massMatBlocks.resize(2); // [k, omega]
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateSizes()
{
    Base:: updateSizes();

    for (short_t i = 0; i < numTMvars; i++)
    {
        m_ddof[i].setZero(m_dofMappers[i].boundarySize(), 1);

        if (this->getAssemblerOptions().dirStrategy == dirichlet::elimination)
        {
            this->computeDirichletDofs(i, i, m_ddof[i]);
        }
    }
    
    m_solution.setZero(m_dofs, 1);

    //m_currentSolField = constructSolution(m_solution, 0);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_matrix.resize(m_dofs, m_dofs);

    m_baseRhs.setZero(m_dofs, 1);
    m_rhs.setZero(m_dofs, 1);

    m_solution.setZero(m_dofs, 1);

    m_blockLinearK.resize(m_kdofs[0], m_kdofs[0]);
    m_blockLinearO.resize(m_kdofs[1], m_kdofs[1]);
    m_blockTimeIterationK.resize(m_kdofs[0], m_kdofs[0]);
    m_blockTimeIterationO.resize(m_kdofs[1], m_kdofs[1]);
    m_blockNonlinearK.resize(m_kdofs[0], m_kdofs[0]);
    m_blockNonlinearO.resize(m_kdofs[1], m_kdofs[1]);
        
    m_rhsLinearK.setZero(m_kdofs[0], 1);
    m_rhsLinearO.setZero(m_kdofs[1], 1);
    m_rhsTimeIterationK.setZero(m_kdofs[0], 1);
    m_rhsTimeIterationO.setZero(m_kdofs[1], 1);
    m_rhsNonlinearK.setZero(m_kdofs[0], 1);
    m_rhsNonlinearO.setZero(m_kdofs[1], 1);

    m_currentFieldK = constructSolution(m_solution, 0);
    m_currentFieldO = constructSolution(m_solution, 1);
    m_oldTimeFieldK = m_currentFieldK;
    m_oldTimeFieldO = m_currentFieldO;

    //m_blockTimeDiscr.resize(m_pshift, m_pshift);
    //m_rhsTimeDiscr.setZero(m_pshift, 1);
    
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateDofMappers()
{
    m_visitorLinearSST.updateDofMappers(m_dofMappers);
    m_visitorTimeIterationSST_K.updateDofMappers(m_dofMappers);
    m_visitorTimeIterationSST_O.updateDofMappers(m_dofMappers);
    m_visitorNonlinearSST_K.updateDofMappers(m_dofMappers);
    m_visitorNonlinearSST_O.updateDofMappers(m_dofMappers);
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateCurrentSolField(gsMatrix<T>& solVector, bool updateSol)
{
    m_currentFieldK = constructSolution(solVector, 0);
    m_visitorNonlinearSST_K.setCurrentSolution(solVector);
    m_paramsPtr->setKSolution(m_currentFieldK);
    m_currentFieldO = constructSolution(solVector, 1);
    m_visitorNonlinearSST_O.setCurrentSolution(solVector);
    m_paramsPtr->setOmegaSolution(m_currentFieldO);

    if (updateSol)
    {
        m_solution = solVector;
        m_oldTimeFieldK = m_currentFieldK;
        m_oldTimeFieldO = m_currentFieldO;
    }
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::assembleLinearPart()
{
    // matrix and rhs cleaning
    m_blockLinearK.resize(m_kdofs[0], m_kdofs[0]);;
    m_blockLinearO.resize(m_kdofs[0], m_kdofs[0]);;
    m_rhsLinearK.setZero();
    m_rhsLinearO.setZero();
    
    // memory allocation
    m_blockLinearK.reserve(gsVector<index_t>::Constant(m_blockLinearK.outerSize(), m_nnzPerRowTM));
    m_blockLinearO.reserve(gsVector<index_t>::Constant(m_blockLinearK.outerSize(), m_nnzPerRowTM));

    // in this case, linear blocks correspond to time discretization blocks only
    this->assembleBlock(m_visitorLinearSST, 0, m_blockLinearK, m_rhsLinearK);
    this->assembleBlock(m_visitorLinearSST, 1, m_blockLinearO, m_rhsLinearO);
    
    //if(m_paramsPtr->getPde().source()) // if the continuity eqn rhs is given
    //    this->assembleRhs(m_visitorG, 1, m_rhsFG);

    m_rhsLinearK = m_blockLinearK * m_solution.topRows(m_kdofs[0]);
    m_rhsLinearO = m_blockLinearO * m_solution.bottomRows(m_kdofs[1]);

    /*
    // mass matrices for velocity and pressure (needed for preconditioners)
    if ( m_paramsPtr->options().getString("lin.solver") == "iter" )
    {
        gsMatrix<T> dummyRhsU(m_pshift, 1);
        gsMatrix<T> dummyRhsP(m_pdofs, 1);

        gsINSVisitorUUmass<T, MatOrder> visitorUUmass(m_paramsPtr);
        gsINSVisitorPPmass<T, MatOrder> visitorPPmass(m_paramsPtr);
        
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

        //     gsINSVisitorPPlaplace<T, MatOrder> visitorPPlaplace(m_paramsPtr);
        //     visitorPPlaplace.initialize();

        //     m_pcdBlocks[0].resize(m_pdofs, m_pdofs);
        //     m_pcdBlocks[0].reserve(gsVector<index_t>::Constant(m_pdofs, m_nnzPerRowP));

        //     this->assembleBlock(visitorPPlaplace, 0, m_pcdBlocks[0], dummyRhsP);

        //     // prepare for BCs
        //     findPressureBoundaryIDs();
        // }
    }
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::assembleNonlinearPart()
{
    // matrix and rhs cleaning
    m_blockNonlinearK.resize(m_kdofs[0], m_kdofs[0]);
    m_blockNonlinearO.resize(m_kdofs[1], m_kdofs[1]);
    m_rhsNonlinearK.setZero();
    m_rhsNonlinearO.setZero();

    // memory allocation
    m_blockNonlinearK.reserve(gsVector<index_t>::Constant(m_blockNonlinearK.outerSize(), m_nnzPerRowTM));
    m_blockNonlinearO.reserve(gsVector<index_t>::Constant(m_blockNonlinearO.outerSize(), m_nnzPerRowTM));

    this->assembleBlock(m_visitorNonlinearSST_K, 0, m_blockNonlinearK, m_rhsNonlinearK);
    this->assembleBlock(m_visitorNonlinearSST_O, 1, m_blockNonlinearO, m_rhsNonlinearO);

    // linear operators needed for PCD preconditioner
    // if ( m_paramsPtr->options().getString("lin.solver") == "iter" &&  m_paramsRef.options().getString("lin.precType").substr(0, 3) == "PCD" )
    // {
    //     gsMatrix<T> dummyRhsP(m_pdofs, 1);

    //     // Robin BC

    //     // TODO

    //     // pressure convection operator

    //     gsINSVisitorPPconvection<T, MatOrder> visitorPPconv(m_paramsPtr);
    //     visitorPPconv.initialize();
    //     visitorPPconv.setCurrentSolution(m_currentVelField);

    //     m_pcdBlocks[2].resize(m_pdofs, m_pdofs);
    //     m_pcdBlocks[2].reserve(gsVector<index_t>::Constant(m_pdofs, m_nnzPerRowP));

    //     this->assembleBlock(visitorPPconv, 0, m_pcdBlocks[2], dummyRhsP);
    // }
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    if (sourceMat.rows() == m_kdofs[0]) // sourceMat is a block for k
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row(), it.col()) += it.value();
    }
    else if (sourceMat.rows() == m_kdofs[1]) // sourceMat is a block for omega
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row() + m_kdofs[0], it.col() + m_kdofs[0]) += it.value();
    }
    else // something wrong happened
        GISMO_ASSERT((sourceMat.rows() == m_kdofs[0]) || (sourceMat.rows() == m_kdofs[1]), "Wrong matrix size for k or omega component in SST model.");
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillBaseSystem() 
{
    gsVector<index_t> nonZerosVector(m_dofs);
    
    nonZerosVector.topRows(m_kdofs[0]) = getNnzVectorPerOuter(m_blockLinearK);
    nonZerosVector.bottomRows(m_kdofs[1]) = getNnzVectorPerOuter(m_blockLinearO);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_baseMatrix.reserve(nonZerosVector);

    fillGlobalMat(m_baseMatrix, m_blockLinearK);
    fillGlobalMat(m_baseMatrix, m_blockLinearO);
    
    m_baseRhs.topRows(m_kdofs[0]) = m_rhsLinearK;
    m_baseRhs.bottomRows(m_kdofs[1]) = m_rhsLinearO;

    m_isBaseReady = true;
    m_isSystemReady = false;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillSystem()
{
    m_matrix = m_baseMatrix;

    fillGlobalMat(m_matrix, m_blockNonlinearK);
    fillGlobalMat(m_matrix, m_blockNonlinearO);

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_rhs = m_baseRhs;
    m_rhs.topRows(m_kdofs[0]) += m_rhsNonlinearK;
    m_rhs.bottomRows(m_kdofs[1]) += m_rhsNonlinearO;

    this->fillGlobalMat(m_matrix, m_blockTimeIterationK);
    this->fillGlobalMat(m_matrix, m_blockTimeIterationO);
    m_rhs.topRows(m_kdofs[0]) += m_rhsTimeIterationK;
    m_rhs.bottomRows(m_kdofs[1]) += m_rhsTimeIterationO;

    m_isSystemReady = true;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::initialize()
{
    gsFlowAssemblerBase<T, MatOrder>::initialize();  

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
    
    m_bInitialized = true;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if(updateSol)
    {
        gsInfo << "Updating time iteration terms in TM solver ... ";
        // matrix and rhs cleaning
        m_blockTimeIterationK.resize(m_kdofs[0], m_kdofs[0]);
        m_blockTimeIterationO.resize(m_kdofs[1], m_kdofs[1]);
        m_rhsTimeIterationK.setZero();
        m_rhsTimeIterationO.setZero();

        // memory allocation
        m_blockTimeIterationK.reserve(gsVector<index_t>::Constant(m_blockTimeIterationK.outerSize(), m_nnzPerRowTM));
        m_blockTimeIterationO.reserve(gsVector<index_t>::Constant(m_blockTimeIterationO.outerSize(), m_nnzPerRowTM));

        this->assembleBlock(m_visitorTimeIterationSST_K, 0, m_blockTimeIterationK, m_rhsTimeIterationK);
        this->assembleBlock(m_visitorTimeIterationSST_O, 1, m_blockTimeIterationO, m_rhsTimeIterationO);

        m_rhsLinearK = m_blockLinearK * m_solution.topRows(m_kdofs[0]);
        m_rhsLinearO = m_blockLinearO * m_solution.bottomRows(m_kdofs[1]);
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();

}

/*
template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const index_t unk)
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
*/

/*
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
*/

/*
template<class T, int MatOrder>
gsField<T> gsTMAssemblerSST<T, MatOrder>::constructSolution(const gsMatrix<T>& solVector, index_t unk) const
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
*/

/*
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
*/

template<class T, int MatOrder>
gsSparseMatrix<T, MatOrder> gsTMAssemblerSST<T, MatOrder>::getBlockKK(bool linPartOnly)
{
    if (m_isSystemReady && !linPartOnly)
    {
        return m_matrix.topLeftCorner(m_kdofs[0], m_kdofs[0]);
    }
    else if (m_isBaseReady && linPartOnly)
    {
        return m_baseMatrix.topLeftCorner(m_kdofs[0], m_kdofs[0]);
    }
    else
    {
        gsSparseMatrix<T, MatOrder> blockK = m_blockLinearK;

        if(!linPartOnly)
            blockK += m_blockNonlinearK;

        blockK.makeCompressed();

        return blockK;
    }
}


template<class T, int MatOrder>
gsSparseMatrix<T, MatOrder> gsTMAssemblerSST<T, MatOrder>::getBlockOO(bool linPartOnly)
{
    if (m_isSystemReady && !linPartOnly)
    {
        return m_matrix.bottomRightCorner(m_kdofs[1], m_kdofs[1]);
    }
    else if (m_isBaseReady && linPartOnly)
    {
        return m_baseMatrix.bottomRightCorner(m_kdofs[1], m_kdofs[1]);
    }
    else
    {
        gsSparseMatrix<T, MatOrder> blockO = m_blockLinearO;

        if(!linPartOnly)
            blockO += m_blockNonlinearO;

        blockO.makeCompressed();

        return blockO;
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsTMAssemblerSST<T, MatOrder>::getRhsK() const
{ 
    if (m_isSystemReady)
        return m_rhs.topRows(m_kdofs[0]);
    else
    {
        return (m_rhsLinearK + m_rhsNonlinearK + m_rhsTimeIterationK);
    }
}


template<class T, int MatOrder>
gsMatrix<T> gsTMAssemblerSST<T, MatOrder>::getRhsO() const
{ 
    if (m_isSystemReady)
        return m_rhs.bottomRows(m_kdofs[1]);
    else
    {
        return (m_rhsLinearO + m_rhsNonlinearO + m_rhsTimeIterationO);
    }
}

// --- PCD functions ---

// template<class T, int MatOrder>
// void gsINSAssembler<T, MatOrder>::findPressureBoundaryIDs()
// {
//     findPressureBoundaryPartIDs(m_paramsPtr->getBndIn(), m_presInIDs);
//     findPressureBoundaryPartIDs(m_paramsPtr->getBndOut(), m_presOutIDs);
//     findPressureBoundaryPartIDs(m_paramsPtr->getBndWall(), m_presWallIDs);
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

} // namespace gismo