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
    //Base::initMembers();

    m_bc = m_paramsPtr->getBCs();

    real_t uFreeStream = m_paramsPtr->options().getReal("TM.uFreeStream");
    real_t turbIntensity = m_paramsPtr->options().getReal("TM.turbIntensity");
    real_t viscRatio = m_paramsPtr->options().getReal("TM.viscosityRatio");

    // boundary conditions are automatically set based on boundary conditions for RANS
    m_kin = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    m_kwall = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    gsInfo << "kin = " << util::to_string(m_kin) << ", ";
    gsInfo << "kwall = " << util::to_string(m_kwall) << std::endl;
    gsFunctionExpr<T> Kin(util::to_string(m_kin), m_tarDim);
    gsFunctionExpr<T> Kwall(util::to_string(m_kwall), m_tarDim);
    std::vector<std::pair<int, boxSide> > bndIn = m_paramsPtr->getBndIn();
    std::vector<std::pair<int, boxSide> > bndWall = m_paramsPtr->getBndWall();
    addBCs(m_bc, bndIn, bndWall, Kin, Kwall, 2);

    m_oin = m_kin / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at inlet
    m_owall = m_kwall / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at the wall
    //real_t oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2)); // other way for prescribing this boundary condition
    gsInfo << "oin = " << util::to_string(m_oin) << ", ";
    gsInfo << "owall = " << util::to_string(m_owall) << std::endl;
    gsFunctionExpr<T> Oin(util::to_string(m_oin), m_tarDim);
    gsFunctionExpr<T> Owall(util::to_string(m_owall), m_tarDim);
    addBCs(m_bc, bndIn, bndWall, Oin, Owall, 3);

    m_paramsPtr->setBCs(m_bc);

    /*
    m_dofMappers.resize(numTMvars);
    for (short_t i = 0; i < numTMvars; i++)
    {
        std::vector<gsMultiBasis<T> > bases = getBases();
        gsBoundaryConditions<T> bc = getBCs();
        //m_bases[i].getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, m_bc, m_dofMappers[i], i);
        bases[i].getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, bc, m_dofMappers[i], i);
    }
        */
    std::vector<gsMultiBasis<T> > bases = m_paramsPtr->getBases();
    gsBoundaryConditions<T> bc = m_paramsPtr->getBCs();
    for (size_t i = 0; i < bases.size(); i++)
        bases[i].getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, bc, m_dofMappers[i], i);

    updateSizes();

    //m_SSTPtr = memory::make_shared_not_owned(new SSTModel<T>(m_paramsPtr));
    
    m_visitorLinearSST_K = gsTMVisitorLinearSST<T, MatOrder>(m_paramsPtr, 2);
    m_visitorLinearSST_O = gsTMVisitorLinearSST<T, MatOrder>(m_paramsPtr, 3);
    m_visitorLinearSST_K.initialize();
    m_visitorLinearSST_O.initialize();
    m_visitorLinearSST_K.updateDofMappers(m_dofMappers);
    m_visitorLinearSST_O.updateDofMappers(m_dofMappers);

    m_visitorTimeIterationSST_K = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 2);
    m_visitorTimeIterationSST_O = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 3);
    m_visitorTimeIterationSST_K.initialize();
    m_visitorTimeIterationSST_O.initialize();
    m_visitorTimeIterationSST_K.updateDofMappers(m_dofMappers);
    m_visitorTimeIterationSST_O.updateDofMappers(m_dofMappers);
    gsField<T> velocity = m_paramsPtr->getVelocitySolution();
    m_visitorTimeIterationSST_K.setCurrentSolution(velocity);
    m_visitorTimeIterationSST_O.setCurrentSolution(velocity);

    m_visitorNonlinearSST_K = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 2);
    m_visitorNonlinearSST_O = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 3);
    //m_visitorNonlinearSST_K.setSSTModelEvaluator(m_SSTPtr);
    //m_visitorNonlinearSST_O.setSSTModelEvaluator(m_SSTPtr);
    m_visitorNonlinearSST_K.initialize();
    m_visitorNonlinearSST_O.initialize();
    m_visitorNonlinearSST_K.setCurrentSolution(m_solution);
    m_visitorNonlinearSST_O.setCurrentSolution(m_solution);
    m_visitorNonlinearSST_K.updateDofMappers(m_dofMappers);
    m_visitorNonlinearSST_O.updateDofMappers(m_dofMappers);

    //m_isMassMatReady = false;
    //m_massMatBlocks.resize(2); // [k, omega]
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateSizes()
{
    Base:: updateSizes();

    m_solution.resize(m_dofs, 1);
    m_solution.topRows(m_kdofs[0]).setConstant(m_kwall);
    m_solution.bottomRows(m_kdofs[1]).setConstant(m_owall);
    //gsInfo << m_solution.sum() << std::endl;

    m_currentFieldK = constructSolution(m_solution, 2);
    m_currentFieldO = constructSolution(m_solution, 3);
    m_paramsPtr->setKSolution(m_currentFieldK);
    m_paramsPtr->setOmegaSolution(m_currentFieldO);

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

    //m_currentFieldK = constructSolution(m_solution, 2);
    //m_currentFieldO = constructSolution(m_solution, 3);
    m_oldTimeFieldK = m_currentFieldK;
    m_oldTimeFieldO = m_currentFieldO;

    //m_blockTimeDiscr.resize(m_pshift, m_pshift);
    //m_rhsTimeDiscr.setZero(m_pshift, 1);
    
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateDofMappers()
{
    m_visitorLinearSST_K.updateDofMappers(m_dofMappers);
    m_visitorLinearSST_O.updateDofMappers(m_dofMappers);
    m_visitorTimeIterationSST_K.updateDofMappers(m_dofMappers);
    m_visitorTimeIterationSST_O.updateDofMappers(m_dofMappers);
    m_visitorNonlinearSST_K.updateDofMappers(m_dofMappers);
    m_visitorNonlinearSST_O.updateDofMappers(m_dofMappers);
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateCurrentSolField(const gsMatrix<T>& solVector, bool updateSol)
{
    m_currentFieldK = constructSolution(solVector, 2);
    //m_visitorNonlinearSST_K.setCurrentSolution(solVector);
    m_paramsPtr->setKSolution(m_currentFieldK);
    m_currentFieldO = constructSolution(solVector, 3);
    //m_visitorNonlinearSST_O.setCurrentSolution(solVector);
    m_paramsPtr->setOmegaSolution(m_currentFieldO);

    m_visitorNonlinearSST_K.setCurrentSolution(solVector);
    m_visitorNonlinearSST_O.setCurrentSolution(solVector);
    
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
    m_blockLinearO.resize(m_kdofs[1], m_kdofs[1]);;
    m_rhsLinearK.setZero();
    m_rhsLinearO.setZero();
    
    // memory allocation
    m_blockLinearK.reserve(gsVector<index_t>::Constant(m_blockLinearK.outerSize(), m_nnzPerRowTM));
    m_blockLinearO.reserve(gsVector<index_t>::Constant(m_blockLinearO.outerSize(), m_nnzPerRowTM));

    // in this case, linear blocks correspond to time discretization blocks only
    this->assembleBlock(m_visitorLinearSST_K, 2, m_blockLinearK, m_rhsLinearK);
    this->assembleBlock(m_visitorLinearSST_O, 3, m_blockLinearO, m_rhsLinearO);

    //gsInfo << "LinearK mat: " << m_blockLinearK.sum() << std::endl;
    //gsInfo << "LinearO mat: " << m_blockLinearO.sum() << std::endl;
    
    //if(m_paramsPtr->getPde().source()) // if the continuity eqn rhs is given
    //    this->assembleRhs(m_visitorG, 1, m_rhsFG);

    //m_rhsLinearK = m_blockLinearK * m_solution.topRows(m_kdofs[0]);
    //m_rhsLinearO = m_blockLinearO * m_solution.bottomRows(m_kdofs[1]);

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

    //m_visitorNonlinearSST_K.setConstants();
    //m_visitorNonlinearSST_O.setConstants();
    this->assembleBlock(m_visitorNonlinearSST_K, 2, m_blockNonlinearK, m_rhsNonlinearK);
    this->assembleBlock(m_visitorNonlinearSST_O, 3, m_blockNonlinearO, m_rhsNonlinearO);
    
    //gsInfo << "NonlinearK mat: " << m_blockNonlinearK.sum() << std::endl;
    //gsInfo << "NonlinearO mat: " << m_blockNonlinearO.sum() << std::endl;
    //gsInfo << "NonlinearK rhs: " << m_rhsNonlinearK.sum() << std::endl;
    //gsInfo << "NonlinearO rhs: " << m_rhsNonlinearO.sum() << std::endl;
    //gsInfo << m_blockNonlinearK.sum() << ", " << m_blockNonlinearO.sum() << "; " << m_rhsNonlinearK.sum() << ", " << m_rhsNonlinearO.sum() << std::endl;
    
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
void gsTMAssemblerSST<T, MatOrder>::fillGlobalMat(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat, index_t unk)
{
    if (unk == 2) // sourceMat is a block for k
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row(), it.col()) += it.value();
    }
    else if (unk == 3) // sourceMat is a block for omega
    {
        for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
            for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
                globalMat.coeffRef(it.row() + m_kdofs[0], it.col() + m_kdofs[0]) += it.value();
    }
    else // something wrong happened
        GISMO_ASSERT((unk == 2) || (unk == 3), "Wrong matrix size for k or omega component in SST model.");
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillBaseSystem() 
{
    gsVector<index_t> nonZerosVector(m_dofs);
    
    nonZerosVector.topRows(m_kdofs[0]) = getNnzVectorPerOuter(m_blockLinearK);
    nonZerosVector.bottomRows(m_kdofs[1]) = getNnzVectorPerOuter(m_blockLinearO);

    m_baseMatrix.resize(m_dofs, m_dofs);
    m_baseMatrix.reserve(nonZerosVector);

    fillGlobalMat(m_baseMatrix, m_blockLinearK, 2);
    fillGlobalMat(m_baseMatrix, m_blockLinearO, 3);
    
    //m_baseRhs.topRows(m_kdofs[0]) = m_rhsLinearK;
    //m_baseRhs.bottomRows(m_kdofs[1]) = m_rhsLinearO;
    m_baseRhs.resize(m_dofs, 1);

    m_isBaseReady = true;
    m_isSystemReady = false;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillSystem()
{
    m_matrix = m_baseMatrix;

    fillGlobalMat(m_matrix, m_blockNonlinearK, 2);
    fillGlobalMat(m_matrix, m_blockNonlinearO, 3);

    fillGlobalMat(m_matrix, m_blockTimeIterationK, 2);
    fillGlobalMat(m_matrix, m_blockTimeIterationO, 3);

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_rhs = m_baseRhs;
    m_rhs.topRows(m_kdofs[0]) += m_rhsLinearK;
    m_rhs.bottomRows(m_kdofs[1]) += m_rhsLinearO;
    m_rhs.topRows(m_kdofs[0]) += m_rhsNonlinearK;
    m_rhs.bottomRows(m_kdofs[1]) += m_rhsNonlinearO;
    m_rhs.topRows(m_kdofs[0]) += m_rhsTimeIterationK;
    m_rhs.bottomRows(m_kdofs[1]) += m_rhsTimeIterationO;

    m_isSystemReady = true;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::initialize()
{
    //m_visitorNonlinearSST_K.setSSTModelEvaluator(m_SSTPtr);
    //m_visitorNonlinearSST_O.setSSTModelEvaluator(m_SSTPtr);

    gsFlowAssemblerBase<T, MatOrder>::initialize();  

    // initialization of distance field
    //gsMultiPatch<T> patches = m_paramsPtr->getPde().patches();    // multipatch representing the computational domain
    //gsMultiBasis<T> basis = m_paramsPtr->getBases()[1];           // pressure basis as base basis for distance computation

    //std::vector<std::pair<int, boxSide> > bndIn = m_paramsPtr->getBndIn();
    //std::vector<std::pair<int, boxSide> > bndWall = m_paramsPtr->getBndWall();
    //index_t numRefs = m_paramsPtr->options().getInt("TM.addRefsDF");

    //gsKnotVector<T> kv = patches.patch(0).knots(0);

    //computeDistanceField(patches, basis, numRefs, bndIn, bndWall, distanceField);
    //computeDistanceField(m_paramsPtr, distanceField);
    gsField<T> distfield = computeDistanceField<T>(m_paramsPtr);
    //gsMatrix<T> mat(2, 2);
    //mat << 0.2, 0.4, 0.6, 0.8;
    //gsInfo << distfield.value(mat, 0) << ", " << distfield.value(mat, 1) << ", " << distfield.value(mat, 2) << std::endl;
    //gsInfo << distfield.nPatches() << std::endl;
    //gsMultiPatch<T> mppom = distfield.patches();
    m_paramsPtr->setDistanceField(distfield);
    //if (m_paramsPtr->options().getSwitch("plot"))
        //plotQuantityFromSolution("distance", distanceField, "distancefield", 10000);
    //gsMultiPatch<T> mppom = m_paramsPtr->getDistanceField().patches();
    //gsInfo << m_paramsPtr->getDistanceField().value(mat, 0) << ", " << m_paramsPtr->getDistanceField().value(mat, 1) << ", " << m_paramsPtr->getDistanceField().value(mat, 2) << std::endl;
    //gsInfo << m_paramsPtr->getDistanceField().nPatches() << std::endl;
    //gsField<T> dfield = m_paramsPtr->getDistanceField();
    gsWriteParaview<T>(distfield, "distanceField", 10000);

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
    
    m_isInitialized = true;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    GISMO_ASSERT(m_isInitialized, "Assembler must be initialized first, call initialize()");

    updateCurrentSolField(solVector, updateSol);
    
    //if (m_solution.sum() != 0)
    Base::updateAssembly();

    if(updateSol)
    {
        //gsInfo << "Updating time iteration terms in TM solver ... ";
        // matrix and rhs cleaning
        m_blockTimeIterationK.resize(m_kdofs[0], m_kdofs[0]);
        m_blockTimeIterationO.resize(m_kdofs[1], m_kdofs[1]);
        m_rhsTimeIterationK.setZero();
        m_rhsTimeIterationO.setZero();

        // memory allocation
        m_blockTimeIterationK.reserve(gsVector<index_t>::Constant(m_blockTimeIterationK.outerSize(), m_nnzPerRowTM));
        m_blockTimeIterationO.reserve(gsVector<index_t>::Constant(m_blockTimeIterationO.outerSize(), m_nnzPerRowTM));

        gsField<T> velocity = m_paramsPtr->getVelocitySolution();
        m_visitorTimeIterationSST_K.setCurrentSolution(velocity);
        m_visitorTimeIterationSST_O.setCurrentSolution(velocity);

        this->assembleBlock(m_visitorTimeIterationSST_K, 2, m_blockTimeIterationK, m_rhsTimeIterationK);
        this->assembleBlock(m_visitorTimeIterationSST_O, 3, m_blockTimeIterationO, m_rhsTimeIterationO);

        //gsInfo << "TimeiterationK mat: " << m_blockTimeIterationK.sum() << std::endl;
        //gsInfo << "TimeiterationO mat: " << m_blockTimeIterationO.sum() << std::endl;
        //gsInfo << "TimeiterationK rhs: " << m_rhsTimeIterationK.sum() << std::endl;
        //gsInfo << "TimeiterationO rhs: " << m_rhsTimeIterationO.sum() << std::endl;

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