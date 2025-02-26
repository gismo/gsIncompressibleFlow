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

    numTMvars = m_bases.size();
    m_dofMappers.resize(numTMvars);
    for (short_t i = 0; i < numTMvars; i++)
    {
        m_bases[i].getMapper(getAssemblerOptions().dirStrategy, getAssemblerOptions().intStrategy, getBCs(), m_dofMappers[i], 1);
    }

    updateSizes();

    /*
    m_visitorUUlin = gsINSVisitorUUlin<T, MatOrder>(m_paramsPtr);
    m_visitorUUlin.initialize();

    m_visitorUUnonlin = gsINSVisitorUUnonlin<T, MatOrder>(m_paramsPtr);
    m_visitorUUnonlin.initialize();
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);

    m_visitorUP = gsINSVisitorPU_withUPrhs<T, MatOrder>(m_paramsPtr);
    m_visitorUP.initialize();

    m_visitorF = gsINSVisitorRhsU<T, MatOrder>(m_paramsPtr);
    m_visitorF.initialize();

    m_visitorG = gsINSVisitorRhsP<T, MatOrder>(m_paramsPtr);
    m_visitorG.initialize();
    
    m_visitorTimeDiscr = gsINSVisitorUUtimeDiscr<T, MatOrder>(m_paramsPtr);
    m_visitorTimeDiscr.initialize();

    m_isMassMatReady = false;
    m_massMatBlocks.resize(2); // [velocity, pressure]
    */
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

    // upravit dale
    /*
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

    m_oldTimeVelField = m_currentVelField;

    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_rhsTimeDiscr.setZero(m_pshift, 1);
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateDofMappers()
{
    /*
    m_visitorUUlin.updateDofMappers(m_dofMappers);
    m_visitorUUnonlin.updateDofMappers(m_dofMappers);
    m_visitorUP.updateDofMappers(m_dofMappers);
    m_visitorF.updateDofMappers(m_dofMappers);
    m_visitorG.updateDofMappers(m_dofMappers);
    m_visitorTimeDiscr.updateDofMappers(m_dofMappers);
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateCurrentSolField(const gsMatrix<T> & solVector, bool updateSol)
{
    if (updateSol)
        m_solution = solVector;

    /*
    m_currentVelField = constructSolution(solVector, 0);
    m_visitorUUnonlin.setCurrentSolution(m_currentVelField);


    if (updateSol)
        m_oldTimeVelField = this->constructSolution(solVector, 0);
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::assembleLinearPart()
{
    /*
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

    if(m_paramsPtr->getPde().source()) // if the continuity eqn rhs is given
        this->assembleRhs(m_visitorG, 1, m_rhsFG);

    // matrix cleaning
    m_blockTimeDiscr.resize(m_pshift, m_pshift);
    m_blockTimeDiscr.reserve(gsVector<index_t>::Constant(m_blockTimeDiscr.outerSize(), m_nnzPerRowU));

    gsMatrix<T> dummyRhs;
    dummyRhs.setZero(m_pshift, 1);

    this->assembleBlock(m_visitorTimeDiscr, 0, m_blockTimeDiscr, dummyRhs);

    m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

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
    /*
    // matrix and rhs cleaning
    m_blockUUnonlin_comp.resize(m_udofs, m_udofs);
    m_rhsUnonlin.setZero();

    // memory allocation
    m_blockUUnonlin_comp.reserve(gsVector<index_t>::Constant(m_blockUUnonlin_comp.outerSize(), m_nnzPerRowU));

    this->assembleBlock(m_visitorUUnonlin, 0, m_blockUUnonlin_comp, m_rhsUnonlin);

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
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillGlobalMat_KK(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    /*
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
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillGlobalMat_OO(gsSparseMatrix<T, MatOrder>& globalMat, const gsSparseMatrix<T, MatOrder>& sourceMat)
{
    for (index_t outer = 0; outer < sourceMat.outerSize(); outer++)
        for (typename gsSparseMatrix<T, MatOrder>::InnerIterator it(sourceMat, outer); it; ++it)
            { }//globalMat.coeffRef(it.row() + m_pshift, it.col() + m_pshift) += it.value();
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillBaseSystem() 
{
    /*
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
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::fillSystem()
{
    /*
    m_matrix = m_baseMatrix;

    fillGlobalMat_UU(m_matrix, m_blockUUnonlin_comp);
    fillGlobalMat_UU(m_matrix, m_blockUUnonlin_whole);

    if (!m_matrix.isCompressed())
        m_matrix.makeCompressed();

    m_rhs = m_baseRhs;
    m_rhs.topRows(m_pshift) += m_rhsUnonlin;


    this->fillGlobalMat_UU(m_matrix, m_blockTimeDiscr);
    m_rhs.topRows(m_pshift) += m_rhsTimeDiscr;

    m_isSystemReady = true;
    */
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::initialize()
{
    Base::initialize();

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillBaseSystem();
    
    m_bInitialized = true;
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::update(const gsMatrix<T> & solVector, bool updateSol)
{
    Base::update(solVector, updateSol);

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();

    /*
    gsFlowAssemblerBase<T, MatOrder>::update(solVector, updateSol);

    if(updateSol)
        m_rhsTimeDiscr = m_blockTimeDiscr * m_solution.topRows(m_pshift);

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();
    */
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
    /*
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
    */
}


template<class T, int MatOrder>
gsMatrix<T> gsTMAssemblerSST<T, MatOrder>::getRhsK() const
{ 
    /*
    if (m_isSystemReady)
        return m_rhs.topRows(m_pshift);
    else
    {
        gsMatrix<T> rhsUpart = (m_rhsFG + m_rhsBtB).topRows(m_pshift);
        return (rhsUpart + m_rhsUlin + m_rhsUnonlin);
    }
    */
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