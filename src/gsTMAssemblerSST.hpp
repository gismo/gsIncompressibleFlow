/** @file gsTMAssemblerSST.hpp

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
    
    m_bc = m_paramsPtr->getBCs();

    real_t uFreeStream = m_paramsPtr->options().getReal("TM.uFreeStream");
    real_t turbIntensity = m_paramsPtr->options().getReal("TM.turbIntensity");
    real_t viscRatio = m_paramsPtr->options().getReal("TM.viscosityRatio");

    // boundary conditions are automatically set based on boundary conditions for RANS
    m_kin = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    //m_kwall = 1.5 * math::pow(uFreeStream * turbIntensity, 2);
    m_kwall = 1e-10;
    gsInfo << "kin = " << util::to_string(m_kin) << ", ";
    gsInfo << "kwall = " << util::to_string(m_kwall) << std::endl;
    gsFunctionExpr<T> Kin(util::to_string(m_kin), m_tarDim);
    gsFunctionExpr<T> Kwall(util::to_string(m_kwall), m_tarDim);
    std::vector<std::pair<int, boxSide> > bndIn = m_paramsPtr->getBndIn();
    std::vector<std::pair<int, boxSide> > bndWall = m_paramsPtr->getBndWall();
    addBCs(m_bc, bndIn, bndWall, Kin, Kwall, 2);

    m_oin = m_kin / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at inlet
    real_t inletWidth = 1.0;
    real_t Re = uFreeStream * inletWidth / m_viscosity; 
    int numSamplePts = 50; //number of sample points for which the distance to the boundary is computed
    real_t maxYplus = 2.5; //maximum dimensionless wall distance which is accepted
    real_t wallDistance = computeDimensionlessWallDistance<T>(m_paramsPtr, m_viscosity, Re, uFreeStream, maxYplus, numSamplePts, true, true);
    gsInfo << "\nminimum wallDistance = " << wallDistance << "\n";
    if (wallDistance > 0)
    {
        real_t beta = 0.0708;
        m_owall = 6 * m_viscosity / (beta * math::pow(wallDistance, 2)); // omega on the wall
    }
    else
        m_owall = 500; 
    //m_owall = m_kwall / (m_viscosity * viscRatio); // need to satisfy nu_T / nu approximately at the wall
    //real_t oBlade = 6 * viscosity / (beta * math::pow(wallDistance, 2)); // other way for prescribing this boundary condition
    gsInfo << "oin = " << util::to_string(m_oin) << ", ";
    gsInfo << "owall = " << util::to_string(m_owall) << std::endl;
    gsFunctionExpr<T> Oin(util::to_string(m_oin), m_tarDim);
    gsFunctionExpr<T> Owall(util::to_string(m_owall), m_tarDim);
    addBCs(m_bc, bndIn, bndWall, Oin, Owall, 3);

    m_paramsPtr->setBCs(m_bc);
    m_paramsPtr->updateDofMappers();
    updateSizes();

    m_visitorLinearSST_K = gsTMVisitorLinearSST<T, MatOrder>(m_paramsPtr, 2);
    m_visitorLinearSST_O = gsTMVisitorLinearSST<T, MatOrder>(m_paramsPtr, 3);
    m_visitorLinearSST_K.initialize();
    m_visitorLinearSST_O.initialize();

    m_visitorTimeIterationSST_K = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 2);
    m_visitorTimeIterationSST_O = gsTMVisitorTimeIterationSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 3);
    m_visitorTimeIterationSST_K.initialize();
    m_visitorTimeIterationSST_O.initialize();
    //gsField<T> velocity = m_paramsPtr->getVelocitySolution();
    //m_visitorTimeIterationSST_K.setCurrentSolution(velocity);
    //m_visitorTimeIterationSST_O.setCurrentSolution(velocity);

    m_visitorNonlinearSST_K = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 2);
    m_visitorNonlinearSST_O = gsTMVisitorNonlinearSST<T, MatOrder>(m_paramsPtr, m_TMModelPtr, 3);
    m_visitorNonlinearSST_K.initialize();
    m_visitorNonlinearSST_O.initialize();
    m_visitorNonlinearSST_K.setCurrentSolution(m_solution);
    m_visitorNonlinearSST_O.setCurrentSolution(m_solution);

}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateSizes()
{
    Base:: updateSizes();

    m_solution.resize(m_dofs, 1);
    m_solution.topRows(m_kdofs[0]).setConstant(m_kin);
    m_solution.bottomRows(m_kdofs[1]).setConstant(m_oin);
    
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

    m_oldTimeFieldK = m_currentFieldK;
    m_oldTimeFieldO = m_currentFieldO;
    
}


template<class T, int MatOrder>
void gsTMAssemblerSST<T, MatOrder>::updateCurrentSolField(const gsMatrix<T>& solVector, bool updateSol)
{
    m_currentFieldK = constructSolution(solVector, 2);
    m_paramsPtr->setKSolution(m_currentFieldK);
    m_currentFieldO = constructSolution(solVector, 3);
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

    this->assembleBlock(m_visitorNonlinearSST_K, 2, m_blockNonlinearK, m_rhsNonlinearK);
    this->assembleBlock(m_visitorNonlinearSST_O, 3, m_blockNonlinearO, m_rhsNonlinearO);

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
    gsFlowAssemblerBase<T, MatOrder>::initialize();  

    // initialization of distance field
    gsField<T> distfield = computeDistanceField<T>(m_paramsPtr);
    m_paramsPtr->setDistanceField(distfield);
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
    
    Base::updateAssembly();

    if(updateSol)
    {
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

        m_rhsLinearK = m_blockLinearK * m_solution.topRows(m_kdofs[0]);
        m_rhsLinearO = m_blockLinearO * m_solution.bottomRows(m_kdofs[1]);
    }

    if (m_paramsPtr->options().getSwitch("fillGlobalSyst"))
        fillSystem();

}


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

} // namespace gismo