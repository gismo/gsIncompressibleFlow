/** @file gsFlowSolverBase.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Honnerova
*/

#pragma once
#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <unordered_set>

namespace gismo
{

template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::initMembers()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
    m_iterationNumber = 0;
    m_relNorm = std::numeric_limits<T>::infinity();

    m_fileOutput = m_paramsPtr->options().getSwitch("fileOutput");
    m_dispOutput = !m_paramsPtr->options().getSwitch("quiet");

    if(m_fileOutput)
        createOutputFile();

    m_initAssembT = 0;
    m_assembT = 0;
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::updateSizes()
{
    m_solution.setZero(getAssembler()->numDofs(), 1);
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::createOutputFile()
{
    std::string fileName = m_paramsPtr->options().getString("outFile");

    if (fileName == "")
        fileName = this->getName() + "_output.txt";

    m_outFile.open(fileName);

    std::stringstream output;
    output << "\n" << m_paramsPtr->options() << "\n";
    gsWriteOutput(m_outFile, output.str(), m_fileOutput, false);
}


template<class T, int MatOrder>
real_t gsFlowSolverBase<T, MatOrder>::stopwatchStart()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsPtr->options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        m_clock.restart();

    return 0.0;
}


template<class T, int MatOrder>
real_t gsFlowSolverBase<T, MatOrder>::stopwatchStop()
{

#ifdef GISMO_WITH_PETSC
    if (m_paramsPtr->options().getSwitch("parallel"))
    {
        MPI_Barrier(PETSC_COMM_WORLD);
        return MPI_Wtime();
    }
    else
#endif
        return m_clock.stop();

}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::initialize()
{ 
    if (!getAssembler()->isInitialized())
    {
        real_t time0 = stopwatchStart();
        getAssembler()->initialize();
        real_t time1 = stopwatchStop();

        m_initAssembT += time1 - time0;
    }

    m_linSolverPtr = createLinSolver<T, MatOrder>(m_paramsPtr, getAssembler());
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::nextIteration(const unsigned numberOfIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

    for (unsigned iter = 0; iter < numberOfIterations; iter++)
        nextIteration();
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::solve(const int maxIterations, const T epsilon, const int minIterations)
{
    GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
    int iter = 0;
    m_relNorm = solutionChangeRelNorm();

    while ((iter < minIterations) || ((m_relNorm > epsilon) && (iter < maxIterations)))
    {
        m_outStream.str("");
        m_outStream << "Iteration number " << m_iterationNumber + 1 << "...";
        gsWriteOutput(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);

        nextIteration();
        m_relNorm = solutionChangeRelNorm();

        m_outStream.str("");
        m_outStream << " Solution change relative norm: " << m_relNorm;
        gsWriteOutputLine(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);

        iter++;
    }
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::updateAssembler(const gsMatrix<T>& sol, bool updateSol)
{ 
    real_t time0 = stopwatchStart();
    getAssembler()->update(sol, updateSol);
    real_t time1 = stopwatchStop();

    m_assembT += time1 - time0;
}


template<class T, int MatOrder>
T gsFlowSolverBase<T, MatOrder>::solutionChangeRelNorm() const
{
    T relNorm;

    if (m_iterationNumber)
        relNorm = solutionChangeRelNorm(getAssembler()->getSolution(), m_solution);
    else
        relNorm = std::numeric_limits<T>::infinity();

    return relNorm;
}


template<class T, int MatOrder>
T gsFlowSolverBase<T, MatOrder>::solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
{
    gsMatrix<T> solChangeVector = solOld - solNew;
    T relNorm = solChangeVector.norm() / solNew.norm();

    return relNorm;
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::writeSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew, std::string solstr)
{
    m_outStream.str("");
    m_outStream << "     " << solstr << " solution change relative norm: ";

    for (int i = 0; i < solOld.cols(); i++)
        m_outStream << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

    gsWriteOutputLine(m_outFile, m_outStream.str(), m_fileOutput, m_dispOutput);
}


template<class T, int MatOrder>
void gsFlowSolverBase<T, MatOrder>::markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk)
{
    getAssembler()->markDofsAsEliminatedZeros(boundaryDofs, unk);
    updateSizes();
}


template<class T, int MatOrder>
int gsFlowSolverBase<T, MatOrder>::checkGeoJacobian(int npts, T dist, T tol)
{
    // default values

    if (npts == -1)
        npts = m_paramsPtr->options().getInt("jac.npts");

    if (dist == -1)
        dist = m_paramsPtr->options().getReal("jac.dist");

    if (tol == -1)
        tol = m_paramsPtr->options().getReal("jac.tol");

    short_t dim = m_paramsPtr->getPde().patches().domainDim();
    size_t np = m_paramsPtr->getPde().patches().nPatches();

    npts = math::pow(math::abs(npts), dim-1); // number of pts on one patch side
    dist = math::abs(dist);
    tol = math::abs(tol);

    for (size_t p = 0; p < np; p++)
    {
        const gsGeometry<T>* patch = &m_paramsPtr->getPde().patches().patch(p);

        gsMatrix<T> parRange = patch->support();
        GISMO_ASSERT(parRange.rows() == dim, "checkGeoJacobian: something went wrong, parRange.rows() != dim.");

        std::unordered_set<gsVector<T>, gsVectorHash<T> > uniquePts; // set of unique instances of generated points
        gsMatrix<T> pts(dim, npts);

        for (short_t i = 0; i < dim; i++) // direction i fixed
        {
            gsMatrix<T> reducedRange(dim-1, parRange.cols()); // parameter range without direction i
            
            if (i > 0)
                reducedRange.topRows(i) = parRange.topRows(i);

            if (i < dim - 1)
                reducedRange.bottomRows(dim - i - 1) = parRange.bottomRows(dim - i - 1);

            // first and last point shifted by dist from walls
            reducedRange.col(0).array() += dist;
            reducedRange.col(1).array() -= dist;

            gsMatrix<T> reducedPts = gsPointGrid(reducedRange, npts); // uniformly distributed points in the reduced range

            if (i > 0)
                pts.topRows(i) = reducedPts.topRows(i);

            if (i < dim - 1)
                pts.bottomRows(dim - i - 1) = reducedPts.bottomRows(dim - i - 1);

            pts.row(i) = gsVector<T>::Constant(pts.cols(), parRange(i,0) + dist); // fixed coordinate in direction i (lower+dist)

            for(int j = 0; j < pts.cols(); j++)
                uniquePts.insert(pts.col(j));

            pts.row(i) = gsVector<T>::Constant(pts.cols(), parRange(i,1) - dist); // fixed coordinate in direction i (upper-dist)

            for(int j = 0; j < pts.cols(); j++)
                uniquePts.insert(pts.col(j));

        }

        for (const gsVector<T>& pt : uniquePts)
        {
            real_t jacDet = patch->jacobian(pt).determinant();

            if (jacDet < 0)
            {
                gsWarn << "Negative geometry jacobian!\n";
                return -1;
            }

            if (jacDet < tol)
            {
                gsWarn << "Geometry jacobian close to zero!\n";
                return 1;
            }
        }

        gsInfo << "\n";
    }

    return 0;
}


} // namespace gismo