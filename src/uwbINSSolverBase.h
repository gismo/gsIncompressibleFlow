/** @file uwbINSSolverBase.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbINSAssemblerBase.h"
#include "uwbINSSolverParams.h"
#include "uwbObjectiveFunctionEvaluator.h"

namespace gismo
{

template<class T>
class uwbINSSolverBase
{
public:
    uwbINSSolverBase()
    {
        m_pAssembler = NULL;
    }

    virtual ~uwbINSSolverBase()
    {
        if (m_pAssembler)
        {
            delete m_pAssembler;
            m_pAssembler = NULL;
        }

        m_ofileRelNorm.close();
    }

protected:
    virtual void initMembers()
    {
        m_solution.setZero(getAssembler()->numDofs(), 1);
        m_iterationNumber = 0;
        m_objFunValue = 0.;
        m_objError = std::numeric_limits<T>::infinity();
        m_relNorm = std::numeric_limits<T>::infinity();
        m_assembT = 0;
        m_solsetupT = 0;
        m_solveT = 0;
        m_ofileRelNorm.open("solChangeRelNorm.txt");

        #ifdef GISMO_WITH_PARDISO
            pardisoSetup(m_solver);
        #endif
    }

    virtual void reinitMembers() { initMembers(); }

    void openCloseFileObj(bool open = true)
    {
        if (open)
            m_ofile.open("objFunValues.txt");
        else
            m_ofile.close();
    }
public:
    virtual void initialize() 
    { 
        if (!getAssembler()->isInitialized())
        {
            m_clock.restart();
            getAssembler()->initialize();
            m_assembT += m_clock.stop();
        }
    }

    static void pardisoSetup(typename gsSparseSolver<T>::PardisoLU& solver)
    {
        solver.setParam(7, 15);
        solver.setParam(9, 13);
        solver.setParam(12, 0);
    }

    virtual void setSolution(const gsMatrix<T> & solVector) { m_solution = solVector; }

    virtual void setStokesSolution()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        gsInfo << "Setting Stokes initial condition...\n";

        #ifdef GISMO_WITH_PARDISO
        typename gsSparseSolver<T>::PardisoLU solver;
        pardisoSetup(solver);
        #else
        typename gsSparseSolver<T>::LU solver;
        #endif 

        gsSparseMatrix<T> stokesMatrix;
        gsMatrix<T> stokesRhs;

        getAssembler()->fillStokesSystem_into(stokesMatrix, stokesRhs);
        
        solver.analyzePattern(stokesMatrix);
        solver.factorize(stokesMatrix);
        m_solution = solver.solve(stokesRhs);

        m_iterationNumber = 0;
    }

    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION }

    virtual void nextIteration(const unsigned numberOfIterations)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        for (unsigned iter = 0; iter < numberOfIterations; iter++)
            nextIteration();
    }

    virtual void updateAssembler()
    {
        m_clock.restart();
        getAssembler()->update(m_solution);
        m_assembT += m_clock.stop();
    }

    virtual void initIteration()
    {
        m_clock.restart();
        m_solver.analyzePattern(getAssembler()->matrix());
        m_solsetupT += m_clock.stop();
    }

    virtual void applySolver(gsMatrix<T>& solution)
    {
        m_clock.restart();
        m_solver.factorize(getAssembler()->matrix());
        m_solsetupT += m_clock.stop();

        m_clock.restart();
        solution = m_solver.solve(getAssembler()->rhs());
        m_solveT += m_clock.stop();
    }

    virtual void applySolver(gsMatrix<T>& solution, real_t alpha_u, real_t alpha_p)
    {
        GISMO_ASSERT(solution.rows() > 0,"The solution in applySolver() is empty!");

        int usize = getAssembler()->getPshift();
        int pdofs = getAssembler()->getPdofs();

        m_clock.restart();
        m_solver.factorize(getAssembler()->matrix());
        m_solsetupT += m_clock.stop();

        m_clock.restart();
        gsMatrix<T> newsol = m_solver.solve(getAssembler()->rhs());
        solution.topRows(usize) = alpha_u * newsol.topRows(usize) + (1 - alpha_u)*solution.topRows(usize);
        solution.bottomRows(pdofs) = alpha_p * newsol.bottomRows(pdofs) + (1 - alpha_p)*solution.bottomRows(pdofs);
        m_solveT += m_clock.stop();
    }

    void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
        int iter = 0;
        m_relNorm = solutionChangeRelNorm();

        while ((iter < minIterations) || ((m_relNorm > epsilon) && (iter < maxIterations)))
        {
            gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

            nextIteration();
            m_relNorm = solutionChangeRelNorm();
            m_ofileRelNorm << m_relNorm << "\n";

            gsInfo << " Solution change relative norm: " << m_relNorm << "\n";

            iter++;
        }
        m_ofileRelNorm.close();
    }

    void solveWithObjective(uwbObjectiveFunctionEvaluator<T>& objFunction, const T tolObjRelVal = 1e-3, const T tolRelNorm = 1e-3, const int maxIterations = 10, const int minIterations = 0)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");
        GISMO_ASSERT(minIterations <= maxIterations, "Minimum number of iterations must be smaller or equal to maximum number of iterations.");

        openCloseFileObj();
        m_objFnVec.setZero(1, maxIterations);

        int iter = 0;
        m_objError = objFunError(objFunction);
        m_relNorm = solutionChangeRelNorm();

        while (((m_objError > tolObjRelVal) || (m_relNorm > tolRelNorm) || (iter < minIterations)) && (iter < maxIterations))
        {
            gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

            nextIteration();
            m_objError = objFunError(objFunction);
            m_ofile << m_objFunValue << "\n";
            m_relNorm = solutionChangeRelNorm();
            m_ofileRelNorm << m_relNorm << "\n";

            if (objFunction.isObjectiveAbsolute())
                gsInfo << " Relative error of the absolute objective function: " << m_objError << "\n";
            else
                gsInfo << " Absolute error of the relative objective function: " << m_objError << "\n";
            gsInfo << " Solution change relative norm: " << m_relNorm << "\n\n";

            iter++;
        }

        m_ofileRelNorm.close();
        openCloseFileObj(false);
    }

    T objFunError(uwbObjectiveFunctionEvaluator<T>& objFunction)
    {
        T objError = 0.;

        if (m_iterationNumber)
        {
            gsField<T> uSol = this->constructSolution(0);
            gsField<T> pSol = this->constructSolution(1);
            objFunction.updateSolution(uSol, pSol);

            T newObjFunValue = objFunction.evaluate();
            gsVector<T> newObjParts = objFunction.getObjectives();
            
            if (objFunction.isObjectiveAbsolute())
            {
                for (int s = 0; s < objFunction.getNumObjParts(); s++)
                {
                    T relValueObjPart = 0.;
                    if (newObjParts[s] < 0.)
                        relValueObjPart = (m_objParts[s] - newObjParts[s]) / math::max(-newObjParts[s], 1e-4);
                    else
                        relValueObjPart = (m_objParts[s] - newObjParts[s]) / math::max(newObjParts[s], 1e-4);
                    
                    if (relValueObjPart < 0.)
                        objError -= m_weights[s] * relValueObjPart;
                    else
                        objError += m_weights[s] * relValueObjPart;
                }
            }
            else
                objError = m_objFunValue - newObjFunValue;

            m_objFnVec(0, m_iterationNumber-1) = newObjFunValue;
            gsInfo << "\n Objective function value: " << newObjFunValue << "\n";

            m_objFunValue = newObjFunValue;
            m_objParts = newObjParts;
        }
        else
        {
            objError = std::numeric_limits<T>::infinity();
            m_objFunValue = objFunction.evaluate();
            m_objParts.setZero(objFunction.getNumObjParts());
            m_weights.setZero(objFunction.getNumObjParts());
            m_objParts = objFunction.getObjectives();
            m_weights = objFunction.getWeights();
        }
        if (objError < 0)
            objError = -objError;
        return objError;
    }

    void plotObjFunValue(std::string const & fileName)
    {
        gsMatrix<T> X, Y;
        X.setZero(1, m_iterationNumber);
        Y.setZero(1, m_iterationNumber);
        Y.middleCols(0, m_iterationNumber) = m_objFnVec.middleCols(0, m_iterationNumber);
        for (int i = 0; i < (int)m_iterationNumber; i++)
        {
            X.coeffRef(0, i) = i+1;
        }

        gsWriteParaviewPoints<>(X, Y, fileName);
    }

    T solutionChangeRelNorm() const
    {
        T relNorm;

        if (m_iterationNumber)
        {
            gsMatrix<T> solChangeVector = getAssembler()->getSolution() - m_solution;
            relNorm = solChangeVector.norm() / m_solution.norm();

        }
        else
        {
            relNorm = std::numeric_limits<T>::infinity();
        }

        return relNorm;
    }

    T solutionChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsMatrix<T> solChangeVector = solOld - solNew;
        T relNorm = solChangeVector.norm() / solNew.norm();

        return relNorm;
    }

    virtual void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     [u, p] solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

        gsInfo << "\n";
    }

    virtual T residualRelNorm() const { GISMO_NO_IMPLEMENTATION }

    void addPressureOutletCondition(int patch, boxSide side)
    {
        getAssembler()->addPressureOutletCondition(patch, side);
        reinitMembers();
    }

    virtual void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    {
        getAssembler()->markDofsAsEliminatedZeros(boundaryDofs, unk);
        reinitMembers();
    }

    gsField<T> constructSolution(int unk, bool relative = false) const
    {
        return getAssembler()->constructSolution(m_solution, unk, relative);
    }

    gsField<T> constructSolution(int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    {
        return getAssembler()->constructSolution(m_solution, unk, patchNumbers, sides, relative);
    }

    gsField<T> constructSolutionCombined(gsVector<size_t> relPatches) const
    {
        return getAssembler()->constructSolutionCombined(m_solution, relPatches);
    }

    T computeDimensionlessWallDistance(gsVector<int> distancePatches, std::vector<boxSide> distanceSides, T viscosity, T reynoldsNumber, T uFreeStream, T maxYplus = 1.0, unsigned npts = 20, bool print = false, bool estimate = true) const
    {
        return getAssembler()->computeDimensionlessWallDistance(m_solution, distancePatches, distanceSides, viscosity, reynoldsNumber, uFreeStream, maxYplus, npts, print, estimate);
    }

    T computeAspectRatio(bool minAR = false) { return getAssembler()->computeAspectRatio(minAR); }

    void plot2DVorticity(std::string const & fn, unsigned npts = 10000)
    {
        gsInfo << "Plotting vorticity from 2D velocity field in Paraview...\n";
        getAssembler()->plot2DVorticity(fn, m_solution, npts);
    }

    void plot2DVorticity(std::string const & fn, gsMatrix<T>& solution, unsigned npts = 10000)
    {
        gsInfo << "Plotting vorticity from 2D velocity field in Paraview...\n";
        getAssembler()->plot2DVorticity(fn, solution, npts);
    }

    void plotPressureCoefficient(std::string const & fn,
                                   index_t referencePatchIndex, gsVector<T>& referencePoint, T freeStreamVelocity,
                                   T density, unsigned npts = 10000)
    {
        gsInfo << "Plotting pressure coefficient from 2D pressure and velocity fields in Paraview...\n";
        getAssembler()->plotPressureCoefficient(fn, m_solution, referencePatchIndex, referencePoint, freeStreamVelocity, density, npts);
    }

    void plotPressureCoefficient(std::string const & fn, gsMatrix<T>& solution,
                                   index_t referencePatchIndex, gsVector<T>& referencePoint, T freeStreamVelocity,
                                   T density, unsigned npts = 10000)
    {
        gsInfo << "Plotting vorticity from 2D pressure and velocity fields in Paraview...\n";
        getAssembler()->plotPressureCoefficient(fn, solution, referencePatchIndex, referencePoint, freeStreamVelocity, density, npts);
    }

    virtual void evalResiduum(const gsMatrix<T> & solVector, std::vector<T> & residuum)
    { GISMO_NO_IMPLEMENTATION }

    virtual void evalResiduum(std::vector<T> & residuum)
    { GISMO_NO_IMPLEMENTATION }

    virtual uwbINSAssemblerBase<T>*    getAssembler() const { return m_pAssembler; }
    const gsMatrix<T> &        getSolution() const { return m_solution; }
    unsigned             getIterationNumber() const { return m_iterationNumber; }
    T                    getObjError() const { return m_objError; }
    T                    getRelNorm() const { return m_relNorm; }
    int                  numDofs() const { return getAssembler()->numDofs(); }
    virtual const T getAssemblyTime() const { return m_assembT; }
    virtual const T getSolverSetupTime() const { return m_solsetupT; }
    virtual const T getSolveTime() const { return m_solveT; }

protected:
    uwbINSAssemblerBase<T>*    m_pAssembler;
    gsMatrix<T>                m_solution;
    unsigned                   m_iterationNumber;
    bool                       m_bPatternAnalyzed;
    gsMatrix<T>                m_objFnVec;
    T                          m_objFunValue;
    gsVector<T>                m_objParts;
    gsVector<T>                m_weights;
    T                          m_objError;
    T                          m_relNorm;

#ifdef GISMO_WITH_PARDISO
                typename gsSparseSolver<T>::PardisoLU m_solver;
#else
                typename gsSparseSolver<T>::LU m_solver;
#endif

    std::ofstream              m_ofile;
    std::ofstream              m_ofileRelNorm;

    gsStopwatch m_clock;
    T m_assembT, m_solsetupT, m_solveT;
    
}; //uwbINSSolverBase

} //namespace gismo
