/** @file uwbADRSolverBase.h

Author(s): E. Turnerova
*/

#pragma once

#include "uwbADRAssemblerBase.h"
#include "uwbADRSolverParams.h"

namespace gismo
{

template<class T>
class uwbADRSolverBase
{
public:
    uwbADRSolverBase()
    {
        m_pAssembler = NULL;
    }

    virtual ~uwbADRSolverBase()
    {
        if (m_pAssembler)
        {
            delete m_pAssembler;
            m_pAssembler = NULL;
            m_ofileRelNorm.close();
        }
    }

protected:
    virtual void initMembers()
    {
        m_solution.setZero(getAssembler()->numDofs(), 1);
        m_iterationNumber = 0;
        m_relNorm = std::numeric_limits<T>::infinity();
        m_ofileRelNorm.open("solChangeRelNorm.txt");

        #ifdef GISMO_WITH_PARDISO
            pardisoSetup(m_solver);
        #endif

    }

    virtual void reinitMembers() { initMembers(); }

    void pardisoSetup(typename gsSparseSolver<T>::PardisoLU& solver)
    {
        solver.setParam(7, 10);
        solver.setParam(12, 0);
    }

public:
    virtual void initialize() { getAssembler()->initialize(); }
    virtual void setSolution(const gsMatrix<T> & solVector) { m_solution = solVector; }
    virtual void nextIteration() { GISMO_NO_IMPLEMENTATION }

    virtual void nextIteration(const unsigned numberOfIterations)
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        for (unsigned iter = 0; iter < numberOfIterations; iter++)
            nextIteration();
    }

    virtual void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
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

    void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     Solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

        gsInfo << "\n";
    }

    virtual T residualRelNorm() const { GISMO_NO_IMPLEMENTATION }

    gsField<T> constructSolution() const
    {
        return getAssembler()->constructSolution(m_solution);
    }

    virtual uwbADRAssemblerSteady<T>*    getAssembler() const { return m_pAssembler; }
    const gsMatrix<T> &        getSolution() const { return m_solution; }
    unsigned             getIterationNumber() const { return m_iterationNumber; }
    T                    getRelNorm() const { return m_relNorm; }
    int                  numDofs() const { return getAssembler()->numDofs(); }

protected:
    uwbADRAssemblerSteady<T>*    m_pAssembler;
    gsMatrix<T>                m_solution;
    unsigned                   m_iterationNumber;
    T                          m_relNorm;
    std::ofstream              m_ofileRelNorm;

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_solver;
    #else
    typename gsSparseSolver<T>::LU m_solver;
    #endif
    
}; //uwbADRSolverBase

} //namespace gismo
