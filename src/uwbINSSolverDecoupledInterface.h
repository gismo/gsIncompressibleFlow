/** @file uwbINSSolverDecoupledInterface.h

Author(s): H. Hornikova
*/

#pragma once

#include "uwbINSSolverDecoupled1.h"
#include "uwbINSSolverDecoupled1PeriodicIterative.h"
#include "uwbINSSolverDecoupled1PeriodicCoupled.h"
#include "uwbINSSolverDecoupled2.h"
#include "uwbINSSolverDecoupled2PeriodicIterative.h"
#include "uwbINSSolverDecoupled2PeriodicCoupled.h"

namespace gismo
{

template<class T>
class uwbINSSolverDecoupledInterface
{

public:
    uwbINSSolverDecoupledInterface(uwbINSSolverParams<T>& params)
    {
        decoupled::projection proj = params.settings().getProjVersion();
        decoupled::method method = params.settings().getDecoupledMethod();

        if (!(params.getBCs().numPeriodic()))
        {
            switch (proj)
            {
            case decoupled::proj1:
                m_pSolver = new uwbINSSolverDecoupled1<T>(params);
                break;
            case decoupled::proj2:
                m_pSolver = new uwbINSSolverDecoupled2<T>(params);
                break;
            default:
                gsInfo << "Projection version invalid, choosing version 1.\n";
                m_pSolver = new uwbINSSolverDecoupled1<T>(params);
                break;
            }
        }
        else
        {
            switch (method)
            {
            case decoupled::iterative:
                gsInfo << "Using iterative method for rotation/periodic conditions.\n";
                switch (proj)
                {
                case decoupled::proj1:
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicIterative<T>(params);
                    break;
                case decoupled::proj2:
                    m_pSolver = new uwbINSSolverDecoupled2PeriodicIterative<T>(params);
                    break;
                default:
                    gsInfo << "Projection version invalid, choosing version 1.\n";
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicIterative<T>(params);
                    break;
                }
                break;
            case decoupled::coupled:
                gsInfo << "Using coupled method for rotation/periodic conditions.\n";
                switch (proj)
                {
                case decoupled::proj1:
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicCoupled<T>(params);
                    break;
                case decoupled::proj2:
                    m_pSolver = new uwbINSSolverDecoupled2PeriodicCoupled<T>(params);
                    break;
                default:
                    gsInfo << "Projection version invalid, choosing version 1.\n";
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicCoupled<T>(params);
                    break;
                }
                break;
            default:
                gsInfo << "Method for rotation/periodic conditions set to 'none' or invalid, choosing iterative with 1 inner iteration. \n";
                switch (proj)
                {
                case decoupled::proj1:
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicIterative<T>(params);
                    break;
                case decoupled::proj2:
                    m_pSolver = new uwbINSSolverDecoupled2PeriodicIterative<T>(params);
                    break;
                default:
                    gsInfo << "Projection version invalid, choosing version 1.\n";
                    m_pSolver = new uwbINSSolverDecoupled1PeriodicIterative<T>(params);
                    break;
                }
                m_pSolver->setInnerIterations(0);
                break;
            }
        }
    }

    ~uwbINSSolverDecoupledInterface()
    {
        delete m_pSolver;
        m_pSolver = NULL;
    }

    void initialize()
    { m_pSolver->initialize(); }

    void setSolution(const gsMatrix<T> & solVector)
    { m_pSolver->setSolution(solVector); }

    void setStokesSolution()
    { m_pSolver->setStokesSolution(); }

    void nextIteration()
    { m_pSolver->nextIteration(); }

    void nextIteration(const unsigned numberOfIterations)
    { m_pSolver->nextIteration(numberOfIterations); }

    void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
    { m_pSolver->solve(maxIterations, epsilon, minIterations); }

    T solutionChangeRelNorm() const
    { return m_pSolver->solutionChangeRelNorm(); }

    void addPressureOutletCondition(int patch, boxSide side)
    { m_pSolver->addPressureOutletCondition(patch, side); }

    void markDofsAsEliminatedZeros(const std::vector< gsMatrix< index_t > > & boundaryDofs, const int unk = 0)
    { m_pSolver->markDofsAsEliminatedZeros(boundaryDofs, unk); }

    gsField<T> constructSolution(const int unk, bool relative = false) const
    { return m_pSolver->constructSolution(unk, relative); }

    gsField<T> constructSolution(int unk, gsVector<index_t> patchNumbers, std::vector<boxSide> sides, bool relative = false) const
    { return m_pSolver->constructSolution(unk, patchNumbers, sides, relative); }

    gsField<T> constructSolutionCombined(gsVector<size_t> relPatches) const
    { return m_pSolver->constructSolutionCombined(relPatches); }

    uwbINSAssemblerBase<T>* getAssembler() const
    { return m_pSolver->getAssembler(); }

    const gsMatrix<T> & getSolution() const
    { return m_pSolver->getSolution(); }

    unsigned getIterationNumber() const
    { return m_pSolver->getIterationNumber(); }

    int numDofs() const
    { return m_pSolver->numDofs(); }

    uwbINSSolverDecoupledBase<T> * getSolver()
    { return m_pSolver; }

protected:
    uwbINSSolverDecoupledBase<T> * m_pSolver;

}; //uwbINSSolverDecoupledInterface

} //namespace gismo