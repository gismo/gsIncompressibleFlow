/** @file uwbTMSolverKOmega.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMSolverBase.h"
#include "uwbTMAssemblerKOmega.h"
#include "uwbTMSUPGAssemblerKOmega.h"
#include "uwbTMCrosswindAssemblerKOmega.h"
#include "uwbTMADAssemblerKOmega.h"
#include "uwbTMisoADAssemblerKOmega.h"
#include "uwbTMAssemblerKOmega_AFC.h"
#include "uwbTMtanhCSDassemblerKOmega.h"

namespace gismo
{

template<class T>
class uwbTMSolverKOmega : public uwbTMSolverBaseUnsteady<T>
{

public:
    typedef uwbTMSolverBaseUnsteady<T> Base;

public:
    uwbTMSolverKOmega(uwbINSSolverParams<T>& params) : Base(params)
    {
        //create assembler
        if (params.settings().get(constantsINS::TMafc) || params.settings().get(constantsINS::TMafcHO))
        {
            params.getAssemblerOptions().dirStrategy = dirichlet::none;
            m_pAssembler = new uwbTMAssemblerKOmega_AFC<T>(params);
        }
        else if (params.settings().get(constantsINS::TMcrosswind))
            m_pAssembler = new uwbTMCrosswindAssemblerKOmega<T>(params);
        else if (params.settings().get(constantsINS::TMsupg))
            m_pAssembler = new uwbTMSUPGAssemblerKOmega<T>(params);
        else if (params.settings().get(constantsINS::TMad))
            m_pAssembler = new uwbTMADAssemblerKOmega<T>(params);
        else if (params.settings().get(constantsINS::TMisoAD))
            m_pAssembler = new uwbTMisoADAssemblerKOmega<T>(params);
        else if (params.settings().get(constantsINS::tanhCSD))
            m_pAssembler = new uwbTMtanhCSDassemblerKOmega<T>(params);
        else
            m_pAssembler = new uwbTMAssemblerKOmega<T>(params);

        #ifdef GISMO_WITH_PARDISO
        uwbTMSolverBase<T>::pardisoSetup(m_kSolver);
        uwbTMSolverBase<T>::pardisoSetup(m_oSolver);
        #endif
    }

    virtual ~uwbTMSolverKOmega()
    {
    }

public:

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        //if(Base::isTMSUPG())
        //    getAssembler()->updateTurbCoeffSolField(m_solution, m_uSolField);

        m_clock.restart();
        getAssembler()->update(m_solution, m_uSolField);
        m_assembT += m_clock.stop();

        // in the first time step the pattern will get analyzed
        if (!m_iterationNumber)
        {
            m_clock.restart();
            m_kSolver.analyzePattern(getAssembler()->matrix(0));
            m_oSolver.analyzePattern(getAssembler()->matrix(1));
            m_solsetupT += m_clock.stop();
        }

        gsMatrix<T> tmpSolution(getAssembler()->numVarDofs(), getAssembler()->getNumVar());

        m_clock.restart();
        m_kSolver.factorize(getAssembler()->matrix(0));
        m_oSolver.factorize(getAssembler()->matrix(1));
        m_solsetupT += m_clock.stop();

        m_clock.restart();
        tmpSolution.col(0) = m_kSolver.solve(getAssembler()->rhs(0));
        tmpSolution.col(1) = m_oSolver.solve(getAssembler()->rhs(1));
        m_solveT += m_clock.stop();

        //------- positive turbulent variables -----
        tmpSolution = tmpSolution.cwiseMax(1e-10);
        //------------------------------------------

        uwbTMSolverBase<T>::dispSolChangeRelNorm(m_solution, tmpSolution);

        int iter = 0;
        T kRelNorm = this->solutionChangeRelNorm(m_solution.col(0), tmpSolution.col(0));
        T oRelNorm = this->solutionChangeRelNorm(m_solution.col(1), tmpSolution.col(1));

        int maxPicardIter = m_innerIter;
        if (!m_iterationNumber)
            maxPicardIter = m_innerFirstIter;
        gsInfo << "        TM Picard's iterations...\n";
        while ((math::max(kRelNorm, oRelNorm) > m_innerTol) && (iter < maxPicardIter))
        {
            gsInfo << "         ";
            //if(Base::isTMSUPG())
            //    getAssembler()->updateTurbCoeffSolField(tmpSolution, m_uSolField);

            m_clock.restart();
            getAssembler()->updatePicard(tmpSolution);
            m_assembT += m_clock.stop();

            gsMatrix<T> oldSol = tmpSolution;

            m_clock.restart();
            m_kSolver.factorize(getAssembler()->matrix(0));
            m_oSolver.factorize(getAssembler()->matrix(1));
            m_solsetupT += m_clock.stop();

            m_clock.restart();
            tmpSolution.col(0) = m_kSolver.solve(getAssembler()->rhs(0));
            tmpSolution.col(1) = m_oSolver.solve(getAssembler()->rhs(1));
            m_solveT += m_clock.stop();
       
            //------- positive turbulent variables -----
            tmpSolution = tmpSolution.cwiseMax(1e-10);
            //------------------------------------------
       
            this->dispSolChangeRelNorm(oldSol, tmpSolution);

            kRelNorm = this->solutionChangeRelNorm(oldSol.col(0), tmpSolution.col(0));
            oRelNorm = this->solutionChangeRelNorm(oldSol.col(1), tmpSolution.col(1));
            iter++;
        }

        m_solution = tmpSolution;

        m_iterationNumber++;
    }

    void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
    {
        GISMO_ASSERT(m_pAssembler->isInitialized(), "Assembler must be initialized first, call initialize()");
        int iter = 0;
        T kRelNorm = this->solutionChangeRelNorm(0);
        T oRelNorm = this->solutionChangeRelNorm(1);

        while ((iter < minIterations) || ((kRelNorm > epsilon || oRelNorm > epsilon) && (iter < maxIterations)))
        {
            gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

            nextIteration();

            kRelNorm = this->solutionChangeRelNorm(0);
            oRelNorm = this->solutionChangeRelNorm(1);

            gsInfo << " Solution change relative norm for k: " << kRelNorm << ", for omega: " << oRelNorm << "\n";

            iter++;
        }
    }

    void dispSolChangeRelNorm(gsMatrix<T> solOld, gsMatrix<T> solNew) const
    {
        gsInfo << "     TM Picard's solution change relative norm: ";

        for (int i = 0; i < solOld.cols(); i++)
            gsInfo << this->solutionChangeRelNorm(solOld.col(i), solNew.col(i)) << ", ";

        gsInfo << "\n";
    }

public:
    virtual uwbTMAssemblerKOmega<T>* getAssembler() const
    {
        uwbTMAssemblerKOmega_AFC<T>* pAssembler = dynamic_cast<uwbTMAssemblerKOmega_AFC<T>*>(m_pAssembler);
        uwbTMSUPGAssemblerKOmega<T>* pSUPGassembler = dynamic_cast<uwbTMSUPGAssemblerKOmega<T>*>(m_pAssembler);
        uwbTMCrosswindAssemblerKOmega<T>* pCrosswindAssembler = dynamic_cast<uwbTMCrosswindAssemblerKOmega<T>*>(m_pAssembler);
        uwbTMADAssemblerKOmega<T>* pADassembler = dynamic_cast<uwbTMADAssemblerKOmega<T>*>(m_pAssembler);
        uwbTMisoADAssemblerKOmega<T>* pIsoADassembler = dynamic_cast<uwbTMisoADAssemblerKOmega<T>*>(m_pAssembler);
        uwbTMtanhCSDassemblerKOmega<T>* pTanhCSDassembler = dynamic_cast<uwbTMtanhCSDassemblerKOmega<T>*>(m_pAssembler);

        if (pAssembler != NULL)
            return pAssembler;
        else if (pCrosswindAssembler != NULL)
            return dynamic_cast<uwbTMCrosswindAssemblerKOmega<T>*>(m_pAssembler);
        else if (pSUPGassembler != NULL)
            return dynamic_cast<uwbTMSUPGAssemblerKOmega<T>*>(m_pAssembler);
        else if (pADassembler != NULL)
            return dynamic_cast<uwbTMADAssemblerKOmega<T>*>(m_pAssembler);
        else if (pIsoADassembler != NULL)
            return dynamic_cast<uwbTMisoADAssemblerKOmega<T>*>(m_pAssembler);
        else if (pTanhCSDassembler != NULL)
            return dynamic_cast<uwbTMtanhCSDassemblerKOmega<T>*>(m_pAssembler);
        else
            return dynamic_cast<uwbTMAssemblerKOmega<T>*>(m_pAssembler);
    }

protected:

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_kSolver, m_oSolver;
    #else
    typename gsSparseSolver<T>::LU m_kSolver, m_oSolver;
    #endif 

    uwbTMEvaluatorKOmega<T>* m_pEvaluator;

    // members from uwbTMSolverBaseUnsteady
    using Base::m_innerIter;
    using Base::m_innerFirstIter;
    using Base::m_innerTol;

    // members from uwbTMSolverBase
    using uwbTMSolverBase<T>::m_pAssembler;
    using uwbTMSolverBase<T>::m_solution;
    using uwbTMSolverBase<T>::m_uSolField;
    using uwbTMSolverBase<T>::m_iterationNumber;
    using uwbTMSolverBase<T>::m_clock;
    using uwbTMSolverBase<T>::m_assembT;
    using uwbTMSolverBase<T>::m_solsetupT;
    using uwbTMSolverBase<T>::m_solveT;

}; // class uwbTMSolverKOmega

} // namespace gismo
