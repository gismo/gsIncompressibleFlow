/** @file uwbTMSolverKOmegaLinSteady.h

Author(s): H. Hornikova, E. Turnerova
*/

#pragma once

#include "uwbTMSolverBase.h"
#include "uwbTMAssemblerKOmegaLinSteady.h"
#include "uwbTMSUPGAssemblerKOmegaLinSteady.h"
#include "uwbTMCrosswindAssemblerKOmegaLinSteady.h"
#include "uwbTMADAssemblerKOmegaLinSteady.h"
#include "uwbTMisoADAssemblerKOmegaLinSteady.h"

namespace gismo
{

template<class T>
class uwbTMSolverKOmegaLinSteady : public uwbTMSolverBase<T>
{

public:
    typedef uwbTMSolverBase<T> Base;

public:
    uwbTMSolverKOmegaLinSteady(uwbINSSolverParams<T>& params) : Base()
    {
        //create assembler
        if (params.settings().get(constantsINS::TMcrosswind))
            m_pAssembler = new uwbTMCrosswindAssemblerKOmegaLinSteady<T>(params);
        else if (params.settings().get(constantsINS::TMsupg))
            m_pAssembler = new uwbTMSUPGAssemblerKOmegaLinSteady<T>(params);
        else if (params.settings().get(constantsINS::TMad))
            m_pAssembler = new uwbTMADAssemblerKOmegaLinSteady<T>(params);
        else if (params.settings().get(constantsINS::TMisoAD))
            m_pAssembler = new uwbTMisoADAssemblerKOmegaLinSteady<T>(params);
        else
            m_pAssembler = new uwbTMAssemblerKOmegaLinSteady<T>(params);

        #ifdef GISMO_WITH_PARDISO
        Base::pardisoSetup(m_kSolver);
        Base::pardisoSetup(m_oSolver);
        #endif
    }

    virtual ~uwbTMSolverKOmegaLinSteady()
    {
    }

public:

    virtual void nextIteration()
    {
        GISMO_ASSERT(getAssembler()->isInitialized(), "Assembler must be initialized first, call initialize()");

        getAssembler()->update(m_solution);

        // in the first time step the pattern will get analyzed
        if (!m_iterationNumber)
        {
            m_kSolver.analyzePattern(getAssembler()->matrix(0));
            m_oSolver.analyzePattern(getAssembler()->matrix(1));
        }

        m_kSolver.factorize(getAssembler()->matrix(0));
        m_solution.col(0) = m_kSolver.solve(getAssembler()->rhs(0));

        m_oSolver.factorize(getAssembler()->matrix(1));
        m_solution.col(1) = m_oSolver.solve(getAssembler()->rhs(1));
        //------- positive turbulent variables -----
        m_solution = m_solution.cwiseMax(1e-10);
        //------------------------------------------

        m_iterationNumber++;
    }

    void solve(const int maxIterations = 10, const T epsilon = 1e-3, const int minIterations = 0)
    {
        GISMO_ASSERT(m_pAssembler->isInitialized(), "Assembler must be initialized first, call initialize()");
        int iter = 0;
        T kRelNorm = solutionChangeRelNorm(0);
        T oRelNorm = solutionChangeRelNorm(1);

        while ((iter < minIterations) || ((kRelNorm > epsilon || oRelNorm > epsilon) && (iter < maxIterations)))
        {
            gsInfo << "Iteration number " << m_iterationNumber + 1 << "...";

            nextIteration();

            kRelNorm = solutionChangeRelNorm(0);
            oRelNorm = solutionChangeRelNorm(1);

            gsInfo << " Solution change relative norm for k: " << kRelNorm << ", for omega: " << oRelNorm << "\n";

            iter++;
        }
    }

    virtual uwbTMAssemblerKOmegaLinSteady<T>* getAssembler() const
    {
        uwbTMSUPGAssemblerKOmegaLinSteady<T>* pSUPGassembler = dynamic_cast<uwbTMSUPGAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        uwbTMCrosswindAssemblerKOmegaLinSteady<T>* pCrosswindAssembler = dynamic_cast<uwbTMCrosswindAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        uwbTMADAssemblerKOmegaLinSteady<T>* pADassembler = dynamic_cast<uwbTMADAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        uwbTMisoADAssemblerKOmegaLinSteady<T>* pIsoADassembler = dynamic_cast<uwbTMisoADAssemblerKOmegaLinSteady<T>*>(m_pAssembler);

        if (pCrosswindAssembler != NULL)
            return dynamic_cast<uwbTMCrosswindAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        else if (pSUPGassembler != NULL)
            return dynamic_cast<uwbTMSUPGAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        else if (pADassembler != NULL)
            return dynamic_cast<uwbTMADAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        else if (pIsoADassembler != NULL)
            return dynamic_cast<uwbTMisoADAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
        else
            return dynamic_cast<uwbTMAssemblerKOmegaLinSteady<T>*>(m_pAssembler);
    }

protected:

    #ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU m_kSolver, m_oSolver;
    #else
    typename gsSparseSolver<T>::LU m_kSolver, m_oSolver;
    #endif 

    // members from uwbTMSolverBase
    using Base::m_pAssembler;
    using Base::m_solution;
    using Base::m_iterationNumber;

    // member functions from uwbTMSolverBase
    using Base::solutionChangeRelNorm;

}; // class uwbTMSolverKOmegaLinSteady

} // namespace gismo
