/** @file uwbADRSolverParams.h

Author(s): E. Turnerova
*/

#pragma once

#include <gismo.h>
#include "uwbADRProblemSettings.h"
#include "uwbADRPde.h"


namespace gismo
{

template<class T>
class uwbADRSolverParams
{

public:
    uwbADRSolverParams(uwbADRPde<T>  &      pde,
        gsMultiBasis<T> &                   bases,
        gsAssemblerOptions                  opt,
        int                                 numThreads = 0) :
        m_pde(pde), m_bases(bases), m_assemblerOptions(opt)
    {
        if (numThreads)
        {
            m_numThreads = numThreads;
        }
        else
        {
            m_numThreads =
            #ifdef _OPENMP 
                std::round(omp_get_max_threads() / 2);
            #else
                1;
            #endif
        }
    }

    uwbADRSolverParams(uwbADRPde<T>  &      pde,
        gsMultiBasis<T> &                   bases,
        gsAssemblerOptions                  opt,
        uwbADRProblemSettings<T>&           probSettings,
        int                                 numThreads = 0) :
        uwbADRSolverParams(pde, bases, opt, numThreads)
    {      
        m_problemSettings = probSettings;
    }

    ~uwbADRSolverParams() { }

public:
    
    const uwbADRPde<T>& getPde() const { return m_pde; }

    const gsBoundaryConditions<T>& getBCs() const { return m_pde.bc(); }

    const gsMultiBasis<T> & getBases() const { return m_bases; }
    gsMultiBasis<T> &       getBases() { return m_bases; }

    const gsAssemblerOptions& getAssemblerOptions() const { return m_assemblerOptions; }
    gsAssemblerOptions& getAssemblerOptions() { return m_assemblerOptions; }
    
    const uwbADRProblemSettings<T>& settings() const { return m_problemSettings; }
    uwbADRProblemSettings<T>&       settings() { return m_problemSettings; }

    int getNumThreads() const { return m_numThreads; }
    
    void setNumThreads(int numThreads)
    {
        m_numThreads = numThreads;
    }

protected:
    //members
    uwbADRPde<T>                        m_pde;
    gsMultiBasis<T>                     m_bases;
    gsAssemblerOptions                  m_assemblerOptions;
    uwbADRProblemSettings<T>            m_problemSettings;
    int                                 m_numThreads;
}; // class uwbADRSolverParams

} // namespace gismo
