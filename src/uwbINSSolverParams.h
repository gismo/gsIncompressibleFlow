/** @file uwbINSSolverParams.h

Author(s): H. Hornikova
*/

#pragma once

#include <gismo.h>
#include "uwbINSProblemSettings.h"
#include "uwbINSPde.h"
#include "uwbPreconditioners.h"

namespace gismo
{

template<class T>
class uwbINSSolverParams
{

public:
    uwbINSSolverParams(uwbINSPde<T>  &        pde,
        std::vector<gsMultiBasis<T> >&       bases,
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

        m_precOpt = uwbINSPreconditioner<real_t>::defaultOptions();
    }

    uwbINSSolverParams(uwbINSPde<T>  &        pde,
        std::vector<gsMultiBasis<T> >&       bases,
        gsAssemblerOptions                  opt,
        uwbINSProblemSettings<T>&           probSettings,
        int                                 numThreads = 0) :
        uwbINSSolverParams(pde, bases, opt, numThreads)
    {      
        m_problemSettings = probSettings;
    }

    ~uwbINSSolverParams()
    {
    }

public:
    
    const uwbINSPde<T>& getPde() const { return m_pde; }

    const gsBoundaryConditions<T>& getBCs() const { return m_pde.bc(); }

    const std::vector<gsMultiBasis<T> > & getBases() const { return m_bases; }
    std::vector<gsMultiBasis<T> > &       getBases() { return m_bases; }

    const gsAssemblerOptions& getAssemblerOptions() const { return m_assemblerOptions; }
    gsAssemblerOptions& getAssemblerOptions() { return m_assemblerOptions; }
    
    const uwbINSProblemSettings<T>& settings() const { return m_problemSettings; }
    uwbINSProblemSettings<T>&       settings() { return m_problemSettings; }

    int getNumThreads() const { return m_numThreads; }
    
    void setNumThreads(int numThreads)
    {
        m_numThreads = numThreads;
    }

    void setPrecOptions(gsOptionList& opt) { m_precOpt = opt; }
    const gsOptionList& getPrecOptions() const { return m_precOpt; }
    gsOptionList& getPrecOptions() { return m_precOpt; }

protected:
    //members
    uwbINSPde<T>                        m_pde;
    std::vector<gsMultiBasis<T> >       m_bases;
    gsAssemblerOptions                  m_assemblerOptions;
    uwbINSProblemSettings<T>            m_problemSettings;
    gsOptionList                        m_precOpt;
    int                                 m_numThreads;
}; // class uwbINSSolverParams

} // namespace gismo
