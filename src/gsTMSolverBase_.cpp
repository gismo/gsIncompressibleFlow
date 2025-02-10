#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMSolverBase.h>
//#include <gsIncompressibleFlow/src/gsTMSolverBase.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMSolverBase<real_t, RowMajor>;
    
    CLASS_TEMPLATE_INST gsTMSolverBase<real_t, ColMajor>;
    
} // namespace gismo