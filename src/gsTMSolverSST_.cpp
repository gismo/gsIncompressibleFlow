#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMSolverSST.h>
#include <gsIncompressibleFlow/src/gsTMSolverSST.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMSolverSST<real_t, RowMajor>;
    
    CLASS_TEMPLATE_INST gsTMSolverSST<real_t, ColMajor>;
    
} // namespace gismo