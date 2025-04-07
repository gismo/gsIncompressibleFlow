#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsTMAssemblerBase.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMAssemblerBase<real_t, RowMajor>;
    
    CLASS_TEMPLATE_INST gsTMAssemblerBase<real_t, ColMajor>;
    
} // namespace gismo