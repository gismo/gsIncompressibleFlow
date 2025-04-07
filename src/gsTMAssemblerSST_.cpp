#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMAssemblerSST.h>
#include <gsIncompressibleFlow/src/gsTMAssemblerSST.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMAssemblerSST<real_t, RowMajor>;
    
    CLASS_TEMPLATE_INST gsTMAssemblerSST<real_t, ColMajor>;
    
} // namespace gismo