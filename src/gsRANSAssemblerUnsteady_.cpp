#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.h>
#include <gsIncompressibleFlow/src/gsRANSAssemblerUnsteady.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsRANSAssemblerUnsteady<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsRANSAssemblerUnsteady<real_t, ColMajor>;
    
} // namespace gismo