#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsRANSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsRANSSolverUnsteady.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsRANSSolverUnsteady<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsRANSSolverUnsteady<real_t, ColMajor>;
    
} // namespace gismo