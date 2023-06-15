#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.h>
#include <gsIncompressibleFlow/src/gsINSSolverUnsteady.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSSolverUnsteady<real_t>;
    CLASS_TEMPLATE_INST gsINSSolverUnsteadyIter<real_t, gsGMRes<real_t> >;
} 
