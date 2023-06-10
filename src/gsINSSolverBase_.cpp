#include <gsCore/gsTemplateTools.h>

#include <gsSolver/gsGMRes.h>
#include <gsIncompressibleFlow/src/gsINSSolverBase.h>
#include <gsIncompressibleFlow/src/gsINSSolverBase.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSSolverBase<real_t>;
    CLASS_TEMPLATE_INST gsINSSolverBaseIter<real_t, gsGMRes<real_t> >;
} 
