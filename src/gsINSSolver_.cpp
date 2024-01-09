#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolver.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSSolver<real_t>;
    CLASS_TEMPLATE_INST gsINSSolverDirect<real_t>;
    CLASS_TEMPLATE_INST gsINSSolverDirectSteady<real_t>;
    CLASS_TEMPLATE_INST gsINSSolverDirectUnsteady<real_t>;


} // namespace gismo
