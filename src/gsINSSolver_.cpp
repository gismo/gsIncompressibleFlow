#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolver.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSSolver<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSSolverSteady<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSSolverUnsteady<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSSolver<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSSolverSteady<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSSolverUnsteady<real_t, ColMajor>;
    
} // namespace gismo
