#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowLinSystSolver<real_t>;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_direct<real_t>;

} // namespace gismo