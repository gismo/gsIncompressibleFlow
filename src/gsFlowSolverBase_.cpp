#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowSolverBase.h>
#include <gsIncompressibleFlow/src/gsFlowSolverBase.hpp>

namespace gismo
{
    
    CLASS_TEMPLATE_INST gsFlowSolverBase<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowSolverBase<real_t, ColMajor>;

} // namespace gismo