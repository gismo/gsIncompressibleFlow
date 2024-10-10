#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.hpp>

namespace gismo
{
    
    CLASS_TEMPLATE_INST gsFlowAssemblerBase<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowAssemblerBase<real_t, ColMajor>;

} // namespace gismo
