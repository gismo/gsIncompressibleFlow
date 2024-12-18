#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowVisitors.h>
#include <gsIncompressibleFlow/src/gsFlowVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowVisitor<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowVisitorVectorValued<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsFlowVisitor<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsFlowVisitorVectorValued<real_t, ColMajor>;

} // namespace gismo