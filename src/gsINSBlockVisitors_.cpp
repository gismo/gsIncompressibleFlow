#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSBlockVisitors.h>
#include <gsIncompressibleFlow/src/gsINSBlockVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSVisitorBase<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockAVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlocksBVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlocksCVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockNVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockMVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockApVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockNpVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockMpVisitor<real_t>;
    CLASS_TEMPLATE_INST gsINSRhsVisitor<real_t>;

} // namespace gismo

