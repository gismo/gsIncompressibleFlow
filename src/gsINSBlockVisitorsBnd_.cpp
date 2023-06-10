#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSBlockVisitorsBnd.h>
#include <gsIncompressibleFlow/src/gsINSBlockVisitorsBnd.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSBlockVisitorBnd<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockVisitorRobinPCD<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockVisitorNitsche<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockAVisitorNitsche<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockBVisitorNitsche<real_t>;

} // namespace gismo
