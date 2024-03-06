#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSAssembler<real_t>;
    CLASS_TEMPLATE_INST gsINSAssemblerSteady<real_t>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady<real_t>;

} // namespace gismo
