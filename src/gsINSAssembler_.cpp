#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSAssemblerBase1<real_t>;
    CLASS_TEMPLATE_INST gsINSAssembler<real_t>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady1<real_t>;

} // namespace gismo
