#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSAssembler<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerSteady<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSAssembler<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerSteady<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady<real_t, ColMajor>;

} // namespace gismo
