#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMModels.h>
#include <gsIncompressibleFlow/src/gsTMModels.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMModelData<real_t>;
    CLASS_TEMPLATE_INST gsTMModelData_SST<real_t>;
    
} // namespace gismo