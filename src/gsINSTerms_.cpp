#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSTerm_PvalUdiv<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_UdivPval<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_UsolGradVal<real_t>;
    
} // namespace gismo