#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSTermPvalUdiv<real_t>;
    CLASS_TEMPLATE_INST gsINSTermUdivPval<real_t>;
    CLASS_TEMPLATE_INST gsINSTermUsolGradVal<real_t>;
    
} // namespace gismo