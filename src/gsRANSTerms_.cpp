#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsRANSTerms.h>
#include <gsIncompressibleFlow/src/gsRANSTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsRANSTerm_SymmetricGradient<real_t>;
    CLASS_TEMPLATE_INST gsRANSTerm_SymmetricGradient_full<real_t>;
    //CLASS_TEMPLATE_INST gsRANSTerm_SymmetricGradientDiag<real_t>;
    //CLASS_TEMPLATE_INST gsRANSTerm_SymmetricGradientOffdiag<real_t>;

} // namespace gismo