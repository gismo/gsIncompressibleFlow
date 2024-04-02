#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSTerm<real_t>;
    CLASS_TEMPLATE_INST gsINSTermNonlin<real_t>;
    CLASS_TEMPLATE_INST gsINSTermValVal<real_t>;
    CLASS_TEMPLATE_INST gsINSTermTimeDiscr<real_t>;
    CLASS_TEMPLATE_INST gsINSTermGradGrad<real_t>;
    CLASS_TEMPLATE_INST gsINSTermDiffusion<real_t>;
    CLASS_TEMPLATE_INST gsINSTermPvalUdiv<real_t>;
    CLASS_TEMPLATE_INST gsINSTermUdivPval<real_t>;
    CLASS_TEMPLATE_INST gsINSTermUsolGradVal<real_t>;
    CLASS_TEMPLATE_INST gsINSTermRhs<real_t>;

} // namespace gismo