#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMTerms.h>
#include <gsIncompressibleFlow/src/gsTMTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMTerm_CoeffGradGrad<real_t>;
    CLASS_TEMPLATE_INST gsTMTerm_VecCoeffGradVal<real_t>;
    CLASS_TEMPLATE_INST gsTMTerm_CoeffValVal<real_t>;

} // namespace gismo