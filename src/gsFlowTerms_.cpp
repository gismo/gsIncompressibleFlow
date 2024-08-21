#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsFlowTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowTerm<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermNonlin<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermValVal<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermTimeDiscr<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermGradGrad<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermDiffusion<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermRhs<real_t>;

} // namespace gismo