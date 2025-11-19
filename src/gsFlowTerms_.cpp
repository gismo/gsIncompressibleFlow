#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowTerms.h>
#include <gsIncompressibleFlow/src/gsFlowTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowTerm<real_t>;
    CLASS_TEMPLATE_INST gsFlowTermNonlin<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_ValVal<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_TimeDiscr<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_GradGrad<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_Diffusion<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_TCSDStabilization_time<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_TCSDStabilization_advection<real_t>;
    CLASS_TEMPLATE_INST gsFlowTerm_rhs<real_t>;

} // namespace gismo