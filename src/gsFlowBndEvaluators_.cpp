#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.h>
#include <gsIncompressibleFlow/src/gsFlowBndEvaluators.hpp>

namespace gismo
{
    
    CLASS_TEMPLATE_INST gsFlowBndEvaluator<real_t>;
    CLASS_TEMPLATE_INST gsFlowBndEvaluator_flowRate<real_t>;

} // namespace gismo