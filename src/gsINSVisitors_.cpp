#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSVisitorUU<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlin<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlin<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorUUtimeDiscr<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPU<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPU_withUPrhs<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorUP<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPP<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlin<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPPnonlin<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPPmass<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlaplace<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorPPconvection<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsU<real_t>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsP<real_t>;

} // namespace gismo