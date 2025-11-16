#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsRANSVisitors.h>
#include <gsIncompressibleFlow/src/gsRANSVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsRANSVisitorUU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorUU_full<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorTCSDStabilization_time<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorTCSDStabilization_advection<real_t, RowMajor>;
    //CLASS_TEMPLATE_INST gsRANSVisitorUUSymmetricGradientDiag<real_t, RowMajor>;
    //CLASS_TEMPLATE_INST gsRANSVisitorUUSymmetricGradientOffdiag<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsRANSVisitorUU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorUU_full<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorTCSDStabilization_time<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsRANSVisitorTCSDStabilization_advection<real_t, ColMajor>;
    //CLASS_TEMPLATE_INST gsRANSVisitorUUSymmetricGradientDiag<real_t, ColMajor>;
    //CLASS_TEMPLATE_INST gsRANSVisitorUUSymmetricGradientOffdiag<real_t, ColMajor>;
    
} // namespace gismo