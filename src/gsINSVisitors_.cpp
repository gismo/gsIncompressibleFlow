#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSVisitorUU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUtimeDiscr<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU_withUPrhs<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUP<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPP<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPnonlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPmass<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlaplace<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPconvection<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsP<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSVisitorUU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUtimeDiscr<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU_withUPrhs<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUP<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPP<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPnonlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPmass<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlaplace<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPconvection<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsP<real_t, ColMajor>;


} // namespace gismo