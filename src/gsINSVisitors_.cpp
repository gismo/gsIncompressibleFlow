#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSVisitors.h>
#include <gsIncompressibleFlow/src/gsINSVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSVisitorUU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUmass<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUrotation<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUtimeDiscr<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUU_TCSD_time<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU_withUPrhs<real_t, RowMajor>;
    //CLASS_TEMPLATE_INST gsINSVisitorPU_SUPG_presssure<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUP<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPP<real_t, RowMajor>;
    // CLASS_TEMPLATE_INST gsINSVisitorPPlin<real_t, RowMajor>;
    // CLASS_TEMPLATE_INST gsINSVisitorPPnonlin<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPmass<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlaplace<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPconvection<real_t, RowMajor>;
    //CLASS_TEMPLATE_INST gsINSVisitorPP_ResidualStabilization_continuity<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsU<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsP<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlinWeakDirichlet<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlinWeakDirichlet<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPUWeakDirichlet<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUPWeakDirichlet<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPWeakDirichlet<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSVisitorUU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUmass<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUrotation<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUtimeDiscr<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUU_TCSD_time<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPU_withUPrhs<real_t, ColMajor>;
    //CLASS_TEMPLATE_INST gsINSVisitorPU_SUPG_presssure<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUP<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPP<real_t, ColMajor>;
    // CLASS_TEMPLATE_INST gsINSVisitorPPlin<real_t, ColMajor>;
    // CLASS_TEMPLATE_INST gsINSVisitorPPnonlin<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPmass<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPlaplace<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPconvection<real_t, ColMajor>;
    //CLASS_TEMPLATE_INST gsINSVisitorPP_ResidualStabilization_continuity<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsU<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorRhsP<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUnonlinWeakDirichlet<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUUlinWeakDirichlet<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPUWeakDirichlet<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorUPWeakDirichlet<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSVisitorPPWeakDirichlet<real_t, ColMajor>;


} // namespace gismo