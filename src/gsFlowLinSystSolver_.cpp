#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.h>
#include <gsIncompressibleFlow/src/gsFlowLinSystSolver.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowLinSystSolver<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_direct<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_iter<real_t, RowMajor, gsGMRes<real_t> >;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_iterSP<real_t, RowMajor, gsGMRes<real_t> >;

    CLASS_TEMPLATE_INST gsFlowLinSystSolver<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_direct<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_iter<real_t, ColMajor, gsGMRes<real_t> >;
    CLASS_TEMPLATE_INST gsFlowLinSystSolver_iterSP<real_t, ColMajor, gsGMRes<real_t> >;

} // namespace gismo