#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSPrecondBlocks.h>
#include <gsIncompressibleFlow/src/gsINSPrecondBlocks.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSPrecondBlock<real_t>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockMod<real_t>;

    CLASS_TEMPLATE_INST gsINSPrecondBlockF<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFwhole<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFdiag<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFmod<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockBt<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurLSC<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurPCD<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurPCDmod<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurAL<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurSIMPLE<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurMSIMPLER<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurStokes<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSPrecondBlockF<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFwhole<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFdiag<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockFmod<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondBlockBt<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurLSC<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurPCD<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurPCDmod<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurAL<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurSIMPLE<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurMSIMPLER<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSPrecondSchurStokes<real_t, ColMajor>;
} 
