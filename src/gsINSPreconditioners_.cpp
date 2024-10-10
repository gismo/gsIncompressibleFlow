#include <gsCore/gsTemplateTools.h>

#include "gsINSPreconditioners.h"
#include "gsINSPreconditioners.hpp"

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSPreconditioner<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondBase<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondLSC<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCD<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCDmod<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondAL<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLE<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLER<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondMSIMPLER<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokes<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokesTriang<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSPreconditioner<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondBase<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondLSC<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCD<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCDmod<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondAL<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLE<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLER<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondMSIMPLER<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokes<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokesTriang<real_t, ColMajor>;
} 
