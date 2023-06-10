#include <gsCore/gsTemplateTools.h>

#include "gsINSPreconditioners.h"
#include "gsINSPreconditioners.hpp"

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSPreconditioner<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondBase<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondLSC<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCD<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondPCDmod<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondAL<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLE<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondSIMPLER<real_t>;
    CLASS_TEMPLATE_INST gsINSBlockPrecondMSIMPLER<real_t>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokes<real_t>;
    CLASS_TEMPLATE_INST gsBlockPrecondStokesTriang<real_t>;
} 
