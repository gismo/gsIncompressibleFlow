#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSTerms.h>
#include <gsIncompressibleFlow/src/gsINSTerms.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSTerm_PvalUdiv<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_UdivPval<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_UsolGradVal<real_t>;

    CLASS_TEMPLATE_INST gsINSTerm_CoeffUvalUval_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsUVal_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_CoeffUvalUvalPenalty_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsUValPenalty_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_CoeffUvalUdiv_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsUdiv_WeakDirichlet<real_t>;
    //CLASS_TEMPLATE_INST gsINSTerm_UvalPval_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsPvalU_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_PvalUval_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsUvalP_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_PvalPval_WeakDirichlet<real_t>;
    CLASS_TEMPLATE_INST gsINSTerm_RhsPvalP_WeakDirichlet<real_t>;
    
} // namespace gismo