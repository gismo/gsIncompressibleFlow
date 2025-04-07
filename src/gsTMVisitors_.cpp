#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsTMVisitors.h>
#include <gsIncompressibleFlow/src/gsTMVisitors.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsTMVisitorLinearSST<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsTMVisitorTimeIterationSST<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsTMVisitorNonlinearSST<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsTMVisitorLinearSST<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsTMVisitorTimeIterationSST<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsTMVisitorNonlinearSST<real_t, ColMajor>;
    
} // namespace gismo