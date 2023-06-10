#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSBlockAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSBlockAssemblerBase.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSBlockAssemblerBase<real_t>;
    
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockAVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlocksBVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockNVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockMVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockMpVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockApVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleBlock<gsINSBlockNpVisitor<real_t> >(gsSparseMatrix<real_t> & matrixBlock, gsMatrix<real_t> & rhs, const int nonZerosPerCol);
    
    TEMPLATE_INST void gsINSBlockAssemblerBase<real_t>::assembleRhs<gsINSRhsVisitor<real_t> >(gsMatrix<real_t> & rhs);

} // namespace gismo
