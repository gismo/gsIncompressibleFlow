#include <gsCore/gsTemplateTools.h>
// #include <gsIncompressibleFlow/gsIncompressibleFlow.h>
#include <gsIncompressibleFlow/src/gsNavStokesPde.h>
#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSSolver.h>

namespace gismo
{

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsIncompressibleFlow(py::module &m)
  {
    gismo::pybind11_init_gsNavStokesPde( m );
    gismo::pybind11_init_gsFlowSolverParams( m );
    gismo::pybind11_init_gsFlowAssemblerBase( m );
    gismo::pybind11_init_gsINSAssembler( m );
    gismo::pybind11_init_gsINSAssemblerSteady( m );
    gismo::pybind11_init_gsINSAssemblerUnsteady( m );
    gismo::pybind11_init_gsINSSolver( m );
    gismo::pybind11_init_gsINSSolverSteady( m );
    gismo::pybind11_init_gsINSSolverUnsteady( m );
    // gismo::pybind11_init_gsRANSSolverUnsteady( m );
  }

#endif

}
