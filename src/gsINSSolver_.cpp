#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolver.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSSolver<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSSolverSteady<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSSolverUnsteady<real_t, RowMajor>;

    CLASS_TEMPLATE_INST gsINSSolver<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSSolverSteady<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSSolverUnsteady<real_t, ColMajor>;
    
#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsINSSolver(py::module &m)
    {
        using Class = gsINSSolver<real_t, RowMajor>;
        py::class_<Class>(m, "gsINSSolver")
            .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())
            .def("getName", &Class::getName, "Returns the name of the class")
            .def("solveStokes", &Class::solveStokes, "Compute the Stokes problem and save solution")
            .def("getAssembler", &Class::getAssembler, "Returns a pointer to the assembler", py::return_value_policy::reference_internal)
            ;
    }

    void pybind11_init_gsINSSolverSteady(py::module &m)
    {
        using BaseClass = gsINSSolver<real_t, RowMajor>;
        using Class = gsINSSolverSteady<real_t, RowMajor>;
        py::class_<Class, BaseClass>(m, "gsINSSolverSteady")
            .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())
            .def("nextIteration", &Class::nextIteration, "Perform next iteration step")
            .def("getAssembler", &Class::getAssembler, "Returns a pointer to the assembler", py::return_value_policy::reference_internal)
            ;
    }

    void pybind11_init_gsINSSolverUnsteady(py::module &m)
    {
        using BaseClass = gsINSSolver<real_t, RowMajor>;
        using Class = gsINSSolverUnsteady<real_t, RowMajor>;
        py::class_<Class, BaseClass>(m, "gsINSSolverUnsteady")
            .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())
            .def(py::init<gsFlowSolverParams<real_t>&, bool>(), py::keep_alive<1, 2>())
            .def("nextIteration", &Class::nextIteration, "Perform next iteration step")
            .def("getAssembler", &Class::getAssembler, "Returns a pointer to the assembler", py::return_value_policy::reference_internal)
            .def("getSimulationTime", &Class::getSimulationTime, "Returns the elapsed simulation time")
            .def("getAvgPicardIterations", &Class::getAvgPicardIterations, "Returns the average number of Picard iterations per time step")
            ;
    }

#endif

} // namespace gismo
