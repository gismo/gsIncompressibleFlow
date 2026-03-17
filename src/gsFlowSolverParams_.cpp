#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowSolverParams.h>

namespace gismo
{
    CLASS_TEMPLATE_INST gsFlowSolverParams<real_t>;
#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsFlowSolverParams(py::module &m)
  {
    using Class = gsFlowSolverParams<real_t>;
    py::class_<Class>(m, "gsFlowSolverParams")

      // Constructors
      .def(py::init<const gsNavStokesPde<real_t> &,
                      const std::vector<gsMultiBasis<real_t> > &
                      >())
      .def(py::init<typename gsNavStokesPde<real_t>::Ptr,
                      const std::vector<gsMultiBasis<real_t> > &
                      >(), py::keep_alive<1, 2>())

      // Simple getters
      .def("getPde", &Class::getPde, "Returns a const reference to the PDE")
      .def("getBCs", &Class::getBCs, "Returns a const reference to the boundary conditions")
      .def("hasPeriodicBC", &Class::hasPeriodicBC, "Returns true if BCs contain periodic conditions")
      .def("isRotation", &Class::isRotation, "Returns true if angular velocity is non-zero")
      .def("updateDofMappers", &Class::updateDofMappers, "Update the stored DOF mappers", py::arg("finalize") = true)

      // Options
      .def("options", py::overload_cast<>(&Class::options, py::const_), "Returns a const reference to the options container", py::return_value_policy::reference_internal)
      .def("options", py::overload_cast<>(&Class::options), "Returns a reference to the options container")
      .def("precOptions", py::overload_cast<>(&Class::precOptions, py::const_), "Returns a const reference to the preconditioner options container", py::return_value_policy::reference_internal)
      .def("precOptions", py::overload_cast<>(&Class::precOptions), "Returns a reference to the preconditioner options container")
      .def("assemblerOptions", py::overload_cast<>(&Class::assemblerOptions, py::const_), "Returns a const reference to the assembly options container", py::return_value_policy::reference_internal)
      .def("assemblerOptions", py::overload_cast<>(&Class::assemblerOptions), "Returns a reference to the assembly options container")
      ;
  }

#endif

} // namespace gismo
