#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.h>
#include <gsIncompressibleFlow/src/gsFlowAssemblerBase.hpp>

#ifdef GISMO_WITH_PYBIND11
#include <pybind11/pybind11.h>
#endif

namespace gismo
{
    
    CLASS_TEMPLATE_INST gsFlowAssemblerBase<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsFlowAssemblerBase<real_t, ColMajor>;

#ifdef GISMO_WITH_PYBIND11

    namespace py = pybind11;

    void pybind11_init_gsFlowAssemblerBase(py::module &m)
    {
        constexpr int Major = RowMajor;
        using Class = gsFlowAssemblerBase<real_t, Major>;
        py::class_<Class>(m, "gsFlowAssemblerBase")

        // Constructor with keep_alive to prevent Python object GC
        .def(py::init<gsFlowSolverParams<real_t>::Ptr>(), py::keep_alive<1, 2>())

        // Member functions
        .def("initialize", &Class::initialize, "Initialize the assembler")
        .def("update", &Class::update, "Update assembler with solution vector", py::arg("solVector"), py::arg("updateSol") = true)
        .def("constructSolution", &Class::constructSolution, "Construct solution field from solution vector", py::arg("solVector"), py::arg("unk"), py::arg("customSwitch") = false)
        .def("numDofs", &Class::numDofs, "Returns number of DOFs")
        .def("getTarDim", &Class::getTarDim, "Returns target dimension")
        .def("isInitialized", &Class::isInitialized, "Returns true if assembler is initialized")

        // Getters
        .def("getDirichletDofs", &Class::getDirichletDofs, "Returns coefficients at Dirichlet boundaries")
        .def("getSolution", &Class::getSolution, "Returns current computed solution")
        .def("getPatches", &Class::getPatches, "Returns multipatch domain")
        .def("getBasis", py::overload_cast<index_t>(&Class::getBasis), "Returns basis for unknown", py::arg("unk"))
        .def("getBCs", &Class::getBCs, "Returns boundary conditions")
        .def("getAssemblerOptions", &Class::getAssemblerOptions, "Returns assembler options")
        .def("options", &Class::options, "Returns flow solver option list")

        ;
    }

#endif

} // namespace gismo
