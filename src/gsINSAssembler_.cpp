#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsINSAssembler.h>
#include <gsIncompressibleFlow/src/gsINSAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsINSAssembler<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerSteady<real_t, RowMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady<real_t, RowMajor>;
    
    CLASS_TEMPLATE_INST gsINSAssembler<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerSteady<real_t, ColMajor>;
    CLASS_TEMPLATE_INST gsINSAssemblerUnsteady<real_t, ColMajor>;
    
    
#ifdef GISMO_WITH_PYBIND11
    
    namespace py = pybind11;
    
    void pybind11_init_gsINSAssembler(py::module &m)
    {
        constexpr int Major = RowMajor;
        using Base = gsFlowAssemblerBase<real_t, Major>;
        using Class = gsINSAssembler<real_t, Major>;
        py::class_<Class,Base>(m, "gsINSAssembler")
        
        // Constructor with keep_alive to prevent Python object GC
        .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())

        // initialization / update
        .def("initialize", &Class::initialize, "Initialize the assembler")
        .def("update", &Class::update, "Update assembler with solution vector", py::arg("solVector"), py::arg("updateSol") = true)

        // stokes / solution helpers
        .def("fillStokesSystem", &Class::fillStokesSystem, "Fill the matrix and right-hand side for the Stokes problem", py::arg("stokesMat"), py::arg("stokesRhs"))
        .def("constructSolution", &Class::constructSolution, "Construct solution field from solution vector", py::arg("solVector"), py::arg("unk"), py::arg("customSwitch") = false)
        .def("computeFlowRate", &Class::computeFlowRate, "Compute flow rate through given patch side", py::arg("patch"), py::arg("side"), py::arg("solution"))

        // getters
        .def("getViscosity", &Class::getViscosity, "Return viscosity value")
        .def("matrix", &Class::matrix, "Returns the assembled matrix", py::return_value_policy::reference_internal)
        .def("rhs", &Class::rhs, "Returns the assembled right-hand side", py::return_value_policy::reference_internal)
        .def("getSolution", &Class::getSolution, "Returns the current computed solution", py::return_value_policy::reference_internal)
        .def("numDofsUnk", &Class::numDofsUnk, "Returns number of DOFs for the i-th unknown", py::arg("i"))
        .def("getUdofs", &Class::getUdofs, "Returns the number of velocity DOFs")
        .def("getPdofs", &Class::getPdofs, "Returns the number of pressure DOFs")
        .def("getPshift", &Class::getPshift, "Returns the DOF shift of pressure")

        .def("getBlockUU", &Class::getBlockUU, "Returns the velocity-velocity block of the linear system", py::arg("linPartOnly") = false)
        .def("getBlockUUcompDiag", &Class::getBlockUUcompDiag, "Returns the diagonal block of velocity-velocity block for i-th component", py::arg("i") = 0, py::arg("linPartOnly") = false)
        .def("getBlockUP", &Class::getBlockUP, "Returns the velocity-pressure block of the linear system", py::return_value_policy::reference_internal)
        .def("getBlockPU", &Class::getBlockPU, "Returns the pressure-velocity block of the linear system")
        .def("getBlockUPcomp", &Class::getBlockUPcomp, "Returns the part of velocity-pressure block for i-th velocity component", py::arg("i"))
        .def("getBlockPUcomp", &Class::getBlockPUcomp, "Returns part of pressure-velocity block for i-th velocity component", py::arg("i"))

        // mass matrix access (non-const via lambda to preserve reference semantics)
        .def("getMassMatrix", static_cast<gsSparseMatrix<double, Major>& (Class::*)(index_t)>(&Class::getMassMatrix), "Returns the mass matrix for unknown with index unk", py::arg("unk"), py::return_value_policy::reference_internal)
        .def("getMassMatrix", static_cast<const gsSparseMatrix<double, Major>& (Class::*)(index_t) const>(&Class::getMassMatrix), "Returns the mass matrix for unknown with index unk", py::arg("unk"), py::return_value_policy::reference_internal)


        .def("getRhsU", &Class::getRhsU, "Returns the velocity part of the right-hand side", py::arg("linPartOnly") = false)
        .def("getRhsUcomp", &Class::getRhsUcomp, "Returns part of the right-hand side for i-th velocity component", py::arg("i"), py::arg("linPartOnly") = false)
        .def("getRhsP", &Class::getRhsP, "Returns the pressure part of the right-hand side")
            ;
    }
    
    void pybind11_init_gsINSAssemblerSteady(py::module &m)
    {
        constexpr int Major = RowMajor;
        using BaseClass = gsINSAssembler<real_t, Major>;
        using Class = gsINSAssemblerSteady<real_t, Major>;
        py::class_<Class, BaseClass>(m, "gsINSAssemblerSteady")
            .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())
            ;
    }

    void pybind11_init_gsINSAssemblerUnsteady(py::module &m)
    {
        constexpr int Major = RowMajor;
        using BaseClass = gsINSAssembler<real_t, Major>;
        using Class = gsINSAssemblerUnsteady<real_t, Major>;
        py::class_<Class, BaseClass>(m, "gsINSAssemblerUnsteady")
            .def(py::init<gsFlowSolverParams<real_t>&>(), py::keep_alive<1, 2>())
            .def("update", &Class::update, "Update assembler with solution vector", py::arg("solVector"), py::arg("updateSol") = true)
            ;
    }
        
#endif
    
    
} // namespace gismo
