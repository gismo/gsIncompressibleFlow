#include <gsCore/gsTemplateTools.h>

#include <gsIncompressibleFlow/src/gsNavStokesPde.h>

namespace gismo
{
    CLASS_TEMPLATE_INST gsNavStokesPde<real_t>;

#ifdef GISMO_WITH_PYBIND11

  namespace py = pybind11;

  void pybind11_init_gsNavStokesPde(py::module &m)
  {
    using Class = gsNavStokesPde<real_t>;
    py::class_<Class>(m, "gsNavStokesPde")

    // Constructors
    .def(py::init<const gsMultiPatch<real_t> &, 
                    const gsBoundaryConditions<real_t> &,
                    const gsFunction<real_t> *,
                    const real_t
                    >())
    .def(py::init<const gsMultiPatch<real_t> &, 
                    const gsBoundaryConditions<real_t> &,
                    const gsFunction<real_t> *,
                    const gsFunction<real_t> *,
                    const real_t
                    >())

    .def("__str__", [] (Class & self)
     {
         std::ostringstream os;
         self.print(os);
         return os.str();
     },
     "Returns a string with information about the object.")
    ;     
  }

#endif

} // namespace gismo
