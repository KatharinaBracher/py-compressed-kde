#include "pybind.hpp"

#include "component.hpp"

void pybind_component(py::module &m) {
    
    py::class_<Component>(m, "Component",
    R"pbdoc(
        Data container for single mixture component.
    )pbdoc")
    .def_property_readonly("scale_factor", [](Component & obj) { return obj.scale_factor; },
     R"pbdoc(
        Component scale factor.
    )pbdoc")
    .def_property_readonly("scale_factor_log", [](Component & obj) { return obj.scale_factor_log; },
    R"pbdoc(
        Component log scale factor.
    )pbdoc")
    .def_property_readonly("location", [](Component & obj) { return obj.location; },
    R"pbdoc(
        Component location.
    )pbdoc")
    .def_property_readonly("bandwidth", [](Component & obj) { return obj.bandwidth; },
    R"pbdoc(
        Component bandwidth.
    )pbdoc");
    
}
