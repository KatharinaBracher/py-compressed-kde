#include "pybind.hpp"

#include "grid.hpp"

void pybind_grid(py::module &m) {
    
    // GRID CLASS
    py::class_<Grid>(m,"Grid",
    R"pbdoc(
        Base class for evaluation grids.
    )pbdoc")
    .def_property_readonly("shape", &Grid::shape,
    R"pbdoc(
        Grid shape (number of grid points for each dimension.)
    )pbdoc")
    .def_property_readonly("size", &Grid::size,
    R"pbdoc(
        Grid size (total number of points in grid).
    )pbdoc")
    .def_property_readonly("ndim", &Grid::ndim,
    R"pbdoc(
        Dimensionality of grid.
    )pbdoc")
    .def_property_readonly("klass", &Grid::klass,
    R"pbdoc(
        Grid type.
    )pbdoc")
    
    .def("to_yaml", [](Grid &obj)->std::string {
        YAML::Emitter out;
        YAML::Node node = obj.to_yaml();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;},
    R"pbdoc(
        Represent grid definition as YAML.
        
        Returns
        -------
        string
        
    )pbdoc" )
    
    .def_static("from_yaml", [](std::string s) {
        YAML::Node node = YAML::Load(s);
        return std::unique_ptr<Grid>( grid_from_yaml( node ) ); },
    py::arg("string"),
    R"pbdoc(
        Construct grid definition from YAML
        
        Parameters
        ----------
        string : string
            YAML string grid representation
        
        Returns
        -------
        Grid
        
    )pbdoc" );
        
}
