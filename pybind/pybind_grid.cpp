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
        
    )pbdoc" )
    
    .def("save_to_yaml", [](const Grid& obj, std::string path) { obj.save_to_yaml( path ); },
    py::arg("path"),
    R"pbdoc(
        Save grid definition to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc" )
    
    .def_static("load_from_yaml", [](std::string path) { return std::unique_ptr<Grid>( load_grid_from_yaml(path) ); }, py::arg("path"),
    R"pbdoc(
        Load grid definition from file.
        
        Parameters
        ----------
        path : string
            path to YAML file
        
        Returns
        -------
        Grid
        
    )pbdoc" )
    
    .def("at_index", [](Grid & obj, py::array_t<unsigned int, py::array::c_style | py::array::forcecast> index) {
        unsigned int ndim = obj.ndim();
        unsigned int npoints;
        
        std::vector<long unsigned int> strides;
        
        // check if index is 2d array of size (n,ndim)
        // or 1d array of size (ndim,)
        auto buf = index.request();
        
        if (buf.ndim==1 && buf.shape[0]==ndim) {
            npoints = 1;
            strides.push_back(sizeof(value));
        } else if (buf.ndim==2 && buf.shape[1]==ndim) {
            npoints = buf.shape[0];
            strides.push_back(sizeof(value)*ndim);
            strides.push_back(sizeof(value));
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        // create output array<value> of same shape
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            1,
            buf.shape,
            strides
        ));
        
        auto result_buf = result.request();
        
        // repeatedly call obj.at_index( in, out )
        unsigned int *p_in = (unsigned int *) buf.ptr;
        value *p_out = (value *) result_buf.ptr;
        
        for (unsigned int k=0; k<npoints; ++k) {
            obj.at_index( p_in, p_out );
            p_in += ndim;
            p_out += ndim;
        }
        
        // return output array
        return result;
        
        },
    py::arg("index"),
    R"pbdoc(
        Retrieve grid values at index
        
        Parameters
        ----------
        index : (ndim,) or (n,ndim) array
            Array of indices
        
        Returns
        -------
        array
            grid values at index
        
    )pbdoc" );
        
}
