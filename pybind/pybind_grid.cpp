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
    .def_property("valid",
        [](Grid &grid)-> py::array_t<bool> {
            
            // short cut if empty array
            if (grid.valid().size()==0) {
                // create output array<value> of same shape
                auto result = py::array( py::buffer_info(
                    nullptr,
                    sizeof(bool),
                    py::format_descriptor<bool>::value,
                    1,
                    {0},
                    {sizeof(bool)})
                );

                return result;
            }

            std::vector<long unsigned int> strides(grid.ndim(), sizeof(bool));
            for (int k=grid.ndim()-2; k>=0; --k) {
                strides[k] = strides[k+1] * grid.shape()[k+1];
            }

            // create output array<value> of same shape
            auto result = py::array( py::buffer_info(
                nullptr,
                sizeof(bool),
                py::format_descriptor<bool>::value,
                grid.ndim(),
                grid.shape(),
                strides
            ));

           auto result_ptr = (bool*) result.request().ptr;
            
            for (auto v : grid.valid()) {
                *result_ptr++ = v;
            }
            
            return result;

        },
        [](Grid &grid, py::array_t<bool, py::array::c_style | py::array::forcecast> valid) {
            // check that shape of input equals the grid shape
            auto buf = valid.request();
            auto gridshape = grid.shape();
            auto bufshape = buf.shape;
            if (buf.size>0 &&
                    (buf.ndim!=grid.ndim() ||
                    !std::equal(gridshape.begin(), gridshape.end(), bufshape.begin())
                    )
                ) {
                throw std::runtime_error("Shape of input array does not match grid shape.");
            }

            // flatten input to vector
            auto v = numpy_array_to_vector(valid);

            grid.set_valid(v);
        },
    R"pbdoc(
        Validity of grid points.
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
        to_yaml() -> str

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
        from_yaml(str) -> Grid

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
        save_to_yaml(path) -> None

        Save grid definition to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc" )
    
    .def_static("load_from_yaml", [](std::string path) { return std::unique_ptr<Grid>( load_grid_from_yaml(path) ); }, py::arg("path"),
    R"pbdoc(
        load_from_yaml(path) -> Grid

        Load grid definition from file.
        
        Parameters
        ----------
        path : string
            path to YAML file
        
        Returns
        -------
        Grid
        
    )pbdoc" )
    
    .def(py::pickle(
        &(pickle_get_state<Grid>),
        [](py::tuple t) {
            if (t.size() != 1)
                throw std::runtime_error("Invalid state!");

            char *buffer_pointer;
            ssize_t length;
            PYBIND11_BYTES_AS_STRING_AND_SIZE(
                t[0].cast<py::bytes>().ptr(), &buffer_pointer, &length
            );

            // Get a pointer to the root object inside the buffer.
            auto ptr = flatbuffers::GetRoot<fb_serialize::Grid>((uint8_t*)buffer_pointer);

            return grid_from_flatbuffers(ptr);
        }
    ))

    .def("at_index", [](Grid & obj, py::array_t<unsigned int, py::array::c_style | py::array::forcecast> index) {
        unsigned int ndim = obj.ndim();
        unsigned int npoints;
        
        std::vector<long unsigned int> strides;
        
        // check if index is 2d array of size (n,ndim)
        // or 1d array of size (ndim,)
        // or if ndim=1, 1d array of size (n,)
        auto buf = index.request();
        
        if (ndim==1 && buf.ndim==1) {
            npoints = buf.shape[0];
        } else if (buf.ndim==1 && buf.shape[0]==ndim) {
            npoints = 1;
        } else if (buf.ndim==2 && buf.shape[1]==ndim) {
            npoints = buf.shape[0];
            strides.push_back(sizeof(value)*ndim);
        } else {
            throw std::runtime_error("Expected a (N," + std::to_string(ndim) + ") 2D array of samples.");
        }
        
        strides.push_back(sizeof(value));

        // create output array<value> of same shape
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            buf.ndim,
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
        at_index(index) -> array
        
        Retrieve grid values at index
        
        Parameters
        ----------
        index : (ndim,) or (n,ndim) array
            Array of indices
        
        Returns
        -------
        array
            grid values at index
        
    )pbdoc" )

    .def("points", [](Grid & obj) {
        
        std::vector<unsigned int> indices(obj.size()*obj.ndim(), 0);
        std::vector<long unsigned int> strides(obj.ndim(), 1);
        
        for (unsigned int d=obj.ndim()-1; d>0; --d) {
            strides[d-1] = strides[d]*obj.shape()[d];
        }

        auto cursor = indices.begin();

        for (long unsigned int k=0; k<obj.size(); ++k) {
            for (unsigned int d=0; d<obj.ndim(); ++d) {
                *cursor = (k / strides[d]) % obj.shape()[d];
                cursor++;
            }
        }
        
        // create output array<value>
        auto result = py::array( py::buffer_info(
            nullptr,
            sizeof(value),
            py::format_descriptor<value>::value,
            2,
            {obj.size(), obj.ndim()},
            {sizeof(value)*obj.ndim(), sizeof(value)}
        ));
        
        auto result_buf = result.request();
        
        // repeatedly call obj.at_index( in, out )
        auto p_in = &indices[0];
        value *p_out = (value *) result_buf.ptr;
        
        for (unsigned int k=0; k<obj.size(); ++k) {
            obj.at_index( p_in, p_out );
            p_in += obj.ndim();
            p_out += obj.ndim();
        }
        
        // return output array
        return result;
        
        },
    R"pbdoc(
        points -> array
        
        Retrieve all grid values
        
        Returns
        -------
        array : (n, ndim) array
            grid values
        
    )pbdoc" );
        
}
