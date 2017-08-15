
#include "pybind.hpp"

#include "space.hpp"
#include "grid.hpp"

void pybind_space(py::module &m) {
    
    // SPACE CLASS
    py::class_<Space>(m, "Space",
    R"pbdoc(Base class for space definitions.)pbdoc")
    .def_property_readonly("ndim", &Space::ndim, 
    R"pbdoc(Number of dimensions.)pbdoc")
    
    .def_property_readonly("nbw", &Space::nbw, 
    R"pbdoc(Number of bandwidth values.)pbdoc")
    
    .def("__eq__", [](const Space& a, const Space& b) { return a==b; },
    R"pbdoc(
        Test if spaces are equal.
        
        Returns
        -------
        bool
        
    )pbdoc" )
    
    .def( "issubspace", &Space::issubspace, py::arg("space"),
    R"pbdoc(
        Test if space is subspace.
        
        Parameters
        ----------
        space : Space
        
        Returns
        -------
        bool
    
    )pbdoc" )
    
    .def( "selection", &Space::selection, py::arg("space"),
    R"pbdoc(
        Selection of dimensions that make up subspace.
        
        Parameters
        ----------
        space : Space
            Proper subspace of this Space
        
        Returns
        -------
        [bools]
            For each dimension if it is part of the subspace
        
    )pbdoc" )
    
    .def_property_readonly( "default_kernel", &Space::default_kernel,
    R"pbdoc(Default kernel for this space.)pbdoc")
    
    .def("to_yaml", [](Space &k)->std::string {
        YAML::Emitter out;
        YAML::Node node = k.toYAML();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;},
    R"pbdoc(
        Represent space definition as YAML.
        
        Returns
        -------
        string
        
    )pbdoc" )
    
    .def("save", [](const Space& obj, std::string path) { obj.save( path ); },
    py::arg("path"),
    R"pbdoc(
        Save space definition to YAML file.
        
        Parameters
        ----------
        path : string
            path tho YAML file
        
    )pbdoc" )
    
    .def_static("load", [](std::string path) { return std::unique_ptr<Space>( load_space(path) ); }, py::arg("path"),
    R"pbdoc(
        Load space definition from file.
        
        Parameters
        ----------
        path : string
            path to YAML file
        
        Returns
        -------
        Space
        
    )pbdoc" )
    
    .def_static("from_yaml", [](std::string s) {
        YAML::Node node = YAML::Load( s );
        return std::unique_ptr<Space>( space_from_YAML( node ) );
    }, py::arg("string"),
    R"pbdoc(
        Construct space definition from YAML
        
        Parameters
        ----------
        string : string
            YAML string space representation
        
        Returns
        -------
        Space
        
    )pbdoc");
    
    
    // EUCLIDEAN SPACE CLASS
    py::class_<EuclideanSpace, Space>(m, "EuclideanSpace",
    R"pbdoc(
        N-dimensional euclidean space definition.
        
        Parameters
        ----------
        labels for dimensions
        kernel : Kernel
            Gaussian, Epanechnikov or Box kernel
        bandwidth : 1d array
            bandwidths for default kernel
            
    )pbdoc")

    .def( "__init__", [](EuclideanSpace & obj, std::vector<std::string> labels, const Kernel & kernel, py::array_t<value, py::array::c_style | py::array::forcecast> bandwidth){
        auto bw = numpy_array_to_vector( bandwidth );
        new (&obj) EuclideanSpace( labels, kernel, bw );
    },
    py::arg("labels"), py::arg("kernel") = GaussianKernel(), py::arg("bandwidth") = std::vector<value>(0), 
    R"pbdoc(
        Constructs Euclidean Space.
        
        Parameters
        ----------
        labels : [strings]
            labels for dimensions
        kernel : Kernel
            Gaussian, Epanechnikov or Box kernel
        bandwidth : 1d array
            bandwidths for default kernel
    )pbdoc")
    
    .def( "grid", [](const EuclideanSpace& obj, std::vector<py::array_t<value, py::array::c_style | py::array::forcecast>> & v, py::array_t<bool, py::array::c_style | py::array::forcecast> & valid, py::array_t<bool, py::array::c_style | py::array::forcecast> & selection ) { 
        std::vector<std::vector<value>> vec( v.size() );
        for (unsigned int k=0; k<v.size(); ++k) {
            vec[k] = numpy_array_to_vector( v[k] );
        }
        
        auto vec_valid = numpy_array_to_vector( valid );
        auto vec_select = numpy_array_to_vector( selection );
        
        // construct grid
        return std::unique_ptr<Grid>( obj.grid(vec, vec_valid, vec_select) );
        
    }, py::arg("vectors"), py::arg("valid") = std::vector<bool>(0), py::arg("selection") = std::vector<bool>(0),
    R"pbdoc(
        Constructs grid from vectors for select dimensions.
        
        Parameters
        ----------
        vectors: [ 1d arrays ]
            A list with for each dimension a list of grid points
        valid : 1d boolean array
            For each grid point in the n-dimensional grid if it is a valid point or not
        selection : (ndim,) boolean array
            For each dimension if a grid vector is specified
        
        Returns
        -------
        Grid
        
    )pbdoc" );
    
    
    // CATEGORICAL SPACE CLASS
    py::class_<CategoricalSpace, Space>(m, "CategoricalSpace",
    R"pbdoc(
        1-dimensional categorical space definition.
        
        Parameters
        ----------
        label : string
            Label for categorical dimension
        categories : [ strings ]
            Category names
        index : int
            Default category index
            
    )pbdoc")
    .def( py::init<std::string, std::vector<std::string>, unsigned int>(), py::arg("label"), py::arg("categories"), py::arg("index")=0,
    R"pbdoc(
        Constructs Categorical Space.
        
        Parameters
        ----------
        label : string
            Label for categorical dimension
        categories : [ strings ]
            Category names
        index : int
            Default category index
        
    )pbdoc" )
    
    .def( "grid", [](const CategoricalSpace& obj) { return std::unique_ptr<Grid>( obj.grid() ); },
    R"pbdoc(
        Constructs grid.
        
        Returns
        -------
        Grid
        
    )pbdoc" );
    
    
    // MULTIPLICATIVE SPACE CLASS
    py::class_<MultiSpace, Space>(m, "MultiSpace",
    R"pbdoc(
        Multiplicative space definition.
        
        Parameters
        ----------
        spaces : [ Spaces ]
        
    )pbdoc")
    .def( py::init<std::vector<const Space*>>(), py::arg("spaces"),
    R"pbdoc(
        Constructs Multiplicative Space.
        
        Parameters
        ----------
        spaces : [ Spaces ]
        
    )pbdoc"  )
    
    .def( "grid", [](const MultiSpace& obj, std::vector<Grid*> & g) { return std::unique_ptr<Grid>( obj.grid( g ) ); }, py::arg("grids") ,
    R"pbdoc(
        Construct grid.
        
        Parameters
        ----------
        grids : [ Grids ]
            Grids for subspaces
        
        Returns
        -------
        Grid
        
    )pbdoc" );
    
    
    // CIRCULAR SPACE CLASS
    py::class_<CircularSpace, Space>(m, "CircularSpace",
    R"pbdoc(
        1-dimensional circular space definition.
        
        Parameters
        ----------
        label : string
            label for circular dimension
        kappa : scalar
            Kappa for default Von Mises kernel
        mu : scalar
            Mu for default Von Mises kernel
        
    )pbdoc")
    .def( py::init<std::string, value, value>(), py::arg("label"), py::arg("kappa")=DEFAULT_KAPPA, py::arg("mu")=DEFAULT_MU,
    R"pbdoc(
        Constructs Circular Space.
        
        Parameters
        ----------
        label : string
            label for circular dimension
        kappa : scalar
            Kappa for default Von Mises kernel
        mu : scalar
            Mu for default Von Mises kernel
        
    )pbdoc"  )
    
    .def( "grid", [](const CircularSpace& obj, unsigned int n) { return std::unique_ptr<Grid>( obj.grid(n) ); }, py::arg("n")=DEFAULT_CIRCULAR_GRID_SIZE,
    R"pbdoc(
        Construct grid.
        
        Parameters
        ----------
        n : int
            number of points in circular grid
        
        Returns
        -------
        Grid
        
    )pbdoc" );
    
    
    // ENCODED SPACE CLASS
    py::class_<EncodedSpace, Space>(m, "EncodedSpace",
    R"pbdoc(
        1-dimensional encoded space definition.
        
        Parameters
        ----------
        label : string
            label for encoded dimension
        distances : (n,n) array
            matrix of squared distances
        bandwidth : scalar
            Bandwidth for default gaussian kernel
        
    )pbdoc")
    .def( "__init__", [](EncodedSpace & obj, std::string label, py::array_t<value, py::array::c_style | py::array::forcecast> distances, value bandwidth) {
        auto vec_distances = numpy_array_to_vector( distances );
        new (&obj) EncodedSpace( label, vec_distances, bandwidth );
    },
    py::arg("label"), py::arg("distances"), py::arg("bandwidth")=DEFAULT_ENCODED_BANDWIDTH,
    R"pbdoc(
        Constructs Encoded Space.
        
        Parameters
        ----------
        label : string
            label for encoded dimension
        distances : (n,n) array
            matrix of squared distances
        bandwidth : scalar
            Bandwidth for default gaussian kernel
        
    )pbdoc"  )
    
    .def( "grid", [](const EncodedSpace& obj, unsigned int delta) { return std::unique_ptr<Grid>( obj.grid(delta) ); },
    py::arg("delta")=DEFAULT_ENCODED_GRID_DELTA,
    R"pbdoc(
        Constructs grid.
        
        Parameters
        ----------
        delta : int
            Sampling interval for grid.
        
        Returns
        -------
        Grid
        
    )pbdoc"  );
    
}
