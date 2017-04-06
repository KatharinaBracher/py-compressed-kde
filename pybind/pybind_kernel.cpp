
#include "pybind.hpp"

#include "kernel.hpp"

void pybind_kernel(py::module &m) {
    
    // KERNEL CLASS
    py::class_<Kernel>(m,"Kernel",
    R"pbdoc(
        Base class for kernel functions.
    )pbdoc")
    .def("to_yaml", [](Kernel &k)->std::string {
        YAML::Emitter out;
        YAML::Node node = k.toYAML();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;},
    R"pbdoc(
        Represent kernel definition as YAML.
        
        Returns
        -------
        string
        
    )pbdoc" )
        
    .def_static("from_yaml", [](std::string s) {
        YAML::Node node = YAML::Load( s );
        return std::unique_ptr<Kernel>( kernel_from_YAML( node ) ); },
    py::arg("string"),
    R"pbdoc(
        Construct kernel definition from YAML
        
        Parameters
        ----------
        string : string
            YAML string kernel representation
        
        Returns
        -------
        Kernel
        
    )pbdoc" );
    
    py::class_<GaussianKernel>(m,"GaussianKernel",py::base<Kernel>(),
    R"pbdoc(
        Gaussian kernel function.
        
        Parameters
        ----------
        cutoff : scalar
            The standard deviation of the gaussian kernel beyond which
            the probability is set to zero.
            
    )pbdoc")
    .def(py::init<value>(), py::arg("cutoff")=DEFAULT_GAUSSIAN_CUTOFF,
    R"pbdoc(
        Constructs Gaussian kernel.
        
        Parameters
        ----------
        cutoff : scalar
            The standard deviation of the gaussian kernel beyond which
            the probability is set to zero.
        
    )pbdoc"
    );
    
    py::class_<EpanechnikovKernel>(m,"EpanechnikovKernel",py::base<Kernel>(),
    R"pbdoc(
        Epanechnikov kernel function.
    )pbdoc")
    .def(py::init<>());
    
    py::class_<BoxKernel>(m,"BoxKernel",py::base<Kernel>(),
    R"pbdoc(
        Box kernel function.
    )pbdoc")
    .def(py::init<>());
    
}
