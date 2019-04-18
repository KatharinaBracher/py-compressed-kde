
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
        YAML::Node node = k.to_yaml();
        out << YAML::Flow;
        out << node;
        std::string s = out.c_str();
        return s;},
    R"pbdoc(
        to_yaml() -> str

        Represent kernel definition as YAML.
        
        Returns
        -------
        string
        
    )pbdoc" )
        
    .def_static("from_yaml", [](std::string s) {
        YAML::Node node = YAML::Load( s );
        return std::unique_ptr<Kernel>( kernel_from_yaml( node ) ); },
    py::arg("string"),
    R"pbdoc(
        from_yaml(str) -> Kernel

        Construct kernel definition from YAML
        
        Parameters
        ----------
        string : string
            YAML string kernel representation
        
        Returns
        -------
        Kernel
        
    )pbdoc" );
    
    py::class_<GaussianKernel, Kernel>(m,"GaussianKernel",
    R"pbdoc(
        Gaussian kernel function.
        
        Parameters
        ----------
        cutoff : scalar
            The standard deviation of the gaussian kernel beyond which
            the probability is set to zero.
            
    )pbdoc")
    .def(py::init<value>(), py::arg("cutoff")=DEFAULT_GAUSSIAN_CUTOFF);
    
    py::class_<EpanechnikovKernel, Kernel>(m,"EpanechnikovKernel",
    R"pbdoc(
        Epanechnikov kernel function.
    )pbdoc")
    .def(py::init<>());
    
    py::class_<BoxKernel, Kernel>(m,"BoxKernel",
    R"pbdoc(
        Box kernel function.
    )pbdoc")
    .def(py::init<>());
    
}
