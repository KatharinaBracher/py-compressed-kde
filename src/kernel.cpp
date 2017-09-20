#include "kernel.hpp"
#include <stdexcept>

// yaml
std::unique_ptr<Kernel> kernel_from_YAML( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["type"] ) {
        throw std::runtime_error("Not a valid YAML description of kernel.");
    }
    
    KernelType ktype = kerneltype_fromstring( node["type"].as<std::string>( "unknown" ) );
    
    if (ktype == KernelType::Gaussian) {
        return GaussianKernel::from_yaml( node["info"] );
    } else if (ktype == KernelType::Epanechnikov) {
        return EpanechnikovKernel::from_yaml( node["info"] );
    } else if (ktype == KernelType::Box) {
        return BoxKernel::from_yaml( node["info"] );
    } else {
        throw std::runtime_error("Unknown kernel type.");
    }

}


// hdf5
std::unique_ptr<Kernel> kernel_from_hdf5( const HighFive::Group & group ) {
    
    std::string type_str;
    HighFive::DataSet ds_type = group.getDataSet("type");
    ds_type.read(type_str);
    
    KernelType ktype = kerneltype_fromstring( type_str );
    
    if (ktype == KernelType::Gaussian) {
        return GaussianKernel::from_hdf5( group.getGroup("info") );
    } else if (ktype == KernelType::Epanechnikov) {
        return EpanechnikovKernel::from_hdf5( group.getGroup("info") );
    } else if (ktype == KernelType::Box) {
        return BoxKernel::from_hdf5( group.getGroup("info") );
    } else {
        throw std::runtime_error("Unknown kernel type.");
    }
    
}




