#include "kernel.hpp"
#include <stdexcept>


Kernel * kernel_from_YAML( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["type"] ) {
        throw std::runtime_error("Not a valid YAML description of kernel.");
    }
    
    Kernel * k = nullptr;
    
    KernelType ktype = kerneltype_fromstring( node["type"].as<std::string>( "unknown" ) );
    
    if (ktype == KernelType::Gaussian) {
        k = GaussianKernel::fromYAML( node["info"] );
    } else if (ktype == KernelType::Epanechnikov) {
        k = EpanechnikovKernel::fromYAML( node["info"] );
    } else {
        k = BoxKernel::fromYAML( node["info"] );
    }
    
    return k;
}






