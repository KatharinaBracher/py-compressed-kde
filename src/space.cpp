#include "space.hpp"
#include <fstream>


// yaml
std::unique_ptr<Space> space_from_yaml( const YAML::Node & node ) {
    
    std::unique_ptr<Space> ptr;
    
    if (!node.IsMap() || !node["class"] ) {
        throw std::runtime_error("Not a valid YAML description of space.");
    }
    
    std::string klass = node["class"].as<std::string>( "unknown" );
    
    if (klass=="multi") {
        ptr = MultiSpace::from_yaml( node["space"] );
    } else if (klass=="euclidean") {
        ptr = EuclideanSpace::from_yaml( node["space"] );
    } else if (klass=="categorical") {
        ptr = CategoricalSpace::from_yaml( node["space"] );
    } else if (klass=="circular") {
        ptr = CircularSpace::from_yaml( node["space"] );
    } else if (klass=="encoded") {
        ptr = EncodedSpace::from_yaml( node["space"] );
    } else {
        throw std::runtime_error("Unknown space.");
    }
    
    auto default_kernel = Component::from_yaml(node["kernel"]);
    
    ptr->set_default_kernel(*default_kernel);
    
    return ptr;
    
}


std::unique_ptr<Space> load_space_from_yaml( std::string path ) {
    
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return space_from_yaml( node );
    
}

// hdf5
std::unique_ptr<Space> space_from_hdf5(const HighFive::Group & group) {
    
    std::unique_ptr<Space> ptr;
    
    std::string klass;
    HighFive::Attribute attr_klass = group.getAttribute("class");
    attr_klass.read(klass);
    
    if (klass=="multi") {
        ptr = MultiSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="euclidean") {
        ptr = EuclideanSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="categorical") {
        ptr = CategoricalSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="circular") {
        ptr = CircularSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="encoded") {
        ptr = EncodedSpace::from_hdf5( group.getGroup("space") );
    } else {
        throw std::runtime_error("Unknown space.");
    }
    
    auto default_kernel = Component::from_hdf5(group.getGroup("kernel"));
    
    ptr->set_default_kernel(*default_kernel);
    
    return ptr;
}
