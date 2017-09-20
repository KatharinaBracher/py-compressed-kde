#include "space.hpp"
#include <fstream>


// yaml
std::unique_ptr<Space> space_from_yaml( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["class"] ) {
        throw std::runtime_error("Not a valid YAML description of space.");
    }
    
    std::string klass = node["class"].as<std::string>( "unknown" );
    
    if (klass=="multi") {
        return MultiSpace::from_yaml( node["space"] );
    } else if (klass=="euclidean") {
        return EuclideanSpace::from_yaml( node["space"] );
    } else if (klass=="categorical") {
        return CategoricalSpace::from_yaml( node["space"] );
    } else if (klass=="circular") {
        return CircularSpace::from_yaml( node["space"] );
    } else if (klass=="encoded") {
        return EncodedSpace::from_yaml( node["space"] );
    } else {
        throw std::runtime_error("Unknown space.");
    }
    
}


std::unique_ptr<Space> load_space_from_yaml( std::string path ) {
    
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return space_from_yaml( node );
    
}

// hdf5
std::unique_ptr<Space> space_from_hdf5(const HighFive::Group & group) {
    
    std::string klass;
    HighFive::Attribute attr_klass = group.getAttribute("class");
    attr_klass.read(klass);
    
    if (klass=="multi") {
        return MultiSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="euclidean") {
        return EuclideanSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="categorical") {
        return CategoricalSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="circular") {
        return CircularSpace::from_hdf5( group.getGroup("space") );
    } else if (klass=="encoded") {
        return EncodedSpace::from_hdf5( group.getGroup("space") );
    } else {
        throw std::runtime_error("Unknown space.");
    }
}
