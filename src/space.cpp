#include "space.hpp"
#include <fstream>

Space * space_from_YAML( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["class"] ) {
        throw std::runtime_error("Not a valid YAML description of space.");
    }
    
    Space * k = nullptr;
    
    std::string klass = node["class"].as<std::string>( "unknown" );
    
    if (klass=="multi") {
        k = MultiSpace::fromYAML( node["space"] );
    } else if (klass=="euclidean") {
        k = EuclideanSpace::fromYAML( node["space"] );
    } else if (klass=="categorical") {
        k = CategoricalSpace::fromYAML( node["space"] );
    } else if (klass=="circular") {
        k = CircularSpace::fromYAML( node["space"] );
    } else {
        throw std::runtime_error("Unknown space.");
    }
    
    return k;
    
}


Space * load_space( std::string path ) {
    
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return space_from_YAML( node );
    
}
