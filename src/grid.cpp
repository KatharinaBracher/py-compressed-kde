#include "grid.hpp"
#include "spacespec.hpp"


Grid * grid_from_YAML( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["class"] || !node["space"]) {
        throw std::runtime_error("Not a valid YAML description of grid.");
    }
    
    SpaceSpecification space = SpaceSpecification::fromYAML( node["space"] );
    std::vector<bool> valid = node["valid"].as<std::vector<bool>>( std::vector<bool>({}) );
    std::vector<long unsigned int> shape = node["shape"].as<std::vector<long unsigned int>>( std::vector<long unsigned int>({}) );
    
    Grid* k = nullptr;
    
    std::string klass = node["class"].as<std::string>( "unknown" );
    
    if (klass=="multi") {
        k = MultiGrid::fromYAML( node["grid"], space, valid );
    } else if (klass=="vector") {
        k = VectorGrid::fromYAML( node["grid"], space, valid );
    } else if (klass=="array") {
        k = ArrayGrid::fromYAML( node["grid"], space, valid, shape );
    } else {
        throw std::runtime_error("Unknown grid.");
    }
    
    return k;
    
}
