#include "component.hpp"

YAML::Node Component::toYAML() const {
    YAML::Node node;
    node["loc"] = location;
    node["bw"] = bandwidth;
    return node;
}

Component* Component::fromYAML( const YAML::Node & node ) {
    
    Component * k = new Component();
    
    k->location = node["loc"].as<std::vector<value>>();
    k->bandwidth = node["bw"].as<std::vector<value>>();
    
    return k;
    
}
