#include "component.hpp"

YAML::Node Component::to_yaml() const {
    YAML::Node node;
    node["loc"] = location;
    node["bw"] = bandwidth;
    return node;
}

std::unique_ptr<Component> Component::from_yaml( const YAML::Node & node ) {
    
    auto k = std::make_unique<Component>();
    
    k->location = node["loc"].as<std::vector<value>>();
    k->bandwidth = node["bw"].as<std::vector<value>>();
    
    return k;
    
}

void Component::to_hdf5(HighFive::Group & group) const {
    
    HighFive::DataSet loc = group.createDataSet<value>(
        "loc", HighFive::DataSpace::From(location));
    
    loc.write(location);
    
    HighFive::DataSet bw = group.createDataSet<value>(
        "bw", HighFive::DataSpace::From(bandwidth));
    
    bw.write(bandwidth);
}

std::unique_ptr<Component> Component::from_hdf5(const HighFive::Group & group) {
    
    auto k = std::make_unique<Component>();
    
    HighFive::DataSet loc = group.getDataSet("loc");
    loc.read(k->location);
    
    HighFive::DataSet bw = group.getDataSet("bw");
    bw.read(k->bandwidth);
    
    return k;
}
