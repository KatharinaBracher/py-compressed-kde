// ---------------------------------------------------------------------
// This file is part of the compressed decoder library.
//
// Copyright (C) 2020 - now Neuro-Electronics Research Flanders
//
// The compressed decoder library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The compressed decoder library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with falcon-core. If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------
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
