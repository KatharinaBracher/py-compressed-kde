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
#include "grid.hpp"
#include "spacespec.hpp"
#include <fstream>

// yaml
std::unique_ptr<Grid> grid_from_yaml( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["class"] || !node["space"]) {
        throw std::runtime_error("Not a valid YAML description of grid.");
    }
    
    SpaceSpecification space = SpaceSpecification::from_yaml( node["space"] );
    
    std::vector<bool> valid = node["valid"].as<std::vector<bool>>(
        std::vector<bool>({}));
        
    std::vector<long unsigned int> shape =
        node["shape"].as<std::vector<long unsigned int>>(
            std::vector<long unsigned int>({}));
    
    std::string klass = node["class"].as<std::string>( "unknown" );
    
    if (klass=="multi") {
        return MultiGrid::from_yaml( node["grid"], space, valid );
    } else if (klass=="vector") {
        return VectorGrid::from_yaml( node["grid"], space, valid );
    } else if (klass=="array") {
        return ArrayGrid::from_yaml( node["grid"], space, valid, shape );
    } else {
        throw std::runtime_error("Unknown grid.");
    }
    
}

std::unique_ptr<Grid> load_grid_from_yaml( std::string path ) {
        
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return grid_from_yaml( node );
    
}

// hdf5
std::unique_ptr<Grid> grid_from_hdf5( const HighFive::Group & group ) {
    
    std::string klass;
    
    HighFive::Attribute attr_klass = group.getAttribute("class");
    attr_klass.read(klass);
    
    auto space = SpaceSpecification::from_hdf5(group.getGroup("space"));
    
    std::vector<unsigned char> tmp;
    HighFive::DataSet ds_valid = group.getDataSet("valid");
    ds_valid.read(tmp);
    
    std::vector<bool> valid(tmp.begin(), tmp.end());
    
    if (klass=="multi") {
        return MultiGrid::from_hdf5( group.getGroup("grid"), space, valid );
    } else if (klass=="vector") {
        return VectorGrid::from_hdf5( group.getGroup("grid"), space, valid );
    } else if (klass=="array") {
        
        std::vector<long unsigned int> shape;
        HighFive::DataSet ds_shape = group.getDataSet("shape");
        ds_shape.read(shape);
        
        return ArrayGrid::from_hdf5( group.getGroup("grid"), space, valid, shape );
    } else {
        throw std::runtime_error("Unknown grid.");
    }
    
}
