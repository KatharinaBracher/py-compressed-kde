#pragma once

#include "grid_base.hpp"
#include "grid_vector.hpp"
#include "grid_array.hpp"
#include "grid_multi.hpp"

// yaml
std::unique_ptr<Grid> grid_from_yaml( const YAML::Node & node );
std::unique_ptr<Grid> load_grid_from_yaml( std::string path );

// hdf5
std::unique_ptr<Grid> grid_from_hdf5( const HighFive::Group & group );




