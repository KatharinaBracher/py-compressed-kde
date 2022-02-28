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
#pragma once

#include "grid_base.hpp"
#include "grid_vector.hpp"
#include "grid_array.hpp"
#include "grid_multi.hpp"

// yaml
std::unique_ptr<Grid> grid_from_yaml( const YAML::Node & node );
std::unique_ptr<Grid> load_grid_from_yaml( std::string path );

// flatbuffers
std::unique_ptr<Grid> grid_from_flatbuffers(const fb_serialize::Grid * grid );

// hdf5
std::unique_ptr<Grid> grid_from_hdf5( const HighFive::Group & group );




