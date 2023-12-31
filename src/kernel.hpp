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

#include "kernel_base.hpp"
#include "kernel_gaussian.hpp"
#include "kernel_epanechnikov.hpp"
#include "kernel_box.hpp"
#include "kernel_vonmises.hpp"

#include "yaml-cpp/yaml.h"

#include "schema_generated.h"

// yaml
std::unique_ptr<Kernel> kernel_from_yaml( const YAML::Node & node );

// flatbuffers
std::unique_ptr<Kernel> kernel_from_flatbuffers( const fb_serialize::Kernel * kernel);

// hdf5
std::unique_ptr<Kernel> kernel_from_hdf5( const HighFive::Group & group );












