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

#include "common.hpp"
#include "yaml-cpp/yaml.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
// without next include, there a compile error:
// invalid use of incomplete type ‘class HighFive::Group’
// possibly due to dependency on the order of inclusion?
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include "datatype_generated.h"

#include <vector>

struct Component {
    
    std::vector<value> location;
    std::vector<value> bandwidth;
    
    value scale_factor;
    value scale_factor_log;
    
    YAML::Node to_yaml() const;
    static std::unique_ptr<Component> from_yaml( const YAML::Node & node );
    
    void to_hdf5(HighFive::Group & group) const;
    static std::unique_ptr<Component> from_hdf5(const HighFive::Group & group);
};

std::vector<std::unique_ptr<Component>> components_from_flatbuffers(
    const fb_serialize::Kernels * kernels
);