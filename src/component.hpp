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
