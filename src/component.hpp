#pragma once

#include "common.hpp"
#include "yaml-cpp/yaml.h"
#include <vector>

struct Component {
    
    std::vector<value> location;
    std::vector<value> bandwidth;
    
    value scale_factor;
    value scale_factor_log;
    
    YAML::Node toYAML() const;
    static Component* fromYAML( const YAML::Node & node );
};
