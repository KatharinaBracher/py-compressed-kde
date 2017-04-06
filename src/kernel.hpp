#pragma once

#include "kernel_base.hpp"
#include "kernel_gaussian.hpp"
#include "kernel_epanechnikov.hpp"
#include "kernel_box.hpp"
#include "kernel_vonmises.hpp"

#include "yaml-cpp/yaml.h"

Kernel * kernel_from_YAML( const YAML::Node & node );


//#include "yaml-cpp/yaml.h"

//#include <vector>
//#include <cmath>
//#include <string>
//#include <stdint.h>
//#include <algorithm>











