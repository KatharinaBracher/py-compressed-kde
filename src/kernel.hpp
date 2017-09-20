#pragma once

#include "kernel_base.hpp"
#include "kernel_gaussian.hpp"
#include "kernel_epanechnikov.hpp"
#include "kernel_box.hpp"
#include "kernel_vonmises.hpp"

#include "yaml-cpp/yaml.h"

// yaml
std::unique_ptr<Kernel> kernel_from_yaml( const YAML::Node & node );

// hdf5
std::unique_ptr<Kernel> kernel_from_hdf5( const HighFive::Group & group );












