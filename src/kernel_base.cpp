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
#include "kernel_base.hpp"

std::string kerneltype_tostring( KernelType k ) {
    std::string s;
    switch(k) {
        case KernelType::Gaussian : s="gaussian"; break;
        case KernelType::Epanechnikov : s="epanechnikov"; break;
        case KernelType::Box : s="box"; break;
    }
    return s;
}

KernelType kerneltype_fromstring( std::string s ) {
    KernelType k;
    
     if (s=="gaussian") {
        k = KernelType::Gaussian;
    } else if (s=="epanechnikov") {
        k = KernelType::Epanechnikov;
    } else if (s=="box") {
        k = KernelType::Box;
    } else {
        throw std::runtime_error("Unknown kernel type.");
    }
    
    return k;
}

// constructor
Kernel::Kernel( KernelType k ) : type_(k) {}

KernelType Kernel::type() const { return type_; }

// methods
value Kernel::scale_factor( unsigned int n, value * bw, bool log ) const {
    return 1.;
}
value Kernel::scale_factor( unsigned int n, value * bw, bool log, 
    std::vector<bool>::const_iterator selection ) const {
    return 1;
}

value Kernel::probability( unsigned int n, const value * loc, const value * bw, 
    const value * point ) const { 
    return 0.;
}
value Kernel::probability( value dsquared ) const { return 0.; }

value Kernel::log_probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
    return -std::numeric_limits<value>::infinity();
}
value Kernel::log_probability( value dsquared ) const { 
    return -std::numeric_limits<value>::infinity();
}

value Kernel::partial_logp( unsigned int n, const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection) const { 
    return -std::numeric_limits<value>::infinity();
}


// yaml
YAML::Node Kernel::to_yaml() const {
    YAML::Node node;
    node["type"] = kerneltype_tostring(  type_ );
    node["info"] = to_yaml_impl();
    return node;
}

YAML::Node Kernel::to_yaml_impl() const {
    throw std::runtime_error("Not implemented.");
}


// hdf5
void Kernel::to_hdf5(HighFive::Group & group) const {
    std::string type_str = kerneltype_tostring(type_);
    HighFive::DataSet ds_type = group.createDataSet<std::string>("type", 
        HighFive::DataSpace::From(type_str));
    ds_type.write(type_str);
    
    HighFive::Group subgroup = group.createGroup("info");
    this->to_hdf5_impl(subgroup);
}

void Kernel::to_hdf5_impl(HighFive::Group & group) const {
    throw std::runtime_error("Not implemented.");
}

