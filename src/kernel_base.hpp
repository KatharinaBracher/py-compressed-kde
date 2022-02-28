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

#include <string>
#include <stdexcept>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include "yaml-cpp/yaml.h"

enum class KernelType { Gaussian, Epanechnikov, Box };

std::string kerneltype_tostring( KernelType k );
KernelType kerneltype_fromstring( std::string k );

class Kernel {
public:
    // constructor
    Kernel( KernelType k );
    
    // clone
    virtual Kernel * clone() const { throw std::runtime_error("not implemented"); }
    virtual ~Kernel(){};

    KernelType type() const;
    virtual std::string to_string() const { return kerneltype_tostring(type()); }
    
    // methods
    virtual value scale_factor( unsigned int n, value * bw, bool log ) const;
    virtual value scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const;
    
    virtual value probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value probability( value dsquared ) const;
    
    virtual value log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value log_probability( value dsquared ) const;
    
    virtual value partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const;
    
    // yaml
    YAML::Node to_yaml() const;
    virtual YAML::Node to_yaml_impl() const;
    
    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    

    
protected:
    KernelType type_;
};

template <typename T>
class KernelBase : public Kernel {
public:
    
    using Kernel::Kernel;

    virtual Kernel* clone() const override {
        return new T(static_cast<T const&>(*this));
    }
};    

