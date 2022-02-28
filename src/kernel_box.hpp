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

static const double BOX_KERNEL_FACTOR = 1.7400570569722662;

value box_scale_factor( unsigned ndim, value det, bool log );

class BoxKernel : public KernelBase<BoxKernel> {
public:
    // constructor
    BoxKernel();
    
    // metods    
    virtual value scale_factor( unsigned int n, value * bw, bool log ) const;
    virtual value scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const;
    
    virtual value probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value probability( value dsquared ) const;
    
    virtual value log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value log_probability( value dsquared ) const;
    
    virtual value partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const;
    
    // yaml
    virtual YAML::Node to_yaml_impl() const;
    static std::unique_ptr<BoxKernel> from_yaml( const YAML::Node & node );
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::Kernel> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const override;
    static std::unique_ptr<BoxKernel> from_flatbuffers(const fb_serialize::BoxKernel * kernel);

    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<BoxKernel> from_hdf5(const HighFive::Group & group);

};
