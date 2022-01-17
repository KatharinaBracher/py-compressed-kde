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

#include "space_base.hpp"

class CategoricalSpace : public SpaceBase<CategoricalSpace> {
public:
    // constructor
    CategoricalSpace( std::string name, std::vector<std::string> labels, unsigned int category=0 );
    
    SpaceSpecification make_spec( std::string, std::vector<std::string> labels );
    Component make_kernel(unsigned int category) const;
    
    // grid construction
    Grid * grid() const;
    
    // methods
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log=false ) const override;
    
    virtual value mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const override;
    
    virtual void merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const override;
    
    virtual value probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const override;
    virtual value log_probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const override;
    
    virtual value partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const;
    
    // yaml
    static std::unique_ptr<CategoricalSpace> from_yaml( const YAML::Node & node );
    virtual YAML::Node to_yaml_impl() const;
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::SpaceData> to_flatbuffers_impl(flatbuffers::FlatBufferBuilder & builder) const;
    static std::unique_ptr<CategoricalSpace> from_flatbuffers(const fb_serialize::Space * space);

    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<CategoricalSpace> from_hdf5(const HighFive::Group & group);
    
    virtual void distance( const value * x, const value * y, value * result ) const {
        if (static_cast<unsigned int>(*x) == static_cast<unsigned int>(*y)) { 
            *result = 0.;
        } else {
            *result = std::numeric_limits<value>::infinity();
        }
    }
    
protected:
    std::vector<std::string> labels_;
};
