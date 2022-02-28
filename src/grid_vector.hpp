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

#include "grid_base.hpp"

std::vector<long unsigned int> shape_from_vectors( const std::vector<std::vector<value>> & vectors );
//SpaceSpecification space_from_grids( const std::vector<Grid*> & grids );

class VectorGrid : public GridBase<VectorGrid> {
public:
    // constructor
    VectorGrid( const std::vector<std::vector<value>> & vectors, 
        const SpaceSpecification & space, const std::vector<bool> & valid );
    
    // methods to compute probability
    virtual void probability( const CategoricalSpace & space, value weight, 
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const CircularSpace & space, value weight, 
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const EncodedSpace & space, value weight, 
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const EuclideanSpace & space, value weight, 
        const value * loc, const value * bw, value * result ) override;
    
    // methods to compute partial log probability
    virtual void partial_logp( const CategoricalSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const CircularSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EncodedSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EuclideanSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const MultiSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) override;
    
    // yaml
    virtual YAML::Node to_yaml_impl() const;
    static std::unique_ptr<Grid> from_yaml( const YAML::Node & node, 
        const SpaceSpecification & space, const std::vector<bool> & valid );
    
    // flatbuffers
    virtual std::vector<flatbuffers::Offset<fb_serialize::FloatArray>> to_flatbuffers_data(flatbuffers::FlatBufferBuilder & builder) const;
    static std::unique_ptr<VectorGrid> from_flatbuffers(const fb_serialize::Grid * grid);

    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<Grid> from_hdf5(const HighFive::Group & group, const SpaceSpecification & space, const std::vector<bool> & valid);
    
    virtual void at_index(const unsigned int * index, value * result) const {
        for (auto & v : vectors_) {
            if (*index >= v.size()) {
                *result = std::numeric_limits<value>::quiet_NaN();
            } else {
                *result = v[*index];
            }
            ++index;
            ++result;
        }
    }
    
protected:
    std::vector<std::vector<value>> vectors_;
    std::vector<std::vector<value>> ptemp_;
};
