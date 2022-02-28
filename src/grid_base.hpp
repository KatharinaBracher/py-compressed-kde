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
#include "component.hpp"
#include "spacespec.hpp"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include "schema_generated.h"

#include <string>
#include <vector>
#include <fstream>

// forward declarations
class Space;
class CategoricalSpace;
class CircularSpace;
class EncodedSpace;
class EuclideanSpace;
class MultiSpace;

class Grid {
public:
    // constructor
    Grid( std::string klass, const SpaceSpecification & space, std::vector<long unsigned int> shape, const std::vector<bool> & valid );
    
    // clone
    virtual Grid* clone() const;
    
    // properties
    std::string klass() const;
    
    const std::vector<long unsigned int> & shape() const;
    unsigned int size() const;
    unsigned int ndim() const;
    
    const std::vector<bool> & valid() const;
    void set_valid(const std::vector<bool> & valid);
    unsigned int ninvalid() const;
    unsigned int nvalid() const;
    
    // space
    const SpaceSpecification & specification() const;
    
    // comparison
    friend bool operator==(const Grid& lhs, const Grid& rhs) {
        return lhs.shape()==rhs.shape() && lhs.specification()==rhs.specification();
        // NOTE: here we do not check if the values in the grid are the same!
    }
    
    //void probability( const Space & space, value weight, const Component & k, value * result );
    
    // methods to compute probability
    virtual void probability( const Space & space, value weight, 
        const value * loc, const value * bw, value * result );
    virtual void probability( const CategoricalSpace & space, value weight, 
        const value * loc, const value * bw, value * result );
    virtual void probability( const CircularSpace & space, value weight, 
        const value * loc, const value * bw, value * result );
    virtual void probability( const EncodedSpace & space, value weight, 
        const value * loc, const value * bw, value * result );
    virtual void probability( const EuclideanSpace & space, value weight, 
        const value * loc, const value * bw, value * result );
    virtual void probability( const MultiSpace & space, value weight, 
        const value * loc, const value * bw, value * result );
    
    //void partial_logp( const Space & space, const Component & k, const std::vector<bool> & selection, value * result );
    
    // methods to compute partial log probability
    virtual void partial_logp( const Space & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    virtual void partial_logp( const CategoricalSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    virtual void partial_logp( const CircularSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    virtual void partial_logp( const EncodedSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    virtual void partial_logp( const EuclideanSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    virtual void partial_logp( const MultiSpace & space, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result );
    
    //virtual void marginal( const Space & space, const Component & k, const std::vector<bool> & selection, value * result );
    
    // yaml
    YAML::Node to_yaml() const;
    virtual YAML::Node to_yaml_impl() const;
    void save_to_yaml( std::ostream & stream, bool flow=false ) const;
    void save_to_yaml( std::string path, bool flow=false ) const;
    
    // flatbuffers
    virtual flatbuffers::Offset<fb_serialize::Grid> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    virtual std::vector<flatbuffers::Offset<fb_serialize::FloatArray>> to_flatbuffers_data(flatbuffers::FlatBufferBuilder & builder) const;
    virtual std::vector<flatbuffers::Offset<fb_serialize::Grid>> to_flatbuffers_grids(flatbuffers::FlatBufferBuilder & builder) const;

    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    
    virtual void at_index(const unsigned int * index, value * result) const;
    
protected:
    std::string klass_;
    SpaceSpecification spec_;
    std::vector<long unsigned int> shape_;
    std::vector<bool> valid_;
    unsigned int ninvalid_;
};

template <typename T>
class GridBase : public Grid {
public:
    
    using Grid::Grid;
    
    virtual Grid* clone() const override {
        return new T(static_cast<T const&>(*this));
    }
};
