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
#include "spacespec.hpp"
#include "component.hpp"
#include "grid.hpp"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include "yaml-cpp/yaml.h"

#include "schema_generated.h"

#include <memory>
#include <string>
#include <fstream>

class Space {
public:
    // constructor
    Space( std::string klass, const SpaceSpecification & spec, const Component k );
    // copy constructor
    Space( const Space & other );
    // clone
    virtual Space* clone() const;
    
    // comparison
    friend bool operator==(const Space& lhs, const Space& rhs) {
        return lhs.specification()==rhs.specification();
    }

    virtual ~Space(){};

    // properties
    std::vector<bool> selection( const Space& space ) const;
    bool issubspace( const Space& space ) const;
    
    std::string klass() const;
    unsigned int ndim() const;
    unsigned int nbw() const;
    //const std::vector<std::string> dimension_names() const;
    //const std::vector<std::string> dimension_types() const;
    //const std::vector<std::string> dimension_details() const;
    
    const SpaceSpecification & specification() const;
    const Component & default_kernel() const;
    
    void set_default_kernel(const Component & k) {
        // check if dimensionalities are OK
        if (k.location.size()!=ndim() || k.bandwidth.size()!=nbw()) {
            throw std::runtime_error("Cannot set default kernel. Dimensions do not match.");
        }
        
        default_kernel_ = k;
    }
    
    std::unique_ptr<Component> kernel() const;
    
    // kernel construction
    template <typename ITER>
    std::unique_ptr<Component> kernel( ITER it ) const {
        auto k =  std::unique_ptr<Component>( new Component( default_kernel_ ) );
        std::copy( it, it + ndim(), std::begin(k->location) );
        return k;
    }
    
    // methods
    void update_scale_factor (Component & k) const;
    value compute_scale_factor( Component & k, bool log = false ) const;
    value compute_scale_factor( Component & k, const std::vector<bool> & selection, 
        bool log = false ) const;
    
    virtual value compute_scale_factor( value * bw, bool log = false ) const { 
        return 1.;
    }
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, 
        value * bw, bool log=false ) const { 
        return 1.;
    }
    
    value mahalanobis_distance_squared( const Component & reference, 
        const Component & target, value threshold) const;
    virtual value mahalanobis_distance_squared( const value * refloc, 
        const value * refbw, const value * targetloc, value threshold) const;
    
    void merge( value w1, Component & first, value w2, const Component & second ) const;
    virtual void merge( value w1, value * loc1, value * bw1, value w2, 
        const value * loc2, const value * bw2 ) const;
    
    value probability( const Component & k, const value * point ) const;
    virtual value probability( const value * loc, const value * bw, 
        const value * point ) const;
    virtual void probability( const value * loc, const value * bw, 
        const value * points, unsigned int n, value * result ) const;
    
    value log_probability( const Component & k, const value * point ) const;
    virtual value log_probability( const value * loc, const value * bw,
        const value * point ) const;
    virtual void log_probability( const value * loc, const value * bw, 
        const value * points, unsigned int n, value * result ) const;
    
    void probability( Grid & grid, value weight, const Component & k, value * result ) {
        probability( grid, weight, k.location.data(), k.bandwidth.data(), result );
    }
    
    virtual void probability( Grid & grid, value weight, const value * loc, 
        const value * bw, value * result ) const {
        throw std::runtime_error("Space::probability(Grid,...) not implemented.");
    }
    
    value partial_logp( const Component & k, const value * point, 
        const std::vector<bool> & selection ) const;
    virtual value partial_logp( const value * loc, const value * bw, 
        const value * point, std::vector<bool>::const_iterator selection ) const;
    
    void partial_logp( Grid & grid, std::vector<bool> & selection, 
        value factor, const Component & k, value * result ) {
        partial_logp( grid, selection.cbegin(), factor, k.location.data(), k.bandwidth.data(), result );
    }
    virtual void partial_logp( Grid & grid, 
        std::vector<bool>::const_iterator selection, value factor, 
            const value * loc, const value * bw, value * result ) const {
        throw std::runtime_error("Space::partial_logp(Grid,...) not implemented.");
    }
    
    // yaml
    YAML::Node to_yaml() const;
    virtual YAML::Node to_yaml_impl() const;
    void save_to_yaml( std::ostream & stream, bool flow=false ) const;
    void save_to_yaml( std::string path, bool flow=false ) const;
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::Space> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    virtual flatbuffers::Offset<fb_serialize::SpaceData> to_flatbuffers_impl(flatbuffers::FlatBufferBuilder & builder) const;

    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    
    virtual void distance( const value * x, const value * y, value * result ) const {
        throw std::runtime_error("Not implemented");
    }
    
protected:
    std::string klass_;
    SpaceSpecification spec_;
    Component default_kernel_;
    size_t hash_;
};

template <typename T>
class SpaceBase : public Space {
public:
    
    using Space::Space;

    virtual Space* clone() const override {
        return new T(static_cast<T const&>(*this));
    }
    
    virtual void probability( Grid & grid, value weight, const value * loc, 
        const value * bw, value * result ) const override final{
        grid.probability( *(static_cast<const T*>(this)), weight, loc, bw, result );
    }
    
    virtual void partial_logp( Grid & grid, 
        std::vector<bool>::const_iterator selection, value factor, 
        const value * loc, const value * bw, value * result ) const override final {
        
        grid.partial_logp( *(static_cast<const T*>(this)), selection, 
            factor, loc, bw, result );
    }
};    


