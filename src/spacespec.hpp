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

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include "datatype_generated.h"

#include "yaml-cpp/yaml.h"
#include <string>


class DimSpecification {
public:
    // constructor
    DimSpecification( std::string name, std::string type, std::string extra );
    
    // properties
    std::string detail() const;
    
    std::string name() const;
    std::string type() const;
    std::string extra() const;
    
    size_t hash() const;
    
    // comparison
    friend bool operator==(const DimSpecification& lhs, 
        const DimSpecification& rhs) {
        return lhs.hash()==rhs.hash();
    }
    
    // yaml
    YAML::Node to_yaml() const;
    static DimSpecification from_yaml( const YAML::Node & node );
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::Dimension> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    static DimSpecification from_flatbuffers(const fb_serialize::Dimension * spec);

    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    static DimSpecification from_hdf5(const HighFive::Group & group);
    
    
protected:    
    std::string name_;
    std::string type_;
    std::string extra_; // or vector?
    size_t hash_;
};

class SpaceSpecification {
public:
    // constructors
    SpaceSpecification( const std::vector<DimSpecification> & dims );
    SpaceSpecification( DimSpecification dim );
    SpaceSpecification();
    
    // properties
    size_t hash() const;
    
    unsigned int ndim() const;
    DimSpecification dim(unsigned int index=0) const;
    const std::vector<DimSpecification> & dims() const;
    
    std::vector<std::string> names() const;
    std::vector<std::string> types() const;
    std::vector<std::string> details() const;
    
    // methods
    void append( DimSpecification dim );
    void append( const std::vector<DimSpecification> & dims );
    void append( const SpaceSpecification & dims );
    
    void prepend( DimSpecification dim );
    void prepend( const std::vector<DimSpecification> & dims );
    void prepend( const SpaceSpecification & dims );
    
    std::vector<bool> selection( const SpaceSpecification & other ) const;
    bool issubspace( const SpaceSpecification & other ) const;
    
    SpaceSpecification select( const std::vector<bool> & selection ) const;
    
    // comparison
    friend bool operator==(const SpaceSpecification& lhs, const SpaceSpecification& rhs) {
        return lhs.hash()==rhs.hash();
    }
    
    // yaml
    YAML::Node to_yaml() const;
    static SpaceSpecification from_yaml( const YAML::Node & node );
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::SpaceSpecification> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    static SpaceSpecification from_flatbuffers(const fb_serialize::SpaceSpecification * spec);

    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    static SpaceSpecification from_hdf5(const HighFive::Group & group);
    
protected:
    void update_hash();
    
protected:
    std::vector<DimSpecification> dims_;
    size_t hash_;
};
