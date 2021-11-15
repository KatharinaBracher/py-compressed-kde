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
#include "spacespec.hpp"
#include "common.hpp"
#include <algorithm>

// constructors
DimSpecification::DimSpecification( std::string name, std::string type, 
    std::string extra )
    : name_(name), type_(type), extra_(extra) {
    
    hash_ = std::hash<std::string>{}( detail() );
}

// properties
std::string DimSpecification::detail() const {
    return name_ + "(" + type_ + ")[" + extra_ + "]";
}

std::string DimSpecification::name() const { return name_; }
std::string DimSpecification::type() const { return type_; }
std::string DimSpecification::extra() const { return extra_; }

size_t DimSpecification::hash() const { return hash_; }

// yaml
YAML::Node DimSpecification::to_yaml() const {
    YAML::Node node;
    node["name"] = name_;
    node["type"] = type_;
    node["extra"] = extra_;
    return node;
}

DimSpecification DimSpecification::from_yaml( const YAML::Node & node ) {
    
    return DimSpecification( 
        node["name"].as<std::string>(),
        node["type"].as<std::string>(),
        node["extra"].as<std::string>() );
}    

// hdf5
void DimSpecification::to_hdf5(HighFive::Group & group) const {
    
    HighFive::DataSet name = group.createDataSet<std::string>(
        "name", HighFive::DataSpace::From(name_));
    name.write(name_);
    
    HighFive::DataSet type = group.createDataSet<std::string>(
        "type", HighFive::DataSpace::From(type_));
    type.write(type_);
    
    HighFive::DataSet extra = group.createDataSet<std::string>(
        "extra", HighFive::DataSpace::From(extra_));
    extra.write(extra_);
    
}

DimSpecification DimSpecification::from_hdf5(const HighFive::Group & group) {
    
    std::string name, type, extra;
    
    HighFive::DataSet ds_name = group.getDataSet("name");
    ds_name.read(name);
    
    HighFive::DataSet ds_type = group.getDataSet("type");
    ds_type.read(type);
    
    HighFive::DataSet ds_extra = group.getDataSet("extra");
    ds_extra.read(extra);
    
    return DimSpecification(name, type, extra);
    
}


// constructors
SpaceSpecification::SpaceSpecification( const std::vector<DimSpecification> & dims )
    : dims_(dims) {
    
    if (!isunique( names() )) {
        throw std::runtime_error("Non-unique dimension names.");
    }
    
    update_hash();
}
    
SpaceSpecification::SpaceSpecification( DimSpecification dim )
    : SpaceSpecification( std::vector<DimSpecification>( {dim} ) ) {}

SpaceSpecification::SpaceSpecification() : hash_(0) {}

// properties
size_t SpaceSpecification::hash() const {
    return hash_;
}

void SpaceSpecification::update_hash() {
    hash_ = 0;
    for (auto & k : dims_) {
        hash_ = hash_combine( hash_, k.hash() );
    }
}

unsigned int SpaceSpecification::ndim() const {
    return dims_.size();
}

DimSpecification SpaceSpecification::dim(unsigned int index) const {
    if (index>=ndim()) {
        throw std::runtime_error("Dimension index out of bound.");
    }
    
    return dims_[index];
}

const std::vector<DimSpecification> & SpaceSpecification::dims() const {
    return dims_;
}

std::vector<std::string> SpaceSpecification::names() const {
    std::vector<std::string> names;
    
    for (auto & k : dims_) {
        names.push_back( k.name() );
    }
    
    return names;
}

std::vector<std::string> SpaceSpecification::types() const {
    std::vector<std::string> types;
    
    for (auto & k : dims_) {
        types.push_back( k.type() );
    }
    
    return types;
}

std::vector<std::string> SpaceSpecification::details() const {
    std::vector<std::string> details;
    
    for (auto & k : dims_) {
        details.push_back( k.detail() );
    }
    
    return details;
}

// methods
void SpaceSpecification::append( DimSpecification dim ) {

    append( std::vector<DimSpecification>( { dim } ) );
}

void SpaceSpecification::append( const std::vector<DimSpecification> & dims ) {
    
    unsigned int n = ndim();
    
    dims_.insert( dims_.end(), dims.cbegin(), dims.cend() );
    
    if (!isunique( names() )) {
        dims_.erase( dims_.begin() + n, dims_.end() );
        throw std::runtime_error("Non-unique dimension names.");
    }
    
    update_hash();
    
}

void SpaceSpecification::append( const SpaceSpecification & space ) {
    
    append( space.dims() );
}

void SpaceSpecification::prepend( DimSpecification dim ) {
    prepend( std::vector<DimSpecification>({ dim }) );
}

void SpaceSpecification::prepend( const std::vector<DimSpecification> & dims ) {
    
    unsigned int n = ndim();
    
    dims_.insert( dims_.begin(), dims.cbegin(), dims.cend() );
    
    if (!isunique( names() )) {
        dims_.erase( dims_.begin() + n, dims_.end() );
        throw std::runtime_error("Non-unique dimension names.");
    }
    
    update_hash();
    
}

void SpaceSpecification::prepend( const SpaceSpecification & space ) {
    prepend( space.dims() );
}

std::vector<bool> SpaceSpecification::selection( const SpaceSpecification & other ) const {
    
    std::vector<bool> sel( ndim(), false);
    
    bool b;
    
    unsigned int m = 0;
    
    for (unsigned int k=0; k<other.ndim(); ++k) {
        
        b = false;
        
        for ( ; m<ndim(); ++m ) {
            if (dim(m)==other.dim(k)) {
                sel[m] = true;
                m+=1;
                b = true;
                break;
            }
        }
        
        if (!b) { throw std::runtime_error("Not a proper subspace."); }
        
    }
    
    return sel;
    
}

bool SpaceSpecification::issubspace( const SpaceSpecification & other ) const {
    
    try {
        this->selection( other );
    } catch (std::exception & me) {
        return false;
    }
    
    return true;
    
}

SpaceSpecification SpaceSpecification::select( const std::vector<bool> & selection ) const {
    
    if (selection.size()!=ndim()) {
        throw std::runtime_error("Incorrect selection vector size.");
    }
    
    std::vector<DimSpecification> spec;
    
    for (unsigned int k=0; k<ndim(); ++k) {
        if (selection[k]) {
            spec.push_back( dims_[k] );
        }
    }
    
    return SpaceSpecification( spec );
    
}


// yaml
YAML::Node SpaceSpecification::to_yaml() const {
    YAML::Node node;
    for (auto & k : dims_) {
        node["dimensions"].push_back( k.to_yaml() );
    }
    return node;
}

SpaceSpecification SpaceSpecification::from_yaml( const YAML::Node & node ) {
    
    std::vector<DimSpecification> dspec;
    
    if (!node["dimensions"] || !node["dimensions"].IsSequence()) {
        throw std::runtime_error("Invalid speace specification.");
    }
    
    for (auto & k : node["dimensions"]) {
        dspec.push_back( DimSpecification::from_yaml( k ) );
    }
    
    return SpaceSpecification( dspec );
}


// hdf5
void SpaceSpecification::to_hdf5(HighFive::Group & group) const {
    
    unsigned int dim = 0;
    
    HighFive::Attribute attr = group.createAttribute<unsigned int>(
            "ndim", HighFive::DataSpace::From(dim));
    attr.write(dims_.size());
    
    for (auto & k : dims_) {
        
        HighFive::Group subgroup = group.createGroup("dim" + std::to_string(dim));
        
        k.to_hdf5(subgroup);
        
        ++dim;
    }
}

SpaceSpecification SpaceSpecification::from_hdf5(const HighFive::Group & group) {
    
    std::vector<DimSpecification> dspec;
    
    unsigned int ndim;
    
    HighFive::Attribute attr = group.getAttribute("ndim");
    attr.read(ndim);
    
    for (unsigned int k=0; k<ndim; ++k) {
        dspec.push_back(DimSpecification::from_hdf5(
            group.getGroup("dim" + std::to_string(k))));
    }
    
    return SpaceSpecification(dspec);
}

