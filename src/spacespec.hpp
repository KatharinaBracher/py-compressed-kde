#pragma once

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

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
    
    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    static SpaceSpecification from_hdf5(const HighFive::Group & group);
    
protected:
    void update_hash();
    
protected:
    std::vector<DimSpecification> dims_;
    size_t hash_;
};
