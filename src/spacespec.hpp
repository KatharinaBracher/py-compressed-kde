#pragma once

#include "yaml-cpp/yaml.h"
#include <string>


class DimSpecification {
public:
    DimSpecification( std::string name, std::string type, std::string extra );
    
    YAML::Node toYAML() const;
    static DimSpecification fromYAML( const YAML::Node & node );
    
    std::string detail() const;
    
    std::string name() const;
    std::string type() const;
    std::string extra() const;
    
    size_t hash() const;
    
    friend bool operator==(const DimSpecification& lhs, const DimSpecification& rhs) {
        return lhs.hash()==rhs.hash();
    }
    
protected:    
    std::string name_;
    std::string type_;
    std::string extra_; // or vector?
    size_t hash_;
};

class SpaceSpecification {
public:
    SpaceSpecification( const std::vector<DimSpecification> & dims );
    SpaceSpecification( DimSpecification dim );
    SpaceSpecification();
    
    YAML::Node toYAML() const;
    static SpaceSpecification fromYAML( const YAML::Node & node );
    
    size_t hash() const;
    
    unsigned int ndim() const;
    DimSpecification dim(unsigned int index=0) const;
    const std::vector<DimSpecification> & dims() const;
    
    std::vector<std::string> names() const;
    std::vector<std::string> types() const;
    std::vector<std::string> details() const;
    
    void append( DimSpecification dim );
    void append( const std::vector<DimSpecification> & dims );
    void append( const SpaceSpecification & dims );
    
    void prepend( DimSpecification dim );
    void prepend( const std::vector<DimSpecification> & dims );
    void prepend( const SpaceSpecification & dims );
    
    std::vector<bool> selection( const SpaceSpecification & other ) const;
    bool issubspace( const SpaceSpecification & other ) const;
    
    SpaceSpecification select( const std::vector<bool> & selection ) const;
    
    friend bool operator==(const SpaceSpecification& lhs, const SpaceSpecification& rhs) {
        return lhs.hash()==rhs.hash();
    }

protected:
    void update_hash();
    
protected:
    std::vector<DimSpecification> dims_;
    size_t hash_;
};
