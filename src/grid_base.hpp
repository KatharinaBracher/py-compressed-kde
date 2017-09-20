#pragma once

#include "common.hpp"
#include "component.hpp"
#include "spacespec.hpp"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include "highfive/H5Group.hpp"

#include <string>
#include <vector>

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
    
    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    
    
protected:
    std::string klass_;
    SpaceSpecification spec_;
    std::vector<long unsigned int> shape_;
    std::vector<bool> valid_;
    unsigned int nvalid_;
};

template <typename T>
class GridBase : public Grid {
public:
    
    using Grid::Grid;
    
    virtual Grid* clone() const override {
        return new T(static_cast<T const&>(*this));
    }
};
