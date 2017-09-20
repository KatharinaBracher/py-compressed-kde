#pragma once

#include "grid_base.hpp"

std::vector<long unsigned int> shape_from_array_args( const std::vector<long unsigned int> & shape, unsigned int array_size, unsigned int ndim );

class ArrayGrid : public GridBase<ArrayGrid> {
public:
    // constructor
    ArrayGrid( const std::vector<value> & array, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
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
    static std::unique_ptr<Grid> from_yaml( const YAML::Node & node, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<Grid> from_hdf5(const HighFive::Group & group, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
    
protected:
    std::vector<value> array_;
};
