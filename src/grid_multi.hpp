#pragma once

#include "grid_base.hpp"

SpaceSpecification space_from_grids( const std::vector<Grid*> & grids );
std::vector<long unsigned int> shape_from_grids( const std::vector<Grid*> & grids );
const std::vector<bool> valid_from_grids(const std::vector<Grid*> & grids,
    const std::vector<bool> & valid);

class MultiGrid : public GridBase<MultiGrid> {
public:
    // constructor
    MultiGrid( const std::vector<Grid*> & grids, const std::vector<bool> & valid );
    // copy constructor
    MultiGrid( const MultiGrid & other );
    
    //properties
    unsigned int ngrids() const;
    const Grid & subgrid(unsigned int index=0);
    
    // method to compute probability
    virtual void probability( const MultiSpace & space, value weight, 
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
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<Grid> from_hdf5(const HighFive::Group & group, 
        const SpaceSpecification & space, const std::vector<bool> & valid);
    
    
    virtual void at_index(const unsigned int * index, value * result) const {
        for (unsigned int k=0; k<grids_.size(); ++k) {
            grids_[k]->at_index(index, result);
            index  += grids_[k]->ndim();
            result += grids_[k]->ndim();
        }
    }
    
protected:
    std::vector<std::unique_ptr<Grid>> grids_;
};

