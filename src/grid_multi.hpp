#pragma once

#include "grid_base.hpp"

std::vector<long unsigned int> shape_from_grids( const std::vector<Grid*> & grids );

class MultiGrid : public GridBase<MultiGrid> {
public:
    MultiGrid( const std::vector<Grid*> & grids, const std::vector<bool> & valid );
    
    MultiGrid( const MultiGrid & other );
    
    const Grid & subgrid(unsigned int index=0);
    
    virtual YAML::Node asYAML() const;
    static Grid * fromYAML( const YAML::Node & node, const SpaceSpecification & space, const std::vector<bool> & valid );
    
    virtual void probability( const MultiSpace & space, value weight, const value * loc, const value * bw, value * result ) override;
    
    virtual void partial_logp( const CategoricalSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const CircularSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EncodedSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EuclideanSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const MultiSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) override;
    
protected:
    std::vector<std::unique_ptr<Grid>> grids_;
};

