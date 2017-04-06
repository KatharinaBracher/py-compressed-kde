#pragma once

#include "space_base.hpp"

class MultiSpace : public SpaceBase<MultiSpace> {
public:
    MultiSpace( std::vector<const Space*> spaces );
    
    MultiSpace( const MultiSpace & other );
    
    SpaceSpecification make_spec( std::vector<const Space*> spaces ) const;
    
    unsigned int nchildren() const { return spaces_.size(); }
    
    const Space & child(unsigned int index) const {
        if (index>=nchildren()) {
            throw std::runtime_error("Invalid child space index,");
        }
        return *spaces_[index];
    }
    
    Component make_kernel( std::vector<const Space*> spaces ) const;
    
    static MultiSpace * fromYAML( const YAML::Node & node );
    virtual YAML::Node asYAML() const;
    
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log=false ) const override;
    
    virtual value mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const override;
    
    virtual void merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const override;
    
    virtual value probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const;
    virtual value log_probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const;
    
    virtual value partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const;
    
    Grid * grid( const std::vector<Grid*> & grids, const std::vector<bool> & valid = {}) const;
    
protected:
    std::vector<std::unique_ptr<Space>> spaces_;
};
