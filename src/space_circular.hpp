#pragma once

#include "space.hpp"

static const value DEFAULT_KAPPA = 5.;
static const value DEFAULT_MU = 0.;
static const unsigned int DEFAULT_CIRCULAR_GRID_SIZE = 24;
static const value DEFAULT_CIRCULAR_GRID_OFFSET = 0.;

class CircularSpace : public SpaceBase<CircularSpace> {
public:
    // constructor
    CircularSpace( std::string name, value kappa = 5., value mu = 0. );
    
    SpaceSpecification make_spec( std::string );
    Component make_kernel( value kappa, value mu) const;
    
    // grid construction
    Grid * grid(unsigned int n=DEFAULT_CIRCULAR_GRID_SIZE, value offset=DEFAULT_CIRCULAR_GRID_OFFSET) const;
    
    // methods
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log=false ) const override;
    
    virtual value mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const override;
    
    virtual void merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const override;
    
    virtual value probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const;
    virtual value log_probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const;
    
    virtual value partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const;
    
    // yaml
    static std::unique_ptr<CircularSpace> from_yaml( const YAML::Node & node );
    virtual YAML::Node to_yaml_impl() const;
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<CircularSpace> from_hdf5(const HighFive::Group & group);
    
    virtual void distance( const value * x, const value * y, value * result ) const {
        *result = circular_difference(*y, *x);
    }
    
};
