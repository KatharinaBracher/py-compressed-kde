#pragma once

#include "space_base.hpp"

class CategoricalSpace : public SpaceBase<CategoricalSpace> {
public:
    // constructor
    CategoricalSpace( std::string name, std::vector<std::string> labels, unsigned int category=0 );
    
    SpaceSpecification make_spec( std::string, std::vector<std::string> labels );
    Component make_kernel(unsigned int category) const;
    
    // grid construction
    Grid * grid() const;
    
    // methods
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log=false ) const override;
    
    virtual value mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const override;
    
    virtual void merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const override;
    
    virtual value probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const override;
    virtual value log_probability( const value * loc, const value * bw, const value * point ) const override;
    virtual void log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const override;
    
    virtual value partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const;
    
    // yaml
    static std::unique_ptr<CategoricalSpace> from_yaml( const YAML::Node & node );
    virtual YAML::Node to_yaml_impl() const;
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<CategoricalSpace> from_hdf5(const HighFive::Group & group);
    
    
    
protected:
    std::vector<std::string> labels_;
};
