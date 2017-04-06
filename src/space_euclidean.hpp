#pragma once

#include "space_base.hpp"
#include "kernel.hpp"

static const value DEFAULT_EUCLIDEAN_BANDWIDTH = 1.;
static const value DEFAULT_EUCLIDEAN_LOCATION = 0.;

class EuclideanSpace : public SpaceBase<EuclideanSpace> {
public:
    EuclideanSpace( std::vector<std::string> names, std::vector<value> bandwith = {}, std::vector<value> location = {});
    EuclideanSpace( std::vector<std::string> names, const Kernel & k, std::vector<value> bandwith = {}, std::vector<value> location = {});
    
    EuclideanSpace( const EuclideanSpace & other ) : SpaceBase<EuclideanSpace>(other), names_(other.names_), kernel_(other.kernel_->clone()) {}
    
    SpaceSpecification make_spec( std::vector<std::string> names, const Kernel & k );
    Component make_kernel(unsigned int n, std::vector<value> bw, std::vector<value> loc, const Kernel & k) const;
    
    static EuclideanSpace * fromYAML( const YAML::Node & node );
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
    
    VectorGrid * grid( const std::vector<std::vector<value>> & vectors, const std::vector<bool> & valid = {}, const std::vector<bool> & selection = {}) const;
    
protected:
    std::vector<std::string> names_;
    std::unique_ptr<Kernel> kernel_;
};
