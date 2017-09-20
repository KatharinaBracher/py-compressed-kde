#pragma once

#include "space_base.hpp"
#include "kernel.hpp"

static const value DEFAULT_ENCODED_BANDWIDTH = 1.;
static const unsigned int DEFAULT_ENCODED_INDEX = 0;
static const unsigned int DEFAULT_ENCODED_GRID_DELTA = 1;

class EncodedSpace : public SpaceBase<EncodedSpace> {
public:
    // constructors
    EncodedSpace( std::string name, std::vector<value> & lut, 
        value bandwidth = DEFAULT_ENCODED_BANDWIDTH, 
        unsigned int index = DEFAULT_ENCODED_INDEX);
    EncodedSpace( std::string name, std::vector<value> & lut, const Kernel & k, 
        value bandwidth = DEFAULT_ENCODED_BANDWIDTH, 
        unsigned int index = DEFAULT_ENCODED_INDEX);
    
    // copy constructor
    EncodedSpace( const EncodedSpace & other )
        : SpaceBase<EncodedSpace>(other), nlut_(other.nlut_), lut_(other.lut_), 
          kernel_(other.kernel_->clone()) {}
    
    SpaceSpecification make_spec( std::string name, std::vector<value> & lut, 
        const Kernel & k );
    Component make_kernel( value bw, unsigned int idx, const Kernel & k) const;
    
    // grid construction
    Grid * grid(unsigned int delta=DEFAULT_ENCODED_GRID_DELTA) const;
    
    // methods
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, 
        value * bw, bool log=false ) const override;
    
    virtual value mahalanobis_distance_squared( const value * refloc, 
        const value * refbw, const value * targetloc, value threshold) const override;
    
    virtual void merge( value w1, value * loc1, value * bw1, value w2, 
        const value * loc2, const value * bw2 ) const override;
    
    virtual value probability( const value * loc, const value * bw, 
        const value * point ) const override;
    virtual void probability( const value * loc, const value * bw, 
        const value * points, unsigned int n, value * result ) const;
    virtual value log_probability( const value * loc, const value * bw, 
        const value * point ) const override;
    virtual void log_probability( const value * loc, const value * bw, 
        const value * points, unsigned int n, value * result ) const;
    
    virtual value partial_logp( const value * loc, const value * bw, 
        const value * point, std::vector<bool>::const_iterator selection ) const;
    
    // yaml
    static std::unique_ptr<EncodedSpace> from_yaml( const YAML::Node & node );
    virtual YAML::Node to_yaml_impl() const;
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<EncodedSpace> from_hdf5(const HighFive::Group & group);
    
protected:
    unsigned int nlut_;
    std::shared_ptr<std::vector<value>> lut_; // set once, read by many
    std::unique_ptr<Kernel> kernel_;
};
