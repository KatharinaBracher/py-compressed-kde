#pragma once

#include "space_base.hpp"
#include "kernel.hpp"

static const value DEFAULT_ENCODED_BANDWIDTH = 1.;
static const unsigned int DEFAULT_ENCODED_INDEX = 0;
static const unsigned int DEFAULT_ENCODED_GRID_DELTA = 1;

value nearest(const std::vector<value> & v, value x);
unsigned int nearest_index(const std::vector<value> & v, value x);

class EncodedSpace : public SpaceBase<EncodedSpace> {
public:
    // constructors
    EncodedSpace( std::string name, const std::vector<value> & lut, 
        value bandwidth = DEFAULT_ENCODED_BANDWIDTH, 
        unsigned int index = DEFAULT_ENCODED_INDEX);
    EncodedSpace( std::string name, const std::vector<value> & points,
        const std::vector<value> & lut, 
        value bandwidth = DEFAULT_ENCODED_BANDWIDTH, 
        unsigned int index = DEFAULT_ENCODED_INDEX);
    EncodedSpace( std::string name, const std::vector<value> & points,
        const std::vector<value> & lut,
        const Kernel & k, value bandwidth = DEFAULT_ENCODED_BANDWIDTH, 
        unsigned int index = DEFAULT_ENCODED_INDEX);
    
    // copy constructor
    EncodedSpace( const EncodedSpace & other )
        : SpaceBase<EncodedSpace>(other), use_index_(other.use_index_),
          nlut_(other.nlut_), points_(other.points_),
          lut_(other.lut_), kernel_(other.kernel_->clone()) {}
    
    SpaceSpecification make_spec( std::string name, const std::vector<value> & lut, const Kernel & k );
    Component make_kernel( value bw, const std::vector<value> & points, unsigned int idx, const Kernel & k) const;
    
    // grid construction
    Grid * grid(unsigned int delta=DEFAULT_ENCODED_GRID_DELTA) const;
    Grid * grid(const std::vector<value> & v, const std::vector<bool> & valid = {}) const;
        
    // methods
    virtual value compute_scale_factor( value * bw, bool log = false ) const override;
    virtual value compute_scale_factor( std::vector<bool>::const_iterator selection, 
        value * bw, bool log=false ) const override;
    
    unsigned int get_index(value x) const;
    
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
    
    bool use_index() const { return use_index_; }
    
    virtual void distance( const value * x, const value * y, value * result ) const {
        *result = std::sqrt( (*lut_)[get_index(*x) + get_index(*y)*nlut_] );
    }
    
protected:
    bool use_index_;
    unsigned int nlut_;
    std::shared_ptr<std::vector<value>> points_; // set once, read by many
    std::shared_ptr<std::vector<value>> lut_; // set once, read by many
    std::unique_ptr<Kernel> kernel_;
};
