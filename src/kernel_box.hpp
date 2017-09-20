#pragma once 

#include "kernel_base.hpp"

static const double BOX_KERNEL_FACTOR = 1.7400570569722662;

value box_scale_factor( unsigned ndim, value det, bool log );

class BoxKernel : public KernelBase<BoxKernel> {
public:
    // constructor
    BoxKernel();
    
    // metods    
    virtual value scale_factor( unsigned int n, value * bw, bool log ) const;
    virtual value scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const;
    
    virtual value probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value probability( value dsquared ) const;
    
    virtual value log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value log_probability( value dsquared ) const;
    
    virtual value partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const;
    
    // yaml
    virtual YAML::Node to_yaml_impl() const;
    static std::unique_ptr<BoxKernel> from_yaml( const YAML::Node & node );
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<BoxKernel> from_hdf5(const HighFive::Group & group);

};
