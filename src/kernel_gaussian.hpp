#pragma once

#include "kernel_base.hpp"

static const double DEFAULT_GAUSSIAN_CUTOFF = 3.;

value gaussian_scale_factor( unsigned int ndim, value det, value cutoff, bool log );

class GaussianKernel : public KernelBase<GaussianKernel> {
public:
    GaussianKernel( value cutoff = DEFAULT_GAUSSIAN_CUTOFF );
    
    value cutoff() const;
    void set_cutoff( value v );
    
    virtual std::string to_string() const { return kerneltype_tostring(type()) + "(cutoff=" + std::to_string(cutoff_) + ")"; }
    
    virtual YAML::Node asYAML() const;
    static GaussianKernel * fromYAML( const YAML::Node & node );
    
    virtual value scale_factor( unsigned int n, value * bw, bool log ) const;
    virtual value scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const;
    virtual value probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value probability( value dsquared ) const;
    virtual value log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const;
    virtual value log_probability( value dsquared ) const;
    virtual value partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const;
    
    
protected:
    value cutoff_;
    value cutoff_squared_;
};
