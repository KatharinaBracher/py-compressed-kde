#include "kernel_gaussian.hpp"

#include <cmath>
#include <limits>

value gaussian_scale_factor( unsigned int ndim, value det, value cutoff, bool log ) {
    
    value s;
    
    if (log) {
        s = -fastlog( det * std::pow(M_PI*2, 0.5*ndim));
        s -= ndim * fastlog( 1 - std::erfc( cutoff/SQRT2 ) );
    } else {
        s = 1./( det * std::pow(M_PI*2,0.5*ndim) );
        s /= std::pow(1 - std::erfc( cutoff/SQRT2 ), ndim );
    }
    
    return s;
}

GaussianKernel::GaussianKernel( value cutoff ) : KernelBase<GaussianKernel>(KernelType::Gaussian), cutoff_(cutoff), cutoff_squared_(cutoff*cutoff) {}

value GaussianKernel::cutoff() const { return cutoff_; }
void GaussianKernel::set_cutoff( value v ) { cutoff_=v; cutoff_squared_=v*v; }

YAML::Node GaussianKernel::asYAML() const {
    YAML::Node node;
    node["cutoff"] = cutoff_;
    return node;
}
GaussianKernel * GaussianKernel::fromYAML( const YAML::Node & node ) {
    return new GaussianKernel( node["cutoff"].as<value>( DEFAULT_GAUSSIAN_CUTOFF ) );
}

value GaussianKernel::scale_factor( unsigned int n, value * bw, bool log ) const {
    value det = std::accumulate( bw, bw+n, 1., std::multiplies<value>() );
    return gaussian_scale_factor( n, det, cutoff_, log );
}
value GaussianKernel::scale_factor( unsigned int n, value * bw, bool log, std::vector<bool>::const_iterator selection ) const {
    unsigned int ndim = 0;
    value det = 1.;
    
    for (unsigned int d=0; d<n; ++d) {
        if (*selection) {
            det *= bw[d];
            ++ndim;
        }
        ++selection;
    }
    return gaussian_scale_factor( ndim, det, cutoff_, log );
}
value GaussianKernel::probability( unsigned int n, const value * loc, const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / (*bw++);
        d += (tmp*tmp);
        
        if (d>=cutoff_squared_) { return 0.; }
        
    }
    
    return fastexp( -0.5*d );
}
value GaussianKernel::probability( value dsquared ) const {
    if (dsquared>=cutoff_squared_) { return 0.; }
    else { return fastexp( -0.5*dsquared ); }
}
value GaussianKernel::log_probability( unsigned int n, const value * loc, const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / (*bw++);
        d += (tmp*tmp);
        
        if (d>=cutoff_squared_) { return -std::numeric_limits<value>::infinity(); }
        
    }
    
    return -0.5*d;
}
value GaussianKernel::log_probability( value dsquared ) const {
    if (dsquared>=cutoff_squared_) { return -std::numeric_limits<value>::infinity(); }
    else { return -0.5*dsquared; }
}
value GaussianKernel::partial_logp( unsigned int n, const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection) const {
    
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        if (*selection++) {
            tmp = (*point - *loc) / *bw;
            d += (tmp*tmp);
            if (d>=cutoff_squared_) { return -std::numeric_limits<value>::infinity(); }
        }
        ++point;
        ++loc;
        ++bw;
    }
    
    return -0.5*d ;
    
}
