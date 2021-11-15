// ---------------------------------------------------------------------
// This file is part of the compressed decoder library.
//
// Copyright (C) 2020 - now Neuro-Electronics Research Flanders
//
// The compressed decoder library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The compressed decoder library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with falcon-core. If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------
#include "kernel_box.hpp"

#include <cmath>
#include <limits>

value box_scale_factor( unsigned ndim, value det, bool log ) {
    value s;
    
    s = std::pow(M_PI,0.5*ndim) / std::tgamma( 0.5*ndim + 1 );
    s = 1. / s;
    s /= det;


    if (log) {
        s = fastlog(s);
    }
    
    return s;
}

// constructor
BoxKernel::BoxKernel() : KernelBase<BoxKernel>(KernelType::Box) {}

// methods
value BoxKernel::scale_factor( unsigned int n, value * bw, bool log ) const {
    value det = std::accumulate( bw, bw+n, std::pow(BOX_KERNEL_FACTOR,n), std::multiplies<value>() );
    return box_scale_factor( n, det, log );
}
value BoxKernel::scale_factor( unsigned int n, value * bw, bool log, 
    std::vector<bool>::const_iterator selection ) const {
        
    unsigned int ndim = 0;
    value det = 1.;
    
    for (unsigned int d=0; d<n; ++d) {
        if (*selection) {
            det *= bw[d] * BOX_KERNEL_FACTOR;
            ++ndim;
        }
        ++selection;
    }
    
    return box_scale_factor( ndim, det, log );
}

value BoxKernel::probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / ((*bw++)*BOX_KERNEL_FACTOR);
        d += (tmp*tmp);
        
        if (d>1.) { return 0.; }
     }
     
     return 1.;
}
value BoxKernel::probability( value dsquared ) const {
    if (dsquared>=1.) { return 0.; }
    else { return 1.; }
}

value BoxKernel::log_probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / ((*bw++)*BOX_KERNEL_FACTOR);
        d += (tmp*tmp);
        
        if (d>1.) { return -std::numeric_limits<value>::infinity(); }
     }
     
     return 0.;
}
value BoxKernel::log_probability( value dsquared ) const {
    if (dsquared>=1.) { return -std::numeric_limits<value>::infinity();; }
    else { return 0.; }
}

value BoxKernel::partial_logp( unsigned int n, const value * loc, 
    const value * bw, const value * point, 
    std::vector<bool>::const_iterator selection) const {
    
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        if (*selection++) {
            tmp = (point[k] - loc[k])/(bw[k]*BOX_KERNEL_FACTOR);
            d += (tmp*tmp);
            if (d>=1.) { return -std::numeric_limits<value>::infinity(); }
        }
    }
    
    return 0.;
}


// yaml
YAML::Node BoxKernel::to_yaml_impl() const {
    return YAML::Node();
}
std::unique_ptr<BoxKernel> BoxKernel::from_yaml( const YAML::Node & node ) {
    return std::make_unique<BoxKernel>();
}

// hdf5
void BoxKernel::to_hdf5_impl(HighFive::Group & group) const {}

std::unique_ptr<BoxKernel> BoxKernel::from_hdf5(const HighFive::Group & group) {
    return std::make_unique<BoxKernel>();
}
