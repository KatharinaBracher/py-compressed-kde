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

// constructor
GaussianKernel::GaussianKernel( value cutoff )
    : KernelBase<GaussianKernel>(KernelType::Gaussian), cutoff_(cutoff), 
        cutoff_squared_(cutoff*cutoff) {}

// properties
value GaussianKernel::cutoff() const { return cutoff_; }
void GaussianKernel::set_cutoff( value v ) { cutoff_=v; cutoff_squared_=v*v; }

// methods
value GaussianKernel::scale_factor( unsigned int n, value * bw, bool log ) const {
    value det = std::accumulate( bw, bw+n, 1., std::multiplies<value>() );
    return gaussian_scale_factor( n, det, cutoff_, log );
}
value GaussianKernel::scale_factor( unsigned int n, value * bw, bool log, 
    std::vector<bool>::const_iterator selection ) const {
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

value GaussianKernel::probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
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

value GaussianKernel::log_probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
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

value GaussianKernel::partial_logp( unsigned int n, const value * loc, 
    const value * bw, const value * point, std::vector<bool>::const_iterator selection) const {
    
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


// yaml
YAML::Node GaussianKernel::to_yaml_impl() const {
    YAML::Node node;
    node["cutoff"] = cutoff_;
    return node;
}

std::unique_ptr<GaussianKernel> GaussianKernel::from_yaml(
    const YAML::Node & node) {
    
    return std::make_unique<GaussianKernel>(
        node["cutoff"].as<value>(DEFAULT_GAUSSIAN_CUTOFF));
}


// hdf5
void GaussianKernel::to_hdf5_impl(HighFive::Group & group) const {
    HighFive::DataSet ds = group.createDataSet<value>(
        "cutoff", HighFive::DataSpace::From(cutoff_));
    ds.write(cutoff_);
}

std::unique_ptr<GaussianKernel> GaussianKernel::from_hdf5(
    const HighFive::Group & group) {
    
    value cutoff;
    HighFive::DataSet ds = group.getDataSet("cutoff");
    ds.read(cutoff);
    
    return std::make_unique<GaussianKernel>(cutoff);
}
