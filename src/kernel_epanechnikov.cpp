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
#include "kernel_epanechnikov.hpp"

#include <cmath>
#include <limits>

value epanechnikov_scale_factor( unsigned ndim, value det, bool log ) {
    value s;
    
    s = std::pow(M_PI,0.5*ndim) / std::tgamma( 0.5*ndim + 1 );
    s = (0.5*ndim + 1) / s;
    s /= det;


    if (log) {
        s = fastlog(s);
    }
    
    return s;
    
}

// constructor
EpanechnikovKernel::EpanechnikovKernel() : KernelBase<EpanechnikovKernel>(KernelType::Epanechnikov) {}

// methods
value EpanechnikovKernel::scale_factor(unsigned int n, value * bw, bool log) const {
    value det = std::accumulate( bw, bw+n, std::pow(EPA_KERNEL_FACTOR,n), std::multiplies<value>() );
    return epanechnikov_scale_factor( n, det, log );
}
value EpanechnikovKernel::scale_factor( unsigned int n, value * bw, bool log, 
    std::vector<bool>::const_iterator selection ) const {
    
    unsigned int ndim = 0;
    value det = 1.;
    
    for (unsigned int d=0; d<n; ++d) {
        if (*selection) {
            det *= bw[d] * EPA_KERNEL_FACTOR;
            ++ndim;
        }
        ++selection;
    }
    
    return epanechnikov_scale_factor( ndim, det, log );
}

value EpanechnikovKernel::probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / ( (*bw++) * EPA_KERNEL_FACTOR );
        d += (tmp*tmp);
        
        if (d>=1.) { return 0.; }
     }
     
     return (1-d);
}
value EpanechnikovKernel::probability( value dsquared ) const {
    if (dsquared>=1.) { return 0.; }
    else { return (1-dsquared); }
}

value EpanechnikovKernel::log_probability( unsigned int n, const value * loc, 
    const value * bw, const value * point ) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        tmp = (*point++ - *loc++) / ( (*bw++) * EPA_KERNEL_FACTOR );
        d += (tmp*tmp);
        
        if (d>=1.) { return -std::numeric_limits<value>::infinity(); }
     }
     
     return fastlog(1-d);
}
value EpanechnikovKernel::log_probability( value dsquared ) const {
    if (dsquared>=1.) { return -std::numeric_limits<value>::infinity(); }
    else { return fastlog(1-dsquared); }
}

value EpanechnikovKernel::partial_logp( unsigned int n, const value * loc, 
    const value * bw, const value * point, std::vector<bool>::const_iterator selection) const {
    value tmp, d=0.;
    
    for (unsigned int k=0; k<n; ++k) {
        if (*selection++) {
            tmp = (point[k] - loc[k])/ (bw[k]*EPA_KERNEL_FACTOR);
            d += (tmp*tmp);
            if (d>=1.) { return -std::numeric_limits<value>::infinity(); }
        }
    }
    
    return fastlog(1-d);
}

// yaml
YAML::Node EpanechnikovKernel::to_yaml_impl() const {
    return YAML::Node();
}
std::unique_ptr<EpanechnikovKernel> EpanechnikovKernel::from_yaml(
    const YAML::Node & node ) {
    return std::make_unique<EpanechnikovKernel>();
}


// flatbuffers
flatbuffers::Offset<fb_serialize::Kernel> EpanechnikovKernel::to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const {

    return fb_serialize::CreateKernel(
        builder,
        builder.CreateString(kerneltype_tostring(type())),
        fb_serialize::KernelType_EpanechnikovKernel,
        fb_serialize::CreateEpanechnikovKernel(builder).Union()
    );
}

std::unique_ptr<EpanechnikovKernel> EpanechnikovKernel::from_flatbuffers(const fb_serialize::EpanechnikovKernel * kernel) {
    return std::make_unique<EpanechnikovKernel>();
}

// hdf5
void EpanechnikovKernel::to_hdf5_impl(HighFive::Group & group) const {}

std::unique_ptr<EpanechnikovKernel> EpanechnikovKernel::from_hdf5(
    const HighFive::Group & group) {
    return std::make_unique<EpanechnikovKernel>();
}
