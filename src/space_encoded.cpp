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
#include "space_encoded.hpp"

// constructors
EncodedSpace::EncodedSpace( std::string name, const std::vector<value> & lut,
    value bandwidth, unsigned int index )
    : EncodedSpace::EncodedSpace(name, {}, lut, GaussianKernel(), bandwidth, index) {}
    
EncodedSpace::EncodedSpace( std::string name, const std::vector<value> & points,
    const std::vector<value> & lut, value bandwidth, unsigned int index )
    : EncodedSpace::EncodedSpace(name, points, lut, GaussianKernel(), bandwidth, index) {}

EncodedSpace::EncodedSpace( std::string name, const std::vector<value> & points,
    const std::vector<value> & lut, const Kernel & k, value bandwidth, unsigned int index)
    : SpaceBase<EncodedSpace>( "encoded", make_spec( name, lut, k),
        make_kernel( bandwidth, points, index, k ) ),
      lut_( new std::vector<value>(lut) ), kernel_(k.clone()) {
    
    nlut_ = static_cast<unsigned int>( std::sqrt( lut_->size() ) );
    if (nlut_*nlut_ != lut_->size()) {
        throw std::runtime_error("Squared distance look-up table needs to be a square matrix.");
    }
    
    //std::cout << "nlut, size: " << std::to_string(nlut_) << ", " << std::to_string(lut_->size()) << std::endl;
    
    if (points.empty()) {
        
        use_index_ = true;
        points_.reset(new std::vector<value>());
        
        if (index>=nlut_) {
            throw std::runtime_error("Index is out of range.");
        }
        
    } else {
        use_index_ = false;
        points_.reset( new std::vector<value>(points) );
        
        if (nlut_!=points_->size()) {
            throw std::runtime_error("Sizes of point vector and look-up table do not match.");
        }
        
        // points have to be sorted
        if (!std::is_sorted(points_->begin(), points_->end())) {
            throw std::runtime_error("Points vector needs to be sorted.");
        }
    }
}

SpaceSpecification EncodedSpace::make_spec( std::string name, const std::vector<value> & lut,
    const Kernel & ktype ) {
    
    std::string extra = "kernel=" + ktype.to_string();
    extra += ", N=" + std::to_string( (unsigned int) std::sqrt( lut.size() ) );
    
    return SpaceSpecification( DimSpecification(name, "encoded", extra) );
}

Component EncodedSpace::make_kernel(value bw, const std::vector<value> & points, unsigned int loc, 
    const Kernel & ktype) const {
    
    Component k;
    
    if (points.empty()) {
        // use index
        k.location = { static_cast<value>(loc) };
    } else {
        if (loc>=points.size()) {
            throw std::runtime_error("Index out of range.");
        }
        k.location = { points[loc] };
    }
    
    k.bandwidth = { bw };
    
    k.scale_factor = ktype.scale_factor( 1, &bw, false );
    k.scale_factor_log = fastlog( k.scale_factor );
    
    return k;
}

// grid construction
Grid * EncodedSpace::grid(unsigned int delta) const {
    
    unsigned int n = (nlut_-1)/delta + 1;
    std::vector<value> v( n );
    
    if (use_index_) {
        for (unsigned int k=0; k<n; ++k) {
            v[k] = delta * k;
        }
    } else {
        for (unsigned int k=0; k<n; ++k) {
            v[k] = (*points_)[delta * k];
        }
    }
    
    return new VectorGrid( { v }, specification(), {} );
}

value nearest(const std::vector<value> & v, value x) {
    
    return v[nearest_index(v,x)];
}

unsigned int nearest_index(const std::vector<value> & v, value x) {
    
    unsigned int idx;
    
    auto it = std::lower_bound( v.begin(), v.end(), x );
    
    if (it==v.end()) {
        if (v.size()>0) { idx = v.size()-1; }
        else {throw std::runtime_error("Cannot look for nearest value in empty vector");}
    } else if (it==v.begin()) {
        idx = 0;
    } else {
        if ( (x-*std::prev(it)) < (*it-x) ) {
            it = std::prev(it);
        }
        idx = std::distance(v.begin(), it);
    }
    
    return idx;
}

Grid * EncodedSpace::grid(const std::vector<value> & v, 
    const std::vector<bool> & valid) const {
    
    if (use_index_) {
        // check if all int(values) are in range [0, nlut_-1]
        auto result = std::find_if(v.begin(), v.end(),
            [this](const value & n) { return static_cast<unsigned int>(n)>=this->nlut_; });
        
        if (result!=v.end()) {
            throw std::runtime_error("Found grid values out of range.");
        }
        
        return new VectorGrid( {v}, specification(), valid );
        
    } else {
        
        // map to nearest value in points_
        std::vector<value> mapped_values(v.size());
        std::transform(v.begin(), v.end(), mapped_values.begin(), [this](value y){return nearest(*(this->points_),y);} );
        
        return new VectorGrid( {mapped_values}, specification(), valid );
    }
    
}

// methods
value EncodedSpace::compute_scale_factor( value * bw, bool log) const {
    
    return kernel_->scale_factor( 1, bw, log );
}

value EncodedSpace::compute_scale_factor( std::vector<bool>::const_iterator selection, 
    value * bw, bool log) const {
    
    return kernel_->scale_factor( 1, bw, log, selection );
}

unsigned int EncodedSpace::get_index(value x) const {
    unsigned int result;
    if (use_index_) {
        result = static_cast<unsigned int>(x);
    } else {
        result = nearest_index(*points_, x);
    }
    
    if (result>=nlut_) {
        std::cout << "get_index: out of range " << std::to_string(result) << " " << std::to_string(nlut_) << std::endl;
        std::cout << "x, points size, use_index: " << std::to_string(x) << " " << std::to_string(points_->size()) << " " << std::to_string(use_index_) << std::endl;
        throw std::runtime_error("get index: out of range.");
    }
    
    //std::cout << "x, index, points[index] " << std::to_string(x) << ", " << std::to_string(result) << ", " << std::to_string( (*points_)[result] ) << std::endl;
    
    return result;
}

value EncodedSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    
    unsigned int index1;
    unsigned int index2;
    
    //std::cout << "encoded space: mahalanobis distance " << std::to_string(*refloc) << " " << std::to_string(*targetloc) << std::endl;
    
    try {
        index1 = get_index(*refloc);
        index2 = get_index(*targetloc);
    } catch (...) {
        return threshold;
    }
    
    //if (index1>=nlut_ || index2>=nlut_) { return threshold; }
    //std::cout << "encoded space: mahalanobis distance returns " << std::to_string(index1) << " " << std::to_string(index2) << std::endl;
    return (*lut_)[index1 + nlut_*index2] / (*refbw * *refbw);
}

void EncodedSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {
    
    //std::cout << "encoded space: merge" << std::endl;
    
    // find lut entry k that minimizes w1*lut[index1,k] + w2*lut[index2,k]
    auto index1 = get_index(*loc1);
    auto index2 = get_index(*loc2);
    
    unsigned int k=0;
    value min_distance = std::numeric_limits<value>::infinity();
    value tmp;
    
    for (unsigned int n=0; n<nlut_; ++n) {
        tmp = w1 * (*lut_)[n+index1*nlut_] + w2 * (*lut_)[n+index2*nlut_];
        if (tmp<min_distance) {
            min_distance=tmp;
            k = n;
        }
    }
    
    if (use_index_) {
        *loc1 = static_cast<value>(k);
    } else {
        *loc1 = (*points_)[k];
    }
    
    value w = w1+w2;
    
    *bw1 = std::sqrt(w1*(*bw1 * *bw1)/w + w2*(*bw2 * *bw2)/w + 
                    w1*w2*(*lut_)[index1+index2*nlut_]/(w*w));
}

value EncodedSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    
    unsigned int idx;
    
    //std::cout << "encoded space: probability" << std::endl;
    
    try {
        idx = get_index( *point );
    } catch (...) {
        return 0.;
    }
    
    //if (idx>=nlut_) { return 0.; }
    
    value d = (*lut_)[idx + nlut_*get_index(*loc)] / (*bw * *bw);
    
    return kernel_->probability( d );
}

void EncodedSpace::probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = probability( loc, bw, points++ );
    }
}

value EncodedSpace::log_probability( const value * loc, const value * bw, 
    const value * point ) const {
    
    unsigned int idx;
    
    //std::cout << "encoded space: log probability" << std::endl;
    
    try {
        idx = get_index( *point );
    } catch (...) {
        return -std::numeric_limits<value>::infinity();
    }
    
    //if (idx>=nlut_) { return -std::numeric_limits<value>::infinity(); }
    
    value d = (*lut_)[idx + nlut_*get_index(*loc)] / (*bw * *bw);
    
    return kernel_->log_probability( d );
}

void EncodedSpace::log_probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = log_probability( loc, bw, points++ );
    }
}

value EncodedSpace::partial_logp( const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection ) const {
    
    value p=0.;
    unsigned int idx;
    
    //std::cout << "encoded space: partial logp " << std::to_string(*loc) << " " << std::to_string(*point) << std::endl;
    //std::cout << "selection: " << std::to_string(*selection) << std::endl;
    //std::cout << "bandwidth: " << std::to_string(*bw) << std::endl;
    
    if (*selection) {
        //unsigned int idx = static_cast<unsigned int>( *point );
        try {
            idx = get_index(*point);
        } catch (...) {
            p = -std::numeric_limits<value>::infinity();
            return p;
        }
        //if (idx>=nlut_) { p=-std::numeric_limits<value>::infinity(); }
        //else {
            p = (*lut_)[idx + nlut_*get_index(*loc)] / (*bw * *bw);
            
            //if (p>=cutoff_squared_) { p=-std::numeric_limits<value>::infinity(); } else { p = -0.5*p; }
            //std::cout << "fast log: " << std::to_string(p) << std::endl;
            p = fastlog( kernel_->probability( p ) );
        //}
    }
    
    return p;
}


// yaml
std::unique_ptr<EncodedSpace> EncodedSpace::from_yaml( const YAML::Node & node ) {
    
    if (!node["name"] || !node["lut"]) {
        throw std::runtime_error("Ill-formed encoded space definition.");
    }
    
    std::string name = node["name"].as<std::string>();
    
    std::unique_ptr<Kernel> k;
    
    if (node["kernel"]) {
        k = kernel_from_yaml( node["kernel"] );
    } else {
        k.reset( new GaussianKernel() );
    }
    
    std::vector<value> lut = node["lut"].as<std::vector<value>>();
    
    bool use_index = node["use_index"].as<bool>();
    
    if (use_index) {
        return std::make_unique<EncodedSpace>(name, std::vector<value>(), lut, *k);
    } else {
        std::vector<value> points = node["points"].as<std::vector<value>>();
        return std::make_unique<EncodedSpace>(name, points, lut, *k);
    }
}

YAML::Node EncodedSpace::to_yaml_impl() const {
    YAML::Node node;
    node["name"] = specification().dim(0).name();;
    node["kernel"] = kernel_->to_yaml();
    node["lut"] = *lut_;
    node["use_index"] = use_index_;
    if (!use_index_) {
        node["points"] = *points_;
    }
    return node;
}


// hdf5
void EncodedSpace::to_hdf5_impl(HighFive::Group & group) const {
    std::string name = specification().dim(0).name();
    
    HighFive::DataSet ds_name = group.createDataSet<std::string>(
        "name", HighFive::DataSpace::From(name));
    ds_name.write(name);
    
    HighFive::Group kernel = group.createGroup("kernel");
    kernel_->to_hdf5(kernel);
    
    HighFive::DataSet ds_lut = group.createDataSet<value>(
        "lut", HighFive::DataSpace::From(*lut_));
    ds_lut.write(*lut_);
    
    if (!use_index_) {
        HighFive::DataSet ds_points = group.createDataSet<value>(
            "points", HighFive::DataSpace::From(*points_));
        ds_points.write(*points_);
    }
    
}

std::unique_ptr<EncodedSpace> EncodedSpace::from_hdf5(const HighFive::Group & group) {
    std::string name;
    
    HighFive::DataSet ds_name = group.getDataSet("name");
    ds_name.read(name);
    
    std::vector<value> lut;
    HighFive::DataSet ds_lut = group.getDataSet("lut");
    ds_lut.read(lut);
    
    std::unique_ptr<Kernel> k = kernel_from_hdf5(group.getGroup("kernel"));
    
    if (group.exist("points")) {
        
        std::vector<value> points;
        HighFive::DataSet ds_points = group.getDataSet("points");
        ds_points.read(points);
    
        return std::make_unique<EncodedSpace>(name, points, lut, *k);
        
    } else {
        return std::make_unique<EncodedSpace>(name, std::vector<value>(), lut, *k);
    }
}
