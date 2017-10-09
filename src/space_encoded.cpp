#include "space_encoded.hpp"

// constructors
EncodedSpace::EncodedSpace( std::string name, std::vector<value> & lut,
    value bandwidth, unsigned int index )
    : EncodedSpace::EncodedSpace(name, lut, GaussianKernel(), bandwidth, index) {}

EncodedSpace::EncodedSpace( std::string name, std::vector<value> & lut, 
    const Kernel & k, value bandwidth, unsigned int index)
    : SpaceBase<EncodedSpace>( "encoded", make_spec( name, lut, k),
        make_kernel( bandwidth, index, k ) ), 
      lut_( new std::vector<value>(lut) ), kernel_(k.clone()) {
    nlut_ = static_cast<unsigned int>( std::sqrt( lut_->size() ) );
    if (nlut_*nlut_ != lut_->size()) {
        throw std::runtime_error("Squared distance look-up table needs to be a square matrix.");
    }
    if (index>=nlut_) {
        throw std::runtime_error("Index is out of range.");
    }
}

SpaceSpecification EncodedSpace::make_spec( std::string name, std::vector<value> & lut, 
    const Kernel & ktype ) {
    
    std::string extra = "kernel=" + ktype.to_string();
    extra += ", N=" + std::to_string( (unsigned int) std::sqrt( lut.size() ) );
    
    return SpaceSpecification( DimSpecification(name, "encoded", extra) );
}

Component EncodedSpace::make_kernel(value bw, unsigned int loc, 
    const Kernel & ktype) const {
    
    Component k;
    
    k.location = { static_cast<value>(loc) };
    k.bandwidth = { bw };
    
    k.scale_factor = ktype.scale_factor( 1, &bw, false );
    k.scale_factor_log = fastlog( k.scale_factor );
    
    return k;
}

// grid construction
Grid * EncodedSpace::grid(unsigned int delta) const {
    
    unsigned int n = (nlut_-1)/delta + 1;
    std::vector<value> v( n );
    
    for (unsigned int k=0; k<n; ++k) {
        v[k] = delta * k;
    }
    
    return new VectorGrid( { v }, specification(), {} );
}

Grid * EncodedSpace::grid(const std::vector<value> & v, 
    const std::vector<bool> & valid) const {
    
    // check if all int(values) are in range [0, nlut_-1]
    auto result = std::find_if(v.begin(), v.end(),
        [this](const value & n) { return static_cast<unsigned int>(n)>=this->nlut_; });
    
    if (result!=v.end()) {
        throw std::runtime_error("Found grid values out of range.");
    }
    
    return new VectorGrid( {v}, specification(), valid );
        
}

// methods
value EncodedSpace::compute_scale_factor( value * bw, bool log) const {
    
    return kernel_->scale_factor( 1, bw, log );
}

value EncodedSpace::compute_scale_factor( std::vector<bool>::const_iterator selection, 
    value * bw, bool log) const {
    
    return kernel_->scale_factor( 1, bw, log, selection );
}


value EncodedSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    
    unsigned int index1 = static_cast<unsigned int>(*refloc);
    unsigned int index2 = static_cast<unsigned int>(*targetloc);
    
    if (index1>=nlut_ || index2>=nlut_) { return threshold; }
    
    return (*lut_)[index1 + nlut_*index2] / (*refbw * *refbw);
}

void EncodedSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {
    
    // find lut entry k that minimizes w1*lut[index1,k] + w2*lut[index2,k]

    unsigned int k=0;
    value min_distance = std::numeric_limits<value>::infinity();
    value tmp;
    
    for (unsigned int n=0; n<nlut_; ++n) {
        tmp = w1 * (*lut_)[n+static_cast<unsigned int>(*loc1)*nlut_] + w2 * (*lut_)[n+static_cast<unsigned int>(*loc2)*nlut_];
        if (tmp<min_distance) {
            min_distance=tmp;
            k = n;
        }
    }
    
    *loc1 = static_cast<value>(k);
    
    value w = w1+w2;
    
    *bw1 = std::sqrt(w1*(*bw1 * *bw1)/w + w2*(*bw2 * *bw2)/w + 
                    w1*w2*(*lut_)[static_cast<unsigned int>(*loc1)+
                                  static_cast<unsigned int>(*loc2)*nlut_]/(w*w));
}

value EncodedSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    
    unsigned int idx = static_cast<unsigned int>( *point );
    
    if (idx>=nlut_) { return 0.; }
    
    value d = (*lut_)[idx + nlut_*static_cast<unsigned int>(*loc)] / (*bw * *bw);
    
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
    
    unsigned int idx = static_cast<unsigned int>( *point );
    
    if (idx>=nlut_) { return -std::numeric_limits<value>::infinity(); }
    
    value d = (*lut_)[idx + nlut_*static_cast<unsigned int>(*loc)] / (*bw * *bw);
    
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
    
    if (*selection) {
        unsigned int idx = static_cast<unsigned int>( *point );
        
        if (idx>=nlut_) { p=-std::numeric_limits<value>::infinity(); }
        else {
            p = (*lut_)[idx + nlut_*static_cast<unsigned int>(*loc)] / (*bw * *bw);
            
            //if (p>=cutoff_squared_) { p=-std::numeric_limits<value>::infinity(); } else { p = -0.5*p; }
            
            p = fastlog( kernel_->probability( p ) );
        }
    }
        
    return p;
}


// yaml
std::unique_ptr<EncodedSpace> EncodedSpace::from_yaml( const YAML::Node & node ) {
    
    if (!node["name"] || !node["lut"]) {
        throw std::runtime_error("Ill-formed encoded space definition.");
    }
    
    std::string name = node["names"].as<std::string>();
    
    std::unique_ptr<Kernel> k;
    
    if (node["kernel"]) {
        k = kernel_from_yaml( node["kernel"] );
    } else {
        k.reset( new GaussianKernel() );
    }
    
    std::vector<value> lut = node["lut"].as<std::vector<value>>();
    
    return std::make_unique<EncodedSpace>(name, lut, *k);
}

YAML::Node EncodedSpace::to_yaml_impl() const {
    YAML::Node node;
    node["name"] = specification().dim(0).name();;
    node["kernel"] = kernel_->to_yaml();
    node["lut"] = *lut_;
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
}

std::unique_ptr<EncodedSpace> EncodedSpace::from_hdf5(const HighFive::Group & group) {
    std::string name;
    
    HighFive::DataSet ds_name = group.getDataSet("name");
    ds_name.read(name);
    
    std::vector<value> lut;
    HighFive::DataSet ds_lut = group.getDataSet("lut");
    ds_lut.read(lut);
    
    std::unique_ptr<Kernel> k = kernel_from_hdf5(group.getGroup("kernel"));
    
    return std::make_unique<EncodedSpace>(name, lut, *k);
}
