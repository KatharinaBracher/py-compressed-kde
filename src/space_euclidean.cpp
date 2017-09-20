#include "space_euclidean.hpp"

// constructors
EuclideanSpace::EuclideanSpace( std::vector<std::string> names,
    std::vector<value> bandwidth, std::vector<value> location)
    : EuclideanSpace::EuclideanSpace( names, GaussianKernel(), bandwidth, location ) {}

EuclideanSpace::EuclideanSpace( std::vector<std::string> names, const Kernel & k, 
    std::vector<value> bandwidth, std::vector<value> location)
    : SpaceBase<EuclideanSpace>( "euclidean", make_spec( names, k), 
        make_kernel(names.size(), bandwidth, location, k) ), 
    kernel_( k.clone() ) {}

SpaceSpecification EuclideanSpace::make_spec( std::vector<std::string> names, 
    const Kernel & ktype ) {
    
    std::vector<DimSpecification> dspec;
    
    std::string extra = "kernel=" + ktype.to_string();
    
    for (unsigned int k=0; k<names.size(); ++k) {
        dspec.push_back( DimSpecification( names[k], "euclidean", extra ) );
    }
    
    return SpaceSpecification( dspec );
}

Component EuclideanSpace::make_kernel(unsigned int n, std::vector<value> bw, 
    std::vector<value> loc, const Kernel & ktype ) const {
    Component k;
    
    if (bw.size()==0) {
        bw.assign( n, DEFAULT_EUCLIDEAN_BANDWIDTH );
    } else if (bw.size()!=n) {
        throw std::runtime_error("Incorrect bandwidth vector size.");
    }
    
    if (loc.size()==0) {
        loc.assign( n, DEFAULT_EUCLIDEAN_LOCATION );
    } else if (loc.size()!=n) {
        throw std::runtime_error("Incorrect location vector size.");
    }
    
    k.location = loc;
    k.bandwidth = bw;
    
    k.scale_factor = ktype.scale_factor( n, bw.data(), false );
    k.scale_factor_log = fastlog( k.scale_factor );
    
    return k;
}

// grid construction
VectorGrid * EuclideanSpace::grid( const std::vector<std::vector<value>> & vectors, 
    const std::vector<bool> & valid, const std::vector<bool> & selection ) const {

    unsigned int nvectors = vectors.size();
    unsigned int nselection = selection.size();
    unsigned int nselected = std::count(selection.cbegin(), selection.cend(),true);
    
    if (nvectors==ndim()) {
        if (nselection!=0 && nselection!=ndim()) {
            throw std::runtime_error("Incorrect size of selection vector.");
        }
    } else if (nvectors==0 || nvectors>ndim()) {
        throw std::runtime_error("Too many or too few grid vectors specified.");
    } else {
        if (nselection!=ndim()) {
            throw std::runtime_error("Incorrect size of selection vector .");
        }
        if (nselected!=nvectors) {
            throw std::runtime_error("Mismatch between number of grid vectors "
                "and number of selected dimensions");
        }
    }
    
    SpaceSpecification spec;
    
    if (nselection!=0) {
        if (nselected==0) {
            throw std::runtime_error("Select at least one dimension.");
        }
        
        spec = specification().select( selection );
        
        if (nvectors>nselected) {
            
            std::vector<std::vector<value>> v;
            for (unsigned int k=0; k<nvectors; ++k) {
                if (selection[k]) {
                    v.push_back( vectors[k] );
                }
            }
            
            return new VectorGrid( v, spec, valid );
            
        } else {
            return new VectorGrid( vectors, spec, valid );
        }
        
    }
    
    return new VectorGrid( vectors, specification(), valid );
}


// methods
value EuclideanSpace::compute_scale_factor( value * bw, bool log) const {
    
    return kernel_->scale_factor( ndim(), bw, log );
}

value EuclideanSpace::compute_scale_factor(
    std::vector<bool>::const_iterator selection, value * bw, bool log) const {
    
    return kernel_->scale_factor( ndim(), bw, log, selection );
}


value EuclideanSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    
    // bandwidths are stored as Gaussian equivalent
    // so computation of mahalanobis distance is the same for all kernel types
    
    value tmp, d=0;
    
    for (unsigned int k=0; k<ndim(); ++k) {
        tmp = (targetloc[k] - refloc[k])/refbw[k];
        d += (tmp*tmp);
        if (d>=threshold) {break;}
    }
    
    return d;
}

void EuclideanSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {
    
    // bandwidths are stored as Gaussian equivalent
    // so merging is the same for all kernel types
    
    value tmp;
    value w = w1 + w2;
    
    for (unsigned int d=0; d<ndim(); ++d) {
        tmp = w1*(bw1[d]*bw1[d] + loc1[d]*loc1[d]) + w2*(bw2[d]*bw2[d] + loc2[d]*loc2[d]);
        loc1[d] *= w1;
        loc1[d] += loc2[d] * w2;
        loc1[d] /= w;
        bw1[d] = std::sqrt( tmp/w - loc1[d]*loc1[d] );
    }
}

value EuclideanSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    
    return kernel_->probability( ndim(), loc, bw, point );
}

void EuclideanSpace::probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = kernel_->probability(1, loc, bw, points++);
    }
}

value EuclideanSpace::log_probability( const value * loc, const value * bw, 
    const value * point ) const {
    
    return kernel_->log_probability( ndim(), loc, bw, point );
}

void EuclideanSpace::log_probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = kernel_->log_probability(1, loc, bw, points++);
    }
}

value EuclideanSpace::partial_logp( const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection ) const {
    
    return kernel_->partial_logp( ndim(), loc, bw, point, selection );
}


// yaml
std::unique_ptr<EuclideanSpace> EuclideanSpace::from_yaml(const YAML::Node & node) {
    
    if (!node["names"]) {
        throw std::runtime_error("Ill-formed euclidean space definition.");
    }
    
    std::vector<std::string> dims = node["names"].as<std::vector<std::string>>();
    
    std::unique_ptr<Kernel> k;
    
    if (node["kernel"]) {
        k = kernel_from_yaml( node["kernel"] );
    } else {
        k.reset( new GaussianKernel() );
    }
    
    return std::make_unique<EuclideanSpace>(dims, *k);
}

YAML::Node EuclideanSpace::to_yaml_impl() const {
    YAML::Node node;
    node["names"] = specification().names();
    node["kernel"] = kernel_->to_yaml();
    return node;
}


// hdf5
void EuclideanSpace::to_hdf5_impl(HighFive::Group & group) const {
    
    HighFive::DataSet ds_name = group.createDataSet<std::string>(
        "names", HighFive::DataSpace::From(specification().names()));
    ds_name.write(specification().names());
    
    HighFive::Group kernel = group.createGroup("kernel");
    kernel_->to_hdf5(kernel);
}

std::unique_ptr<EuclideanSpace> EuclideanSpace::from_hdf5(
    const HighFive::Group & group) {
    
    std::vector<std::string> names;
    
    HighFive::DataSet ds_name = group.getDataSet("names");
    ds_name.read(names);
    
    std::unique_ptr<Kernel> k = kernel_from_hdf5(group.getGroup("kernel"));
    
    return std::make_unique<EuclideanSpace>(names, *k);
}
