#include "space_categorical.hpp"

// constructor
CategoricalSpace::CategoricalSpace( std::string name, std::vector<std::string> labels, unsigned int category ) :
SpaceBase<CategoricalSpace>( "categorical", make_spec(name, labels), make_kernel(category) ), labels_(labels) {
    if (!isunique(labels)) {
        throw std::runtime_error("Labels are not unique.");
    }
}

SpaceSpecification CategoricalSpace::make_spec( std::string name, std::vector<std::string> labels ) {
    
    std::string extra = "labels=[";
    for (auto & k : labels) {
        extra += k + ",";
    }
    extra += "]";
    
    return SpaceSpecification( DimSpecification( name, "categorical", extra ) );
    
}

Component CategoricalSpace::make_kernel(unsigned int category) const {
    Component k;
    k.location = { static_cast<value>(category) };
    k.scale_factor = 1;
    k.scale_factor_log = 0.;
    return k;
}

Grid * CategoricalSpace::grid() const {
    
    std::vector<value> v(labels_.size());
    std::iota( v.begin(), v.end(), 0 );
    
    return new VectorGrid( { v }, specification(), {} );
    
}
 
// methods
value CategoricalSpace::compute_scale_factor( value * bw, bool log ) const {
    if (log) { return 0.; } else { return 1.; };
}

value CategoricalSpace::compute_scale_factor( 
    std::vector<bool>::const_iterator selection, value * bw, bool log ) const {
    if (log) { return 0.; } else { return 1.; };
}

value CategoricalSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    if (static_cast<unsigned int>(*refloc) == static_cast<unsigned int>(*targetloc)) { 
        return 0.;
    } else {
        return std::numeric_limits<value>::infinity();
    }
}

void CategoricalSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {
    // nothing to do
}

value CategoricalSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    return static_cast<value>( static_cast<unsigned int>(*loc)==static_cast<unsigned int>(*point) );
}

void CategoricalSpace::probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = probability( loc, bw, points++ );
    }
}

value CategoricalSpace::log_probability( const value * loc, const value * bw, 
    const value * point ) const {
    value p = 0;
    if (static_cast<unsigned int>(*loc)!=static_cast<unsigned int>(*point)) {
        p = -std::numeric_limits<value>::infinity();
    }
    return p;
}

void CategoricalSpace::log_probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = log_probability( loc, bw, points++ );
    }
}

value CategoricalSpace::partial_logp( const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection ) const {
    value p=0.;
    
    if (*selection && !(static_cast<unsigned int>(*loc)==static_cast<unsigned int>(*point)) ){
        p = -std::numeric_limits<value>::infinity();
    }
    
    return p;
}


// yaml
std::unique_ptr<CategoricalSpace> CategoricalSpace::from_yaml(
    const YAML::Node & node ) {
    
    if (!node["name"] || !node["labels"]) {
        throw std::runtime_error("Ill-formed categorical space definition.");
    }
    
    return std::make_unique<CategoricalSpace>(node["name"].as<std::string>(), 
        node["labels"].as<std::vector<std::string>>());
    
} 

YAML::Node CategoricalSpace::to_yaml_impl() const {
    YAML::Node node;
    node["name"] = specification().dim(0).name();
    node["labels"] = labels_;
    return node;
}


// hdf5
void CategoricalSpace::to_hdf5_impl(HighFive::Group & group) const {
    
    std::string name = specification().dim(0).name();
    
    HighFive::DataSet ds_name = group.createDataSet<std::string>(
        "name", HighFive::DataSpace::From(name));
    ds_name.write(name);
    
    HighFive::DataSet ds_labels = group.createDataSet<std::string>(
        "labels", HighFive::DataSpace::From(labels_));
    ds_labels.write(labels_);
}

std::unique_ptr<CategoricalSpace> CategoricalSpace::from_hdf5(
    const HighFive::Group & group) {
    
    std::string name;
    std::vector<std::string> labels;
    
    HighFive::DataSet ds_name = group.getDataSet("name");
    ds_name.read(name);
    
    HighFive::DataSet ds_labels = group.getDataSet("labels");
    ds_labels.read(labels);
    
    return std::make_unique<CategoricalSpace>(name, labels);
}
