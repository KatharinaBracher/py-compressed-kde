#include "space_base.hpp"

#include <iostream>
#include <limits>
#include <algorithm>

// constructor
Space::Space( std::string klass, const SpaceSpecification & spec, const Component k ) :
klass_(klass), spec_(spec), default_kernel_(k) {}

// copy constructor
Space::Space( const Space & other ) :
klass_(other.klass_), spec_(other.spec_), default_kernel_(other.default_kernel_) {}

// clone
Space* Space::clone() const { throw std::runtime_error("not implemented"); }

// properties
std::vector<bool> Space::selection( const Space& space ) const {
    
    return specification().selection( space.specification() );
    
}

bool Space::issubspace( const Space& space ) const {
    
    return specification().issubspace( space.specification() );
    
}

std::string Space::klass() const { return klass_; }
unsigned int Space::ndim() const { return spec_.ndim(); }
unsigned int Space::nbw() const { return default_kernel_.bandwidth.size(); }
//const std::vector<std::string> Space::dimension_names() const {
    //std::vector<std::string> names;
    
    //for (auto & k : spec_) {
        //names.push_back( k.name() );
    //}
    
    //return names;
//}

//const std::vector<std::string> Space::dimension_types() const {
    //std::vector<std::string> types;
    
    //for (auto & k : spec_) {
        //types.push_back( k.type() );
    //}
    
    //return types;
//}

//const std::vector<std::string> Space::dimension_details() const {
    //std::vector<std::string> details;
    
    //for (auto & k : spec_) {
        //details.push_back( k.detail() );
    //}
    
    //return details;
//}

const SpaceSpecification & Space::specification() const {
    return spec_;
}
    
const Component & Space::default_kernel() const { return default_kernel_; }

std::unique_ptr<Component> Space::kernel() const {
    return std::unique_ptr<Component>( new Component( default_kernel_ ) );
}

// methods
void Space::update_scale_factor (Component & k) const {
    k.scale_factor = compute_scale_factor( k );
    //k.scale_factor_log = ( k.scale_factor );
}

value Space::compute_scale_factor( Component & k, bool log ) const {
    return compute_scale_factor( k.bandwidth.data() , log );
}
value Space::compute_scale_factor( Component & k, const std::vector<bool> & selection, bool log  ) const {
    return compute_scale_factor( selection.cbegin(), k.bandwidth.data(), log );
}

value Space::mahalanobis_distance_squared( const Component & reference, const Component & target, value threshold) const {
    return mahalanobis_distance_squared( reference.location.data(), reference.bandwidth.data(), target.location.data(), threshold );
}
value Space::mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const {
    return threshold;
}

void Space::merge( value w1, Component & first, value w2, const Component & second ) const {
    merge( w1, first.location.data(), first.bandwidth.data(), w2, second.location.data(), second.bandwidth.data() );
    update_scale_factor( first );
}
void Space::merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const {}

value Space::probability( const Component & k, const value * point ) const {
    return k.scale_factor * probability( k.location.data(), k.bandwidth.data(), point );
}
value Space::probability( const value * loc, const value * bw, const value * point ) const {
    return 0.;
}
void Space::probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    return;
}

value Space::log_probability( const Component & k, const value * point ) const {
    return k.scale_factor_log + log_probability( k.location.data(), k.bandwidth.data(), point );
}
value Space::log_probability( const value * loc, const value * bw, const value * point ) const {
    return -std::numeric_limits<value>::infinity();
}
void Space::log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    return;
}

value Space::partial_logp( const Component & k, const value * point, const std::vector<bool> & selection ) const {
    return partial_logp( k.location.data(), k.bandwidth.data(), point, selection.cbegin() );
}
value Space::partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const {
    return -std::numeric_limits<value>::infinity();
}


// yaml
YAML::Node Space::to_yaml() const {
    YAML::Node node;
    node["class"] = klass();
    node["space"] = this->to_yaml_impl();
    node["kernel"] = this->default_kernel_.to_yaml();
    return node;
}
YAML::Node Space::to_yaml_impl() const {
    throw std::runtime_error("Not implemented.");
}

void Space::save_to_yaml( std::ostream & stream, bool flow ) const {
    auto node = to_yaml(); 
    YAML::Emitter out;
    if (flow) { out << YAML::Flow; }
    out << node; 
    stream << out.c_str();
}
void Space::save_to_yaml( std::string path, bool flow ) const {
    std::ofstream fout(path); save_to_yaml(fout, flow);
}


// hdf5
void Space::to_hdf5(HighFive::Group & group) const {
    
    HighFive::Attribute attr = group.createAttribute<std::string>(
            "class", HighFive::DataSpace::From(klass()));
    attr.write(klass());
    
    HighFive::Group space = group.createGroup("space");
    this->to_hdf5_impl(space);
    
    HighFive::Group kernel = group.createGroup("kernel");
    this->default_kernel_.to_hdf5(kernel);
}
void Space::to_hdf5_impl(HighFive::Group & group) const {
    throw std::runtime_error("Not implemented.");
}
    



