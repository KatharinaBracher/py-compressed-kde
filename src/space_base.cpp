#include "space_base.hpp"

#include <iostream>
#include <limits>
#include <algorithm>

Space::Space( std::string klass, const SpaceSpecification & spec, const Component k ) :
klass_(klass), spec_(spec), default_kernel_(k) {
    
    //if (!isunique( dimension_names() )) {
        //throw std::runtime_error("Dimension names are not unique.");
    //}
    
    //hash_ = 0;
    //for (auto & k : spec_) {
        //hash_ = hash_combine( hash_, k.hash() );
    //}
    
}

Space::Space( const Space & other ) :
klass_(other.klass_), spec_(other.spec_), default_kernel_(other.default_kernel_) {} //, hash_(other.hash_) {}

Space* Space::clone() const { throw std::runtime_error("not implemented"); }

YAML::Node Space::toYAML() const {
    YAML::Node node;
    node["class"] = klass();
    node["space"] = this->asYAML();
    node["kernel"] = this->default_kernel_.toYAML();
    return node;
}
YAML::Node Space::asYAML() const {
    throw std::runtime_error("Not implemented.");
}
void Space::save( std::ostream & stream ) const {
    auto node = toYAML(); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}
void Space::save( std::string path ) const {
    std::ofstream fout(path); save(fout);
}

const SpaceSpecification & Space::specification() const {
    return spec_;
}

//size_t Space::hash() const { return hash_; }

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
    
const Component & Space::default_kernel() const { return default_kernel_; }

std::unique_ptr<Component> Space::kernel() const {
    return std::unique_ptr<Component>( new Component( default_kernel_ ) );
}

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


    



