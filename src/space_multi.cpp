#include "space_multi.hpp"
#include "space.hpp"

MultiSpace::MultiSpace( std::vector<const Space*> spaces ) :
SpaceBase<MultiSpace>( "multi", make_spec(spaces), make_kernel(spaces) ) {
    for (auto & k : spaces) {
        if (k->klass()=="multi") {
            // merge
            auto other = dynamic_cast<const MultiSpace*>(k);
            if (other==nullptr) { throw std::runtime_error("Internal error: cannot cast to MultiSpace"); }
            for (unsigned int c=0; c<other->nchildren(); ++c) {
                spaces_.emplace_back( other->child(c).clone() );
            }
        } else {
            spaces_.emplace_back( k->clone() );
        }
    }
}

MultiSpace::MultiSpace( const MultiSpace & other ) :
SpaceBase<MultiSpace>( other ) {
    for (auto & k : other.spaces_) {
        spaces_.emplace_back( k->clone() );
    }
}

SpaceSpecification MultiSpace::make_spec( std::vector<const Space*> spaces ) const {
    
    SpaceSpecification spec;
    
    for (auto s : spaces) {
        spec.append( s->specification() );
    }
    
    return spec;
}

Component MultiSpace::make_kernel( std::vector<const Space*> spaces ) const {
    Component k;
    
    k.scale_factor = 1.;
    k.scale_factor_log = 0.;
    
    for (auto s : spaces) {
        k.bandwidth.insert( std::end(k.bandwidth), std::begin(s->default_kernel().bandwidth), std::end(s->default_kernel().bandwidth) );
        k.location.insert( std::end(k.location), std::begin(s->default_kernel().location), std::end(s->default_kernel().location) );
        k.scale_factor *= s->default_kernel().scale_factor;
        k.scale_factor_log += s->default_kernel().scale_factor_log;
    }
    
    return k;
}

MultiSpace * MultiSpace::fromYAML( const YAML::Node & node ) {
    
    std::vector<std::unique_ptr<Space>> spaces;
    
    if (!node["spaces"] || !node["spaces"].IsSequence()) {
        throw std::runtime_error("Ill-formed multiplicative space definition.");
    }
    
    unsigned int nspaces = node["spaces"].size();
    
    for (unsigned int k=0; k<nspaces; ++k) {
        spaces.emplace_back( space_from_YAML( node["spaces"][k] ) );
    }
    
    std::vector<const Space*> pspaces;
    for (auto & k : spaces) {
        pspaces.push_back( k.get() );
    }
    
    MultiSpace * p = new MultiSpace( pspaces );
    
    return p;
}

YAML::Node MultiSpace::asYAML() const {
    YAML::Node node;
    for (auto & k : spaces_) {
        node["spaces"].push_back( k->toYAML() );
    }
    return node;
}



value MultiSpace::compute_scale_factor( value * bw, bool log ) const {
   value s;
    if (log) {
        s = 0.;
        for (auto & k : spaces_ ) {
            s += k->compute_scale_factor(bw, log);
            bw += k->nbw();
        }
    } else {
        s = 1.;
        for (auto & k : spaces_ ) {
            s *= k->compute_scale_factor(bw, log);
            bw += k->nbw();
        }
    }
    return s;
}

value MultiSpace::compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log ) const {
   value s;
    if (log) {
        s = 0.;
        for (auto & k : spaces_ ) {
            s += k->compute_scale_factor(selection, bw, log);
            selection += k->ndim();
            bw += k->nbw();
        }
    } else {
        s = 1.;
        for (auto & k : spaces_ ) {
            s *= k->compute_scale_factor(bw, log);
            selection += k->ndim();
            bw += k->nbw();
        }
    }
    return s;
}

value MultiSpace::mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const {
    
    value d = 0.;
    for (unsigned int k=0; k<spaces_.size(); ++k) {
        
        d+=spaces_[k]->mahalanobis_distance_squared( refloc, refbw, targetloc, threshold );
        
        if (d>=threshold) {break;}
        
        refloc+=spaces_[k]->ndim();
        refbw+=spaces_[k]->nbw();
        targetloc+=spaces_[k]->ndim();
        
    }
    
    return d;
}

void MultiSpace::merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const {

    for (unsigned int k=0; k<spaces_.size(); ++k) {
        
        spaces_[k]->merge( w1, loc1, bw1, w2, loc2, bw2 );
        
        loc1 += spaces_[k]->ndim();
        bw1 += spaces_[k]->nbw();
        loc2 += spaces_[k]->ndim();
        bw2 += spaces_[k]->nbw();
        
    }
        
}

value MultiSpace::probability( const value * loc, const value * bw, const value * point ) const {
    
    value p = 1.;
    
    for (auto & k : spaces_) {
        p *= k->probability(loc, bw, point);
        if (p==0.) {return p;}
        loc += k->ndim();
        bw += k->nbw();
        point += k->ndim();
    }
    return p;
    
}

void MultiSpace::probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    throw std::runtime_error("Not implemented: MultiSpace::probability");
}

value MultiSpace::log_probability( const value * loc, const value * bw, const value * point ) const {
    
    value p = 0.;
    
    for (auto & k : spaces_) {
        p += k->log_probability(loc, bw, point);
        if (std::isinf(p)) {return p;}
        loc += k->ndim();
        bw += k->nbw();
        point += k->ndim();
    }
    return p;
    
}

void MultiSpace::log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    throw std::runtime_error("Not implemented: MultiSpace::probability");
}

value MultiSpace::partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const {
    value p = 0;
    for (auto & k : spaces_) {
        p += k->partial_logp( loc, bw, point, selection ); // do not propagate scale
        if (std::isinf(p)) {return p;}
        loc += k->ndim();
        bw += k->nbw();
        point += k->ndim();
        selection += k->ndim();
    }
    return p;
}

   
Grid * MultiSpace::grid( const std::vector<Grid*> & grids, const std::vector<bool> & valid ) const {
    
    std::unique_ptr<MultiGrid> g = std::unique_ptr<MultiGrid>( new MultiGrid( grids, valid ) );
    
    if (!specification().issubspace( g->specification() )) {
        throw std::runtime_error("Grid space is not proper subspace.");
    }
    
    return g.release();
    
}
