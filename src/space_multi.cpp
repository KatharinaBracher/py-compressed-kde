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
#include "space_multi.hpp"
#include "space.hpp"

// constructor
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

// copy constructor
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

// grid construction
Grid * MultiSpace::grid( const std::vector<Grid*> & grids, 
    const std::vector<bool> & valid ) const {
    
    std::unique_ptr<MultiGrid> g = std::unique_ptr<MultiGrid>(
        new MultiGrid( grids, valid ) );
    
    if (!specification().issubspace( g->specification() )) {
        throw std::runtime_error("Grid space is not proper subspace.");
    }
    
    return g.release();
    
}    

// methods
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

value MultiSpace::compute_scale_factor( std::vector<bool>::const_iterator selection, 
    value * bw, bool log ) const {
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

value MultiSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    
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

void MultiSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {

    for (unsigned int k=0; k<spaces_.size(); ++k) {
        
        spaces_[k]->merge( w1, loc1, bw1, w2, loc2, bw2 );
        
        loc1 += spaces_[k]->ndim();
        bw1 += spaces_[k]->nbw();
        loc2 += spaces_[k]->ndim();
        bw2 += spaces_[k]->nbw();
        
    }
        
}

value MultiSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    
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

void MultiSpace::probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    throw std::runtime_error("Not implemented: MultiSpace::probability");
}

value MultiSpace::log_probability( const value * loc, const value * bw, 
    const value * point ) const {
    
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

void MultiSpace::log_probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    throw std::runtime_error("Not implemented: MultiSpace::probability");
}

value MultiSpace::partial_logp( const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection ) const {
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


// yaml
std::unique_ptr<MultiSpace> MultiSpace::from_yaml( const YAML::Node & node ) {
    
    std::vector<std::unique_ptr<Space>> spaces;
    
    if (!node["spaces"] || !node["spaces"].IsSequence()) {
        throw std::runtime_error("Ill-formed multiplicative space definition.");
    }
    
    unsigned int nspaces = node["spaces"].size();
    
    for (unsigned int k=0; k<nspaces; ++k) {
        //spaces.emplace_back( space_from_YAML( node["spaces"][k] ) );
        spaces.push_back(std::move(space_from_yaml(node["spaces"][k])));
    }
    
    std::vector<const Space*> pspaces;
    for (auto & k : spaces) {
        pspaces.push_back( k.get() );
    }
    
    return std::make_unique<MultiSpace>(pspaces);
}

YAML::Node MultiSpace::to_yaml_impl() const {
    YAML::Node node;
    for (auto & k : spaces_) {
        node["spaces"].push_back( k->to_yaml() );
    }
    return node;
}


// flatbuffers
flatbuffers::Offset<fb_serialize::SpaceData> MultiSpace::to_flatbuffers_impl(flatbuffers::FlatBufferBuilder & builder) const {

    std::vector<flatbuffers::Offset<fb_serialize::Space>> subspace_vector;

    for (auto & k : spaces_) {
        auto v = k->to_flatbuffers(builder);
        subspace_vector.push_back(v);
    }

    auto subspaces = builder.CreateVector(subspace_vector);

    return fb_serialize::CreateSpaceData(
        builder,
        fb_serialize::SpaceType_NONE,
        0,
        subspaces
    );
}

std::unique_ptr<MultiSpace> MultiSpace::from_flatbuffers(const fb_serialize::Space * space) {

    auto saved_klass = space->klass()->str();

    if (saved_klass!="multi") {
        throw std::runtime_error("Expected multi, but got " + saved_klass);
    }

    auto default_kernel = components_from_flatbuffers(space->default_kernel());

    auto subspaces = space->data()->subspaces();

    unsigned int nspace = subspaces->size();

    std::vector<std::unique_ptr<Space>> s;

    for (unsigned int k=0; k<nspace; ++k) {
        s.push_back(std::move(space_from_flatbuffers(subspaces->Get(k))));
    }

    std::vector<const Space*> s_ptr;
    for (auto & k : s) {
        s_ptr.push_back( k.get() );
    }

    auto ptr = std::make_unique<MultiSpace>(s_ptr);

    ptr->set_default_kernel(*(default_kernel[0]));

    return ptr;
}


// hdf5
void MultiSpace::to_hdf5_impl(HighFive::Group & group) const {
    unsigned int s = 0;
    
    HighFive::Attribute attr = group.createAttribute<unsigned int>(
            "nspace", HighFive::DataSpace::From(s));
    attr.write(spaces_.size());
    
    for (auto & k : spaces_) {
        HighFive::Group subgroup = group.createGroup("space" + std::to_string(s));
        
        k->to_hdf5(subgroup);
        
        ++s;
    }
}

std::unique_ptr<MultiSpace> MultiSpace::from_hdf5(const HighFive::Group & group) {
    
    unsigned int nspace;
    
    HighFive::Attribute attr = group.getAttribute("nspace");
    attr.read(nspace);
    
    std::vector<std::unique_ptr<Space>> s;
    for (unsigned int k=0; k<nspace; ++k) {
        HighFive::Group subgroup = group.getGroup("space" + std::to_string(k));
        s.push_back(std::move(space_from_hdf5(subgroup)));
    }
    
    std::vector<const Space*> ptr;
    for (auto & k : s) {
        ptr.push_back( k.get() );
    }
    
    return std::make_unique<MultiSpace>( ptr );
    
}
