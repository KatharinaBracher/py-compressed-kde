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
#include "grid_multi.hpp"
#include "space.hpp"

SpaceSpecification space_from_grids( const std::vector<Grid*> & grids ) {
    SpaceSpecification spec;
    for (auto & k : grids) {
        spec.append( k->specification() );
    }
    return spec;
}

std::vector<long unsigned int> shape_from_grids( const std::vector<Grid*> & grids ) {
    std::vector<long unsigned int> shape;
    for (auto & k : grids) {
        shape.insert( std::end(shape), k->shape().cbegin(), k->shape().cend() );
    }

    return shape;
}

const std::vector<bool> valid_from_grids(const std::vector<Grid*> & grids,
    const std::vector<bool> & valid) {

    if (valid.size()>0) {
        return valid;
    }

    auto shape = shape_from_grids(grids);
    unsigned int D = grids.size();
    unsigned int N = 1;

    for (auto & n : shape) {
        N *= n;
    }

    std::vector<bool> result(N, true);

    std::vector<bool> has_valid_vector(D, false);
    std::vector<unsigned int> cursor(D,0);
    std::vector<unsigned int> end(D);

    for (unsigned int d=0;d<D;++d) {
        end[d] = grids[d]->size();
        has_valid_vector[d] = (grids[d]->ninvalid()>0);
    }

    for (unsigned int k=0; k<N; ++k) {
        for (unsigned int d=0; d<D; ++d) {
            if (has_valid_vector[d]) {
                result[k] = result[k] && grids[d]->valid()[cursor[d]];
            }
        }
               
        for (int d=D-1; d>=0; --d) {
            
            ++cursor[d];
            if (cursor[d]>=end[d]) {
                cursor[d] = 0;
            } else {
                break;
            }
        }
    }

    return result;
}

// constructor
MultiGrid::MultiGrid( const std::vector<Grid*> & grids, const std::vector<bool> & valid )
    : GridBase<MultiGrid>( "multi", space_from_grids( grids ), 
        shape_from_grids( grids ), valid_from_grids(grids, valid) ) {
    
    for (auto & k : grids) {
        if (k->klass()=="multi") {
            // merge
            auto other = dynamic_cast<MultiGrid*>(k);
            
            if (other==nullptr) {
                throw std::runtime_error("Internal error: cannot cast to MultiGrid");
            }
            
            for (unsigned int c=0; c<other->ngrids(); ++c) {
                grids_.emplace_back( other->subgrid(c).clone() );
            }
        } else {
            grids_.emplace_back( k->clone() );
        }
    }
}

// copy constructor
MultiGrid::MultiGrid( const MultiGrid & other )
    : GridBase<MultiGrid>( other ) {
    
    for (auto & k : other.grids_) {
        grids_.emplace_back( k->clone() );
    }
}

// properties
unsigned int MultiGrid::ngrids() const {
    return grids_.size();
}
const Grid & MultiGrid::subgrid(unsigned int index) {
    if (index>=grids_.size()) {
        throw std::runtime_error("Invalid subgrid index.");
    }
    
    return *grids_[index];
}

// method to compute probability
void MultiGrid::probability( const MultiSpace & space, value weight, const value * loc, const value * bw, value * result ) {
    if (space.nchildren() != grids_.size()) {
        throw std::runtime_error("Invalid number of grids/spaces.");
    }
    
    // loop through all subgrids and subspaces
    // evaluate and collect probabilities
    // combine into result vector
    std::vector<std::vector<value>> tmp( grids_.size() );
    
    for (unsigned int k=0; k<grids_.size(); ++k) {
        tmp[k].resize( grids_[k]->size() );
        space.child(k).probability( *grids_[k], 1., loc, bw, tmp[k].data() );
        loc += space.child(k).ndim();
        bw += space.child(k).nbw();
    }
    
    if (ninvalid()>0) {
        multiply_add_vectors( tmp, size(), weight, result, valid() );
    } else {
        multiply_add_vectors( tmp, size(), weight, result );
    }
    
}


// methods to compute partial log probability
void MultiGrid::partial_logp( const CategoricalSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // single grid with compatible space?
    if (grids_.size()!=1 || !(grids_[0]->specification()==space.specification())) {
        throw std::runtime_error("Incompatible space.");
    }
    
    grids_[0]->partial_logp( space, selection, factor, loc, bw, result );
    //space.partial_logp( *grids_[0], selection, factor, loc, bw, result );
}
void MultiGrid::partial_logp( const CircularSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // single grid with compatible space?
    if (grids_.size()!=1 || !(grids_[0]->specification()==space.specification())) {
        throw std::runtime_error("Incompatible space.");
    }
    
    grids_[0]->partial_logp( space, selection, factor, loc, bw, result );
    //space.partial_logp( *grids_[0], selection, factor, loc, bw, result );
}
void MultiGrid::partial_logp( const EncodedSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // single grid with compatible space?
    if (grids_.size()!=1 || !(grids_[0]->specification()==space.specification())) {
        throw std::runtime_error("Incompatible space.");
    }
    
    grids_[0]->partial_logp( space, selection, factor, loc, bw, result );
    //space.partial_logp( *grids_[0], selection, factor, loc, bw, result );
}
void MultiGrid::partial_logp( const EuclideanSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // single grid with compatible space?
    if (grids_.size()!=1 || !(grids_[0]->specification()==space.specification())) {
        throw std::runtime_error("Incompatible space.");
    }
    
    grids_[0]->partial_logp( space, selection, factor, loc, bw, result );
    //space.partial_logp( *grids_[0], selection, factor, loc, bw, result );
}
void MultiGrid::partial_logp( const MultiSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {

    // loop through all matching subspaces
    // if (std::count( selection, selection+subspace.ndim(), true)>0)
    // then subspace.log_probability( *grids_[index], 0., loc, bw, tmp[index].data() )
    // ++index
    // }
    // loc += subspace.ndim()
    // bw += subspace.nbw()
    // selection += subspace.ndim()
    
    // evaluate and collect probabilities
    // combine into result vector
    std::vector<std::vector<value>> tmp( grids_.size() );
    unsigned int index = 0;
    
    for (unsigned int k=0; k<space.nchildren(); ++k) {
        
        if (std::count( selection, selection+space.child(k).ndim(), true)>0) {
            
            tmp[index].resize( grids_[index]->size() );
            space.child(k).partial_logp( *grids_[index], selection, 0., loc, bw, tmp[index].data() );
            ++index;
            
        }
        loc += space.child(k).ndim();
        bw += space.child(k).nbw();
        selection += space.child(k).ndim();
        
    }
    
    if (ninvalid()>0) {
        add_assign_vectors( tmp, size(), factor, result, valid() );
    } else {
        add_assign_vectors( tmp, size(), factor, result );
    }
}


// yaml
YAML::Node MultiGrid::to_yaml_impl() const {
    YAML::Node node;
    for (auto & k : grids_) {
        node["grids"].push_back( k->to_yaml() );
    }
    return node;
}

std::unique_ptr<Grid> MultiGrid::from_yaml( const YAML::Node & node, 
    const SpaceSpecification & space, const std::vector<bool> & valid ) {
    std::vector<std::unique_ptr<Grid>> g;
    for (auto & k : node["grids"]) {
        g.push_back(std::move(grid_from_yaml( k )));
    }
    std::vector<Grid*> ptr;
    for (auto & k : g) {
        ptr.push_back( k.get() );
    }
    
    return std::make_unique<MultiGrid>( ptr, valid );
}


// flatbuffers
std::vector<flatbuffers::Offset<fb_serialize::Grid>> MultiGrid::to_flatbuffers_grids(flatbuffers::FlatBufferBuilder & builder) const {

    std::vector<flatbuffers::Offset<fb_serialize::Grid>> grid_vector;

    for (auto & k : grids_) {
        grid_vector.push_back(
            k->to_flatbuffers(builder)
        );
    }

    return grid_vector;
}

std::unique_ptr<MultiGrid> MultiGrid::from_flatbuffers(const fb_serialize::Grid * grid) {

    auto tmp = grid->valid();
    std::vector<bool> valid(tmp->begin(), tmp->end());

    unsigned int ndim = grid->grids()->size();

    std::vector<std::unique_ptr<Grid>> g;

    for (unsigned int k=0; k<ndim; ++k) {
        g.push_back(std::move(grid_from_flatbuffers(grid->grids()->Get(k))));
    }

    std::vector<Grid*> ptr;
    for (auto & k : g) {
        ptr.push_back( k.get() );
    }

    return std::make_unique<MultiGrid>( ptr, valid );
}


// hdf5
void MultiGrid::to_hdf5_impl(HighFive::Group & group) const {
    
    unsigned int g = 0;
    
    HighFive::Attribute attr = group.createAttribute<unsigned int>(
            "ndim", HighFive::DataSpace::From(g));
    attr.write(grids_.size());
    
    for (auto & k : grids_) {
        HighFive::Group subgroup = group.createGroup("grid" + std::to_string(g));
        
        k->to_hdf5(subgroup);
        
        ++g;
    }
}

std::unique_ptr<Grid> MultiGrid::from_hdf5(const HighFive::Group & group, 
    const SpaceSpecification & space, const std::vector<bool> & valid) {
    
    unsigned int ndim;
    
    HighFive::Attribute attr = group.getAttribute("ndim");
    attr.read(ndim);
    
    std::vector<std::unique_ptr<Grid>> g;
    for (unsigned int k=0; k<ndim; ++k) {
        HighFive::Group subgroup = group.getGroup("grid" + std::to_string(k));
        g.push_back(std::move(grid_from_hdf5(subgroup)));
    }
    
    std::vector<Grid*> ptr;
    for (auto & k : g) {
        ptr.push_back( k.get() );
    }
    
    return std::make_unique<MultiGrid>( ptr, valid );
    
}
    
