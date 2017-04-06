#include "grid_vector.hpp"
#include "space.hpp"

#include <iostream>

std::vector<long unsigned int> shape_from_vectors( const std::vector<std::vector<value>> & vectors ) {
    std::vector<long unsigned int> shape( vectors.size() );
    for (unsigned int k=0; k<vectors.size(); ++k) {
        shape[k] = vectors[k].size();
    }
    return shape;
}

VectorGrid::VectorGrid( const std::vector<std::vector<value>> & vectors, const SpaceSpecification & space, const std::vector<bool> & valid ) : 
GridBase<VectorGrid>( "vector", space, shape_from_vectors(vectors), valid ), vectors_(vectors) {
    
    if (shape().size() != vectors_.size() || ndim()!=vectors_.size()) {
        throw std::runtime_error("Incompatible number of vectors.");
    }
    
    ptemp_ = vectors_;
      
}

YAML::Node VectorGrid::asYAML() const {
    YAML::Node node;
    node["vectors"] = vectors_;
    return node;
}

Grid * VectorGrid::fromYAML( const YAML::Node & node, const SpaceSpecification & space, const std::vector<bool> & valid ) {
    std::vector<std::vector<value>> v = node["vectors"].as<std::vector<std::vector<value>>>();
    return new VectorGrid( v, space, valid );
}

void VectorGrid::probability( const CategoricalSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for categorical variables
    value * ptr = vectors_[0].data();
    for (unsigned int n=0; n<vectors_[0].size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void VectorGrid::probability( const CircularSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for circular variables
    value * ptr = vectors_[0].data();
    for (unsigned int n=0; n<vectors_[0].size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void VectorGrid::probability( const EncodedSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for encoded variables
    value * ptr = vectors_[0].data();
    for (unsigned int n=0; n<vectors_[0].size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void VectorGrid::probability( const EuclideanSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    for (unsigned int k=0; k<vectors_.size(); ++k) {
        space.probability( loc++, bw++, vectors_[k].data(), vectors_[k].size(), ptemp_[k].data() );
    }
    
    if (nvalid()>0) {
        multiply_add_vectors( ptemp_, size(), weight, result, valid() );
    } else {
        multiply_add_vectors( ptemp_, size(), weight, result );
    }
}


void VectorGrid::partial_logp( const CategoricalSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    if (*selection) {
        for (auto & k : vectors_[0]) {
            if (static_cast<unsigned int>(*loc)!=static_cast<unsigned int>(k)) {
                *result++ = -std::numeric_limits<value>::infinity();
            } else {
                *result++ = factor;
            }
        }
    }
}
void VectorGrid::partial_logp( const CircularSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    
    value * ptr = vectors_[0].data();
    for (unsigned int k=0; k<vectors_[0].size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
    
    // alternative: first test for selection, then compute log probability
}
void VectorGrid::partial_logp( const EncodedSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    
    value * ptr = vectors_[0].data();
    for (unsigned int k=0; k<vectors_[0].size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
}
void VectorGrid::partial_logp( const EuclideanSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // for each selected grid vector compute log probability
    // add selected vectors and assign to result
    
    unsigned int index = 0;
    
    for (unsigned int k=0; k<space.ndim(); ++k) {
        if (*selection++) {
            space.log_probability( loc, bw, vectors_[index].data(), vectors_[index].size(), ptemp_[index].data() );
            ++index;
        }
        ++loc;
        ++bw;
    }
    
    if (nvalid()>0) {
        add_assign_vectors( ptemp_, size(), factor, result, valid() );
    } else {
        add_assign_vectors( ptemp_, size(), factor, result );
    }
}

void VectorGrid::partial_logp( const MultiSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    // search for child space that has same specification
    unsigned int index;
    for (index=0; index<space.nchildren(); ++index) {
        if (space.child(index).specification()==specification()) { break; }
        selection+=space.child(index).ndim();
        loc+=space.child(index).ndim();
        bw+=space.child(index).nbw();
    }
    
    if (index>=space.nchildren()) {
        throw std::runtime_error("Incompatible space.");
    }
    
    space.child(index).partial_logp( *this, selection, factor, loc, bw, result );
    
}
