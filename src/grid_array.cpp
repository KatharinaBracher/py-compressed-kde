#include "grid_array.hpp"
#include "space.hpp"

std::vector<long unsigned int> shape_from_array_args( const std::vector<long unsigned int> & shape, unsigned int array_size, unsigned int ndim ) {
    
    std::vector<long unsigned int> s;
    
    long unsigned int npoints = array_size / ndim;
    
    if (npoints * ndim != array_size) {
        std::runtime_error("Array does not contain multiple of ndim samples.");
    }
    
    if (shape.size()==0) {
        s.push_back( npoints );
    } else if (shape.size()!=ndim) {
        std::runtime_error("Incompatible array shape.");
    } else {
        s = shape;
    }
    
    return s;
}

ArrayGrid::ArrayGrid( const std::vector<value> & array, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape ) :
GridBase<ArrayGrid>( "array", space, shape_from_array_args( shape, array.size(), space.ndim() ), valid ), array_(array) {
    
    if (this->valid().size()>0 && array_.size()!=ndim()*nvalid()) {
        throw std::runtime_error("Incorrect number of points in array.");
    }
    
    unsigned int npoints = array_.size() / ndim();
    
    if ( (this->valid().size()==0 && npoints!=size()) || (this->valid().size()>0 && npoints!=nvalid())) {
        throw std::runtime_error("Number of points in array incompatible with validity vector.");
    }
    
}

YAML::Node ArrayGrid::asYAML() const {
    YAML::Node node;
    node["array"] = array_;
    return node;
}

Grid * ArrayGrid::fromYAML( const YAML::Node & node, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape ) {
    std::vector<value> a = node["array"].as<std::vector<value>>();
    return new ArrayGrid( a, space, valid, shape );
}

void ArrayGrid::probability( const CategoricalSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for categorical variables
    value * ptr = array_.data();
    for (unsigned int n=0; n<array_.size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void ArrayGrid::probability( const CircularSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for circular variables
    value * ptr = array_.data();
    for (unsigned int n=0; n<array_.size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void ArrayGrid::probability( const EncodedSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    // ignore valid vector for encoded variables
    value * ptr = array_.data();
    for (unsigned int n=0; n<array_.size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void ArrayGrid::probability( const EuclideanSpace & space, value weight, const value * loc, const value * bw, value * result )  {
    value * ptr = array_.data();
    
    if (nvalid()>0) {
        auto vptr = valid().cbegin();
        // loop through all grid points
        for (unsigned int k=0; k<size(); ++k) {
            if (*vptr++==true) {
                *result = space.probability( loc, bw, ptr );
                ptr += ndim();
            }
            ++result;
        }
    } else {
        for (unsigned int k=0; k<array_.size()/ndim(); ++k) {
            *result++ += space.probability( loc, bw, ptr );
            ptr += ndim();
        }
    }

}

void ArrayGrid::partial_logp( const CategoricalSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    if (*selection) {
        for (auto & k : array_) {
            if (static_cast<unsigned int>(*loc)!=static_cast<unsigned int>(k)) {
                *result++ = -std::numeric_limits<value>::infinity();
            } else {
                *result++ = factor;
            }
        }
    }
}
void ArrayGrid::partial_logp( const CircularSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    
    value * ptr = array_.data();
    for (unsigned int k=0; k<array_.size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
    
    // alternative: first test for selection, then compute log probability
}
void ArrayGrid::partial_logp( const EncodedSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    
    value * ptr = array_.data();
    for (unsigned int k=0; k<array_.size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
}
void ArrayGrid::partial_logp( const EuclideanSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
    
    // for each point in grid array
    // compute partial_logp + factor and assign to result
    unsigned int npoints = array_.size() / ndim();
    value * ptr = array_.data();
    
    for (unsigned int k=0; k<npoints; ++k) {
        *result++ = factor + space.partial_logp(loc,bw,ptr,selection);
        ptr += ndim();
    }
    
}
void ArrayGrid::partial_logp( const MultiSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
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

