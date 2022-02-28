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
    
    long unsigned int product = std::accumulate(s.begin(), s.end(), 1, std::multiplies<long unsigned int>());
    
    if (product!=npoints) {
        std::runtime_error("Shape does not match number of points.");
    }
    
    return s;
}

// constructor
ArrayGrid::ArrayGrid( const std::vector<value> & array,
    const SpaceSpecification & space, const std::vector<bool> & valid,
    std::vector<long unsigned int> shape )
    : GridBase<ArrayGrid>( "array", space,
        shape_from_array_args(shape, array.size(), space.ndim()), valid),
        array_(array) {
    
    if (this->valid().size()>0 && array_.size()!=ndim()*nvalid()) {
        throw std::runtime_error("Expecting the same number of grid points "
            "in array as the number of valid points.");
    }
    
    unsigned int npoints = array_.size() / ndim();
    
    if ( (this->valid().size()==0 && npoints!=size()) || 
         (this->valid().size()>0 && npoints!=nvalid())) {
        throw std::runtime_error("Number of points in array incompatible "
            "with validity vector.");
    }
    
}


// methods to compute probability
void ArrayGrid::probability( const CategoricalSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
    // ignore valid vector for categorical variables
    value * ptr = array_.data();
    for (unsigned int n=0; n<array_.size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void ArrayGrid::probability( const CircularSpace & space, value weight,
    const value * loc, const value * bw, value * result )  {
    // ignore valid vector for circular variables
    value * ptr = array_.data();
    for (unsigned int n=0; n<array_.size(); ++n) {
        *result++ += weight*space.probability( loc, bw, ptr++ );
    }
}
void ArrayGrid::probability( const EncodedSpace & space, value weight,
    const value * loc, const value * bw, value * result )  {
    // ignore valid vector for encoded variables
    value * ptr = array_.data();
    if (ninvalid()>0) {
        auto vptr = valid().cbegin();
        // loop through all grid points
        for (unsigned int k=0; k<array_.size(); ++k) {
            if (*vptr++==true) {
                *result += weight*space.probability( loc, bw, ptr++ );
            }
            ++result;
        }
    } else {
        for (unsigned int n=0; n<array_.size(); ++n) {
            *result++ += weight*space.probability( loc, bw, ptr++ );
        }
    }
}
void ArrayGrid::probability( const EuclideanSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
    value * ptr = array_.data();
    
    if (ninvalid()>0) {
        auto vptr = valid().cbegin();
        // loop through all grid points
        for (unsigned int k=0; k<size(); ++k) {
            if (*vptr++==true) {
                *result += weight*space.probability( loc, bw, ptr );
                ptr += ndim();
            }
            ++result;
        }
    } else {
        for (unsigned int k=0; k<array_.size()/ndim(); ++k) {
            *result++ += weight*space.probability( loc, bw, ptr );
            ptr += ndim();
        }
    }

}

// methods to compute partial log probability
void ArrayGrid::partial_logp( const CategoricalSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, const value * loc, 
    const value * bw, value * result ) {
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
void ArrayGrid::partial_logp( const CircularSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, const value * loc, 
    const value * bw, value * result ) {
    
    value * ptr = array_.data();
    for (unsigned int k=0; k<array_.size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
    
    // alternative: first test for selection, then compute log probability
}
void ArrayGrid::partial_logp( const EncodedSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, const value * loc, 
    const value * bw, value * result ) {
    
    value * ptr = array_.data();
    for (unsigned int k=0; k<array_.size(); ++k) {
        *result++ = factor + space.partial_logp( loc, bw, ptr++, selection );
    }
}
void ArrayGrid::partial_logp( const EuclideanSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, const value * loc, 
    const value * bw, value * result ) {
    
    // for each point in grid array
    // compute partial_logp + factor and assign to result
    unsigned int npoints = array_.size() / ndim();
    value * ptr = array_.data();
    
    for (unsigned int k=0; k<npoints; ++k) {
        *result++ = factor + space.partial_logp(loc,bw,ptr,selection);
        ptr += ndim();
    }
    
}
void ArrayGrid::partial_logp( const MultiSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, const value * loc, 
    const value * bw, value * result ) {
    
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

// yaml
YAML::Node ArrayGrid::to_yaml_impl() const {
    YAML::Node node;
    node["array"] = array_;
    return node;
}

std::unique_ptr<Grid> ArrayGrid::from_yaml(const YAML::Node & node, 
    const SpaceSpecification & space, const std::vector<bool> & valid,
    std::vector<long unsigned int> shape ) {
    
    std::vector<value> a = node["array"].as<std::vector<value>>();
    return std::make_unique<ArrayGrid>( a, space, valid, shape );
}

// flatbuffers
std::vector<flatbuffers::Offset<fb_serialize::FloatArray>> ArrayGrid::to_flatbuffers_data(flatbuffers::FlatBufferBuilder & builder) const {

    std::vector<flatbuffers::Offset<fb_serialize::FloatArray>> floatarray_vector;

    floatarray_vector.push_back(
        fb_serialize::CreateFloatArray(
            builder,
            builder.CreateVector(array_)
        )
    );

    return floatarray_vector;
}

std::unique_ptr<ArrayGrid> ArrayGrid::from_flatbuffers(const fb_serialize::Grid * grid) {

    auto space = SpaceSpecification::from_flatbuffers(grid->space());

    auto tmp = grid->valid();
    std::vector<bool> valid(tmp->begin(), tmp->end());

    auto tmp_shape = grid->shape();
    std::vector<long unsigned int> shape(tmp_shape->begin(), tmp_shape->end());

    auto tmp_data = grid->data()->Get(0)->data();
    std::vector<value> data(tmp_data->begin(), tmp_data->end());

    return std::make_unique<ArrayGrid>(data, space, valid, shape);
}

// hdf5
void ArrayGrid::to_hdf5_impl(HighFive::Group & group) const {
    
    HighFive::DataSet dataset = group.createDataSet<value>(
        "array", HighFive::DataSpace::From(array_));
    
    dataset.write(array_);
    
}

std::unique_ptr<Grid> ArrayGrid::from_hdf5(const HighFive::Group & group, 
    const SpaceSpecification & space, const std::vector<bool> & valid, 
    std::vector<long unsigned int> shape ) {
    
    std::vector<value> a;
    HighFive::DataSet dataset = group.getDataSet("array");
    dataset.read(a);
    
    return std::make_unique<ArrayGrid>( a, space, valid, shape );
}
  
