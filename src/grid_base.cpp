#include "grid_base.hpp"

// constructor
Grid::Grid( std::string klass, const SpaceSpecification & space, 
    std::vector<long unsigned int> shape, const std::vector<bool> & valid )
    : klass_(klass), spec_(space), shape_(shape), valid_(valid) {

    if (valid_.size()!=0 && valid_.size()!=size()) {
        throw std::runtime_error("Incompatible size of valid vector.");
    }

    if (shape.size()!=1 && shape.size()!=ndim()) {
        throw std::runtime_error("Incompatible shape vector.");
    }
    
    nvalid_ = std::count( valid_.cbegin(), valid_.cend(), true );
    
}

// clone
Grid* Grid::clone() const { throw std::runtime_error("not implemented"); }

// properties
std::string Grid::klass() const {
    return klass_;
}

const std::vector<long unsigned int> & Grid::shape() const { return shape_; }
unsigned int Grid::size() const { 
    return std::accumulate(shape_.begin(), shape_.end(), 1.,
        std::multiplies<long unsigned int>() );
}
unsigned int Grid::ndim() const { return spec_.ndim(); }

const std::vector<bool> & Grid::valid() const { return valid_; }
unsigned int Grid::nvalid() const { return nvalid_; }

// space
const SpaceSpecification & Grid::specification() const {
    return spec_;
}

//void Grid::probability( const Space & space, value weight, const Component & k, value * result ) {
    //probability( space, weight, k.location.data(), k.bandwidth.data(), result );
//}

// methods to compute probability
void Grid::probability( const Space & space, value weight, const value * loc, 
    const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
 }
void Grid::probability( const CategoricalSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented CategoricalSpace");
 }
void Grid::probability( const CircularSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented CircularSpace");
 }
void Grid::probability( const EncodedSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented EncodedSpace");
 }
void Grid::probability( const EuclideanSpace & space, value weight, 
    const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented EuclideanSpace");
}
void Grid::probability( const MultiSpace & space, value weight, 
    const value * loc, const value * bw, value * result ) {
    throw std::runtime_error("Not implemented MultiSpace");
}
    

//void Grid::partial_logp( const Space & space, const Component & k, const std::vector<bool> & selection, value * result ) {}

// methods to compute partial log probability
void Grid::partial_logp( const Space & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const CategoricalSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const CircularSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const EncodedSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const EuclideanSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const MultiSpace & space, 
    std::vector<bool>::const_iterator selection, value factor, 
    const value * loc, const value * bw, value * result ) {
    
    throw std::runtime_error("Not implemented: Space");
}

//void Grid::marginal( const Space & space, const Component & k, const std::vector<bool> & selection, value * result ) {}

// yaml
YAML::Node Grid::to_yaml() const {
    YAML::Node node;
    node["class"] = klass();
    node["space"] = spec_.to_yaml();
    node["shape"] = shape_;
    node["valid"] = valid_;
    node["grid"] = this->to_yaml_impl();
    return node;
}

YAML::Node Grid::to_yaml_impl() const {
    throw std::runtime_error("Not implemented.");
}

// hdf5
void Grid::to_hdf5(HighFive::Group & group) const {
    
    HighFive::Attribute attr = group.createAttribute<std::string>(
            "class", HighFive::DataSpace::From(klass()));
    attr.write(klass());
    
    HighFive::Group space = group.createGroup("space");
    spec_.to_hdf5(space);
    
    HighFive::DataSet shape = group.createDataSet<long unsigned int>(
        "shape", HighFive::DataSpace::From(shape_));
    shape.write(shape_);
    
    // no support for writing std::vector<bool>
    // work around: copy to vector
    std::vector<unsigned char> tmp(valid_.begin(), valid_.end());
    
    HighFive::DataSet valid = group.createDataSet<unsigned char>(
        "valid", HighFive::DataSpace::From(tmp));
    valid.write(tmp);
    
    HighFive::Group grid = group.createGroup("grid");
    this->to_hdf5_impl(grid);
    
}

void Grid::to_hdf5_impl(HighFive::Group & group) const {
    throw std::runtime_error("Not implemented.");
}



void Grid::at_index(const unsigned int * index, value * result) const {
    throw std::runtime_error("Not implemented.");
}
