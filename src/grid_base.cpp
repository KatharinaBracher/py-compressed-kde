#include "grid_base.hpp"

Grid::Grid( std::string klass, const SpaceSpecification & space, std::vector<long unsigned int> shape, const std::vector<bool> & valid ) :
klass_(klass), spec_(space), shape_(shape), valid_(valid) {

    if (valid_.size()!=0 && valid_.size()!=size()) {
        throw std::runtime_error("Incompatible size of valid vector.");
    }

    if (shape.size()!=1 && shape.size()!=ndim()) {
        throw std::runtime_error("Incompatible shape vector.");
    }
    
    nvalid_ = std::count( valid_.cbegin(), valid_.cend(), true );
    
}

Grid* Grid::clone() const { throw std::runtime_error("not implemented"); }

std::string Grid::klass() const {
    return klass_;
}

const std::vector<long unsigned int> & Grid::shape() const { return shape_; }
unsigned int Grid::size() const { return std::accumulate( shape_.begin(), shape_.end(), 1., std::multiplies<long unsigned int>() ); }
unsigned int Grid::ndim() const { return spec_.ndim(); }

const std::vector<bool> & Grid::valid() const { return valid_; }
unsigned int Grid::nvalid() const { return nvalid_; }

const SpaceSpecification & Grid::specification() const {
    return spec_;
}

YAML::Node Grid::toYAML() const {
    YAML::Node node;
    node["class"] = klass();
    node["space"] = spec_.toYAML();
    node["shape"] = shape_;
    node["valid"] = valid_;
    node["grid"] = this->asYAML();
    return node;
}
YAML::Node Grid::asYAML() const {
    throw std::runtime_error("Not implemented.");
}


//void Grid::probability( const Space & space, value weight, const Component & k, value * result ) {
    //probability( space, weight, k.location.data(), k.bandwidth.data(), result );
//}

void Grid::probability( const Space & space, value weight, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
 }
void Grid::probability( const CategoricalSpace & space, value weight, const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented CategoricalSpace");
 }
void Grid::probability( const CircularSpace & space, value weight, const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented CircularSpace");
 }
void Grid::probability( const EncodedSpace & space, value weight, const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented EncodedSpace");
 }
void Grid::probability( const EuclideanSpace & space, value weight, const value * loc, const value * bw, value * result )  {
     throw std::runtime_error("Not implemented EuclideanSpace");
}
void Grid::probability( const MultiSpace & space, value weight, const value * loc, const value * bw, value * result ) {
    throw std::runtime_error("Not implemented MultiSpace");
}
    

//void Grid::partial_logp( const Space & space, const Component & k, const std::vector<bool> & selection, value * result ) {}

void Grid::partial_logp( const Space & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const CategoricalSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const CircularSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const EncodedSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const EuclideanSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}
void Grid::partial_logp( const MultiSpace & space, std::vector<bool>::const_iterator selection, value factor, const value * loc, const value * bw, value * result ) {
     throw std::runtime_error("Not implemented: Space");
}

void Grid::marginal( const Space & space, const Component & k, const std::vector<bool> & selection, value * result ) {}
