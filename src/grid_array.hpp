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
#pragma once

#include "grid_base.hpp"

std::vector<long unsigned int> shape_from_array_args( const std::vector<long unsigned int> & shape, unsigned int array_size, unsigned int ndim );

class ArrayGrid : public GridBase<ArrayGrid> {
public:
    // constructor
    ArrayGrid( const std::vector<value> & array, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
    // methods to compute probability
    virtual void probability( const CategoricalSpace & space, value weight, 
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const CircularSpace & space, value weight,
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const EncodedSpace & space, value weight,
        const value * loc, const value * bw, value * result ) override;
    virtual void probability( const EuclideanSpace & space, value weight,
        const value * loc, const value * bw, value * result ) override;
    
    // methods to compute partial log probability
    virtual void partial_logp( const CategoricalSpace & space, 
        std::vector<bool>::const_iterator selection, value factor,
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const CircularSpace & space,
        std::vector<bool>::const_iterator selection, value factor,
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EncodedSpace & space,
        std::vector<bool>::const_iterator selection, value factor,
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const EuclideanSpace & space,
        std::vector<bool>::const_iterator selection, value factor,
        const value * loc, const value * bw, value * result ) override;
    virtual void partial_logp( const MultiSpace & space,
        std::vector<bool>::const_iterator selection, value factor,
        const value * loc, const value * bw, value * result ) override;
    
    // yaml
    virtual YAML::Node to_yaml_impl() const;
    static std::unique_ptr<Grid> from_yaml( const YAML::Node & node, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
    // hdf5
    virtual void to_hdf5_impl(HighFive::Group & group) const;
    static std::unique_ptr<Grid> from_hdf5(const HighFive::Group & group, const SpaceSpecification & space, const std::vector<bool> & valid, std::vector<long unsigned int> shape );
    
    virtual void at_index(const unsigned int * index, value * result) const {
        
        //unsigned int npoints = array_.size()/ndim();
        
        // using shape(), convert index to linear index into array
        // for this conversion, the strides have to be computed
        // next, get the ndim() values at linear index in array
        // if index out of range, set result to std::numeric_limits<value>::quiet_NaN();
        
        std::runtime_error("ArrayGrid::at_index has not been implemented.");
        
    }
    
protected:
    std::vector<value> array_;
};
