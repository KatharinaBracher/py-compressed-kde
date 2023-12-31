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
#include "space_categorical.hpp"

// constructor
CategoricalSpace::CategoricalSpace( std::string name, std::vector<std::string> labels, unsigned int category ) :
SpaceBase<CategoricalSpace>( "categorical", make_spec(name, labels), make_kernel(category) ), labels_(labels) {
    if (!isunique(labels)) {
        throw std::runtime_error("Labels are not unique.");
    }
}

SpaceSpecification CategoricalSpace::make_spec( std::string name, std::vector<std::string> labels ) {
    
    std::string extra = "labels=[";
    for (auto & k : labels) {
        extra += k + ",";
    }
    extra += "]";
    
    return SpaceSpecification( DimSpecification( name, "categorical", extra ) );
    
}

Component CategoricalSpace::make_kernel(unsigned int category) const {
    Component k;
    k.location = { static_cast<value>(category) };
    k.scale_factor = 1;
    k.scale_factor_log = 0.;
    return k;
}

Grid * CategoricalSpace::grid() const {
    
    std::vector<value> v(labels_.size());
    std::iota( v.begin(), v.end(), 0 );
    
    return new VectorGrid( { v }, specification(), {} );
    
}
 
// methods
value CategoricalSpace::compute_scale_factor( value * bw, bool log ) const {
    if (log) { return 0.; } else { return 1.; };
}

value CategoricalSpace::compute_scale_factor( 
    std::vector<bool>::const_iterator selection, value * bw, bool log ) const {
    if (log) { return 0.; } else { return 1.; };
}

value CategoricalSpace::mahalanobis_distance_squared( const value * refloc, 
    const value * refbw, const value * targetloc, value threshold) const {
    if (static_cast<unsigned int>(*refloc) == static_cast<unsigned int>(*targetloc)) { 
        return 0.;
    } else {
        return std::numeric_limits<value>::infinity();
    }
}

void CategoricalSpace::merge( value w1, value * loc1, value * bw1, value w2, 
    const value * loc2, const value * bw2 ) const {
    // nothing to do
}

value CategoricalSpace::probability( const value * loc, const value * bw, 
    const value * point ) const {
    return static_cast<value>( static_cast<unsigned int>(*loc)==static_cast<unsigned int>(*point) );
}

void CategoricalSpace::probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = probability( loc, bw, points++ );
    }
}

value CategoricalSpace::log_probability( const value * loc, const value * bw, 
    const value * point ) const {
    value p = 0;
    if (static_cast<unsigned int>(*loc)!=static_cast<unsigned int>(*point)) {
        p = -std::numeric_limits<value>::infinity();
    }
    return p;
}

void CategoricalSpace::log_probability( const value * loc, const value * bw, 
    const value * points, unsigned int n, value * result ) const {
    for (unsigned int k=0; k<n; ++k) {
        *result++ = log_probability( loc, bw, points++ );
    }
}

value CategoricalSpace::partial_logp( const value * loc, const value * bw, 
    const value * point, std::vector<bool>::const_iterator selection ) const {
    value p=0.;
    
    if (*selection && !(static_cast<unsigned int>(*loc)==static_cast<unsigned int>(*point)) ){
        p = -std::numeric_limits<value>::infinity();
    }
    
    return p;
}


// yaml
std::unique_ptr<CategoricalSpace> CategoricalSpace::from_yaml(
    const YAML::Node & node ) {
    
    if (!node["name"] || !node["labels"]) {
        throw std::runtime_error("Ill-formed categorical space definition.");
    }
    
    return std::make_unique<CategoricalSpace>(node["name"].as<std::string>(), 
        node["labels"].as<std::vector<std::string>>());
    
} 

YAML::Node CategoricalSpace::to_yaml_impl() const {
    YAML::Node node;
    node["name"] = specification().dim(0).name();
    node["labels"] = labels_;
    return node;
}

// flatbuffers
flatbuffers::Offset<fb_serialize::SpaceData> CategoricalSpace::to_flatbuffers_impl(flatbuffers::FlatBufferBuilder & builder) const {

    return fb_serialize::CreateSpaceData(
        builder,
        fb_serialize::SpaceType_CategoricalSpace,
        fb_serialize::CreateCategoricalSpace(
            builder,
            builder.CreateString(specification().dim(0).name()),
            builder.CreateVectorOfStrings(labels_)
        ).Union()
    );
}

std::unique_ptr<CategoricalSpace> CategoricalSpace::from_flatbuffers(const fb_serialize::Space * space) {

    auto saved_klass = space->klass()->str();

    if (saved_klass!="categorical") {
        throw std::runtime_error("Expected categorical, but got " + saved_klass);
    }

    auto default_kernel = components_from_flatbuffers(space->default_kernel());

    auto data = space->data()->value_as_CategoricalSpace();

    std::string name = data->name()->str();

    std::vector<std::string> labels;

    for (auto k : *data->labels()) {
        labels.push_back(k->str());
    }

    auto ptr = std::make_unique<CategoricalSpace>(name, labels);

    ptr->set_default_kernel(*(default_kernel[0]));

    return ptr;
}


// hdf5
void CategoricalSpace::to_hdf5_impl(HighFive::Group & group) const {
    
    std::string name = specification().dim(0).name();
    
    HighFive::DataSet ds_name = group.createDataSet<std::string>(
        "name", HighFive::DataSpace::From(name));
    ds_name.write(name);
    
    HighFive::DataSet ds_labels = group.createDataSet<std::string>(
        "labels", HighFive::DataSpace::From(labels_));
    ds_labels.write(labels_);
}

std::unique_ptr<CategoricalSpace> CategoricalSpace::from_hdf5(
    const HighFive::Group & group) {
    
    std::string name;
    std::vector<std::string> labels;
    
    HighFive::DataSet ds_name = group.getDataSet("name");
    ds_name.read(name);
    
    HighFive::DataSet ds_labels = group.getDataSet("labels");
    ds_labels.read(labels);
    
    return std::make_unique<CategoricalSpace>(name, labels);
}
