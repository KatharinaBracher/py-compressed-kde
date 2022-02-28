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
#include "kernel.hpp"
#include <stdexcept>

// yaml
std::unique_ptr<Kernel> kernel_from_yaml( const YAML::Node & node ) {
    
    if (!node.IsMap() || !node["type"] ) {
        throw std::runtime_error("Not a valid YAML description of kernel.");
    }
    
    KernelType ktype = kerneltype_fromstring( node["type"].as<std::string>( "unknown" ) );
    
    if (ktype == KernelType::Gaussian) {
        return GaussianKernel::from_yaml( node["info"] );
    } else if (ktype == KernelType::Epanechnikov) {
        return EpanechnikovKernel::from_yaml( node["info"] );
    } else if (ktype == KernelType::Box) {
        return BoxKernel::from_yaml( node["info"] );
    } else {
        throw std::runtime_error("Unknown kernel type.");
    }

}

// flatbuffer
std::unique_ptr<Kernel> kernel_from_flatbuffers(const fb_serialize::Kernel * kernel) {

    std::string type_str = kernel->type()->str();
    KernelType ktype = kerneltype_fromstring( type_str );

    if (ktype == KernelType::Gaussian) {
        return GaussianKernel::from_flatbuffers(kernel->kernel_as_GaussianKernel());
    } else if (ktype == KernelType::Epanechnikov) {
        return EpanechnikovKernel::from_flatbuffers(kernel->kernel_as_EpanechnikovKernel());
    } else if (ktype == KernelType::Box) {
        return BoxKernel::from_flatbuffers(kernel->kernel_as_BoxKernel());
    } else {
        throw std::runtime_error("Unknown kernel type: "+type_str);
    }
}


// hdf5
std::unique_ptr<Kernel> kernel_from_hdf5( const HighFive::Group & group ) {
    
    std::string type_str;
    HighFive::DataSet ds_type = group.getDataSet("type");
    ds_type.read(type_str);
    
    KernelType ktype = kerneltype_fromstring( type_str );
    
    if (ktype == KernelType::Gaussian) {
        return GaussianKernel::from_hdf5( group.getGroup("info") );
    } else if (ktype == KernelType::Epanechnikov) {
        return EpanechnikovKernel::from_hdf5( group.getGroup("info") );
    } else if (ktype == KernelType::Box) {
        return BoxKernel::from_hdf5( group.getGroup("info") );
    } else {
        throw std::runtime_error("Unknown kernel type.");
    }
    
}




