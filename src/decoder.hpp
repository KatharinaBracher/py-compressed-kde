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

#include "common.hpp"
#include "likelihood.hpp"
#include "schema_generated.h"

#include <memory>

class Decoder {
public:
    // constructors
    Decoder( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihoods,
        const std::vector<value> & prior = {} );
    Decoder( std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> & likelihoods,
        const std::vector<std::vector<value>> & prior = {} );
    
    // decoding methods
    void decode( std::vector<value*> events, std::vector<unsigned int> nevents, 
        value delta_t, std::vector<value*> result, bool normalize=true );
    
    void decode ( std::vector<value*> events, std::vector<unsigned int> nevents,
        value delta_t, value* result, unsigned int index=0, bool normalize=true );
    
    void decode ( std::vector<std::vector<value>> events, value delta_t,
        std::vector<value*> result, bool normalize=true );
    
    void decode ( std::vector<std::vector<value>> events, value delta_t,
        value* result, unsigned int index=0, bool normalize=true );
    
    // properties
    unsigned int nsources() const;
    bool is_union() const;
    unsigned int n_union() const;
    unsigned int grid_size(unsigned int index=0) const;
    const std::vector<unsigned int> & grid_sizes() const;
    
    std::vector<long unsigned int> grid_shape(unsigned int index=0) const;
    const std::vector<std::vector<long unsigned int>> & grid_shapes() const;
    
    const Grid & grid(unsigned int index=0) const;
    
    std::shared_ptr<StimulusOccupancy> stimulus(unsigned int index=0);
    
    // likelihood getter
    std::shared_ptr<PoissonLikelihood> likelihood( unsigned int source,
        unsigned int index = 0 );
    
    //void add_likelihood( std::shared_ptr<PoissonLikelihood> likelihood );
    //void add_likelihood( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihood );
    //void remove_likelihood( unsigned int source );
    
    // enable/disable sources
    unsigned int nenabled_sources() const;
    void enable_source( unsigned int source );
    void enable_all_sources();
    void enable_one_source( unsigned int source );
    void enable_sources( const std::vector<bool> & state );
    void disable_source( unsigned int source );
    const std::vector<bool> & enabled_sources() const;
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::Decoder> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    static std::unique_ptr<Decoder> from_flatbuffers(const fb_serialize::Decoder * decoder);

    // hdf5 
    void to_hdf5(HighFive::Group & group) const;
    static std::unique_ptr<Decoder> from_hdf5(const HighFive::Group & group);
    
    void save_to_hdf5(std::string filename,
        int flags = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
        std::string path = "");
        
    static std::unique_ptr<Decoder> load_from_hdf5(std::string filename,
        std::string path="");
    
protected:
    // outer vector: sources
    // inner vector: union
    std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> likelihoods_;
    
    // prior for each union member
    std::vector<std::vector<value>> prior_;
    
    // grid size for each union member
    std::vector<unsigned int> grid_sizes_;
    std::vector<std::vector<long unsigned int>> grid_shapes_;
    
    std::vector<bool> likelihood_selection_;
};
