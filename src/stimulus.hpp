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
#include "space.hpp"
#include "grid.hpp"
#include "mixture.hpp"
#include "schema_generated.h"

#include <memory>
#include <mutex>

class StimulusOccupancy {
public:
    // constructor
    StimulusOccupancy( const Space & space, const Grid & grid, double stimulus_duration = 1., value compression = 1.2 );
    
    // properties
    const Space & space() const;
    const Grid & grid() const;
    value compression() const;
    double stimulus_duration() const;
    unsigned int ndim() const;
    
    bool random_insertion() const;
    void set_random_insertion(bool val);
        
    value stimulus_time();
    
    void occupancy( std::vector<value> & out );
    void prob( std::vector<value> & out ) ;
    void logp( std::vector<value> & out );
    
    void occupancy( value * out );
    void prob( value * out );
    void logp( value * out );
    
    // methods
    void add_stimulus( const std::vector<value> & stimuli, unsigned repetitions = 1 );
    void add_stimulus( const value * stimuli, unsigned int n, 
        unsigned int repetitions = 1 );
    
    // yaml
    YAML::Node to_yaml( ) const;
    void save_to_yaml( std::ostream & stream ) const;
    void save_to_yaml( std::string path ) const;
    static std::unique_ptr<StimulusOccupancy> from_yaml(const YAML::Node & node);
    
    // flatbuffers
    flatbuffers::Offset<fb_serialize::StimulusOccupancy> to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const;
    static std::unique_ptr<StimulusOccupancy> from_flatbuffers(const fb_serialize::StimulusOccupancy * stimulus);

    // hdf5
    void to_hdf5(HighFive::Group & group) const;
    static std::unique_ptr<StimulusOccupancy> from_hdf5(const HighFive::Group & group);
    
    void save_to_hdf5(std::string filename,
        int flags = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
        std::string path = "");
        
    static std::unique_ptr<StimulusOccupancy> load_from_hdf5(std::string filename, 
        std::string path="");
    
protected:
    double stimulus_duration_;
    value compression_;
    bool random_insertion_;
    std::unique_ptr<Mixture> stimulus_distribution_;
    std::unique_ptr<Grid> stimulus_grid_;
    std::mutex lock_;
};
