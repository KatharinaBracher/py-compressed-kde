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
#include "stimulus.hpp"

#include <memory>

class PoissonLikelihood {
protected:
    // default constructor
    PoissonLikelihood();
public:
    // constructors
    PoissonLikelihood( Space & stimulus_space, Grid & grid, double stimulus_duration = 1., value compression=1. ); // stimulus duration is fixed
    PoissonLikelihood( Space & event_space, Space & stimulus_space, Grid & grid, double stimulus_duration = 1., value compression=1. ); // stimulus duration is fixed
    PoissonLikelihood( Space & event_space, std::shared_ptr<StimulusOccupancy> stimulus );
    PoissonLikelihood( std::shared_ptr<StimulusOccupancy> stimulus );
    
    // properties
    bool changed() const;
    bool random_insertion() const;
    void set_random_insertion(bool val);
    
    //value rate_offset() const;
    //void set_rate_offset(value val);
    
    value rate_scale() const;
    void set_rate_scale(value val);
    
    unsigned int ndim() const;
    unsigned int ndim_stimulus() const;
    unsigned int ndim_events() const;
    
    const Grid & grid() const;
    const Mixture & event_distribution() const;
    
    std::shared_ptr<StimulusOccupancy> stimulus();
    value mu() const;

    const std::vector<value> & stimulus_logp() const;
    const std::vector<value> & event_rate() const;

    //const std::vector<value> & offset() const;
        
    // methods
    void add_events( const std::vector<value> & events, unsigned int repetitions = 1 );
    void add_events( const value * events, unsigned int n, unsigned int repetitions = 1 );
    
    void precompute();

    void likelihood( value * events, unsigned int n, value delta_t, value * result );
    void logL( value * events, unsigned int n, value delta_t, value * result );
    void event_prob( value * events, unsigned int n, value * result );
    void event_logp( value * events, unsigned int n, value * result );
    
    // yaml
    YAML::Node to_yaml( bool save_stimulus=true ) const;
    
    void save_to_yaml( std::ostream & stream, bool save_stimulus=true ) const;
    void save_to_yaml( std::string path, bool save_stimulus=true ) const;
    
    static std::unique_ptr<PoissonLikelihood> from_yaml(const YAML::Node & node, 
        std::shared_ptr<StimulusOccupancy> stimulus = nullptr);
    
    // hdf5
    void to_hdf5(HighFive::Group & group, bool save_stimulus=true) const;
    
    static std::unique_ptr<PoissonLikelihood> from_hdf5(
        const HighFive::Group & group, std::shared_ptr<StimulusOccupancy> stimulus = nullptr);
    
    void save_to_hdf5(std::string filename, bool save_stimulus=true,
        int flags = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
        std::string path = "");
        
    static std::unique_ptr<PoissonLikelihood> load_from_hdf5(std::string filename, 
        std::string path="", std::shared_ptr<StimulusOccupancy> stimulus = nullptr);
    
protected:
    std::unique_ptr<Mixture> event_distribution_; // full space
    std::shared_ptr<StimulusOccupancy> stimulus_distribution_;
    std::unique_ptr<Grid> stimulus_grid_; // stimulus space
    
    std::vector<value> logp_stimulus_; // pi(x)
    std::vector<value> event_rate_; // p(x)
    std::unique_ptr<PartialMixture> p_event_; // p(a,x) @ x
    //std::vector<value> offset_;
    
    bool changed_;
    bool random_insertion_;
    
    //value rate_offset_;
    value rate_scale_;
};
