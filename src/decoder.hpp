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

#include <memory>

inline std::vector<std::vector<value>> load_prior_from_file(std::string filename, std::string path = ""){
    HighFive::File file(filename, HighFive::File::ReadOnly);
    unsigned int nunion;
    if (path.empty()) {
        path = "/";
    }

    HighFive::Group group = file.getGroup(path);
    // load priors
    group.getDataSet("nunion").read(nunion);
    std::vector<std::vector<value>> priors(nunion);

    HighFive::Group grp_priors = group.getGroup("priors");

    for (unsigned int k=0; k<nunion; ++k) {
        grp_priors.getDataSet("prior" + std::to_string(k)).read(priors[k]);
    }
    return priors;
};
/**
 * @brief compute_posterior based on likelihood result with multiple stimulus spaces (unions)
 * @param result number of unions x grid size
 * @param normalize
 */
void compute_posterior(std::vector<value *> result,
                       std::vector<std::vector<value>> prior,
                       std::vector<unsigned int> grid_sizes, bool normalize);
/**
 * @brief compute_posterior based on likelihood result with 1 stimulus space
 * @param result - grid size
 * @param normalize
 */
void compute_posterior(value * result, std::vector<value> prior, unsigned int grid_size, bool normalize);

class Decoder {
public:
    // constructors
    Decoder( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihoods,
        const std::vector<value> & prior = {} );
    Decoder( std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> & likelihoods,
        const std::vector<std::vector<value>> & prior = {} );
    
    // decoding methods

    /**
     * @brief decode with multiple sources and multiple union
     * @param events  each element of the vector is a pointer to an array of events for one source
     * @param nevents each element of the vector contains the number of event for one source
     * @param delta_t size of the time bin in which events are observed
     * @param result pre-initialized to be the size of the number of stimulus space and then contains in each element an array of grid size
     * @param normalize
     */
    void decode( std::vector<value*> events, std::vector<unsigned int> nevents,
        value delta_t, std::vector<value*> result, bool normalize=true );

    /**
     * @brief decode with multiple sources and 1 stimulus space
     * @param events  each element of the vector contains the event for one source
     * @param nevents each element of the vector contains the number of event for one source
     * @param delta_t size of the time bin in which events are observed
     * @param result pre-initialized array of grid size
     * @param index index of the stimulus space
     * @param normalize
     */
    void decode ( std::vector<value*> events, std::vector<unsigned int> nevents,
        value delta_t, value* result, unsigned int index=0, bool normalize=true );

    /**
     * @brief decode with multiple sources with multiple stimulus spaces (union) - used to reshape the events dimension from a std::vector to an array before calling
     * the decode method upper
     * @param events first dim is the number of sources, second dim the number of events per source
     * @param delta_t size of the time bin in which events are observed
     * @param result pre-initialized to be the size of the number of stimulus space and then contains in each element an array of grid size (?)
     * @param normalize
     */
    void decode ( std::vector<std::vector<value>> events, value delta_t,
        std::vector<value*> result, bool normalize=true );

    /**
     * @brief decode with multiple sources with 1 stimulus space - used to reshape the events dimension from a std::vector to an array before calling
     * the decode method upper
     * @param events first dim is the number of sources, second dim the number of events per source
     * @param delta_t size of the time bin in which events are observed
     * @param result pre-initialized array of grid size
     * @param normalize
     */
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
