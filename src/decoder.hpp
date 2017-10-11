#pragma once

#include "common.hpp"
#include "likelihood.hpp"

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
