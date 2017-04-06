#pragma once

#include "common.hpp"
#include "space.hpp"
#include "grid.hpp"
#include "mixture.hpp"

#include <memory>
#include <mutex>

class StimulusOccupancy {
public:
    StimulusOccupancy( Space & space, Grid & grid, double stimulus_duration = 1., value compression = 1.2 );
    
    const Space & space() const;
    const Grid & grid() const;
    value compression() const;
    double stimulus_duration() const;
    unsigned int ndim() const;
    
    bool random_insertion() const;
    void set_random_insertion(bool val);
    
    void add_stimulus( const std::vector<value> & stimuli, unsigned repetitions = 1 );
    
    void add_stimulus( const value * stimuli, unsigned int n, unsigned int repetitions = 1 );
    
    value stimulus_time();
    
    void occupancy( std::vector<value> & out );
    
    void prob( std::vector<value> & out ) ;
    
    void logp( std::vector<value> & out );
    
    void occupancy( value * out );
    
    void prob( value * out );
    
    void logp( value * out );
    
    YAML::Node toYAML( ) const;
    
    void save( std::ostream & stream ) const;
    
    void save( std::string path ) const;
    
protected:
    double stimulus_duration_;
    value compression_;
    bool random_insertion_;
    std::unique_ptr<Mixture> stimulus_distribution_;
    std::unique_ptr<Grid> stimulus_grid_;
    std::mutex lock_;
};

class PoissonLikelihood {
public:
    PoissonLikelihood( Space & stimulus_space, Grid & grid, double stimulus_duration = 1., value compression=1. ); // stimulus duration is fixed
    PoissonLikelihood( Space & event_space, Space & stimulus_space, Grid & grid, double stimulus_duration = 1., value compression=1. ); // stimulus duration is fixed
    PoissonLikelihood( Space & event_space, std::shared_ptr<StimulusOccupancy> stimulus );
    PoissonLikelihood( std::shared_ptr<StimulusOccupancy> stimulus );
    
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
    
    void add_events( const std::vector<value> & events, unsigned int repetitions = 1 );
    
    void add_events( const value * events, unsigned int n, unsigned int repetitions = 1 );
    
    value mu() const;
    
    const std::vector<value> & stimulus_logp() const;
    
    const std::vector<value> & event_rate() const;
    
    //const std::vector<value> & offset() const;
    
    void precompute();
    
    void likelihood( value * events, unsigned int n, value delta_t, value * result );
    
    void logL( value * events, unsigned int n, value delta_t, value * result );
    
    void event_prob( value * events, unsigned int n, value * result );
    
    void event_logp( value * events, unsigned int n, value * result );
    
    YAML::Node toYAML( bool save_stimulus=true ) const;
    
    void save( std::ostream & stream, bool save_stimulus=true ) const;
    
    void save( std::string path, bool save_stimulus=true ) const;
    
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

class Decoder {
public:
    Decoder( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihoods, const std::vector<value> & prior = {} );
    
    Decoder( std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> & likelihoods, const std::vector<std::vector<value>> & prior = {} );
    
    void decode( std::vector<value*> events, std::vector<unsigned int> nevents, value delta_t, std::vector<value*> result, bool normalize=true );
    
    void decode ( std::vector<value*> events, std::vector<unsigned int> nevents, value delta_t, value* result, unsigned int index=0, bool normalize=true );
    
    void decode ( std::vector<std::vector<value>> events, value delta_t, std::vector<value*> result, bool normalize=true );
    
    void decode ( std::vector<std::vector<value>> events, value delta_t, value* result, unsigned int index=0, bool normalize=true );
    
    
    unsigned int nsources() const;
    bool is_union() const;
    unsigned int n_union() const;
    unsigned int grid_size(unsigned int index=0) const;
    const std::vector<unsigned int> & grid_sizes() const;
    
    std::shared_ptr<PoissonLikelihood> likelihood( unsigned int source, unsigned int index = 0 );
    
    //void add_likelihood( std::shared_ptr<PoissonLikelihood> likelihood );
    //void add_likelihood( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihood );
    //void remove_likelihood( unsigned int source );
    
    unsigned int nenabled_sources() const;
    void enable_source( unsigned int source );
    void enable_all_sources();
    void enable_one_source( unsigned int source );
    void enable_sources( const std::vector<bool> & state );
    void disable_source( unsigned int source );
    const std::vector<bool> & enabled_sources() const;
    
protected:
    std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> likelihoods_;
    std::vector<std::vector<value>> prior_;
    std::vector<unsigned int> grid_sizes_;
    std::vector<bool> likelihood_selection_;
    
};
