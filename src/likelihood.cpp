#include "likelihood.hpp"

StimulusOccupancy::StimulusOccupancy( Space & space, Grid & grid, double stimulus_duration, value compression ) :
stimulus_duration_(stimulus_duration), compression_(compression) {
    
    if (!(space.specification()==grid.specification())) {
        throw std::runtime_error("Grid does not match stimulus space.");
    }
    
    stimulus_distribution_.reset( new Mixture( space, compression ) );
    stimulus_grid_.reset( grid.clone() );

}
    
const Space & StimulusOccupancy::space() const { return stimulus_distribution_->space(); }

const Grid & StimulusOccupancy::grid() const { return *stimulus_grid_; }

value StimulusOccupancy::compression() const { return compression_; }

double StimulusOccupancy::stimulus_duration() const { return stimulus_duration_; }

unsigned int StimulusOccupancy::ndim() const { return stimulus_grid_->ndim(); }

bool StimulusOccupancy::random_insertion() const { return random_insertion_; }
void StimulusOccupancy::set_random_insertion(bool val) { random_insertion_ = val; }
    
void StimulusOccupancy::add_stimulus( const std::vector<value> & stimuli, unsigned repetitions ) {
    
    unsigned int n = stimuli.size() / stimulus_distribution_->space().ndim();
    if (n*stimulus_distribution_->space().ndim() != stimuli.size()) {
        throw std::runtime_error("Not a whole number of samples.");
    }
    
    add_stimulus( stimuli.data(), n, repetitions );
}
    
void StimulusOccupancy::add_stimulus( const value * stimuli, unsigned int n, unsigned int repetitions ) {
    
    lock_.lock();
    stimulus_distribution_->merge_samples( stimuli, n, random_insertion_, static_cast<value>(repetitions) );
    lock_.unlock();
}

value StimulusOccupancy::stimulus_time() {
    
    lock_.lock();
    value t = stimulus_distribution_->sum_of_weights() * stimulus_duration_;
    lock_.unlock();
    return t;
}

void StimulusOccupancy::occupancy( std::vector<value> & out ) {

    out.resize( stimulus_grid_->size() );
    occupancy( out.data() );
    
}
    
void StimulusOccupancy::prob( std::vector<value> & out ) {
    
    out.resize( stimulus_grid_->size() );
    prob( out.data() );
}

void StimulusOccupancy::logp( std::vector<value> & out ) {
    
    out.resize( stimulus_grid_->size() );
    logp( out.data() );
}
    
void StimulusOccupancy::occupancy( value * out ) {
    
    value factor;
    
    lock_.lock();
    stimulus_distribution_->evaluate( *stimulus_grid_, out );
    factor = stimulus_distribution_->sum_of_weights() * stimulus_duration_;
    lock_.unlock();
    
    std::transform( out, out + stimulus_grid_->size(), out, [factor](const value &a) { return a * factor;  } );
    
}
    
void StimulusOccupancy::prob( value * out ) {

    lock_.lock();
    stimulus_distribution_->evaluate( *stimulus_grid_, out );
    lock_.unlock();
}

void StimulusOccupancy::logp( value * out ) {
    
    prob( out );
    std::transform( out, out + stimulus_grid_->size(), out, [](const value &a) { return fastlog(a);  } );
}

YAML::Node StimulusOccupancy::toYAML( ) const {
    
    YAML::Node node;
    
    node["stimulus_duration"] = stimulus_duration_;
    node["compression"] = compression_;
    node["random_insertion"]  = random_insertion_;
    
    node["stimulus_distribution"] = stimulus_distribution_->toYAML();
    node["stimulus_grid"] = stimulus_grid_->toYAML();
    
    return node;
}

void StimulusOccupancy::save( std::ostream & stream ) const {
    
    auto node = toYAML(); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}

void StimulusOccupancy::save( std::string path ) const {
    
    std::ofstream fout(path); save(fout);
}


PoissonLikelihood::PoissonLikelihood( Space & stimulus_space, Grid & grid, double stimulus_duration, value compression ) :
changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    if (!(stimulus_space.specification()==grid.specification())) {
        throw std::runtime_error("Grid does not match stimulus space.");
    }
    
    event_distribution_.reset( new Mixture( stimulus_space, compression ) );
    //stimulus_distribution_.reset( new Mixture( stimulus_space, compression ) );
    stimulus_distribution_.reset( new StimulusOccupancy( stimulus_space, grid, stimulus_duration, compression ) );
    
    stimulus_grid_.reset( grid.clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
    
}

PoissonLikelihood::PoissonLikelihood( Space & event_space, Space & stimulus_space, Grid & grid, double stimulus_duration, value compression ) :
changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    if (!(stimulus_space.specification()==grid.specification())) {
        throw std::runtime_error("Grid does not match stimulus space.");
    }
    
    auto full_space = MultiSpace( { &event_space, &stimulus_space } );
    
    event_distribution_.reset( new Mixture( full_space, compression ) );
    //stimulus_distribution_.reset( new Mixture( stimulus_space, compression ) );
    stimulus_distribution_.reset( new StimulusOccupancy( stimulus_space, grid, stimulus_duration, compression ) );
    //stimulus_grid_.reset( new MultiGrid( { &grid }, {} ) ); // TODO: fix so that wrapping in MultiGrid is unnecessary
    stimulus_grid_.reset( grid.clone() ); // TODO: fix so that wrapping in MultiGrid is unnecessary
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

PoissonLikelihood::PoissonLikelihood( Space & event_space, std::shared_ptr<StimulusOccupancy> stimulus ) :
changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    const Space * ptr = &(stimulus->space());
    
    auto full_space = MultiSpace( {&event_space, ptr} );
    
    event_distribution_.reset( new Mixture( full_space, stimulus->compression() ) );
    
    stimulus_distribution_ = stimulus;
    
    stimulus_grid_.reset( stimulus->grid().clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

PoissonLikelihood::PoissonLikelihood( std::shared_ptr<StimulusOccupancy> stimulus ) :
changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    event_distribution_.reset( new Mixture( stimulus->space(), stimulus->compression() ) );
    
    stimulus_distribution_ = stimulus;
    
    stimulus_grid_.reset( stimulus->grid().clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

bool PoissonLikelihood::changed() const { return changed_; }
bool PoissonLikelihood::random_insertion() const { return random_insertion_; }
void PoissonLikelihood::set_random_insertion(bool val) { random_insertion_ = val; }

//value PoissonLikelihood::rate_offset() const { return rate_offset_; }
//void PoissonLikelihood::set_rate_offset(value val) { rate_offset_ = val; changed_=true; }

value PoissonLikelihood::rate_scale() const { return rate_scale_; }
void PoissonLikelihood::set_rate_scale(value val) { rate_scale_ = val; }

unsigned int PoissonLikelihood::ndim() const { return event_distribution_->space().ndim(); }
unsigned int PoissonLikelihood::ndim_stimulus() const { return stimulus_grid_->ndim(); }
unsigned int PoissonLikelihood::ndim_events() const { return ndim()-ndim_stimulus(); }

const Grid & PoissonLikelihood::grid() const { return *stimulus_grid_; }
const Mixture & PoissonLikelihood::event_distribution() const { return *event_distribution_; }

std::shared_ptr<StimulusOccupancy> PoissonLikelihood::stimulus() { return stimulus_distribution_; }

void PoissonLikelihood::add_events( const std::vector<value> & events, unsigned int repetitions ) {
    
    unsigned int n = events.size() / event_distribution_->space().ndim();
    
    if (n*event_distribution_->space().ndim() != events.size()) {
        throw std::runtime_error("Not a whole number of samples,");
    }
    
    add_events( events.data(), n, repetitions );
}

void PoissonLikelihood::add_events( const value * events, unsigned int n, unsigned int repetitions ) {
    
    if (repetitions==0) { return; }
    
    event_distribution_->merge_samples( events, n, random_insertion_, static_cast<value>(repetitions) );
    
    changed_ = true;
}

value PoissonLikelihood::mu() const {
    return event_distribution_->sum_of_weights() / stimulus_distribution_->stimulus_time();
}

const std::vector<value> & PoissonLikelihood::stimulus_logp() const { 
    return logp_stimulus_;
}

const std::vector<value> & PoissonLikelihood::event_rate() const { 
    return event_rate_;
}

//const std::vector<value> & PoissonLikelihood::offset() const { 
//    return offset_;
//}

void PoissonLikelihood::precompute() { 
    
    stimulus_distribution_->prob( logp_stimulus_ );
    
    //if (rate_offset_>0) {
    //    offset_ = logp_stimulus_;
    //    value factor = rate_offset_ / (rate_scale_*mu());
    //    std::transform( offset_.begin(), offset_.end(), offset_.begin(), [factor](const value & a) { return factor*a; } );
    //}
    
    p_event_.reset( new PartialMixture( event_distribution_.get(), *stimulus_grid_ ) );
    
    event_rate_.resize( stimulus_grid_->size() );
    
    p_event_->marginal( event_rate_.data() );
    
    std::transform( event_rate_.begin(), event_rate_.end(), logp_stimulus_.begin(), event_rate_.begin(), [](const value & a, const value & b) {return a/b;} );
    std::transform( logp_stimulus_.begin(), logp_stimulus_.end(), logp_stimulus_.begin(), [](const value & a) { return fastlog(a); } );
    
    changed_ = false;
    
}

void PoissonLikelihood::likelihood( value * events, unsigned int n, value delta_t, value * result ) {
    
    logL( events, n, delta_t, result );
    std::transform( result, result + stimulus_grid_->size(), result, [](const value & a) { return fastexp(a); } );
}

void PoissonLikelihood::logL( value * events, unsigned int n, value delta_t, value * result ) {
    
    if (changed_) { precompute(); }
    
    event_logp( events, n, result );
    
    value constant =  n*fastlog(delta_t*rate_scale_*mu());
    std::transform( result, result + stimulus_grid_->size(), result, [constant](const value & a) { return a + constant; } );
    
    // subtract n*log(p_stimulus_)
    std::transform( result, result + stimulus_grid_->size(), logp_stimulus_.begin(), result, [n](const value & a, const value & b) { return a - n*b; } );
    
    // subtract delta_t * p_event_stimulus_/p_stimulus_
    constant = delta_t*rate_scale_*mu();
    //value offset = rate_offset_/(rate_scale_*mu());
    std::transform( result, result + stimulus_grid_->size(), event_rate_.begin(), result, [constant](const value & a, const value & b) { return a - constant*b; } );
    
}

void PoissonLikelihood::event_prob( value * events, unsigned int n, value * result ) {
    
    event_logp( events, n, result );
    std::transform( result, result + stimulus_grid_->size(), result, [](const value & a) { return fastexp(a); } );
}

void PoissonLikelihood::event_logp( value * events, unsigned int n, value * result ) {
    
    //if (rate_offset_>0) {
    //    p_event_->complete_multi( events, n, result, offset_.data() );
    //} else {
    p_event_->complete_multi( events, n, result );
    //}
}

YAML::Node PoissonLikelihood::toYAML( bool save_stimulus ) const {
    
    YAML::Node node;
    
    node["rate_scale"] = rate_scale_;
    node["random_insertion"]  = random_insertion_;
    
    node["event_distribution"] = event_distribution_->toYAML();
    
    if (save_stimulus) {
        node["stimulus_distribution"] = stimulus_distribution_->toYAML();
        node["stimulus_grid"] = stimulus_grid_->toYAML();
    }
    
    return node;
}

void PoissonLikelihood::save( std::ostream & stream, bool save_stimulus ) const {
    
    auto node = toYAML(save_stimulus); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}

void PoissonLikelihood::save( std::string path, bool save_stimulus ) const {
    
    std::ofstream fout(path); save(fout, save_stimulus);
}



Decoder::Decoder( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihoods, const std::vector<value> & prior ) :
prior_( {prior} ) {
    
    unsigned int nsources = likelihoods.size();
    if (nsources==0) {
        throw std::runtime_error("Please provide at least one source.");
    }
    
    // check if all sources have same stimulus grid space/size
    for (unsigned int source=1; source<nsources; ++source) {
        if ( ! (likelihoods[0]->grid() == likelihoods[source]->grid()) ) {
            throw std::runtime_error("All sources need to have the same stimulus grid shape and space.");
        }
    }
    
    // check if prior.size()==0 || prior.size()==grid.size()
    if ( prior.size()!=0 && prior.size()!=likelihoods[0]->grid().size() ) {
        throw std::runtime_error("Prior does not have correct number of elements.");
    }
    
    likelihoods_.resize( likelihoods.size() );
    likelihood_selection_.assign( likelihoods.size(), true );
    
    for ( unsigned int k =0 ; k<likelihoods.size(); ++k ) {
        likelihoods_[k].push_back( likelihoods[k] );
    }
    
    grid_sizes_.push_back( likelihoods_[0][0]->grid().size() );
    
}
    
Decoder::Decoder( std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> & likelihoods, const std::vector<std::vector<value>> & prior ) :
likelihoods_(likelihoods), prior_(prior) {
    
    unsigned int nsources = likelihoods.size();
    if (nsources==0) {
        throw std::runtime_error("Please provide at least one source.");
    }
    
    unsigned int nunion = likelihoods[0].size();
    for (unsigned int source=1; source<nsources; ++source) {
        if (likelihoods[source].size()!=nunion) {
            throw std::runtime_error("All sources need to have the same number of likelihoods.");
        }
        for (unsigned int index=0; index<nunion; ++index) {
            if ( ! (likelihoods[source][index]->grid()==likelihoods[0][index]->grid() ) ) {
                throw std::runtime_error("Union likelihoods across sources need to have the same grid size and space.");
            }
        }
    }
    
    for (unsigned int source=0; source<nsources; ++source) {
        
        for (unsigned int index = 1; index<nunion; ++index) {
            // check if likelihoods[source][index] has the same event space as likelihoods[source][0]
        }
        
    }
    
    // test if prior[index].size()==0 || prior[index].size()==likelihoods[0][index]->grid().size()
    for (unsigned int index=0; index<nunion; ++index) {
        grid_sizes_.push_back( likelihoods[0][index]->grid().size() );
        if (prior[index].size()!=0 && prior[index].size()!=grid_sizes_.back()) {
            throw std::runtime_error("Prior does not have correct number of elements.");
        }
    }
    
    likelihood_selection_.assign( nsources, true );
    
}

void Decoder::decode( std::vector<value*> events, std::vector<unsigned int> nevents, value delta_t, std::vector<value*> result, bool normalize ) {
    
    // check events and result vectors
    if ( events.size() != nsources() || nevents.size() != nsources() ) {
        std::runtime_error("Incorrect number of sources.");
    }
    
    if (result.size()!=n_union()) {
        std::runtime_error("Incorrect number of outputs.");
    }
    
    unsigned int n;
    for (unsigned int source=0; source<nsources(); ++source) {
        
        if (!likelihood_selection_[source]) {continue;}
        
        // determine number of spikes for source
        n = nevents[source]/likelihoods_[source][0]->ndim_events();
        if ( n * likelihoods_[source][0]->ndim_events() != nevents[source] ) {
            throw std::runtime_error("Incomplete samples.");
        }
        
        for (unsigned int index=0; index<n_union(); ++index) {
            likelihoods_[source][index]->logL( events[source], n, delta_t, result[index] );
        }
    }
    
    // add log prior
    for (unsigned int index=0; index<n_union(); ++index) {
        if (prior_[index].size()>0) {
            std::transform( result[index], result[index]+prior_[index].size(), prior_[index].data(), result[index], std::plus<value>() );
        }
    }
    
    
    // normalize across union
    if (normalize) {
    
        // find maximum across union
        value max = *std::max_element( result[0], result[0] + grid_sizes_[0] );
        for (unsigned int index=1; index<n_union(); ++index) {
            max = std::max( max, *std::max_element( result[index], result[index] + grid_sizes_[index] ) );
        }
        
        // compute exp( x - max )
        for (unsigned int index=0; index<n_union(); ++index) {
            std::transform( result[index], result[index] + grid_sizes_[index], result[index], [max]( const value & a ) { return fastexp( a - max ); } );
        }
        
        // compute sum across union
        value sum = 0.;
        for (unsigned int index=0; index<n_union(); ++index) {
            sum += std::accumulate( result[index], result[index] + grid_sizes_[index], 0. );
        }
        
        // divide by sum
        for (unsigned int index=0; index<n_union(); ++index) {
            std::transform( result[index], result[index] + grid_sizes_[index], result[index], [sum]( const value & a ) { return a/sum; } );
        }
    
    }
    
}

void Decoder::decode( std::vector<value*> events, std::vector<unsigned int> nevents, value delta_t, value* result, unsigned int index, bool normalize ) {

    if ( events.size() != nsources() || nevents.size() != nsources() ) {
        std::runtime_error("Incorrect number of sources.");
    }
    
    // only for selected part of union
    if (index>=n_union()) {
        throw std::runtime_error("Union index out of bounds.");
    }
    
    unsigned int n;
    
    for (unsigned int source=0; source<nsources(); ++source) {
        
        if (!likelihood_selection_[source]) {continue;}
        
        n = nevents[source]/likelihoods_[source][0]->ndim_events();
        if ( n * likelihoods_[source][0]->ndim_events() != nevents[source] ) {
            throw std::runtime_error("Incomplete samples.");
        }
        
        // sum log likelihoods
        likelihoods_[source][index]->logL( events[source], n, delta_t, result );
    }
    
    // add log prior
    if (prior_[index].size()>0) {
        std::transform( result, result+prior_[index].size(), prior_[index].data(), result, std::plus<value>() );
    }
    
    // normalize
    if (normalize) {
        
        // find maximum
        value max = *std::max_element( result, result + grid_sizes_[0] );
        // compute exp( x - max )
        std::transform( result, result + grid_sizes_[0], result, [max]( const value & a ) { return fastexp( a - max ); } );
        // compute sum
        value sum = std::accumulate( result, result + grid_sizes_[0], 0. );
        // divide by sum
        std::transform( result, result + grid_sizes_[0], result, [sum]( const value & a ) { return a/sum; } );
    }
}

void Decoder::decode ( std::vector<std::vector<value>> events, value delta_t, std::vector<value*> result, bool normalize ) {
    
    std::vector<value*> events_ptr;
    std::vector<unsigned int> events_n;
    
    for (unsigned int k=0; k<events.size(); ++k) {
        events_ptr.push_back( events[k].data() );
        events_n.push_back( events[k].size() );
    }
    
    return decode( events_ptr, events_n, delta_t, result, normalize );
    
}
    
    //// check events and result vectors
    //if ( events.size() != nsources() ) {
        //std::runtime_error("Incorrect number of sources.");
    //}
    
    //if (result.size()!=n_union()) {
        //std::runtime_error("Incorrect number of outputs.");
    //}
    
    //unsigned int n;
    //for (unsigned int source=0; source<nsources(); ++source) {
        
        //// determine number of spikes for source
        //n = events[source].size()/likelihoods_[source][0]->ndim_events();
        //if ( n * likelihoods_[source][0]->ndim_events() != events[source].size() ) {
            //throw std::runtime_error("Incomplete samples.");
        //}
        
        //for (unsigned int index=0; index<n_union(); ++index) {
            //likelihoods_[source][index]->logL( events[source].data(), n, delta_t, result[index] );
        //}
    //}
    
    //// add log prior
    //for (unsigned int index=0; index<n_union(); ++index) {
        //if (prior_[index].size()>0) {
            //std::transform( result[index], result[index]+prior_[index].size(), prior_[index].data(), result[index], std::plus<value>() );
        //}
    //}
    
    
    //// normalize across union
    //{
    
        //// find maximum across union
        //value max = *std::max_element( result[0], result[0] + grid_sizes_[0] );
        //for (unsigned int index=1; index<n_union(); ++index) {
            //max = std::max( max, *std::max_element( result[index], result[index] + grid_sizes_[index] ) );
        //}
        
        //// compute exp( x - max )
        //for (unsigned int index=0; index<n_union(); ++index) {
            //std::transform( result[index], result[index] + grid_sizes_[index], result[index], [max]( const value & a ) { return fastexp( a - max ); } );
        //}
        
        //// compute sum across union
        //value sum = 0.;
        //for (unsigned int index=0; index<n_union(); ++index) {
            //sum += std::accumulate( result[index], result[index] + grid_sizes_[index], 0. );
        //}
        
        //// divide by sum
        //for (unsigned int index=0; index<n_union(); ++index) {
            //std::transform( result[index], result[index] + grid_sizes_[index], result[index], [sum]( const value & a ) { return a/sum; } );
        //}
    
    //}
//}

void Decoder::decode ( std::vector<std::vector<value>> events, value delta_t, value* result, unsigned int index, bool normalize ) {
    
    std::vector<value*> events_ptr;
    std::vector<unsigned int> events_n;
    
    for (unsigned int k=0; k<events.size(); ++k) {
        events_ptr.push_back( events[k].data() );
        events_n.push_back( events[k].size() );
    }
    
    return decode( events_ptr, events_n, delta_t, result, index, normalize );
    
}
    
    //if ( events.size() != nsources() ) {
        //std::runtime_error("Incorrect number of sources.");
    //}
    
    //// only for selected part of union
    //if (index>=n_union()) {
        //throw std::runtime_error("Union index out of bounds.");
    //}
    
    //unsigned int n;
    
    //for (unsigned int source=0; source<nsources(); ++source) {
        
        //n = events[source].size()/likelihoods_[source][0]->ndim_events();
        //if ( n * likelihoods_[source][0]->ndim_events() != events[source].size() ) {
            //throw std::runtime_error("Incomplete samples.");
        //}
        
        //// sum log likelihoods
        //likelihoods_[source][index]->logL( events[source].data(), n, delta_t, result );
    //}
    
    //// add log prior
    //if (prior_[index].size()>0) {
        //std::transform( result, result+prior_[index].size(), prior_[index].data(), result, std::plus<value>() );
    //}
    
    //// normalize
    //{
        //// find maximum
        //value max = *std::max_element( result, result + grid_sizes_[0] );
        //// compute exp( x - max )
        //std::transform( result, result + grid_sizes_[0], result, [max]( const value & a ) { return fastexp( a - max ); } );
        //// compute sum
        //value sum = std::accumulate( result, result + grid_sizes_[0], 0. );
        //// divide by sum
        //std::transform( result, result + grid_sizes_[0], result, [sum]( const value & a ) { return a/sum; } );
    //}
    
//}

unsigned int Decoder::nsources() const { return likelihoods_.size(); }

unsigned int Decoder::nenabled_sources() const { return std::count(likelihood_selection_.begin(), likelihood_selection_.end(), true ); }

bool Decoder::is_union() const { return (!likelihoods_.empty()) && likelihoods_[0].size()>1; }

unsigned int Decoder::n_union() const {
    
    if (likelihoods_.empty()) {
        return 0.;
    } else {
        return likelihoods_[0].size();
    }
}

unsigned int Decoder::grid_size(unsigned int index) const { return grid_sizes_[index]; }

const std::vector<unsigned int> & Decoder::grid_sizes() const { return grid_sizes_; }

std::shared_ptr<PoissonLikelihood> Decoder::likelihood( unsigned int source, unsigned int index ) {
   
    if (source>=nsources() || index>=n_union()) {
        throw std::runtime_error("Source and/or union index out of bounds.");
    }
    
    return likelihoods_[source][index];
}

const std::vector<bool> & Decoder::enabled_sources() const {
    return likelihood_selection_;
}

void Decoder::enable_source( unsigned int source ) {
    
    if (source>=nsources()) {
        throw std::runtime_error("Likelihood index out of range.");
    }
    
    likelihood_selection_[source] = true;
}

void Decoder::enable_all_sources() {
    
    likelihood_selection_.assign( nsources(), true );
}

void Decoder::enable_one_source( unsigned int source ) {
    
    if (source>=nsources()) {
        throw std::runtime_error("Likelihood index out of range.");
    }
    
    likelihood_selection_.assign( nsources(), false );
    likelihood_selection_[source] = true;
}

void Decoder::disable_source( unsigned int source ) {
    
    if (source>=nsources()) {
        throw std::runtime_error("Likelihood index out of range.");
    }
    
    likelihood_selection_[source] = false;
}

void Decoder::enable_sources( const std::vector<bool> & state ) {
    
    if (state.size()!=nsources()) {
        throw std::runtime_error("Invalid vector size, it does not match number of likelihoods.");
    }
    
    likelihood_selection_ = state;
    
}

