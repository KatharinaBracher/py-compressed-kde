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
#include "decoder.hpp"

// constructors
Decoder::Decoder( std::vector<std::shared_ptr<PoissonLikelihood>> & likelihoods, 
    const std::vector<value> & prior )
     : prior_( {prior} ) {
    
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
    grid_shapes_.push_back( likelihoods_[0][0]->grid().shape() );
    
}
    
Decoder::Decoder( std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> & likelihoods, 
    const std::vector<std::vector<value>> & prior )
    : likelihoods_(likelihoods), prior_(prior) {
    
    unsigned int nsources = likelihoods.size();
    if (nsources==0) {
        throw std::runtime_error("Please provide at least one source.");
    }
    
    unsigned int nunion = likelihoods[0].size();
    for (unsigned int source=1; source<nsources; ++source) {
        if (likelihoods[source].size()!=nunion) {
            throw std::runtime_error("All sources need to have the same "
                "number of likelihoods.");
        }
        for (unsigned int index=0; index<nunion; ++index) {
            if ( ! (likelihoods[source][index]->grid()==likelihoods[0][index]->grid() ) ) {
                throw std::runtime_error("Union likelihoods across sources need "
                    "to have the same grid size and space.");
            }

            // to do: check if likelihoods[source][index] has the same event space as likelihoods[source][0]
        }
    }
    
    // test if prior[index].size()==0 || prior[index].size()==likelihoods[0][index]->grid().size()
    for (unsigned int index=0; index<nunion; ++index) {
        grid_sizes_.push_back( likelihoods[0][index]->grid().size() );
        grid_shapes_.push_back( likelihoods[0][index]->grid().shape() );
        if (prior[index].size()!=0 && prior[index].size()!=grid_sizes_.back()) {
            throw std::runtime_error("Prior does not have correct number of elements.");
        }
    }
    
    likelihood_selection_.assign( nsources, true );
}
void compute_posterior(std::vector<value *> result,
                       std::vector<std::vector<value>> prior,
                       std::vector<unsigned int> grid_sizes, bool normalize){

    // add log prior
    for (unsigned int index=0; index<result.size(); ++index) {
        if (prior[index].size()>0) {
            std::transform( result[index], result[index]+prior[index].size(), prior[index].data(), result[index], std::plus<value>() );
        }
    }

    // normalize across union
    if (normalize) {

        // find maximum across union
        value max = *std::max_element( result[0], result[0] + grid_sizes[0] );
        for (unsigned int index=1; index<result.size(); ++index) {
            max = std::max( max, *std::max_element( result[index], result[index] + grid_sizes[index] ) );
        }

        // compute exp( x - max )
        for (unsigned int index=0; index<result.size(); ++index) {
            std::transform( result[index], result[index] + grid_sizes[index], result[index], [max]( const value & a ) { return fastexp( a - max ); } );
        }

        // compute sum across union
        value sum = 0.;
        for (unsigned int index=0; index<result.size(); ++index) {
            sum += std::accumulate( result[index], result[index] + grid_sizes[index], 0. );
        }

        // divide by sum
        for (unsigned int index=0; index<result.size(); ++index) {
            std::transform( result[index], result[index] + grid_sizes[index], result[index], [sum]( const value & a ) { return a/sum; } );
        }

    }

}

// decoding methods
void Decoder::decode( std::vector<value*> events, std::vector<unsigned int> nevents, 
    value delta_t, std::vector<value*> result, bool normalize ) {
    
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
    
    compute_posterior(result, prior_, grid_sizes_, normalize);

}

void compute_posterior(value * result, std::vector<value> prior, unsigned int grid_size, bool normalize){
    // add log prior
    if (prior.size()>0) {
        std::transform( result, result+prior.size(), prior.data(), result, std::plus<value>() );
    }

    // normalize
    if (normalize) {

        // find maximum
        value max = *std::max_element( result, result + grid_size );
        // compute exp( x - max )
        std::transform( result, result + grid_size, result, [max]( const value & a ) { return fastexp( a - max ); } );
        // compute sum
        value sum = std::accumulate( result, result + grid_size, 0. );
        // divide by sum
        std::transform( result, result + grid_size, result, [sum]( const value & a ) { return a/sum; } );
    }
}
void Decoder::decode ( std::vector<std::vector<value>> events, value delta_t,
    std::vector<value*> result, bool normalize ) {

    std::vector<value*> events_ptr;
    std::vector<unsigned int> events_n;

    for (unsigned int k=0; k<events.size(); ++k) {
        events_ptr.push_back( events[k].data() );
        events_n.push_back( events[k].size() );
    }

    return decode( events_ptr, events_n, delta_t, result, normalize );

}

void Decoder::decode( std::vector<value*> events, std::vector<unsigned int> nevents, 
    value delta_t, value* result, unsigned int index, bool normalize ) {

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
    
    compute_posterior(result, prior_[index], grid_sizes_[index], normalize);
}


    
void Decoder::decode ( std::vector<std::vector<value>> events, value delta_t, 
    value* result, unsigned int index, bool normalize ) {
    
    std::vector<value*> events_ptr;
    std::vector<unsigned int> events_n;
    
    for (unsigned int k=0; k<events.size(); ++k) {
        events_ptr.push_back( events[k].data() );
        events_n.push_back( events[k].size() );
    }
    
    return decode( events_ptr, events_n, delta_t, result, index, normalize );
    
}
    
// properties
unsigned int Decoder::nsources() const { return likelihoods_.size(); }

unsigned int Decoder::nenabled_sources() const { 
    return std::count(likelihood_selection_.begin(),
        likelihood_selection_.end(), true );
}

bool Decoder::is_union() const {
    return (!likelihoods_.empty()) && likelihoods_[0].size()>1;
}

unsigned int Decoder::n_union() const {
    
    if (likelihoods_.empty()) {
        return 0.;
    } else {
        return likelihoods_[0].size();
    }
}

unsigned int Decoder::grid_size(unsigned int index) const {
    return grid_sizes_[index];
}

const std::vector<unsigned int> & Decoder::grid_sizes() const {
    return grid_sizes_;
}

std::vector<long unsigned int> Decoder::grid_shape(unsigned int index) const {
    return grid_shapes_[index];
}

const std::vector<std::vector<long unsigned int>> & Decoder::grid_shapes() const {
    return grid_shapes_;
}

const Grid & Decoder::grid(unsigned int index) const {
    
    if (nsources()==0) {
        throw std::runtime_error("No likelihoods.");
    }
    
    if (index>=n_union()) {
        throw std::runtime_error("Index out of bounds.");
    }
    
    return likelihoods_[0][index]->grid();
}

std::shared_ptr<StimulusOccupancy> Decoder::stimulus(unsigned int index) {
    if (nsources()==0) {
        throw std::runtime_error("No likelihoods.");
    }
    
    if (index>=n_union()) {
        throw std::runtime_error("Index out of bounds.");
    }
    
    return likelihoods_[0][index]->stimulus();
}

// likelihood getter
std::shared_ptr<PoissonLikelihood> Decoder::likelihood(
    unsigned int source, unsigned int index) {
   
    if (source>=nsources() || index>=n_union()) {
        throw std::runtime_error("Source and/or union index out of bounds.");
    }
    
    return likelihoods_[source][index];
}

// enable/disable sources
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

// hdf5
void Decoder::to_hdf5(HighFive::Group & group) const {
    
    HighFive::DataSet ds_nsources = group.createDataSet<unsigned int>(
        "nsources", HighFive::DataSpace::From(nsources()));
    ds_nsources.write(nsources());
    
    HighFive::DataSet ds_nunion = group.createDataSet<unsigned int>(
        "nunion", HighFive::DataSpace::From(n_union()));
    ds_nunion.write(n_union());
    
    // save likelihoods
    // for each union element, go over all sources
    // check if the stimulus was already saved
    // if yes: save path to stimulus
    // if no: save stimulus, save path to stimulus
    
    HighFive::Group grp_stim = group.createGroup("stimulus");
    HighFive::Group grp_like = group.createGroup("likelihood");
    
    std::map<std::string, std::shared_ptr<StimulusOccupancy>> stim_map;
    
    bool found;
    std::string key;
    
    for (unsigned int u=0; u<n_union(); ++u) {
        
        for (unsigned int s=0; s<nsources(); ++s) {
            
            key = "stimulus_" + std::to_string(s) + "_" + std::to_string(u);
            found = false;
            
            // is stimulus already saved?
            for (auto & it : stim_map) {
                if (it.second == likelihoods_[s][u]->stimulus()) {
                    key = it.first;
                    found = true;
                    break;
                }
            }
            
            // if no: save stimulus and add to map
            if (!found) {
                stim_map[key] = likelihoods_[s][u]->stimulus();
                
                HighFive::Group g1 = grp_stim.createGroup(key);
                likelihoods_[s][u]->stimulus()->to_hdf5(g1);
            }
            
            // save likelihood (without stimulus)
            HighFive::Group g2 = grp_like.createGroup("likelihood_" + std::to_string(s) + "_" + std::to_string(u));
            likelihoods_[s][u]->to_hdf5(g2, false);
            
            HighFive::Attribute attr = g2.createAttribute<std::string>("stimulus", HighFive::DataSpace::From(key));
            attr.write(key);
        }
    }
    
    // save priors
    HighFive::Group priors = group.createGroup("priors");
    
    unsigned int n = 0;
    
    for (auto & k : prior_) {
        HighFive::DataSet ds_prior = priors.createDataSet<value>(
            "prior" + std::to_string(n), HighFive::DataSpace::From(k));
        
        ds_prior.write(k);
        
        ++n;
    }
    
    // save likelihood selection
    std::vector<unsigned char> tmp(likelihood_selection_.begin(), 
        likelihood_selection_.end());
    
    HighFive::DataSet ds_sel = group.createDataSet<unsigned char>(
        "selection", HighFive::DataSpace::From(tmp));
    ds_sel.write(tmp);
}

std::unique_ptr<Decoder> Decoder::from_hdf5(const HighFive::Group & group) {
    
    unsigned int nsources;
    unsigned int nunion;
    std::string key, stim_key;
    
    group.getDataSet("nsources").read(nsources);
    group.getDataSet("nunion").read(nunion);
    
    // load map of stimuli
    std::map<std::string, std::shared_ptr<StimulusOccupancy>> stim_map;
    
    HighFive::Group grp_stim = group.getGroup("stimulus");
    
    for (unsigned int u=0; u<nunion; ++u) {
        
        for (unsigned int s=0; s<nsources; ++s) {
            
            key = "stimulus_" + std::to_string(s) + "_" + std::to_string(u);
            
            if (!grp_stim.exist(key)) {
                continue;
            }
            
            std::shared_ptr<StimulusOccupancy> stim = StimulusOccupancy::from_hdf5(grp_stim.getGroup(key));
            
            stim_map[key] = stim;
        }
    }
            
    // for each union element, go over all sources
    // load path to stimulus
    // load likelihood (and supply corresponding stimulus)
    std::vector<std::vector<std::shared_ptr<PoissonLikelihood>>> likelihoods(nsources);
    
    HighFive::Group grp_like = group.getGroup("likelihood");
    
    for (unsigned int u=0; u<nunion; ++u) {
        
        for (unsigned int s=0; s<nsources; ++s) {
            
            key = "likelihood_" + std::to_string(s) + "_" + std::to_string(u);
            
            HighFive::Attribute attr = grp_like.getGroup(key).getAttribute("stimulus");
            attr.read(stim_key);
            
            std::shared_ptr<PoissonLikelihood> L= PoissonLikelihood::from_hdf5(grp_like.getGroup(key), stim_map[stim_key]);
            
            likelihoods[s].push_back(L);
            
        }
    }
    
    // load priors
    std::vector<std::vector<value>> priors(nunion);
    
    HighFive::Group grp_priors = group.getGroup("priors");
    
    for (unsigned int k=0; k<nunion; ++k) {
        grp_priors.getDataSet("prior" + std::to_string(k)).read(priors[k]);
    }
    
    // create Decoder object
    auto dec = std::make_unique<Decoder>(likelihoods, priors);
    
    // load likelihood selection
    std::vector<unsigned char> tmp;
    HighFive::DataSet ds_sel = group.getDataSet("selection");
    ds_sel.read(tmp);
    
    std::vector<bool> selection(tmp.begin(), tmp.end());
    
    // set likelihood selection
    dec->enable_sources(selection);
    
    return dec;
}

void Decoder::save_to_hdf5(std::string filename, int flags, std::string path) {
    HighFive::File file(filename, flags);
    
    // to do: create attributes for version, etc.
    
    if (!path.empty()) {
        HighFive::Group group = file.createGroup(path);
        this->to_hdf5(group);
    } else {
        HighFive::Group group = file.getGroup("/");
        this->to_hdf5(group);
    }
    
    file.flush();
}
    
std::unique_ptr<Decoder> Decoder::load_from_hdf5(std::string filename,
    std::string path) {
    
    HighFive::File file(filename, HighFive::File::ReadOnly);
    
    if (!path.empty()) {
        HighFive::Group group = file.getGroup(path);
        return Decoder::from_hdf5(group);
    } else {
        HighFive::Group group = file.getGroup("/");
        return Decoder::from_hdf5(group);
    }    
}
