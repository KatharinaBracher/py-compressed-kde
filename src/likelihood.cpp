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
#include "likelihood.hpp"

// default constructor
PoissonLikelihood::PoissonLikelihood():
changed_(true), random_insertion_(true), rate_scale_(1.) {}

// constructors
PoissonLikelihood::PoissonLikelihood( Space & stimulus_space, Grid & grid, 
    double stimulus_duration, value compression )
    : changed_(true), random_insertion_(true), rate_scale_(1.) {
    
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

PoissonLikelihood::PoissonLikelihood( Space & event_space, Space & stimulus_space, 
    Grid & grid, double stimulus_duration, value compression )
    : changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    if (!(stimulus_space.specification()==grid.specification())) {
        throw std::runtime_error("Grid does not match stimulus space.");
    }
    
    auto full_space = MultiSpace( { &event_space, &stimulus_space } );
    
    event_distribution_.reset( new Mixture( full_space, compression ) );
    stimulus_distribution_.reset( new StimulusOccupancy( stimulus_space, grid, 
        stimulus_duration, compression ) );
    stimulus_grid_.reset( grid.clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

PoissonLikelihood::PoissonLikelihood( Space & event_space, 
    std::shared_ptr<StimulusOccupancy> stimulus )
    : changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    const Space * ptr = &(stimulus->space());
    
    auto full_space = MultiSpace( {&event_space, ptr} );
    
    event_distribution_.reset( new Mixture( full_space, stimulus->compression() ) );
    
    stimulus_distribution_ = stimulus;
    
    stimulus_grid_.reset( stimulus->grid().clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

PoissonLikelihood::PoissonLikelihood( std::shared_ptr<StimulusOccupancy> stimulus )
    : changed_(true), random_insertion_(true), rate_scale_(1.) {
    
    event_distribution_.reset( new Mixture( stimulus->space(), stimulus->compression() ) );
    
    stimulus_distribution_ = stimulus;
    
    stimulus_grid_.reset( stimulus->grid().clone() );
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    event_rate_.assign( stimulus_grid_->size(), 0. );
}

// properties
bool PoissonLikelihood::changed() const { return changed_; }
bool PoissonLikelihood::random_insertion() const { return random_insertion_; }
void PoissonLikelihood::set_random_insertion(bool val) { random_insertion_ = val; }

//value PoissonLikelihood::rate_offset() const { return rate_offset_; }
//void PoissonLikelihood::set_rate_offset(value val) { rate_offset_ = val; changed_=true; }

value PoissonLikelihood::rate_scale() const { return rate_scale_; }
void PoissonLikelihood::set_rate_scale(value val) { rate_scale_ = val; }

unsigned int PoissonLikelihood::ndim() const { 
    return event_distribution_->space().ndim();
}
unsigned int PoissonLikelihood::ndim_stimulus() const { 
    return stimulus_grid_->ndim();
}
unsigned int PoissonLikelihood::ndim_events() const { 
    return ndim()-ndim_stimulus();
}

const Grid & PoissonLikelihood::grid() const { return *stimulus_grid_; }
const Mixture & PoissonLikelihood::event_distribution() const { 
    return *event_distribution_;
}

std::shared_ptr<StimulusOccupancy> PoissonLikelihood::stimulus() { 
    return stimulus_distribution_;
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

// methods
void PoissonLikelihood::add_events( const std::vector<value> & events, 
    unsigned int repetitions ) {
    
    unsigned int n = events.size() / event_distribution_->space().ndim();
    
    if (n*event_distribution_->space().ndim() != events.size()) {
        throw std::runtime_error("Not a whole number of samples,");
    }
    
    add_events( events.data(), n, repetitions );
}

void PoissonLikelihood::add_events( const value * events, unsigned int n, 
    unsigned int repetitions ) {
    
    if (repetitions==0) { return; }
    
    event_distribution_->merge_samples( events, n, random_insertion_, 
        static_cast<value>(repetitions) );
    
    changed_ = true;
}

void PoissonLikelihood::precompute() { 
    
    logp_stimulus_.assign( stimulus_grid_->size(), 0. );
    stimulus_distribution_->prob( logp_stimulus_ );
    
    //if (rate_offset_>0) {
    //    offset_ = logp_stimulus_;
    //    value factor = rate_offset_ / (rate_scale_*mu());
    //    std::transform( offset_.begin(), offset_.end(), offset_.begin(), [factor](const value & a) { return factor*a; } );
    //}
    
    p_event_.reset( new PartialMixture( event_distribution_.get(), *stimulus_grid_ ) );
    
    event_rate_.assign(stimulus_grid_->size(), 0.);
    p_event_->marginal( event_rate_.data() );
    
    std::transform( event_rate_.begin(), event_rate_.end(), logp_stimulus_.begin(), event_rate_.begin(), [](const value & a, const value & b) {return a/b;} );
    std::transform( logp_stimulus_.begin(), logp_stimulus_.end(), logp_stimulus_.begin(), [](const value & a) { return fastlog(a); } );
    
    changed_ = false;
    
}

void PoissonLikelihood::likelihood( value * events, unsigned int n, value delta_t, 
    value * result ) {
    
    logL( events, n, delta_t, result );
    std::transform( result, result + stimulus_grid_->size(), result, 
                    [](const value & a) { return fastexp(a); } );
}

void PoissonLikelihood::logL( value * events, unsigned int n, value delta_t,
    value * result ) {
    
    if (changed_) { precompute(); }
    
    event_logp( events, n, result );
    
    value constant =  n*fastlog(delta_t*rate_scale_*mu());
    std::transform( result, result + stimulus_grid_->size(), result, 
        [constant](const value & a) { return a + constant; } );
    
    // subtract n*log(p_stimulus_)
    std::transform( result, result + stimulus_grid_->size(), logp_stimulus_.begin(), 
        result, [n](const value & a, const value & b) { return a - n*b; } );
    
    // subtract delta_t * p_event_stimulus_/p_stimulus_
    constant = delta_t*rate_scale_*mu();
    //value offset = rate_offset_/(rate_scale_*mu());
    std::transform( result, result + stimulus_grid_->size(), event_rate_.begin(), 
        result, [constant](const value & a, const value & b) { return a - constant*b; } );
    
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


// yaml
YAML::Node PoissonLikelihood::to_yaml( bool save_stimulus ) const {
    
    YAML::Node node;
    
    node["rate_scale"] = rate_scale_;
    node["random_insertion"]  = random_insertion_;
    
    node["event_distribution"] = event_distribution_->to_yaml();
    
    if (save_stimulus) {
        node["stimulus_distribution"] = stimulus_distribution_->to_yaml();
    }
    
    return node;
}

void PoissonLikelihood::save_to_yaml( std::ostream & stream, bool save_stimulus ) const {
    
    auto node = to_yaml(save_stimulus); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}

void PoissonLikelihood::save_to_yaml( std::string path, bool save_stimulus ) const {
    
    std::ofstream fout(path);
    save_to_yaml(fout, save_stimulus);
}

std::unique_ptr<PoissonLikelihood> PoissonLikelihood::from_yaml(
    const YAML::Node & node, std::shared_ptr<StimulusOccupancy> stimulus) {
    
    if (!(node["rate_scale"] && node["random_insertion"] && node["event_distribution"])) {
        throw std::runtime_error( "Cannot retrieve likelihood properties.");
    }
    
    value rate_scale = node["rate_scale"].as<value>();
    bool random_insertion = node["random_insertion"].as<bool>();
    auto event_dist = Mixture::from_yaml(node["event_distribution"]);
    
    // if stimulus distribution was not saved, it should be passed in as shared_ptr
    if (!node["stimulus_distribution"]) {
        if (stimulus==nullptr) {
            throw std::runtime_error("Stimulus distribution was not saved.");
        }
    } else {
        if (stimulus!=nullptr) {
            throw std::runtime_error("Found saved stimulus distribution and "
                "non-null stimulus argument.");
        }
        stimulus = StimulusOccupancy::from_yaml(node["stimulus_distribution"]);
    }
    
    // let create "empty" PoissonLikelihood using protected default constructor
    auto p = std::unique_ptr<PoissonLikelihood>(new PoissonLikelihood());
    
    // and set member variables
    p->event_distribution_ = std::move(event_dist);
    p->stimulus_distribution_ = stimulus;
    p->stimulus_grid_.reset(stimulus->grid().clone());
    p->logp_stimulus_.assign(p->stimulus_grid_->size(), 0.);
    p->event_rate_.assign(p->stimulus_grid_->size(), 0.);
    
    p->rate_scale_ = rate_scale;
    p->random_insertion_ = random_insertion;
    
    return p;
}


// flatbuffers
flatbuffers::Offset<fb_serialize::PoissonLikelihood> PoissonLikelihood::to_flatbuffers(flatbuffers::FlatBufferBuilder &builder, bool save_stimulus) const {

    auto events = event_distribution_->to_flatbuffers(builder);

    flatbuffers::Offset<fb_serialize::StimulusOccupancy> stim;

    if (save_stimulus) {
        stim = stimulus_distribution_->to_flatbuffers(builder);
    }

    fb_serialize::PoissonLikelihoodBuilder likelihood_builder(builder);

    likelihood_builder.add_rate_scale(rate_scale_);
    likelihood_builder.add_random_insertion(random_insertion_);
    likelihood_builder.add_event_distribution(events);

    if (save_stimulus) {
        likelihood_builder.add_stimulus_distribution(stim);
    }

    auto fb_lhood = likelihood_builder.Finish();

    return fb_lhood;
}

std::unique_ptr<PoissonLikelihood> PoissonLikelihood::from_flatbuffers(const fb_serialize::PoissonLikelihood * likelihood, std::shared_ptr<StimulusOccupancy> stimulus) {

    auto rate_scale = likelihood->rate_scale();
    bool random_insertion = likelihood->random_insertion();

    auto event_dist = Mixture::from_flatbuffers(likelihood->event_distribution());

    // if stimulus distribution was not saved, it should be passed in as shared_ptr
    if (likelihood->stimulus_distribution()==nullptr) {
        if (stimulus==nullptr) {
            throw std::runtime_error("Stimulus distribution was not saved and should be provided.");
        }
    } else {
        if (stimulus!=nullptr) {
            throw std::runtime_error("Found both saved stimulus distribution and non-null stimulus argument.");
        }
        stimulus = StimulusOccupancy::from_flatbuffers(likelihood->stimulus_distribution());
    }

    // let create "empty" PoissonLikelihood using protected default constructor
    auto p = std::unique_ptr<PoissonLikelihood>(new PoissonLikelihood());

    // and set member variables
    p->event_distribution_ = std::move(event_dist);
    p->stimulus_distribution_ = stimulus;
    p->stimulus_grid_.reset(stimulus->grid().clone());
    p->logp_stimulus_.assign(p->stimulus_grid_->size(), 0.);
    p->event_rate_.assign(p->stimulus_grid_->size(), 0.);

    p->rate_scale_ = rate_scale;
    p->random_insertion_ = random_insertion;

    return p;
}


// hdf5
void PoissonLikelihood::to_hdf5(HighFive::Group & group, bool save_stimulus) const {
    HighFive::DataSet ds_scale = group.createDataSet<value>("rate_scale", HighFive::DataSpace::From(rate_scale_));
    ds_scale.write(rate_scale_);
    
    HighFive::DataSet ds_rnd = group.createDataSet<bool>("random_insertion", HighFive::DataSpace::From(random_insertion_));
    ds_rnd.write(random_insertion_);
    
    HighFive::Group event_dist = group.createGroup("event_distribution");
    event_distribution_->to_hdf5(event_dist);
    
    if (save_stimulus) {
        HighFive::Group stim_dist = group.createGroup("stimulus_distribution");
        stimulus_distribution_->to_hdf5(stim_dist);
    }
}

std::unique_ptr<PoissonLikelihood> PoissonLikelihood::from_hdf5(
    const HighFive::Group & group, std::shared_ptr<StimulusOccupancy> stimulus) {
    
    value rate_scale;
    group.getDataSet("rate_scale").read(rate_scale);
    
    bool random_insertion;
    group.getDataSet("random_insertion").read(random_insertion);
    
    auto event_dist = Mixture::from_hdf5(group.getGroup("event_distribution"));
    
    // if stimulus distribution was not saved, it should be passed in as shared_ptr
    if (!group.exist("stimulus_distribution")) {
        if (stimulus==nullptr) {
            throw std::runtime_error("Stimulus distribution was not saved and should be provided.");
        }
    } else {
        if (stimulus!=nullptr) {
            throw std::runtime_error("Found both saved stimulus distribution and non-null stimulus argument.");
        }
        stimulus = StimulusOccupancy::from_hdf5(group.getGroup("stimulus_distribution"));
    }
    
    // let create "empty" PoissonLikelihood using protected default constructor
    auto p = std::unique_ptr<PoissonLikelihood>(new PoissonLikelihood());
    
    // and set member variables
    p->event_distribution_ = std::move(event_dist);
    p->stimulus_distribution_ = stimulus;
    p->stimulus_grid_.reset(stimulus->grid().clone());
    p->logp_stimulus_.assign(p->stimulus_grid_->size(), 0.);
    p->event_rate_.assign(p->stimulus_grid_->size(), 0.);
    
    p->rate_scale_ = rate_scale;
    p->random_insertion_ = random_insertion;
    
    return p;

}

void PoissonLikelihood::save_to_hdf5(std::string filename, bool save_stimulus, 
    int flags, std::string path) {
    
    HighFive::File file(filename, flags);
    
    // to do: create attributes for version, etc.
    
    if (!path.empty()) {
        HighFive::Group group = file.createGroup(path);
        this->to_hdf5(group, save_stimulus);
    } else {
        HighFive::Group group = file.getGroup("/");
        this->to_hdf5(group, save_stimulus);
    }
    
    file.flush();
}
    
std::unique_ptr<PoissonLikelihood> PoissonLikelihood::load_from_hdf5(
    std::string filename, std::string path, std::shared_ptr<StimulusOccupancy> stimulus) {
    
    HighFive::File file(filename, HighFive::File::ReadOnly);
    
    if (!path.empty()) {
        HighFive::Group group = file.getGroup(path);
        return PoissonLikelihood::from_hdf5(group, stimulus);
    } else {
        HighFive::Group group = file.getGroup("/");
        return PoissonLikelihood::from_hdf5(group, stimulus);
    }
}


