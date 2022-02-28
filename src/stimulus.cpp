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
#include "stimulus.hpp"

// constructor
StimulusOccupancy::StimulusOccupancy( const Space & space, const Grid & grid, double stimulus_duration, value compression ) :
stimulus_duration_(stimulus_duration), compression_(compression) {
    
    if (!(space.specification()==grid.specification())) {
        throw std::runtime_error("Grid does not match stimulus space.");
    }
    
    stimulus_distribution_.reset( new Mixture( space, compression ) );
    stimulus_grid_.reset( grid.clone() );

}

// properties
const Space & StimulusOccupancy::space() const {
    return stimulus_distribution_->space();
}

const Grid & StimulusOccupancy::grid() const { return *stimulus_grid_; }

value StimulusOccupancy::compression() const { return compression_; }

double StimulusOccupancy::stimulus_duration() const { return stimulus_duration_; }

unsigned int StimulusOccupancy::ndim() const { return stimulus_grid_->ndim(); }

bool StimulusOccupancy::random_insertion() const { return random_insertion_; }
void StimulusOccupancy::set_random_insertion(bool val) { random_insertion_ = val; }

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

// methods
void StimulusOccupancy::add_stimulus( const std::vector<value> & stimuli, 
    unsigned repetitions ) {
    
    unsigned int n = stimuli.size() / stimulus_distribution_->space().ndim();
    if (n*stimulus_distribution_->space().ndim() != stimuli.size()) {
        throw std::runtime_error("Not a whole number of samples.");
    }
    
    add_stimulus( stimuli.data(), n, repetitions );
}
    
void StimulusOccupancy::add_stimulus( const value * stimuli, unsigned int n, 
    unsigned int repetitions ) {
    
    lock_.lock();
    stimulus_distribution_->merge_samples( stimuli, n, random_insertion_, 
        static_cast<value>(repetitions) );
    lock_.unlock();
}

// yaml
YAML::Node StimulusOccupancy::to_yaml( ) const {
    
    YAML::Node node;
    
    node["stimulus_duration"] = stimulus_duration_;
    node["compression"] = compression_;
    node["random_insertion"]  = random_insertion_;
    
    node["stimulus_distribution"] = stimulus_distribution_->to_yaml();
    node["stimulus_grid"] = stimulus_grid_->to_yaml();
    
    return node;
}

void StimulusOccupancy::save_to_yaml( std::ostream & stream ) const {
    
    auto node = to_yaml(); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}

void StimulusOccupancy::save_to_yaml( std::string path ) const {
    
    std::ofstream fout(path);
    save_to_yaml(fout);
}

std::unique_ptr<StimulusOccupancy> StimulusOccupancy::from_yaml(
    const YAML::Node & node ) {
        
    if ( !(node["stimulus_duration"] && node["compression"] && 
           node["random_insertion"]) ) {
        throw std::runtime_error( "Cannot retrieve stimulus properties.");
    }
    
    if ( !(node["stimulus_distribution"] && node["stimulus_grid"]) ) {
        throw std::runtime_error( "Cannot retrieve stimulus distribution or grid.");
    }
    
    double duration = node["stimulus_duration"].as<double>();
    value compression = node["compression"].as<value>();
    bool random_insertion = node["random_insertion"].as<bool>();
    
    auto grid = grid_from_yaml(node["stimulus_grid"]);
    
    auto mix = Mixture::from_yaml(node["stimulus_distribution"]);
    
    auto stim = std::make_unique<StimulusOccupancy>(mix->space(), *grid,
        duration, compression);
    stim->set_random_insertion(random_insertion);
    
    stim->stimulus_distribution_ = std::move(mix);
    
    return stim;
}

// flatbuffers
flatbuffers::Offset<fb_serialize::StimulusOccupancy> StimulusOccupancy::to_flatbuffers(flatbuffers::FlatBufferBuilder &builder) const {

    auto stim = stimulus_distribution_->to_flatbuffers(builder);
    auto grid = stimulus_grid_->to_flatbuffers(builder);

    return fb_serialize::CreateStimulusOccupancy(
        builder,
        stimulus_duration_,
        compression_,
        random_insertion_,
        stim,
        grid
    );
}

std::unique_ptr<StimulusOccupancy> StimulusOccupancy::from_flatbuffers(const fb_serialize::StimulusOccupancy * stimulus) {

    auto duration = stimulus->stimulus_duration();
    auto compression = stimulus->compression();
    auto random_insertion = stimulus->random_insertion();

    auto grid = grid_from_flatbuffers(stimulus->stimulus_grid());

    auto mix = Mixture::from_flatbuffers(stimulus->stimulus_distribution());

    auto stim = std::make_unique<StimulusOccupancy>(mix->space(), *grid, 
        duration, compression);

    stim->set_random_insertion(random_insertion);
    stim->stimulus_distribution_ = std::move(mix);

    return stim;
}


// hdf5
void StimulusOccupancy::to_hdf5(HighFive::Group & group) const {
    
    HighFive::DataSet ds_dur = group.createDataSet<double>("stimulus_duration", HighFive::DataSpace::From(stimulus_duration_));
    ds_dur.write(stimulus_duration_);
    
    HighFive::DataSet ds_comp = group.createDataSet<value>("compression", HighFive::DataSpace::From(compression_));
    ds_comp.write(compression_);
    
    HighFive::DataSet ds_rnd = group.createDataSet<bool>("random_insertion", HighFive::DataSpace::From(random_insertion_));
    ds_rnd.write(random_insertion_);
    
    HighFive::Group stim = group.createGroup("stimulus_distribution");
    stimulus_distribution_->to_hdf5(stim);
    
    HighFive::Group grid = group.createGroup("stimulus_grid");
    stimulus_grid_->to_hdf5(grid);
}

std::unique_ptr<StimulusOccupancy> StimulusOccupancy::from_hdf5(
    const HighFive::Group & group) {
    
    double duration;
    group.getDataSet("stimulus_duration").read(duration);
    
    value compression;
    group.getDataSet("compression").read(compression);
    
    bool random_insertion;
    group.getDataSet("random_insertion").read(random_insertion);
    
    auto grid = grid_from_hdf5(group.getGroup("stimulus_grid"));
    
    auto mix = Mixture::from_hdf5(group.getGroup("stimulus_distribution"));
    
    auto stim = std::make_unique<StimulusOccupancy>(mix->space(), *grid, 
        duration, compression);
    stim->set_random_insertion(random_insertion);
    
    stim->stimulus_distribution_ = std::move(mix);
    
    return stim;
    
}

void StimulusOccupancy::save_to_hdf5(std::string filename, int flags, 
    std::string path) {

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
    
std::unique_ptr<StimulusOccupancy> StimulusOccupancy::load_from_hdf5(
    std::string filename, std::string path) {
        
    HighFive::File file(filename, HighFive::File::ReadOnly);
    
    if (!path.empty()) {
        HighFive::Group group = file.getGroup(path);
        return StimulusOccupancy::from_hdf5(group);
    } else {
        HighFive::Group group = file.getGroup("/");
        return StimulusOccupancy::from_hdf5(group);
    }
}

