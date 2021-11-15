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
#include "mixture.hpp"
#include <random>
#include <algorithm>

#include <iostream>

// constructor
Mixture::Mixture( const Space & space, value threshold ) :
sum_of_weights_(0), sum_of_nsamples_(0),
threshold_(threshold), threshold_squared_(threshold*threshold),
space_(space.clone()) {}

// copy constructor
Mixture::Mixture( const Mixture& other ) :
sum_of_weights_(other.sum_of_weights_), sum_of_nsamples_(other.sum_of_nsamples_),
threshold_(other.threshold_), threshold_squared_(other.threshold_squared_), 
space_(other.space_->clone()), weights_(other.weights_) {
    for (auto & c : other.kernels_ ) {
        kernels_.emplace_back( new Component( *c ) );
    }
}

void Mixture::clear() {
    sum_of_weights_ = 0;
    sum_of_nsamples_ = 0;
    kernels_.clear();
    weights_.clear();
}

// properties
value Mixture::sum_of_weights() const { return sum_of_weights_; }
value Mixture::sum_of_nsamples() const { return sum_of_nsamples_; }
value Mixture::threshold() const { return threshold_; }
unsigned int Mixture::ncomponents() const { return kernels_.size(); }

const std::vector<value> & Mixture::weights() const {
    return weights_;
}

const std::vector<std::unique_ptr<Component>> & Mixture::components() const {
    return kernels_;
}

void Mixture::set_threshold( value v )  {
    if (v<0.) {
        throw std::runtime_error("Threshold should be larger than or equal to 0.");
    }
    
    threshold_ = v;
    threshold_squared_=threshold_*threshold_;
}

// methods
void Mixture::add_samples( const value * samples, unsigned int n, value w, value attenuation ) {
    
    for (unsigned int k=0; k<n; ++k) {
        try {
            kernels_.emplace_back( space_->kernel( samples ) ); // note that sample are not checked and could contain invalid values!
        } catch (std::exception & me) { // ToDo: catch proper exception
            // remove all new components_
            kernels_.resize( kernels_.size() - k );
            throw;
        }
        samples += space_->ndim();
    }
    
    value weight = update_weights_( n, w, attenuation );
    weights_.insert( weights_.end(), n, weight );
    
}

void Mixture::merge_samples( const value * samples, unsigned int n, bool random, value w, value attenuation ) {
    
    if (threshold_==0.) {
        add_samples( samples, n );
        return;
    }
    
    unsigned int index=0;

    //// create vector of new components
    std::vector<std::unique_ptr<Component>> new_kernels;
    
    for (unsigned int k=0; k<n; ++k) {
        new_kernels.emplace_back( space_->kernel( samples ) );
        samples += space_->ndim();
    }
    
    if (random) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle( new_kernels.begin(), new_kernels.end(), g );
    }
    
    value weight = update_weights_( n, w, attenuation );
    
    for (auto & c : new_kernels) {
        if (closest( *c, index )) {
            space_->merge( weights_[index], *kernels_[index], weight, *c );
            weights_[index]+=weight;
        } else { // add
            kernels_.push_back( std::move( c ) );
            weights_.push_back(weight);
        }
    }
    
}

void Mixture::evaluate( const value * points, unsigned int n, value * result ) {
    
    auto weight = weights_.cbegin();
    
    //for (unsigned int k=0; k<n; ++k ) {
        //*result = 0;
        //weight = weights_.cbegin();
        //for (auto & c : kernels_ ) {
            //*result += *weight++ * space_->probability( *c, points );
        //}
        //points += space_->ndim();
        //++result;
    //}
    
    const value * ptr;
    value * res;
    
    std::fill( result, result+n, 0. );
    
    for (auto & c : kernels_) {
        
        ptr = points;
        res = result;
        
        for (unsigned int k=0; k<n; ++k) {
            *res += *weight * space_->probability( *c, ptr );
            ++res;
            ptr += space_->ndim();
        }
        
        ++weight;
        
    }
    
}

void Mixture::evaluate( Grid & grid, value * result ) const {
    
    if (!(grid.specification()==space_->specification())) {
        throw std::runtime_error("Grid does not have the required space specification.");
    }
        
    
    // result array is assumed to have size grid.size()
    
    std::fill( result, result+grid.size(), 0. );
    
    auto weight = weights_.cbegin();
    
    for (auto & c : kernels_) {
        
        space_->probability( grid, *weight * c->scale_factor, *c, result );
        
        ++weight;
        
    }
    
}

void Mixture::partial( const value * points, unsigned int n, const std::vector<bool> & selection, value * result ) const {
    
    if (selection.size()!=space_->ndim()) {
        throw std::runtime_error("Incorrect selection.");
    }
    
    unsigned int ndim = std::count( selection.begin(), selection.end(), true );
    
    value log_scale;
    const value * ptr = points;
    
    for (auto & c : kernels_) {
        
        log_scale = space_->compute_scale_factor( *c, selection, true );
                
        for (unsigned int s=0; s<n; ++s) {
            
            *result = space_->partial_logp( *c, ptr, selection ) + log_scale;
            
            ++result;
            ptr += ndim;
            
        }
        
        ptr = points;
        
    }
}

PartialMixture* Mixture::partial( const value * points, unsigned int n, const std::vector<bool> & selection) const {
    return new PartialMixture(this, selection, points, n);
}

void Mixture::partial( Grid & grid, value * result ) const {
    
    value log_scale;
    
    auto selection = space().specification().selection( grid.specification() );
    
    for (auto & c : kernels_) {
        
        log_scale = space_->compute_scale_factor( *c, selection, true );
        
        space_->partial_logp( grid, selection, log_scale, *c, result );
        
        result += grid.size();
        
    }
    
}

PartialMixture* Mixture::partial( Grid & grid ) const {
    return new PartialMixture(this, grid);
}

void Mixture::marginal( const value * points, unsigned int n, const std::vector<bool> & selection, value * result ) const {
    
    if (selection.size()!=space_->ndim()) {
        throw std::runtime_error("Incorrect selection.");
    }
    
    unsigned int ndim = std::count( selection.begin(), selection.end(), true );
    
    value * out;
    value tmp;
    
    std::vector<value>::const_iterator weight = weights_.cbegin();
    
    value log_scale;
    const value * ptr = points;
    
    for (auto & c : kernels_) {
        
        out = result;
        
        log_scale = space_->compute_scale_factor( *c, selection, true );
        
        for (unsigned int s=0; s<n; ++s) {
            tmp = space_->partial_logp( *c, ptr, selection ) + log_scale;
            
            if (!std::isinf(tmp)) {
                *out += *weight * fastexp( tmp );
            }
            
            ++out;
            ptr += ndim;
        }
        
        ++weight;
        ptr = points;
    }
}

void Mixture::marginal( Grid & grid, value * result ) const {
    
    value log_scale;
    auto selection = space().specification().selection( grid.specification() );
    std::vector<value>::const_iterator weight = weights_.cbegin();
    std::vector<value> tmp(grid.size());
    
    for (auto & c : kernels_) {
        log_scale = space_->compute_scale_factor( *c, selection, true );
        log_scale += fastlog( *weight );
        
        space_->partial_logp( grid, selection, log_scale, *c, tmp.data() );
        
        for (unsigned int k=0; k<grid.size(); ++k) {
            if (!std::isinf(tmp[k])) {
                result[k] += fastexp( tmp[k] );
            }
        }
        ++weight;
    }
}

// protected methods
value Mixture::update_weights_( unsigned int nsamples ) {
    
    sum_of_weights_ += nsamples;
    sum_of_nsamples_ += nsamples;
    
    value mixing_factor_old = (sum_of_weights_-nsamples)/sum_of_weights_;
    value mixing_factor_new = nsamples/sum_of_weights_;
    
    value weight = mixing_factor_new / nsamples;
    
    //adjust weights of existing components
    for (auto & w : weights_) {
        w *= mixing_factor_old;
    }
    
    return weight;
    
}

value Mixture::update_weights_( unsigned int nsamples, value weight, value attenuation ) {
    
    value attenuated_sum_weights = sum_of_weights_ * attenuation;
    value sum_sample_weights = nsamples * weight; //assuming weight per sample of 1
    sum_of_weights_ = attenuated_sum_weights + sum_sample_weights;
    
    value attenuated_sum_n = sum_of_nsamples_ * attenuation;
    sum_of_nsamples_ = attenuated_sum_n + nsamples;
    
    value mixing_factor_old = attenuated_sum_weights/sum_of_weights_;
    value mixing_factor_new = sum_sample_weights/sum_of_weights_;
    
    value w = mixing_factor_new / sum_sample_weights;
    
    //adjust weights of existing components
    for (auto & k : weights_) {
        k *= mixing_factor_old;
    }
    
    return w;
    
}

bool Mixture::closest( const Component & target, unsigned int & index, value threshold_squared) const {
    
    value min_distance = threshold_squared;
    value distance;
    
    for (unsigned int k=0; k<kernels_.size(); ++k) {
        //distance = kernels_[k]->mahalanobis_distance_squared( target, threshold_squared );
        distance = space_->mahalanobis_distance_squared( *kernels_[k], target, threshold_squared );
        if (distance<min_distance) {
            min_distance=distance;
            index = k;
        }
    }
    
    return (min_distance<threshold_squared);
}

bool Mixture::closest( const Component & k, unsigned int & index ) const {
    return closest( k, index, threshold_squared_ );
}


// yaml
YAML::Node Mixture::to_yaml() const {
    YAML::Node node;
    
    node["sum_of_weights"] = sum_of_weights_;
    node["sum_of_nsamples"] = sum_of_nsamples_;
    node["threshold"] = threshold_;
    node["nkernels"] = kernels_.size();
    node["space"] = space_->to_yaml();
    
    node["kernels"] = YAML::Load("[]");
    
    for (auto & k : kernels_) {
        node["kernels"].push_back( k->to_yaml() );
    }
    node["weights"] = weights_;
    
    return node;
    
}

void Mixture::save_to_yaml( std::ostream & stream ) const {
    auto node = to_yaml(); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}
void Mixture::save_to_yaml( std::string path ) const {
    std::ofstream fout(path);
    save_to_yaml(fout);
}

std::unique_ptr<Mixture> Mixture::from_yaml( const YAML::Node & node ) {
    
    if (!node["space"] || !node["weights"].IsSequence() || 
        !node["kernels"].IsSequence() || 
         node["weights"].size()!=node["kernels"].size()) {
        throw std::runtime_error( "Cannot retrieve weights or kernels.");
    }
    
    auto space = space_from_yaml( node["space"] );
    value threshold = node["threshold"].as<value>( THRESHOLD );
    
    auto m = std::make_unique<Mixture>(*space, threshold);
    
    unsigned int nkernels = node["kernels"].size();
    
    m->sum_of_weights_ = node["sum_of_weights"].as<value>( nkernels );
    m->sum_of_nsamples_ = node["sum_of_nsamples"].as<value>( nkernels );
    
    m->weights_ = node["weights"].as<std::vector<value>>();
    
    for (unsigned int k=0; k<nkernels; ++k){
        try {
            m->kernels_.push_back(std::move(Component::from_yaml( node["kernels"][k])));
            
            if (m->kernels_.back()->location.size()!=space->ndim() || 
                m->kernels_.back()->bandwidth.size()!=space->nbw()) {
                throw std::runtime_error("Component vector sizes do not match space.");
            }
            
            space->update_scale_factor( *m->kernels_.back() );
            
        } catch (std::exception & me) { // ToDo: catch proper exception
            throw std::runtime_error("Cannot load kernel data.");
        }
    }
    
    return m;
}
std::unique_ptr<Mixture> Mixture::load_from_yaml( std::string path ) {
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return Mixture::from_yaml( node );
}


// hdf5
void Mixture::to_hdf5(HighFive::Group & group) const {
    
    HighFive::DataSet ds_sow = group.createDataSet<value>("sum_of_weights", HighFive::DataSpace::From(sum_of_weights_));
    ds_sow.write(sum_of_weights_);
    
    HighFive::DataSet ds_son = group.createDataSet<value>("sum_of_nsamples", HighFive::DataSpace::From(sum_of_nsamples_));
    ds_son.write(sum_of_nsamples_);
    
    HighFive::DataSet ds_th = group.createDataSet<value>("threshold", HighFive::DataSpace::From(threshold_));
    ds_th.write(threshold_);
    
    HighFive::DataSet ds_nk = group.createDataSet<unsigned int>("nkernels", HighFive::DataSpace::From(kernels_.size()));
    ds_nk.write(kernels_.size());
    
    HighFive::Group space_group = group.createGroup("space");
    space_->to_hdf5(space_group);
    
    HighFive::DataSet ds_w = group.createDataSet<value>("weights", HighFive::DataSpace::From(weights_));
    ds_w.write(weights_);
    
    HighFive::Group subgroup = group.createGroup("kernels");
    
    HighFive::DataSet ds_loc = subgroup.createDataSet<value>("location",
        HighFive::DataSpace({space_->ndim(),kernels_.size()}));
    
    HighFive::DataSet ds_bw = subgroup.createDataSet<value>("bandwidth",
        HighFive::DataSpace({space_->nbw(), kernels_.size()}));
    
    unsigned int n = 0;
    
    for (auto & k : kernels_) {
        
        ds_loc.select({0,n},{space_->ndim(),1}).write(k->location);
        ds_bw.select({0,n},{space_->nbw(),1}).write(k->bandwidth);
        
        ++n;
    }
    
}

void Mixture::save_to_hdf5( std::string filename, int flags, std::string path) {
    
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

std::unique_ptr<Mixture> Mixture::load_from_hdf5(std::string filename, 
    std::string path) {
        
    HighFive::File file(filename, HighFive::File::ReadOnly);
    
    if (!path.empty()) {
        HighFive::Group group = file.getGroup(path);
        return Mixture::from_hdf5(group);
    } else {
        HighFive::Group group = file.getGroup("/");
        return Mixture::from_hdf5(group);
    }
}

std::unique_ptr<Mixture> Mixture::from_hdf5(const HighFive::Group & group) {
    
    auto space = space_from_hdf5(group.getGroup("space"));
    
    value threshold;
    group.getDataSet("threshold").read(threshold);
    
    auto m = std::make_unique<Mixture>(*space, threshold);
    
    unsigned int nkernels;
    group.getDataSet("nkernels").read(nkernels);
        
    group.getDataSet("sum_of_weights").read(m->sum_of_weights_);
    group.getDataSet("sum_of_nsamples").read(m->sum_of_nsamples_);
    group.getDataSet("weights").read(m->weights_);
    
    HighFive::DataSet loc = group.getGroup("kernels").getDataSet("location");
    HighFive::DataSet bw = group.getGroup("kernels").getDataSet("bandwidth");
    
    for (unsigned int k=0; k<nkernels; ++k){
        try {
            auto comp = std::make_unique<Component>();
            loc.select({0,k},{space->ndim(),1}).read(comp->location);
            bw.select({0,k},{space->nbw(),1}).read(comp->bandwidth);
            
            m->kernels_.push_back(std::move(comp));
            
            space->update_scale_factor( *m->kernels_.back() );
            
        } catch (std::exception & me) { // ToDo: catch proper exception
            throw std::runtime_error("Cannot load kernel data.");
        }
    }
    
    return m;    
}



// constructors
PartialMixture::PartialMixture( const Mixture * source, const std::vector<bool> & selection, const value * points, unsigned int n ) :
mixture_(*source), nsamples_(n), selection_(selection), inverted_selection_(selection) {        
    
    partial_logp_.resize( mixture_.ncomponents() * nsamples_ );
    mixture_.partial( points, nsamples_, selection_, partial_logp_.data() );
    inverted_selection_.flip();
    partial_shape_ = { nsamples_ };
}

PartialMixture::PartialMixture( const Mixture * source, Grid & grid ) :
mixture_(*source), nsamples_(grid.size()), selection_(source->space().specification().selection(grid.specification())), inverted_selection_(selection_) {        
    
    partial_logp_.resize( mixture_.ncomponents() * nsamples_ );
    mixture_.partial( grid, partial_logp_.data() );
    inverted_selection_.flip();
    partial_shape_ = grid.shape();
}

// properties
const Mixture & PartialMixture::mixture() const {
    return mixture_;
}
    
unsigned int PartialMixture::ncomponents() const {
    return partial_logp_.size() / nsamples_;
}

unsigned int PartialMixture::nsamples() const {
    return nsamples_;
}

const std::vector<bool> & PartialMixture::selection() const { 
    return selection_;
}

const std::vector<bool> & PartialMixture::inverse_selection() const { 
    return inverted_selection_;
}

const std::vector<long unsigned int> & PartialMixture::partial_shape() const {
    return partial_shape_;
}
   
const std::vector<value> & PartialMixture::partial_logp() const {
    return partial_logp_;
}

// methods
void PartialMixture::complete ( const value * points, unsigned int n, value * result ) const {
    
    if (ncomponents() != mixture().ncomponents()) {
        throw std::runtime_error("Number of kernels in source mixture has changed.");
    }
    
    value x=0;
    value scale = 0.;
    const value * ptr = points;
    value * presult = result;
    
    unsigned int ndim = std::count( inverted_selection_.begin(), inverted_selection_.end(), true );
    
    auto w = mixture_.weights().cbegin();
    auto it = partial_logp_.cbegin();
    
    for (auto & c : mixture_.components()) {
        
        scale = mixture_.space().compute_scale_factor( *c, inverted_selection_, true );
        
        ptr = points;
        presult = result;
        
        for (unsigned int k=0; k<n; ++k) {
            
            x = mixture_.space().partial_logp( *c, ptr, inverted_selection_) + scale;
            
            ptr += ndim;
            
            if (std::isinf(x)) { presult+=nsamples_; continue; }
            
            for (unsigned int s=0; s<nsamples_; ++s) {
                
                if (std::isinf(it[s])) { ++presult; continue; }
                
                *presult++ += (*w) * fastexp(it[s] + x);
            }
            
        }
        
        ++w;
        it+=nsamples_;
        
    }
    
}

void PartialMixture::complete_multi ( const value * points, unsigned int n, value * result ) const { //, value * offset ) const {
    
    if (ncomponents() != mixture().ncomponents()) {
        throw std::runtime_error("Number of kernels in source mixture has changed.");
    }
    
    value x=0;
    value scale;
    const value * ptr = points;
    
    std::vector<value> tmp(nsamples_);
    
    unsigned int ndim = std::count( inverted_selection_.begin(), inverted_selection_.end(), true );
    
    for (unsigned int k=0; k<n; ++k) {
        
        auto it = partial_logp_.cbegin();
        auto w = mixture_.weights().cbegin();
        
        std::fill( tmp.begin(), tmp.end(), 0. );
        
        for (auto & c : mixture_.components()) {
            
            x = mixture_.space().partial_logp( *c, ptr, inverted_selection_);
            if (std::isinf(x)) { it+=nsamples_; ++w; continue; }
            
            scale = mixture_.space().compute_scale_factor( *c, inverted_selection_, true );
            x += scale;
            
            for (unsigned int s=0; s<nsamples_; ++s) {
                tmp[s] += (*w) * fastexp(*it + x);
                ++it;
            }
            
            ++w;
            
        }
        
        ptr += ndim;
               
        // add offset
        //if (offset!=nullptr) {
        //   std::transform( tmp.begin(), tmp.end(), offset, tmp.begin(), std::plus<value>() );
        //}
        
        // add log of temporary vector to result
        std::transform( tmp.begin(), tmp.end(), result, result, [](const value & a, const value & b) { return fastlog(a) + b; } );
    }
    
}

