
#include "mixture.hpp"
#include <random>
#include <algorithm>

#include <iostream>

Mixture::Mixture( const Space & space, value threshold ) :
sum_of_weights_(0), sum_of_nsamples_(0),
threshold_(threshold), threshold_squared_(threshold*threshold),
space_(space.clone()) {}

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

// setters/getters
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

// WEIGHT UPDATING
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

YAML::Node Mixture::toYAML() const {
    YAML::Node node;
    
    node["sum_of_weights"] = sum_of_weights_;
    node["sum_of_nsamples"] = sum_of_nsamples_;
    node["threshold"] = threshold_;
    node["nkernels"] = kernels_.size();
    node["space"] = space_->toYAML();
    
    node["kernels"] = YAML::Load("[]");
    
    for (auto & k : kernels_) {
        node["kernels"].push_back( k->toYAML() );
    }
    node["weights"] = weights_;
    
    return node;
    
}

Mixture * mixture_from_YAML( const YAML::Node & node ) {
    
    if (!node["space"] || !node["weights"].IsSequence() || !node["kernels"].IsSequence() || node["weights"].size()!=node["kernels"].size()) {
        throw std::runtime_error( "Cannot retrieve weights or kernels.");
    }
    
    auto space = std::unique_ptr<Space>( space_from_YAML( node["space"] ) );
    value threshold = node["threshold"].as<value>( THRESHOLD );
    
    auto m = std::unique_ptr<Mixture>( new Mixture( *space, threshold ) );
    
    unsigned int nkernels = node["kernels"].size();
    
    m->sum_of_weights_ = node["sum_of_weights"].as<value>( nkernels );
    m->sum_of_nsamples_ = node["sum_of_nsamples"].as<value>( nkernels );
    
    m->weights_ = node["weights"].as<std::vector<value>>();
    
    for (unsigned int k=0; k<nkernels; ++k){
        try {
            m->kernels_.push_back( std::unique_ptr<Component>( Component::fromYAML( node["kernels"][k] ) ) );
            
            if (m->kernels_.back()->location.size()!=space->ndim() || m->kernels_.back()->bandwidth.size()!=space->nbw()) {
                throw std::runtime_error("Component vector sizes do not match space.");
            }
            
            space->update_scale_factor( *m->kernels_.back() );
            
            //m->kernels_.back()->update_fromYAML( node["kernels"][k] );
        } catch (std::exception & me) { // ToDo: catch proper exception
            throw std::runtime_error("Cannot load kernel data.");
        }
    }
    
    return m.release();
}

Mixture * load_mixture( std::string path ) {
    std::ifstream ifs(path, std::ifstream::in);
    
    auto node = YAML::Load( ifs );
    
    return mixture_from_YAML( node );
}

void Mixture::save( std::ostream & stream ) const {
    auto node = toYAML(); 
    YAML::Emitter out; 
    out << YAML::Flow; 
    out << node; 
    stream << out.c_str();
}
void Mixture::save( std::string path ) const {
    std::ofstream fout(path); save(fout);
}


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


//void PartialMixture::complete ( const value * points, unsigned int n, value * result ) const {
    
    //if (ncomponents() != mixture().ncomponents()) {
        //throw std::runtime_error("Number of kernels in source mixture has changed.");
    //}
    
    //value x=0;
    //value scale = 0.;
    //const value * ptr = points;
    
    //unsigned int ndim = std::count( inverted_selection_.begin(), inverted_selection_.end(), true );
    
    //for (unsigned int k=0; k<n; ++k) {
        
        //auto it = partial_logp_.cbegin();
        //auto w = mixture_.weights().cbegin();
        
        //for (auto & c : mixture_.components()) {
            
            //scale = mixture_.space().compute_scale_factor( *c, inverted_selection_, true );
            //x = mixture_.space().partial_logp( *c, ptr, inverted_selection_) + scale;
            
            //if (std::isinf(x)) { it+=nsamples_; ++w; continue; }
            
            //for (unsigned int s=0; s<nsamples_; ++s) {
                //result[s] += (*w) * fastexp(*it + x);
                //++it;
            //}
            
            //++w;
            
        //}
        
        //ptr += ndim;
        //result += nsamples_; // TODO: shouldn't this be += nsamples instead of += n ??
    //}
    
//}

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
