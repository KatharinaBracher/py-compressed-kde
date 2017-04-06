#pragma once

#include "space.hpp"

#include <vector>
#include <memory>

static const value THRESHOLD = 1.;

class PartialMixture;

class Mixture {
public:
    Mixture( const Space & space, value threshold = THRESHOLD );
    // copy constructor
    Mixture( const Mixture& other );
    
    friend Mixture * mixture_from_YAML( const YAML::Node & node );
    
    YAML::Node toYAML() const;
    
    void save( std::ostream & stream ) const;
    void save( std::string path ) const;
    
    void clear();
    
    const std::vector<value> & weights() const;
    const std::vector<std::unique_ptr<Component>> & components() const;
        
    // setters/getters
    value sum_of_weights() const;
    value sum_of_nsamples() const;
    value threshold() const;
    unsigned int ncomponents() const;
    
    const Space & space() const { return *space_; }
    
    void set_threshold( value v );
    
    void add_samples( const value * samples, unsigned int n, value w=1., value attenuation=1. );
    void merge_samples( const value * samples, unsigned int n, bool random = true, value w=1., value attenuation=1. );
    
    void evaluate( const value * points, unsigned int n, value * result );
    void evaluate( Grid & grid, value * result ) const;
    
    void partial( const value * points, unsigned int n, const std::vector<bool> & selection, value * result ) const;
    PartialMixture* partial( const value * points, unsigned int n, const std::vector<bool> & selection) const;
    void partial( Grid & grid, value * result ) const;
    PartialMixture* partial( Grid & grid ) const;
    
    void marginal( const value * points, unsigned int n, const std::vector<bool> & selection, value * result ) const;
    void marginal( Grid & grid, value * result ) const;
    
    //void partial( const Grid * grid, const std::vector<bool> & selection, value * result ) const {
        //// grid space has to be subspace of mixture space
        //// if array grid:
        //// partial( grid->full_grid().data(), grid->size(), selection, result )
        //// if (multi)vector grid:
        //// for each kernel:
        ////   evaluate at grid: space_->partial_logp( grid, selection, <tmp> )
        ////   combine probability vectors -> result
    //}
    //void marginal( const Grid * grid, const std::vector<bool> & selection, value * result ) const;
    
protected:
    value update_weights_( unsigned int nsamples );
    value update_weights_( unsigned int nsamples, value weight, value attenuation );
    
    bool closest( const Component & c, unsigned int & index, value threshold_squared) const;
    bool closest( const Component & c, unsigned int & index ) const;
    
protected:
    value sum_of_weights_;
    value sum_of_nsamples_;
    value threshold_;
    value threshold_squared_;
    
    std::unique_ptr<Space> space_;
    std::vector<std::unique_ptr<Component>> kernels_;
    std::vector<value> weights_;
};

Mixture * mixture_from_YAML( const YAML::Node & node );

Mixture * load_mixture( std::string path );

class PartialMixture {
public:
    PartialMixture( const Mixture * source, const std::vector<bool> & selection, const value * points, unsigned int n );
    PartialMixture( const Mixture * source, Grid & grid );
    
    const Mixture & mixture() const;
    
    unsigned int ncomponents() const;
    unsigned int nsamples() const;
    const std::vector<bool> & selection() const;
    const std::vector<bool> & inverse_selection() const;
    
    const std::vector<long unsigned int> & partial_shape() const;
    
    void complete ( const value * points, unsigned int n, value * result ) const;
    void complete_multi ( const value * points, unsigned int n, value * result) const; //, value * offset = nullptr ) const;
        
    template <class result_it>
    void marginal(result_it result) {
        
        auto it = partial_logp_.cbegin();
        
        for (auto & w : mixture_.weights()) {
            
            for (unsigned int s=0; s<nsamples_; ++s) {
                result[s] += w * fastexp(*it);
                ++it;
            }
        }
        
    }
    
    const std::vector<value> & partial_logp() const;
    
    const Space & space() const { return mixture_.space(); }
    //Space & partialspace() const { return mixture_.template_component()->dataspace_selection(selection_.cbegin()); }
    //Space & inversespace() const { return mixture_.template_component()->dataspace_selection(inverted_selection_.cbegin()); }
    
protected:
    Mixture mixture_;
    unsigned int nsamples_;
    std::vector<bool> selection_;
    std::vector<bool> inverted_selection_;
    std::vector<value> partial_logp_;
    std::vector<long unsigned int> partial_shape_;
};
