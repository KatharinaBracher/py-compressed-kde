#include "space_circular.hpp"
#include "kernel_vonmises.hpp"

CircularSpace::CircularSpace( std::string name, value kappa, value mu ) :
SpaceBase<CircularSpace>( "circular", make_spec(name), make_kernel(kappa, mu) ) {}

SpaceSpecification CircularSpace::make_spec( std::string name ) {
    
   return SpaceSpecification( DimSpecification(name, "circular", "") );
}
    
Component CircularSpace::make_kernel(value kappa, value mu) const {
    Component k;
    k.location = { mu };
    k.bandwidth = { kappa };
    k.scale_factor = vonmises_scale_factor(kappa, false);
    k.scale_factor_log = fastlog( k.scale_factor );
    return k;
}

CircularSpace* CircularSpace::fromYAML( const YAML::Node & node ) {
    
    if (!node["name"]) {
        throw std::runtime_error("Ill-formed circular space definition.");
    }
    
    CircularSpace * p = new CircularSpace( node["name"].as<std::string>() );
    
    return p;
    
} 

YAML::Node CircularSpace::asYAML() const {
    YAML::Node node;
    node["name"] = specification().dim(0).name();
    return node;
}

value CircularSpace::compute_scale_factor( value * bw, bool log ) const {
    return vonmises_scale_factor( *bw, log );
}

value CircularSpace::compute_scale_factor( std::vector<bool>::const_iterator selection, value * bw, bool log ) const {
    value s;
    if (*selection) {
        s = vonmises_scale_factor( *bw, log );
    } else {
        if (log) { s = 0.;} 
        else { s = 1.;}
    }
    return s;
}


value CircularSpace::mahalanobis_distance_squared( const value * refloc, const value * refbw, const value * targetloc, value threshold) const {
    // for large kappa, the von mises distribution can be approximated
    // by a Gaussian distribution with variance = 1/kappa
    // we will use this approximation to compute an equivalent
    // mahalanobis distance
    
    // for small kappa, this approximation will lead to an overestimate
    // of the distance
    
    value d = circular_difference(*targetloc, *refloc);
    d = (d*d)*(*refbw);
    
    return d;
}

void CircularSpace::merge( value w1, value * loc1, value * bw1, value w2, const value * loc2, const value * bw2 ) const {
    
    // mean: G( mu1 + (w2/(w1+w2)) * F( mu2 - mu1 ) )
    // where G(x): x+2*PI if x<0
    //             x-2*PI if x>2*PI
    //             x otherwise
    // where F(x): x+2*PI if x<=-PI
    //             x-2*PI if x>PI
    //             x otherwise
    // kappa: 1/K = w1/K1 + w2/K2 - w1*w2*( PI-|PI-|mu1-mu2||)^2
    
    value sum_w = w1 + w2;
    value tmp = *loc2 - *loc1;
    
    *bw1 = w1/(*bw1 * sum_w) + w2/(*bw2 * sum_w) + w1*w2*std::pow( M_PI-std::fabs(M_PI-std::fabs(tmp)), 2 )/(sum_w*sum_w);
    *bw1 = 1./(*bw1);
    
    if (tmp<=-M_PI) {tmp+=2*M_PI;}
    else if (tmp>M_PI) {tmp-=2*M_PI;}
    
    *loc1 = *loc1 + (tmp * w2)/sum_w;
    
    if (*loc1<0) {*loc1+=2*M_PI;}
    else if (*loc1>2*M_PI) {*loc1-=2*M_PI;}
    
}

value CircularSpace::probability( const value * loc, const value * bw, const value * point ) const {
    
    value p;
    if (*bw>KAPPA_GAUSS_APPROX) {
        p = circular_difference( *point, *loc );
        p = -0.5 * (p*p) * (*bw);
    } else {
        p = *bw * std::cos( *point - *loc );
    }
    return fastexp(p);
    
}

void CircularSpace::probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    
    for (unsigned int k=0; k<n; ++k) {
        *result++ = probability( loc, bw, points++ );
    }
    
}

value CircularSpace::log_probability( const value * loc, const value * bw, const value * point ) const {
    
    value p;
    if (*bw>KAPPA_GAUSS_APPROX) {
        p = circular_difference( *point, *loc );
        p = -0.5 * (p*p) * (*bw);
    } else {
        p = *bw * std::cos( *point - *loc );
    }
    return p;
    
}

void CircularSpace::log_probability( const value * loc, const value * bw, const value * points, unsigned int n, value * result ) const {
    
    for (unsigned int k=0; k<n; ++k) {
        *result++ = log_probability( loc, bw, points++ );
    }
    
}

value CircularSpace::partial_logp( const value * loc, const value * bw, const value * point, std::vector<bool>::const_iterator selection ) const {
    value p = 0.;
    
    if (*selection) {
        
        if (*bw>KAPPA_GAUSS_APPROX) {
            p = circular_difference( *point, *loc );
            p = -0.5 * (p*p) * (*bw);
        } else {
            p = *bw * std::cos( *point - *loc );
        }
        
    }
    
    return p;
}
    

Grid * CircularSpace::grid(unsigned int n) const {
    
    std::vector<value> v(n);
    for (unsigned int k=0; k<n; ++k) {
        v[k] = ( (2*M_PI*k) / n );
    }
    
    return new VectorGrid( { v }, specification(), {} );
    
}
